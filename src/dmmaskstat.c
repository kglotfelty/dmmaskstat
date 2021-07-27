/*                                                                
**  Copyright (C) 2004-2008, 2021  Smithsonian Astrophysical Observatory 
*/                                                                

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */

#include <math.h>
#include <float.h>
#include <limits.h>

#include <dslib.h>
#include <dsnan.h>
#include <histlib.h>
#include <cxcregion.h>

#include <dmfilters.h>
#include <dmimgio.h>


/* Input parameters */
typedef struct {
    char infile[DS_SZ_PATHNAME]; // Input file name
    char maskfile[DS_SZ_PATHNAME]; // Mask file name
    char outfile[DS_SZ_PATHNAME]; // Output file name
    short verbose;      // chatter level
    short clobber;      // rm outptu file if it exists?
} Parameters;


/* Information about mask pixel values */
typedef struct {
    long min_mask;
    long max_mask;
    double *vals_in_mask;
} StatsBuffer;

/* The statistics to be computed */
typedef struct {
    long mask_id;       // mask id
    long area;          // number of pixels in mask, may be diff than 'count' if image contains NaN/Null
    long perimeter;     // perimeter around the mask (include holes)
    double compact;     // compactness measure: perimeter^2/area (1=circle)
    double xcen;        // centroid in x-dir
    double ycen;        // centroid in y-dir
    double xavg;        // average x value
    double yavg;        // average y value
    double xmin;        // min x value (lower left x)
    double xmax;        // max x value (upper right x)
    double ymin;        // min y value (lower left y)
    double ymax;        // max y value (upper right y)
    double *more_stats; // stats computed from dmfilters routines, see LIST_OF_STATISTICS
} Stats;


/* Final output statistics */
typedef struct {    
    long num_masks;         // number of masks (not all values between min and max may be present)
    Stats **vals_per_mask;  // arrays of stats values
} MaskStats;



// ------------------------
// Function prototypes

int dmmaskstat();
int convert_coords( Image *image, double x_in, double y_in, double *x_out, double *y_out);
Parameters *get_parameters();
StatsBuffer *get_mask_range(Image *mask);
double *compute_more_stats(double *vals, long nvals);
Stats *get_mask_stats(Image *image, Image *mask, StatsBuffer *sbuff, long mask_id);
MaskStats *loop_over_masks( Image *image, Image *mask, StatsBuffer *sbuff);
int write_output(char *outfile, MaskStats *mask_stats, Image *image);

// ---------

/*
 * These stats come from the dmfilters library.
 */
#define LIST_OF_STATISTICS   {"min", "max", "mean", "count", \
                              "sum", "median", "mode", "mid",  \
                              "sigma", "extreme", "range", "q25", \
                              "q33", "q67", "q75", "olympic", NULL }

// ----------------------
// Functions


/*
 *  Convert coordinates from C-index to Physical coordinates
 */
int convert_coords( Image *image, double x_in, double y_in, double *x_out, double *y_out)
{
    if (image->xdesc == NULL) {
        *x_out = x_in;
        *y_out = y_in;
        return(0);
    }

    double logical[2];
    double physical[2];

    // Convert from 0-based array index to physical coords
    logical[0] = x_in + 1;
    logical[1] = y_in + 1;
    dmCoordCalc_d(image->xdesc, logical, physical);
    if (image->ydesc){
        dmCoordCalc_d(image->ydesc, logical+1, physical+1);
    }
    *x_out = physical[0];
    *y_out = physical[1];

    return(0);
}

/*
 * Load parameters from .par file
 */
Parameters *get_parameters()
{
    
    Parameters *pars;
    if (NULL == (pars = calloc(1,sizeof(Parameters)))) {
        err_msg("ERROR: problem allocating memory");
        return(NULL);
    }

    clgetstr("infile", pars->infile, DS_SZ_PATHNAME);
    clgetstr("maskfile", pars->maskfile, DS_SZ_PATHNAME);
    clgetstr("outfile", pars->outfile, DS_SZ_PATHNAME);
    pars->verbose = clgeti("verbose");
    pars->clobber = clgetb("clobber");

    if (pars->verbose >= 1){
        printf("dmmaskstat\n");
        printf("%15s = %-s\n", "infile", pars->infile);
        printf("%15s = %-s\n", "maskfile", pars->maskfile);
        printf("%15s = %-s\n", "outfile", pars->outfile);
        printf("%15s = %d\n", "verbose", pars->verbose);
        printf("%15s = %-s\n", "clobber", pars->clobber ? "yes" : "no");
    }

    return(pars);        
}


/*
 * Setup the StatsBuff.  We determine the range of mask values
 * in the input image and setup a buffer to be used by the 
 * dmfilters routines.   We set the buffer to the same size
 * as the image -- guaranteed there won't be more pixels than that.
 */
StatsBuffer *get_mask_range(Image *mask)
{

    StatsBuffer *sbuff;
    
    if ( NULL == (sbuff = calloc(1,sizeof(StatsBuffer)))) {
        err_msg("ERROR: mem alloc failed");
        return NULL;
    }

    sbuff->min_mask = LONG_MAX;
    sbuff->max_mask = LONG_MIN;
  
    long xx,yy;

    for(yy=mask->lAxes[1];yy--;) {
        for (xx=mask->lAxes[0];xx--; ) {
            double  mval;
            mval = get_image_value(mask, xx, yy);
            if ( ds_dNAN(mval) ) continue;

            sbuff->min_mask = MIN( sbuff->min_mask, mval );
            sbuff->max_mask = MAX( sbuff->max_mask, mval );
        }
    }


    if (NULL == (sbuff->vals_in_mask = (double*)calloc(mask->lAxes[0]*mask->lAxes[1], sizeof(double)))) {
        err_msg("ERROR: mem alloc failed");
        return NULL;
    }
  
    return sbuff;
}


/*
 * Compute the statistics from the dmfilters library.
 */
double *compute_more_stats(double *vals, long nvals)
{
    // Compute array of statistics for each grid
    char *list_of_stats[] = LIST_OF_STATISTICS;

    // Count the number of statistcs
    int num_stats = 0;
    while(list_of_stats[++num_stats]) {};

    // +1 so has NULL to terminate list
    double *stats = NULL;
    if (NULL == (stats = calloc(num_stats+1,sizeof(double)))) {
        err_msg("ERROR: Problem allocating memory");
        return(NULL);
    }

    int ii;
    for (ii=0;ii<num_stats;ii++) {
        // dmfilter.h
        double (*func)( double *vals, long nvals ) = NULL;
        if (NULL == (func = get_method(list_of_stats[ii]))) {
            err_msg("ERROR: Unknown filter %s", list_of_stats[ii]);
            return(NULL);
        }

        double val;
        val = func(vals, nvals);

        // Special case -- if 0 counts, set to 0 not NaN
        if ((strcmp(list_of_stats[ii], "count") == 0) && (ds_dNAN(val))) {
            val = 0;
        }

        stats[ii] = val;

    } // end for ii


    return(stats);
}

/*
 * Compute the statistics for the current mask_id value.
 * 
 * it computes some stats based on the mask (perimeter,
 */
Stats *get_mask_stats(Image *image, Image *mask, StatsBuffer *sbuff, long mask_id)
{

    Stats *retvals;
    if ( NULL == ( retvals = calloc(1,sizeof(Stats)))) {
        err_msg("ERROR: alloc failed");
        return NULL;
    }

    retvals->mask_id = mask_id;

    long nvals = 0;
    double sum_up = 0;
    retvals->area = 0;
    retvals->xavg = 0;
    retvals->yavg = 0;
    retvals->xcen = 0;
    retvals->ycen = 0;
    retvals->perimeter=0;

    retvals->xmax = retvals->ymax = -1;
    retvals->xmin = image->lAxes[0]+1;
    retvals->ymin = image->lAxes[1]+1;

    long xx,yy;
    for(yy=image->lAxes[1];yy--;) {
      double last_mval;
      last_mval=sbuff->min_mask-1;

      for (xx=image->lAxes[0];xx--; ) {
        double  mval;
        mval = get_image_value(mask, xx, yy);

        if ((mval != mask_id) && (last_mval == mask_id)) retvals->perimeter++;
        if ((mval == mask_id) && (last_mval != mask_id)) retvals->perimeter++;        
        last_mval = mval;
        if (mval != mask_id) {
            continue;
        }

        if ( xx==0 ) retvals->perimeter+=1; /* if at edge then count +1 */
        retvals->area++;
        
        retvals->xavg += xx; 
        retvals->yavg += yy;
        retvals->xmin = MIN(xx,retvals->xmin);
        retvals->xmax = MAX(xx,retvals->xmax);
        retvals->ymin = MIN(yy,retvals->ymin);
        retvals->ymax = MAX(yy,retvals->ymax);

        double val;
        val = get_image_value( image, xx, yy);
        if ( ds_dNAN(val) ) {
            continue;
        }

        sbuff->vals_in_mask[nvals] = val;
        nvals++;

        val = fabs(val);
        retvals->xcen += (xx*val); // use magnitude to compute centroid
        retvals->ycen += (yy*val);
        sum_up += val;

                
      } /* end xx */
    } /* end yy */

    if ( 0 == nvals ) {
        return NULL;
    }

    retvals->more_stats = compute_more_stats( sbuff->vals_in_mask, nvals);

    retvals->xavg /= nvals; 
    retvals->yavg /= nvals;
    retvals->xcen /= sum_up;
    retvals->ycen /= sum_up;
    
    /* Similar to above, but do x then y */
    /* We only need to work inside of min/max range found above */
    for (xx=retvals->xmax;xx>=retvals->xmin;xx--) {
      double last_mval;
      last_mval=sbuff->min_mask-1;

      for(yy=retvals->ymax;yy>=retvals->ymin;yy--) {

        double mval;

        mval = get_image_value(mask, xx, yy);

        if ((mval != mask_id) && (last_mval == mask_id)) retvals->perimeter++;
        if ((mval == mask_id) && (last_mval != mask_id)) retvals->perimeter++;        
        last_mval = mval;
        if (mval != mask_id) {
            continue;
        }
        // If we are at the edge of image and current value is mask
        // then add perimeter
        if ( yy==retvals->ymin) retvals->perimeter+=1;
        
      } /* end xx */
    } /* end yy */

    retvals->compact = ((double)retvals->perimeter*retvals->perimeter/retvals->area);

    return retvals;
    
}

/*
 *  Driver script to loop over mask_id values.
 * 
 *  Keeps track of how many mask_id's have values and stores
 *  values for only those that do.
 */
#define BUFFER_INC 10

MaskStats *loop_over_masks( Image *image, Image *mask, StatsBuffer *sbuff)
{

    MaskStats *retvals;
    
    if (NULL == (retvals = calloc(1,sizeof(MaskStats)))) {
        err_msg("ERROR: mem alloc failed");
        return NULL;
    }
    
    retvals->num_masks = 0;
    retvals->vals_per_mask = NULL;

    long mask_id;    
    for (mask_id = sbuff->min_mask; mask_id <= sbuff->max_mask; mask_id++ ) {

        Stats *mask_stat;
        if (NULL == (mask_stat = get_mask_stats(image, mask, sbuff, mask_id))) {
            continue;
        }


        if ((retvals->num_masks % BUFFER_INC) == 0) {
            if (NULL == retvals->vals_per_mask) {
                retvals->vals_per_mask = (Stats**)calloc(BUFFER_INC,sizeof(Stats*));
            } else {
                retvals->vals_per_mask = (Stats**)realloc(retvals->vals_per_mask, 
                                            (retvals->num_masks+BUFFER_INC)*sizeof(Stats*));
            }
        }
        retvals->vals_per_mask[retvals->num_masks] = mask_stat;
        retvals->num_masks++;        
        
    } // end for mask_id
 
    return retvals;
}

/*
 * Save values to output table.
 * 
 */
int write_output(char *outfile, MaskStats *mask_stats, Image *image)
{

    dmBlock *outBlock;
    if (NULL == (outBlock = dmTableCreate( outfile ) )) {
        err_msg("ERROR: Cannot create output file '%s'\n", outfile );
        return(1);
    }

    char units[100];
    memset( units, 0, sizeof(char)*100);
    dmGetUnit( dmImageGetDataDescriptor( image->block ), units, 99 );

    /* All the header stuff */
    Header_Type *hdr;
    hdr = getHdr( image->block, hdrDM_FILE );
    putHdr( outBlock, hdrDM_FILE, hdr, BASIC_STS, "dmmaskstat" );
    put_param_hist_info(outBlock, "dmmaskstat", NULL, 0);

    if (dmBlockGetNo(outBlock) != 1) {
        dmDataset *outDs = dmBlockGetDataset(outBlock);
        dmBlock *primary = dmDatasetMoveToBlock(outDs, 1);
        putHdr(primary, hdrDM_FILE, hdr, PRIMARY_STS, "dmmaskstat");
    }


    // Write data
    // TODO Vector column and copy WCS from image to x,y cols
    dmDescriptor *mask_id_col = dmColumnCreate(outBlock, "mask", dmLONG, 0, "", "Mask ID number");
    dmDescriptor *area_col = dmColumnCreate(outBlock, "area", dmLONG, 0, "pixel", "Number of pixels in mask");
    dmDescriptor *perimeter_col = dmColumnCreate(outBlock, "perimeter", dmLONG, 0, "pixel", "Perimeter length");
    dmDescriptor *compact_col = dmColumnCreate(outBlock, "compactness", dmDOUBLE, 0, "", "Measure of roundness: perimeter^2/area");
    dmDescriptor *xcen_col = dmColumnCreate(outBlock, "x_centroid", dmDOUBLE, 0, "pixel", "centroid in x");
    dmDescriptor *ycen_col = dmColumnCreate(outBlock, "y_centroid", dmDOUBLE, 0, "pixel", "centroid in y");
    dmDescriptor *xavg_col = dmColumnCreate(outBlock, "x_average", dmDOUBLE, 0, "pixel", "average x value");
    dmDescriptor *yavg_col = dmColumnCreate(outBlock, "y_average", dmDOUBLE, 0, "pixel", "average y value");
    dmDescriptor *xmin_col = dmColumnCreate(outBlock, "x_min", dmDOUBLE, 0, "pixel", "minimum x value");
    dmDescriptor *xmax_col = dmColumnCreate(outBlock, "x_max", dmDOUBLE, 0, "pixel", "maximum x value");
    dmDescriptor *ymin_col = dmColumnCreate(outBlock, "y_min", dmDOUBLE, 0, "pixel", "minimum y value");
    dmDescriptor *ymax_col = dmColumnCreate(outBlock, "y_max", dmDOUBLE, 0, "pixel", "maximum y value");

    // Compute array of statistics for each grid
    char *list_of_stats[] = LIST_OF_STATISTICS;

    // Count the number of statistcs
    int num_stats = 0;
    while(list_of_stats[++num_stats]) {};

    dmDescriptor **more_stats_col = (dmDescriptor **)calloc(num_stats,sizeof(dmDescriptor *));
    long ii;
    for (ii=0;ii<num_stats;ii++) {
        dmDataType dt;
        dt = strcmp(list_of_stats[ii], "count") ? dmDOUBLE : dmLONG;
        more_stats_col[ii] = dmColumnCreate(outBlock, list_of_stats[ii], dt, 0, units, "");
    }

    // Loop over mask ids
    for (ii=0;ii< mask_stats->num_masks;ii++) {
        Stats *at = mask_stats->vals_per_mask[ii];
        double xx,yy;
        long jj;

        dmSetScalar_l( mask_id_col, at->mask_id);
        dmSetScalar_l( area_col, at->area);
        dmSetScalar_l( perimeter_col, at->perimeter);
        dmSetScalar_d( compact_col, at->compact);
        
        convert_coords( image, at->xavg, at->yavg, &xx, &yy);
        dmSetScalar_d( xavg_col, xx);
        dmSetScalar_d( yavg_col, yy);
        
        convert_coords( image, at->xcen, at->ycen, &xx, &yy);
        dmSetScalar_d( xcen_col, xx);
        dmSetScalar_d( ycen_col, yy);

        convert_coords( image, at->xmin, at->ymin, &xx, &yy);
        dmSetScalar_d( xmin_col, xx);
        dmSetScalar_d( ymin_col, yy);

        convert_coords( image, at->xmax, at->ymax, &xx, &yy);
        dmSetScalar_d( xmax_col, xx);
        dmSetScalar_d( ymax_col, yy);
        
        for (jj=0;jj<num_stats; jj++) {
            dmSetScalar_d( more_stats_col[jj], at->more_stats[jj]);
        }

        dmTablePutRow(outBlock, NULL);        
    }

    dmTableClose(outBlock);

    return 0;
}

/*
 * Main routine.
 */
int dmmaskstat()
{

    Parameters *pars;
    if (NULL == (pars = get_parameters())) {
        return(-9);
    }

    if ( ds_clobber( pars->outfile, pars->clobber, NULL ) != 0 ) {
        return(-1);
    }

    Image *image;
    if (NULL == (image = load_image(pars->infile))) {
        return(-1);
    }

    Image *mask;
    if (NULL == (mask = load_image(pars->maskfile))) {
        return(-1);
    }

    if ( ( dmFLOAT == mask->dt) || ( dmDOUBLE==mask->dt)) {
        err_msg("ERROR: Mask file must be integer data type\n");
        return(-1);
    }

    if ( (mask->lAxes[0]!=image->lAxes[0]) || (mask->lAxes[1]!=image->lAxes[1])) {
        err_msg("ERROR: infile and maskfile do not have same image dimensions\n");
        return(-1);
    }
  
    StatsBuffer *sbuff;
    if ( NULL == (sbuff = get_mask_range( mask ) )) {
        return(-1);
    }

    MaskStats *mask_stats;
    if ( NULL == (mask_stats = loop_over_masks(image, mask, sbuff))) {
        return(-1);
    }

    if ( 0 != write_output(pars->outfile, mask_stats, image)) {
        return(-1);
    }

    return(0);

}
