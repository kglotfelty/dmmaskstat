/*                                                                
**  Copyright (C) 2004-2008  Smithsonian Astrophysical Observatory 
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

#include <dslib.h>
#include <dsnan.h>
#include <histlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <cxcregion.h>
#include <dmimgfilt.h>

double get_image_value( void *data, dmDataType dt, long *lAxes, 
                        long xx, long yy, regRegion *dss, 
                        dmDescriptor *xAxis, dmDescriptor *yAxis,
                        long nullval, short has_null )
{

  long npix = xx + (yy * lAxes[0] );
  double retval;

  /* Okay, first get all the data from the different data types.  
     Cast everything to doubles */

  if ( ( xx < 0 ) || ( xx >= lAxes[0] ) ||
       ( yy < 0 ) || ( yy >= lAxes[1] ) ) {
    ds_MAKE_DNAN( retval );
    return(retval);
  }


  switch ( dt ) {
    
  case dmBYTE: {
    unsigned char *img = (unsigned char*)data;
    retval = img[npix];
    break;
  }
    
  case dmSHORT: {
    short *img = (short*)data;
    retval = img[npix];
    break;
  }
    
  case dmUSHORT: {
    unsigned short *img = (unsigned short*)data;
    retval = img[npix];
    break;
  }
    
  case dmLONG: {
    long *img = (long*)data;
    retval = img[npix];
    break;
  }
    
  case dmULONG: {
    unsigned long *img = (unsigned long*)data;
    retval = img[npix];
    break;
  }
    
  case dmFLOAT: {
    float *img = (float*)data;
    retval = img[npix];
    break;
  }
  case dmDOUBLE: {
    double *img = (double*)data;
    retval = img[npix];
    break;
  }
  default:
    ds_MAKE_DNAN( retval );

  }

  
  /* Now ... if it is an integer data type, it could possibly have a
     null value. Check for that */

  if ( has_null && ( retval == nullval ) ) {
    ds_MAKE_DNAN( retval );
    return(retval);
  }
  /* If the image has a data sub space (aka a region filter applied)
     then need to convert coords to physical and check */
  if ( dss && xAxis ) {
    double pos[2];
    double loc[2];
    pos[0]=xx+1;
    pos[1]=yy+1;

    if (yAxis) {  /* If no y axis, then xAxis has 2 components */
      dmCoordCalc_d( xAxis, pos, loc );
      dmCoordCalc_d( yAxis, pos+1, loc+1 );
    } else {
      dmCoordCalc_d( xAxis, pos, loc );
    }
    if ( !regInsideRegion( dss, loc[0], loc[1] ) )
      ds_MAKE_DNAN( retval );
  }

  return(retval);

}



/* Load the data into memory,  check for DSS, null values */
dmDataType get_image_data( dmBlock *inBlock, void **data, long **lAxes,
                           regRegion **dss, long *nullval, short *nullset )
{

  dmDescriptor *imgDesc;
  dmDataType dt;
  dmDescriptor *grp;
  dmDescriptor *imgdss;

  long naxes;
  long npix;
  char ems[1000];

  *nullval = INDEFL;
  *dss = NULL;
  *nullset = 0;
  
  imgDesc = dmImageGetDataDescriptor( inBlock );

  /* Sanity check, only 2D images */
  naxes = dmGetArrayDimensions( imgDesc, lAxes );
  if ( naxes != 2 ) {
    return( dmUNKNOWNTYPE );
  }
  npix = (*lAxes)[0] * (*lAxes)[1];
  dt = dmGetDataType( imgDesc );


  /* Okay, first lets get the image descriptor */
  grp = dmArrayGetAxisGroup( imgDesc, 1 );
  dmGetName( grp, ems, 1000);
  imgdss = dmSubspaceColOpen( inBlock, ems );
  if ( imgdss )
    *dss = dmSubspaceColGetRegion( imgdss);
  
  
  switch ( dt ) 
    {
    case dmBYTE:
      *data = ( void *)calloc( npix, sizeof(char ));
      dmGetArray_ub( imgDesc, (unsigned char*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmSHORT:
      *data = ( void *)calloc( npix, sizeof(short ));
      dmGetArray_s( imgDesc, (short*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmUSHORT:
      *data = ( void *)calloc( npix, sizeof(short ));
      dmGetArray_us( imgDesc, (unsigned short*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmLONG:
      *data = ( void *)calloc( npix, sizeof(long ));
      dmGetArray_l( imgDesc, (long*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmULONG:
      *data = ( void *)calloc( npix, sizeof(long ));
      dmGetArray_ul( imgDesc, (unsigned long*) *data, npix );
      if ( dmDescriptorGetNull_l( imgDesc, nullval) == 0 ) {
        *nullset=0;
      } else
        *nullset=1;
      break;
      
    case dmFLOAT:
      *data = ( void *)calloc( npix, sizeof(float ));
      dmGetArray_f( imgDesc, (float*) *data, npix );
      *nullset = 0;
      break;
      
    case dmDOUBLE:
      *data = ( void *)calloc( npix, sizeof(double ));
      dmGetArray_d( imgDesc, (double*) *data, npix );
      *nullset = 0;
      break;
      
    default:
      return( dmUNKNOWNTYPE );
    }

  return(dt);

}


/* Get the WCS descriptor */
short  get_image_wcs( dmBlock *imgBlock, dmDescriptor **xAxis, 
                      dmDescriptor **yAxis )
{
  

  dmDescriptor *imgData;
  long n_axis_groups;

  imgData = dmImageGetDataDescriptor( imgBlock );
  n_axis_groups = dmArrayGetNoAxisGroups( imgData );
  

  /* This is the usual trick ... can have 1 axis group w/ 
     dimensionality 2 (eg a vector column) or can have
     2 axis groups w/ dimensionaity 1 (eg 2 disjoint columns)*/

  if ( n_axis_groups == 1 ) {
    dmDescriptor *pos = dmArrayGetAxisGroup( imgData, 1 );
    dmDescriptor *xcol;
    long n_components;
    
    n_components = dmGetElementDim( pos );
    if ( n_components != 2 ) {
      err_msg("ERROR: could not find 2D image\n");
      return(-1);
    }
    
    xcol = dmGetCpt( pos, 1 );
    
    *xAxis = pos;
    *yAxis = NULL;
    
  } else if ( n_axis_groups == 2 ) {
    dmDescriptor *xcol;
    dmDescriptor *ycol;
  
    xcol = dmArrayGetAxisGroup( imgData, 1 );
    ycol = dmArrayGetAxisGroup( imgData, 2 );

    *xAxis = xcol;
    *yAxis = ycol;
    
  } else {
    err_msg("Invalid number of axis groups\n");
    *xAxis = NULL;
    *yAxis = NULL;
    return(-1);
  }

  return(0);

}




int dmmaskstat()
{

  char infile[DS_SZ_PATHNAME];
  char maskfile[DS_SZ_PATHNAME];
  char outfile[DS_SZ_PATHNAME];
  short clobber;
  short verbose;

  void *data;
  long *lAxes;
  regRegion *dss;
  long null;
  short has_null;
  dmDataType dt;
  dmBlock *inBlock;
  dmDescriptor *xdesc, *ydesc;
  
  void *mask_data;
  long *mask_lAxes;
  regRegion *mask_dss;
  long mask_null;
  short mask_has_null;
  dmDataType mask_dt;
  dmBlock *mask_inBlock;
  dmDescriptor *mask_xdesc, *mask_ydesc;
  long min_mask, max_mask;
  

  double *vals;
  long nvals;

  dmBlock *outBlock;

  typedef struct {
    long mask_vals;
    double min;
    double max;
    double mean;
    double median;
    double mode;
    double sigma;
    double sum;
    long count;
    long perimeter;
    double compact;
    double range;
    double q25;
    double q33;
    double q67;
    double q75;
    double xcen;
    double ycen;
    double xavg;
    double yavg;


  } Stats; 

  Stats *stats;

  long xx, yy,ii;



  /* Get the parameters */
  clgetstr( "infile", infile, DS_SZ_FNAME );
  clgetstr( "maskfile", maskfile, DS_SZ_FNAME );
  clgetstr( "outfile", outfile, DS_SZ_FNAME );
  clobber = clgetb( "clobber" );
  verbose = clgeti( "verbose" );
  


  if ( ds_clobber( outfile, clobber, NULL ) != 0 ) {
    return(-1);
  }


  if ( NULL == ( inBlock = dmImageOpen( infile) ) ) {
    err_msg("ERROR: Cannot open image '%s'\n", infile );
    return(-1);
  }

  if ( dmUNKNOWNTYPE == ( dt = get_image_data( inBlock, &data,  &lAxes, 
                                               &dss, &null, &has_null ) ) ) {
    err_msg("ERROR: Cannot get image data or unknown image data-type for "
            "file '%s'\n", infile);
    return(-1);
  }

  if ( 0 != get_image_wcs( inBlock, &xdesc, &ydesc ) ) {
    err_msg("ERROR: Cannot load WCS for file '%s'\n", infile );
    return(-1);
  }


  /* cut-n-paste of above, for mask file */
  if ( NULL == ( mask_inBlock = dmImageOpen( maskfile) ) ) {
    err_msg("ERROR: Cannot open image '%s'\n", maskfile );
    return(-1);
  }

  if (dmUNKNOWNTYPE==(mask_dt=get_image_data(mask_inBlock, &mask_data,  
					     &mask_lAxes, &mask_dss, 
					     &mask_null, &mask_has_null ))) {
    err_msg("ERROR: Cannot get image data or unknown image data-type for "
            "file '%s'\n", maskfile);
    return(-1);
  }

  if ( ( dmFLOAT == mask_dt) || ( dmDOUBLE==mask_dt)) {
    err_msg("ERROR: Mask file must be integer data type\n");
    return(-1);
  }
  
  if ( (mask_lAxes[0]!=lAxes[0]) || (mask_lAxes[1]!=lAxes[1])) {
    err_msg("ERROR: infile and maskfile do not have same image dimensions\n");
    return(-1);
  }
  

  if ( 0 != get_image_wcs( mask_inBlock, &mask_xdesc, &mask_ydesc ) ) {
    err_msg("ERROR: Cannot load WCS for file '%s'\n", maskfile );
    return(-1);
  }


  min_mask = LONG_MAX;
  max_mask = LONG_MIN;

  

  for(yy=lAxes[1];yy--;) {
    for (xx=lAxes[0];xx--; ) {
      double  mval;
      mval = get_image_value( mask_data, mask_dt, mask_lAxes, xx, yy,
			      mask_dss, mask_xdesc, mask_ydesc, mask_null, 
			      mask_has_null );
      if ( ds_dNAN(mval) ) continue;

      min_mask = MIN( min_mask, mval );
      max_mask = MAX( max_mask, mval );

    }
  }

  vals = (double*)calloc(lAxes[0]*lAxes[1], sizeof(double));
  
  stats=(Stats*)calloc((max_mask-min_mask+1), sizeof(Stats));


  
  for (ii=min_mask; ii<=max_mask;ii++) {
    double xavg, yavg;
    double xcen, ycen;
    long perimeter;
    nvals = 0;
    xavg = 0;
    yavg = 0;
    xcen = 0;
    ycen = 0;
    perimeter=0;

    for(yy=lAxes[1];yy--;) {
      double last_mval;
      last_mval=min_mask-1;

      for (xx=lAxes[0];xx--; ) {
	double val;
	double  mval;
	mval = get_image_value( mask_data, mask_dt, mask_lAxes, xx, yy,
				mask_dss, mask_xdesc, mask_ydesc, mask_null, 
				mask_has_null );
	if ( ds_dNAN(mval) ) { 
	  if ( last_mval == ii ) perimeter++;
	  last_mval=mval; 
	  continue; 
	}
	if ( mval != ii ) { 
	  if ( last_mval == ii ) perimeter++;
	  last_mval=mval; 
	  continue; 
	}
	
	val = get_image_value( data, dt, lAxes, xx, yy, dss, xdesc,
			       ydesc, null, has_null );
	if ( last_mval != mval ) {
	  perimeter+=1;
	}
	if ( xx==0 ) perimeter+=1; /* if at edge then count +1 */
	
	if ( ds_dNAN(val) ) { last_mval=mval; continue;}

	xavg += xx; 
	yavg += yy;
	xcen += (xx*val);
	ycen += (yy*val);
	vals[nvals] = val;
	nvals++;
	last_mval = mval;
		
      } /* end xx */
    } /* end yy */

    
    /* Similar to above, but do x then y */
    for (xx=lAxes[0];xx--; ) {
      double last_mval;
      last_mval=min_mask-1;

      for(yy=lAxes[1];yy--;) {

	double val;
	double  mval;

	mval = get_image_value( mask_data, mask_dt, mask_lAxes, xx, yy,
				mask_dss, mask_xdesc, mask_ydesc, mask_null, 
				mask_has_null );
	if ( ds_dNAN(mval) ) { 
	  if ( last_mval == ii ) perimeter++;
	  last_mval=mval; 
	  continue; 
	}
	if ( mval != ii ) { 
	  if ( last_mval == ii ) perimeter++;
	  last_mval=mval; 
	  continue; 
	}
	
	val = get_image_value( data, dt, lAxes, xx, yy, dss, xdesc,
			       ydesc, null, has_null );
	if ( last_mval != mval ) {
	  perimeter+=1;
	}
	if ( yy==0 ) perimeter+=1;
	
	if ( ds_dNAN(val) ) { last_mval=mval; continue;}

	last_mval = mval;

      } /* end xx */
    } /* end yy */


    if ( 0 == nvals ) {
      ds_MAKE_DNAN( xavg );
      ds_MAKE_DNAN( yavg );
      ds_MAKE_DNAN( xcen );
      ds_MAKE_DNAN( ycen );
    }
    
    stats[ii-min_mask].mask_vals = ii;
    stats[ii-min_mask].min = _filtMIN( vals, nvals );
    stats[ii-min_mask].max = _filtMAX( vals, nvals );
    stats[ii-min_mask].mean = _filtMEAN( vals, nvals );
    stats[ii-min_mask].median = _filtMEDIAN( vals, nvals );
    stats[ii-min_mask].mode = _filtMODE( vals, nvals );
    stats[ii-min_mask].sigma = _filtSIG( vals, nvals );
    stats[ii-min_mask].sum = _filtSUM( vals, nvals );
    stats[ii-min_mask].count = _filtCOUNT( vals, nvals );
    stats[ii-min_mask].range = _filtRANGE( vals, nvals );
    stats[ii-min_mask].q25 = _filtQUANTILE_25( vals, nvals );
    stats[ii-min_mask].q33 = _filtQUANTILE_33( vals, nvals );
    stats[ii-min_mask].q67 = _filtQUANTILE_67( vals, nvals );
    stats[ii-min_mask].q75 = _filtQUANTILE_75( vals, nvals );
    stats[ii-min_mask].xcen = xcen / stats[ii-min_mask].sum ;
    stats[ii-min_mask].ycen = ycen / stats[ii-min_mask].sum ;
    stats[ii-min_mask].xavg = xavg / stats[ii-min_mask].count ;
    stats[ii-min_mask].yavg = yavg / stats[ii-min_mask].count ;
    stats[ii-min_mask].perimeter = perimeter;
    stats[ii-min_mask].compact = ((double)perimeter) * perimeter / 
      ((double)nvals);
    

  } /* end for ii, loop over mask values */



  if ( NULL == ( outBlock = dmTableCreate( outfile ))){
    err_msg("ERROR: Cannot create output table '%s'\n", outfile );
    return(-1);
  } else {
    dmDescriptor *mask_col, *min_col, *max_col, *mean_col, *median_col;
    dmDescriptor *mode_col, *sig_col, *sum_col, *count_col, *range_col;
    dmDescriptor *q25_col,*q33_col,*q67_col,*q75_col;
    dmDescriptor *avg_col, *cen_col, *per_col, *com_col;

    char *avg[] = { "X_AVGERAGE", "Y_AVGERAGE" };
    char *cen[] = { "X_CENTROID", "Y_CENTROID"};

    Header_Type *hdr;


    mask_col=dmColumnCreate( outBlock, "MASK", dmLONG, 0, NULL, "mask value");
    min_col=dmColumnCreate( outBlock, "MIN", dmDOUBLE, 0, NULL, "min value");
    max_col=dmColumnCreate( outBlock, "MAX", dmDOUBLE, 0, NULL, "max value");
    mean_col=dmColumnCreate( outBlock, "MEAN", dmDOUBLE, 0, NULL, "mean value");
    median_col=dmColumnCreate( outBlock, "MEDIAN", dmDOUBLE, 0, NULL, "median value");
    mode_col=dmColumnCreate( outBlock, "MODE", dmDOUBLE, 0, NULL, "mode value");
    sig_col=dmColumnCreate( outBlock, "SIGMA", dmDOUBLE, 0, NULL, "sigma value");
    sum_col=dmColumnCreate( outBlock, "SUM", dmDOUBLE, 0, NULL, "sum of pixels");
    count_col=dmColumnCreate( outBlock, "COUNT", dmLONG, 0, NULL, "no. pixels");
    range_col=dmColumnCreate( outBlock, "RANGE", dmDOUBLE, 0, NULL, "max-min");
    per_col=dmColumnCreate( outBlock, "PERIMETER", dmLONG, 0, NULL, "length of perimeter");
    com_col=dmColumnCreate( outBlock, "COMPACTNESS", dmDOUBLE, 0, NULL,
			    "perimeter**2/area" );

    q25_col=dmColumnCreate( outBlock, "QUANT_25", dmDOUBLE, 0, NULL, "25%");
    q33_col=dmColumnCreate( outBlock, "QUANT_33", dmDOUBLE, 0, NULL, "33%");
    q67_col=dmColumnCreate( outBlock, "QUANT_67", dmDOUBLE, 0, NULL, "67%");
    q75_col=dmColumnCreate( outBlock, "QUANT_75", dmDOUBLE, 0, NULL, "75%");


    avg_col=dmColumnCreateVector( outBlock, "AVERAGE", dmDOUBLE, 0, NULL, 
				  "average location", avg, 2);
    cen_col=dmColumnCreateVector( outBlock, "CENTROID", dmDOUBLE, 0, NULL, 
				  "centroid location", cen, 2);

    

    hdr = getHdr( inBlock, hdrDM_FILE );
    putHdr( outBlock, hdrDM_FILE, hdr, BASIC_STS, "dmmaskstat");

    put_param_hist_info( outBlock, "dmmaskstat", NULL, 0 );

    for (ii=min_mask; ii<=max_mask;ii++) {
      double dvals[2];

      dmSetScalar_l(mask_col, stats[ii-min_mask].mask_vals);
      dmSetScalar_d( min_col,stats[ii-min_mask].min);
      dmSetScalar_d( max_col,stats[ii-min_mask].max);
      dmSetScalar_d( mean_col,stats[ii-min_mask].mean);
      dmSetScalar_d( median_col,stats[ii-min_mask].median);
      dmSetScalar_d( mode_col,stats[ii-min_mask].mode);
      dmSetScalar_d( sig_col,stats[ii-min_mask].sigma);
      dmSetScalar_d( sum_col,stats[ii-min_mask].sum);
      dmSetScalar_d( count_col,stats[ii-min_mask].count);
      dmSetScalar_d( range_col,stats[ii-min_mask].range);
      dmSetScalar_d( q25_col,stats[ii-min_mask].q25);
      dmSetScalar_d( q33_col,stats[ii-min_mask].q33);
      dmSetScalar_d( q67_col, stats[ii-min_mask].q67);
      dmSetScalar_d( q75_col, stats[ii-min_mask].q75);
      dmSetScalar_l( per_col, stats[ii-min_mask].perimeter);
      dmSetScalar_d( com_col, stats[ii-min_mask].compact);


      dvals[0] = stats[ii-min_mask].xavg+1;
      dvals[1] = stats[ii-min_mask].yavg+1;
      if ( xdesc ) {
	double loc[2];
	if (ydesc) {  /* If no y desc, then xDesc has 2 components */
	  dmCoordCalc_d( xdesc, dvals, loc );
	  dmCoordCalc_d( ydesc, dvals+1, loc+1 );
	} else {
	  dmCoordCalc_d( xdesc, dvals, loc );
	}
	dvals[0] =loc[0];
	dvals[1] = loc[1];
      }
      dmSetVector_d( avg_col, dvals, 2 );

      dvals[0] = stats[ii-min_mask].xcen+1;
      dvals[1] = stats[ii-min_mask].ycen+1;
      if ( xdesc ) {
	double loc[2];
	if (ydesc) {  /* If no y desc, then xDesc has 2 components */
	  dmCoordCalc_d( xdesc, dvals, loc );
	  dmCoordCalc_d( ydesc, dvals+1, loc+1 );
	} else {
	  dmCoordCalc_d( xdesc, dvals, loc );
	}
	dvals[0] =loc[0];
	dvals[1] = loc[1];
      }
      dmSetVector_d( cen_col, dvals, 2 );


      /*
	dmSetScalar_d( foo_col, stats[ii-min_mask].xcen);
	dmSetScalar_d( foo_col, stats[ii-min_mask].ycen);
	dmSetScalar_d( foo_col, stats[ii-min_mask].xavg);
	dmSetScalar_d( foo_col, stats[ii-min_mask].yavg);
      */

      dmTablePutRow( outBlock, NULL);
    }
    
    dmTableClose(outBlock );

  }
  dmImageClose( inBlock );
  return(0);


}
