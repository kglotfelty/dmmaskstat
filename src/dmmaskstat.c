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
#include <stdlib.h>
#include <string.h>

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

typedef struct {
    double convex_area;         // area enclosed by convex hull
    double convex_perimeter;    // perimeter of convex hull
    double feret_diameter;      // maximum caliper distance
    double min_feret_diameter;  // minimum caliper distance
    double max_inscribed_circle_diameter;
    double eq_circle_diameter;
    double mbr_long_side;
    double mbr_short_side;
    double aspect_ratio;
    double area_perimeter_ratio;
    double circularity;
    double elongation;
    double convexity;
    double solidity;
    long num_holes;
    double thinness_ratio;
    double contour_temperature;
    double orientation;
    double major_axis;
    double minor_axis;
    double fractal_box_dimension;
    long holes_area;
    double porosity;
} ShapeStats;

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
    ShapeStats shape;   // shape statistics
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
int compute_shape_stats(Image *mask, long mask_id, Stats *stats);
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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

typedef struct {
    double x;
    double y;
} Point;

typedef struct {
    double x1;
    double y1;
    double x2;
    double y2;
} Segment;

static double invalid_double(void)
{
    return sqrt(-1.0);
}

static int compare_points(const void *aa, const void *bb)
{
    const Point *a = (const Point *)aa;
    const Point *b = (const Point *)bb;

    if (a->x < b->x) return -1;
    if (a->x > b->x) return 1;
    if (a->y < b->y) return -1;
    if (a->y > b->y) return 1;
    return 0;
}

static double cross(Point oo, Point aa, Point bb)
{
    return (aa.x - oo.x) * (bb.y - oo.y) - (aa.y - oo.y) * (bb.x - oo.x);
}

static long convex_hull(Point *points, long npoints, Point *hull)
{
    long ii, kk = 0;

    if (npoints <= 1) {
        if (npoints == 1) hull[0] = points[0];
        return npoints;
    }

    qsort(points, npoints, sizeof(Point), compare_points);

    for (ii = 0; ii < npoints; ii++) {
        if ((ii > 0) &&
            (points[ii].x == points[ii-1].x) &&
            (points[ii].y == points[ii-1].y)) {
            continue;
        }

        while ((kk >= 2) && (cross(hull[kk-2], hull[kk-1], points[ii]) <= 0)) kk--;
        hull[kk++] = points[ii];
    }

    long lower = kk;
    for (ii = npoints - 2; ii >= 0; ii--) {
        if ((points[ii].x == points[ii+1].x) &&
            (points[ii].y == points[ii+1].y)) {
            if (ii == 0) break;
            continue;
        }

        while ((kk > lower) && (cross(hull[kk-2], hull[kk-1], points[ii]) <= 0)) kk--;
        hull[kk++] = points[ii];
        if (ii == 0) break;
    }

    if (kk > 1) kk--;
    return kk;
}

static double distance_points(Point aa, Point bb)
{
    double dx = aa.x - bb.x;
    double dy = aa.y - bb.y;
    return sqrt(dx*dx + dy*dy);
}

static double polygon_area(Point *poly, long npoly)
{
    double area = 0.0;
    long ii;

    if (npoly < 3) return 0.0;

    for (ii = 0; ii < npoly; ii++) {
        Point a = poly[ii];
        Point b = poly[(ii + 1) % npoly];
        area += (a.x * b.y) - (b.x * a.y);
    }

    return fabs(area) / 2.0;
}

static double polygon_perimeter(Point *poly, long npoly)
{
    double perimeter = 0.0;
    long ii;

    if (npoly < 2) return 0.0;

    for (ii = 0; ii < npoly; ii++) {
        perimeter += distance_points(poly[ii], poly[(ii + 1) % npoly]);
    }

    return perimeter;
}

static void hull_calipers(Point *hull, long nhull, ShapeStats *shape)
{
    double max_distance = 0.0;
    double min_width = DBL_MAX;
    double best_long = 0.0;
    double best_short = 0.0;
    double best_area = DBL_MAX;
    long ii, jj;

    if (nhull == 0) {
        shape->feret_diameter = invalid_double();
        shape->min_feret_diameter = invalid_double();
        shape->mbr_long_side = invalid_double();
        shape->mbr_short_side = invalid_double();
        return;
    }

    for (ii = 0; ii < nhull; ii++) {
        for (jj = ii + 1; jj < nhull; jj++) {
            max_distance = MAX(max_distance, distance_points(hull[ii], hull[jj]));
        }
    }

    if (nhull == 1) {
        shape->feret_diameter = 0.0;
        shape->min_feret_diameter = 0.0;
        shape->mbr_long_side = 0.0;
        shape->mbr_short_side = 0.0;
        return;
    }

    for (ii = 0; ii < nhull; ii++) {
        Point a = hull[ii];
        Point b = hull[(ii + 1) % nhull];
        double ex = b.x - a.x;
        double ey = b.y - a.y;
        double elen = sqrt(ex*ex + ey*ey);
        double ux, uy, vx, vy;
        double min_u = DBL_MAX, max_u = -DBL_MAX;
        double min_v = DBL_MAX, max_v = -DBL_MAX;

        if (elen <= 0.0) continue;

        ux = ex / elen;
        uy = ey / elen;
        vx = -uy;
        vy = ux;

        for (jj = 0; jj < nhull; jj++) {
            double proj_u = (hull[jj].x * ux) + (hull[jj].y * uy);
            double proj_v = (hull[jj].x * vx) + (hull[jj].y * vy);
            min_u = MIN(min_u, proj_u);
            max_u = MAX(max_u, proj_u);
            min_v = MIN(min_v, proj_v);
            max_v = MAX(max_v, proj_v);
        }

        double width = max_v - min_v;
        double side_u = max_u - min_u;
        double side_v = width;
        double area = side_u * side_v;

        min_width = MIN(min_width, width);

        if (area < best_area) {
            best_area = area;
            best_long = MAX(side_u, side_v);
            best_short = MIN(side_u, side_v);
        }
    }

    shape->feret_diameter = max_distance;
    shape->min_feret_diameter = (min_width == DBL_MAX) ? 0.0 : min_width;
    shape->mbr_long_side = best_long;
    shape->mbr_short_side = best_short;
}

static double point_segment_distance(double px, double py, Segment seg)
{
    double vx = seg.x2 - seg.x1;
    double vy = seg.y2 - seg.y1;
    double wx = px - seg.x1;
    double wy = py - seg.y1;
    double len2 = vx*vx + vy*vy;
    double tt;

    if (len2 <= 0.0) {
        double dx = px - seg.x1;
        double dy = py - seg.y1;
        return sqrt(dx*dx + dy*dy);
    }

    tt = ((wx * vx) + (wy * vy)) / len2;
    if (tt < 0.0) tt = 0.0;
    if (tt > 1.0) tt = 1.0;

    double qx = seg.x1 + (tt * vx);
    double qy = seg.y1 + (tt * vy);
    double dx = px - qx;
    double dy = py - qy;
    return sqrt(dx*dx + dy*dy);
}

static int is_mask_pixel(Image *mask, long mask_id, long xx, long yy)
{
    if ((xx < 0) || (yy < 0) || (xx >= mask->lAxes[0]) || (yy >= mask->lAxes[1])) {
        return 0;
    }

    return (get_image_value(mask, xx, yy) == mask_id);
}

static int add_boundary_segment(Image *mask, long mask_id, long xx, long yy, int side, Segment *segments, long *nsegments)
{
    int exposed = 0;
    Segment seg;

    switch (side) {
        case 0:
            exposed = !is_mask_pixel(mask, mask_id, xx, yy - 1);
            seg.x1 = xx - 0.5; seg.y1 = yy - 0.5;
            seg.x2 = xx + 0.5; seg.y2 = yy - 0.5;
            break;
        case 1:
            exposed = !is_mask_pixel(mask, mask_id, xx + 1, yy);
            seg.x1 = xx + 0.5; seg.y1 = yy - 0.5;
            seg.x2 = xx + 0.5; seg.y2 = yy + 0.5;
            break;
        case 2:
            exposed = !is_mask_pixel(mask, mask_id, xx, yy + 1);
            seg.x1 = xx + 0.5; seg.y1 = yy + 0.5;
            seg.x2 = xx - 0.5; seg.y2 = yy + 0.5;
            break;
        default:
            exposed = !is_mask_pixel(mask, mask_id, xx - 1, yy);
            seg.x1 = xx - 0.5; seg.y1 = yy + 0.5;
            seg.x2 = xx - 0.5; seg.y2 = yy - 0.5;
            break;
    }

    if (exposed) {
        segments[*nsegments] = seg;
        (*nsegments)++;
    }

    return 0;
}

static void count_holes(Image *mask, long mask_id, Stats *stats)
{
    long xmin = stats->xmin;
    long xmax = stats->xmax;
    long ymin = stats->ymin;
    long ymax = stats->ymax;
    long width = (xmax - xmin + 3);
    long height = (ymax - ymin + 3);
    long total = width * height;
    unsigned char *blocked = NULL;
    unsigned char *visited = NULL;
    long *queue = NULL;
    long holes = 0;
    long hole_area = 0;
    long xx, yy;

    blocked = (unsigned char *)calloc(total, sizeof(unsigned char));
    visited = (unsigned char *)calloc(total, sizeof(unsigned char));
    queue = (long *)calloc(total, sizeof(long));

    if ((blocked == NULL) || (visited == NULL) || (queue == NULL)) {
        free(blocked);
        free(visited);
        free(queue);
        stats->shape.num_holes = 0;
        stats->shape.holes_area = 0;
        return;
    }

    for (yy = ymin; yy <= ymax; yy++) {
        for (xx = xmin; xx <= xmax; xx++) {
            if (is_mask_pixel(mask, mask_id, xx, yy)) {
                long gx = xx - xmin + 1;
                long gy = yy - ymin + 1;
                blocked[(gy * width) + gx] = 1;
            }
        }
    }

    long head = 0, tail = 0;
    queue[tail++] = 0;
    visited[0] = 1;

    while (head < tail) {
        long idx = queue[head++];
        long gx = idx % width;
        long gy = idx / width;
        long dirs[4][2] = {{1,0},{-1,0},{0,1},{0,-1}};
        int dd;

        for (dd = 0; dd < 4; dd++) {
            long nx = gx + dirs[dd][0];
            long ny = gy + dirs[dd][1];
            long nidx;

            if ((nx < 0) || (ny < 0) || (nx >= width) || (ny >= height)) continue;
            nidx = (ny * width) + nx;
            if (visited[nidx] || blocked[nidx]) continue;

            visited[nidx] = 1;
            queue[tail++] = nidx;
        }
    }

    for (yy = 0; yy < height; yy++) {
        for (xx = 0; xx < width; xx++) {
            long idx = (yy * width) + xx;
            long component_area = 0;
            head = tail = 0;

            if (blocked[idx] || visited[idx]) continue;

            holes++;
            visited[idx] = 1;
            queue[tail++] = idx;

            while (head < tail) {
                long cur = queue[head++];
                long gx = cur % width;
                long gy = cur / width;
                long dirs[4][2] = {{1,0},{-1,0},{0,1},{0,-1}};
                int dd;

                if ((gx > 0) && (gx < width - 1) && (gy > 0) && (gy < height - 1)) {
                    component_area++;
                }

                for (dd = 0; dd < 4; dd++) {
                    long nx = gx + dirs[dd][0];
                    long ny = gy + dirs[dd][1];
                    long nidx;

                    if ((nx < 0) || (ny < 0) || (nx >= width) || (ny >= height)) continue;
                    nidx = (ny * width) + nx;
                    if (visited[nidx] || blocked[nidx]) continue;

                    visited[nidx] = 1;
                    queue[tail++] = nidx;
                }
            }

            hole_area += component_area;
        }
    }

    stats->shape.num_holes = holes;
    stats->shape.holes_area = hole_area;

    free(blocked);
    free(visited);
    free(queue);
}

static double fractal_box_dimension(Point *centers, long ncenters, long xmin, long ymin)
{
    int box_sizes[] = {2, 3, 4, 6, 8, 12, 16, 32, 64};
    int nbins = sizeof(box_sizes) / sizeof(box_sizes[0]);
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    int nfit = 0;
    int bb;

    for (bb = 0; bb < nbins; bb++) {
        int size = box_sizes[bb];
        long max_bx = 0;
        long max_by = 0;
        long ii;
        unsigned char *occupied = NULL;
        long boxes_x, boxes_y, count = 0;

        for (ii = 0; ii < ncenters; ii++) {
            long bx = (long)((centers[ii].x - xmin) / size);
            long by = (long)((centers[ii].y - ymin) / size);
            max_bx = MAX(max_bx, bx);
            max_by = MAX(max_by, by);
        }

        boxes_x = max_bx + 1;
        boxes_y = max_by + 1;
        occupied = (unsigned char *)calloc(boxes_x * boxes_y, sizeof(unsigned char));
        if (occupied == NULL) continue;

        for (ii = 0; ii < ncenters; ii++) {
            long bx = (long)((centers[ii].x - xmin) / size);
            long by = (long)((centers[ii].y - ymin) / size);
            occupied[(by * boxes_x) + bx] = 1;
        }

        for (ii = 0; ii < boxes_x * boxes_y; ii++) {
            if (occupied[ii]) count++;
        }

        free(occupied);

        if (count > 0) {
            double xval = log(1.0 / (double)size);
            double yval = log((double)count);
            sx += xval;
            sy += yval;
            sxx += xval * xval;
            sxy += xval * yval;
            nfit++;
        }
    }

    double denom = (nfit * sxx) - (sx * sx);
    if ((nfit < 2) || (fabs(denom) <= 0.0)) return invalid_double();

    return ((nfit * sxy) - (sx * sy)) / denom;
}

int compute_shape_stats(Image *mask, long mask_id, Stats *stats)
{
    long area = stats->area;
    long xx, yy, ii, jj;
    Point *corners = NULL;
    Point *centers = NULL;
    Point *hull = NULL;
    Segment *segments = NULL;
    long ncorners = 0;
    long ncenters = 0;
    long nhull = 0;
    long nsegments = 0;
    double max_radius = 0.0;
    double mu20 = 0.0, mu02 = 0.0, mu11 = 0.0;
    double shape_xavg = 0.0, shape_yavg = 0.0;
    double common, lambda_major, lambda_minor;

    if (area <= 0) return 1;

    corners = (Point *)calloc(area * 4, sizeof(Point));
    centers = (Point *)calloc(area, sizeof(Point));
    hull = (Point *)calloc(area * 4, sizeof(Point));
    segments = (Segment *)calloc(area * 4, sizeof(Segment));

    if ((corners == NULL) || (centers == NULL) || (hull == NULL) || (segments == NULL)) {
        free(corners);
        free(centers);
        free(hull);
        free(segments);
        return 1;
    }

    for (yy = stats->ymin; yy <= stats->ymax; yy++) {
        for (xx = stats->xmin; xx <= stats->xmax; xx++) {
            if (!is_mask_pixel(mask, mask_id, xx, yy)) continue;

            centers[ncenters].x = xx;
            centers[ncenters].y = yy;
            shape_xavg += xx;
            shape_yavg += yy;
            ncenters++;

            corners[ncorners].x = xx - 0.5; corners[ncorners].y = yy - 0.5; ncorners++;
            corners[ncorners].x = xx + 0.5; corners[ncorners].y = yy - 0.5; ncorners++;
            corners[ncorners].x = xx + 0.5; corners[ncorners].y = yy + 0.5; ncorners++;
            corners[ncorners].x = xx - 0.5; corners[ncorners].y = yy + 0.5; ncorners++;

            add_boundary_segment(mask, mask_id, xx, yy, 0, segments, &nsegments);
            add_boundary_segment(mask, mask_id, xx, yy, 1, segments, &nsegments);
            add_boundary_segment(mask, mask_id, xx, yy, 2, segments, &nsegments);
            add_boundary_segment(mask, mask_id, xx, yy, 3, segments, &nsegments);
        }
    }

    nhull = convex_hull(corners, ncorners, hull);

    stats->shape.convex_area = polygon_area(hull, nhull);
    stats->shape.convex_perimeter = polygon_perimeter(hull, nhull);
    hull_calipers(hull, nhull, &stats->shape);

    for (ii = 0; ii < ncenters; ii++) {
        double nearest = DBL_MAX;
        for (jj = 0; jj < nsegments; jj++) {
            nearest = MIN(nearest, point_segment_distance(centers[ii].x, centers[ii].y, segments[jj]));
        }
        if (nearest != DBL_MAX) max_radius = MAX(max_radius, nearest);
    }

    stats->shape.max_inscribed_circle_diameter = 2.0 * max_radius;
    stats->shape.eq_circle_diameter = sqrt((4.0 * (double)area) / M_PI);
    stats->shape.aspect_ratio = (stats->shape.mbr_short_side > 0.0) ?
        stats->shape.mbr_long_side / stats->shape.mbr_short_side : invalid_double();
    stats->shape.area_perimeter_ratio = (stats->perimeter > 0) ?
        (double)area / (double)stats->perimeter : invalid_double();
    stats->shape.circularity = (area > 0) ?
        ((double)stats->perimeter * (double)stats->perimeter) / (double)area : invalid_double();
    stats->shape.elongation = (stats->shape.mbr_long_side > 0.0) ?
        1.0 - (stats->shape.mbr_short_side / stats->shape.mbr_long_side) : invalid_double();
    stats->shape.convexity = (stats->perimeter > 0) ?
        stats->shape.convex_perimeter / (double)stats->perimeter : invalid_double();
    stats->shape.solidity = (stats->shape.convex_area > 0.0) ?
        (double)area / stats->shape.convex_area : invalid_double();
    stats->shape.thinness_ratio = (stats->perimeter > 0) ?
        (4.0 * M_PI * (double)area) / ((double)stats->perimeter * (double)stats->perimeter) :
        invalid_double();
    stats->shape.contour_temperature = ((stats->perimeter > 0) &&
                                        ((double)stats->perimeter > stats->shape.convex_perimeter)) ?
        (log((2.0 * (double)stats->perimeter) /
             ((double)stats->perimeter - stats->shape.convex_perimeter)) / log(2.0)) - 1.0 :
        invalid_double();
    stats->shape.fractal_box_dimension = fractal_box_dimension(centers, ncenters, stats->xmin, stats->ymin);

    shape_xavg /= (double)ncenters;
    shape_yavg /= (double)ncenters;

    for (ii = 0; ii < ncenters; ii++) {
        double dx = centers[ii].x - shape_xavg;
        double dy = centers[ii].y - shape_yavg;
        mu20 += dx * dx;
        mu02 += dy * dy;
        mu11 += dx * dy;
    }
    mu20 /= (double)area;
    mu02 /= (double)area;
    mu11 /= (double)area;

    common = sqrt(((mu20 - mu02) * (mu20 - mu02)) + (4.0 * mu11 * mu11));
    lambda_major = (mu20 + mu02 + common) / 2.0;
    lambda_minor = (mu20 + mu02 - common) / 2.0;

    stats->shape.orientation = atan2(2.0 * mu11, mu20 - mu02) * 90.0 / M_PI;
    if (stats->shape.orientation < 0.0) stats->shape.orientation += 180.0;
    stats->shape.major_axis = 4.0 * sqrt(MAX(lambda_major, 0.0));
    stats->shape.minor_axis = 4.0 * sqrt(MAX(lambda_minor, 0.0));

    count_holes(mask, mask_id, stats);
    stats->shape.porosity = ((area + stats->shape.holes_area) > 0) ?
        (double)stats->shape.holes_area / (double)(area + stats->shape.holes_area) : invalid_double();

    free(corners);
    free(centers);
    free(hull);
    free(segments);

    return 0;
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

    if (0 != compute_shape_stats(mask, mask_id, retvals)) {
        err_msg("ERROR: Could not compute shape statistics for mask %ld", mask_id);
        return NULL;
    }

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
    dmDescriptor *convex_area_col = dmColumnCreate(outBlock, "convex_area", dmDOUBLE, 0, "pixel**2", "Area enclosed by the convex hull");
    dmDescriptor *convex_perimeter_col = dmColumnCreate(outBlock, "convex_perimeter", dmDOUBLE, 0, "pixel", "Perimeter of the convex hull");
    dmDescriptor *feret_diameter_col = dmColumnCreate(outBlock, "feret_diameter", dmDOUBLE, 0, "pixel", "Maximum Feret diameter");
    dmDescriptor *min_feret_diameter_col = dmColumnCreate(outBlock, "min_feret_diameter", dmDOUBLE, 0, "pixel", "Minimum Feret diameter");
    dmDescriptor *max_inscribed_circle_diameter_col = dmColumnCreate(outBlock, "max_inscribed_circle_diameter", dmDOUBLE, 0, "pixel", "Diameter of the maximum inscribed circle");
    dmDescriptor *eq_circle_diameter_col = dmColumnCreate(outBlock, "eq_circle_diameter", dmDOUBLE, 0, "pixel", "Diameter of a circle with the same area");
    dmDescriptor *mbr_long_side_col = dmColumnCreate(outBlock, "mbr_long_side", dmDOUBLE, 0, "pixel", "Long side of the minimum bounding rectangle");
    dmDescriptor *mbr_short_side_col = dmColumnCreate(outBlock, "mbr_short_side", dmDOUBLE, 0, "pixel", "Short side of the minimum bounding rectangle");
    dmDescriptor *aspect_ratio_col = dmColumnCreate(outBlock, "aspect_ratio", dmDOUBLE, 0, "", "Minimum bounding rectangle long side divided by short side");
    dmDescriptor *area_perimeter_ratio_col = dmColumnCreate(outBlock, "area_perimeter_ratio", dmDOUBLE, 0, "pixel", "Area divided by perimeter");
    dmDescriptor *circularity_col = dmColumnCreate(outBlock, "circularity", dmDOUBLE, 0, "", "Perimeter squared divided by area");
    dmDescriptor *elongation_col = dmColumnCreate(outBlock, "elongation", dmDOUBLE, 0, "", "One minus short side divided by long side");
    dmDescriptor *convexity_col = dmColumnCreate(outBlock, "convexity", dmDOUBLE, 0, "", "Convex hull perimeter divided by perimeter");
    dmDescriptor *solidity_col = dmColumnCreate(outBlock, "solidity", dmDOUBLE, 0, "", "Area divided by convex hull area");
    dmDescriptor *num_holes_col = dmColumnCreate(outBlock, "num_holes", dmLONG, 0, "", "Number of holes inside the mask");
    dmDescriptor *thinness_ratio_col = dmColumnCreate(outBlock, "thinness_ratio", dmDOUBLE, 0, "", "Four pi area divided by perimeter squared");
    dmDescriptor *contour_temperature_col = dmColumnCreate(outBlock, "contour_temperature", dmDOUBLE, 0, "", "Contour temperature shape measure");
    dmDescriptor *orientation_col = dmColumnCreate(outBlock, "orientation", dmDOUBLE, 0, "degree", "Major-axis orientation counter-clockwise from positive x");
    dmDescriptor *major_axis_col = dmColumnCreate(outBlock, "major_axis", dmDOUBLE, 0, "pixel", "Major axis from second moments");
    dmDescriptor *minor_axis_col = dmColumnCreate(outBlock, "minor_axis", dmDOUBLE, 0, "pixel", "Minor axis from second moments");
    dmDescriptor *fractal_box_dimension_col = dmColumnCreate(outBlock, "fractal_box_dimension", dmDOUBLE, 0, "", "Fractal dimension from box counting");
    dmDescriptor *holes_area_col = dmColumnCreate(outBlock, "holes_area", dmLONG, 0, "pixel", "Area of holes inside the mask");
    dmDescriptor *porosity_col = dmColumnCreate(outBlock, "porosity", dmDOUBLE, 0, "", "Hole area divided by object plus hole area");
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
        dmSetScalar_d( convex_area_col, at->shape.convex_area);
        dmSetScalar_d( convex_perimeter_col, at->shape.convex_perimeter);
        dmSetScalar_d( feret_diameter_col, at->shape.feret_diameter);
        dmSetScalar_d( min_feret_diameter_col, at->shape.min_feret_diameter);
        dmSetScalar_d( max_inscribed_circle_diameter_col, at->shape.max_inscribed_circle_diameter);
        dmSetScalar_d( eq_circle_diameter_col, at->shape.eq_circle_diameter);
        dmSetScalar_d( mbr_long_side_col, at->shape.mbr_long_side);
        dmSetScalar_d( mbr_short_side_col, at->shape.mbr_short_side);
        dmSetScalar_d( aspect_ratio_col, at->shape.aspect_ratio);
        dmSetScalar_d( area_perimeter_ratio_col, at->shape.area_perimeter_ratio);
        dmSetScalar_d( circularity_col, at->shape.circularity);
        dmSetScalar_d( elongation_col, at->shape.elongation);
        dmSetScalar_d( convexity_col, at->shape.convexity);
        dmSetScalar_d( solidity_col, at->shape.solidity);
        dmSetScalar_l( num_holes_col, at->shape.num_holes);
        dmSetScalar_d( thinness_ratio_col, at->shape.thinness_ratio);
        dmSetScalar_d( contour_temperature_col, at->shape.contour_temperature);
        dmSetScalar_d( orientation_col, at->shape.orientation);
        dmSetScalar_d( major_axis_col, at->shape.major_axis);
        dmSetScalar_d( minor_axis_col, at->shape.minor_axis);
        dmSetScalar_d( fractal_box_dimension_col, at->shape.fractal_box_dimension);
        dmSetScalar_l( holes_area_col, at->shape.holes_area);
        dmSetScalar_d( porosity_col, at->shape.porosity);
        
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
