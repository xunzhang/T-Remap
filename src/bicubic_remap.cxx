#include "bicubic_remap.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "grid.h"
#include "grid_func.h"
#include "t_get_nbr.h"
#define L 0.2*PI

bicubic_remap_2D_operator::bicubic_remap_2D_operator(char *grid1_name, 
                                                       char *grid2_name, 
                                                       char *alg_name)
: Common_remap(grid1_name, grid2_name, alg_name)
{
	grid *grid_src, *grid_dst;
	int i,j;
    
	grid_src = cpl_grids->get_grid(grid_name_src);
	grid_dst = cpl_grids->get_grid(grid_name_dst);
    
	nlon_each_lat_src=grid_src->get_nlon_each_lat();
	nlon_each_lat_dst=grid_dst->get_nlon_each_lat();
    
	lat_begindx_src = new int [num_lats_src];
	lat_begindx_dst = new int [num_lats_dst];
	num_wgts = npts_dst * 16;
	wgt_indx_src = new int [num_wgts/4];
	wgt_indx_dst = NULL;
	wgt_values = new double [num_wgts];
    coefficients = new double [num_wgts];
    derivative_of_data_src = new double [npts_src*4];
	lat_begindx_src[0] = 0;
	for(i = 1; i < num_lats_src; i ++) 
		lat_begindx_src[i] = lat_begindx_src[i - 1] + nlon_each_lat_src[i - 1];
    
	lat_begindx_dst[0] = 0;
	for(i = 1; i < num_lats_dst; i ++)
		lat_begindx_dst[i] = lat_begindx_dst[i - 1] + nlon_each_lat_dst[i - 1];
}

bicubic_remap_2D_operator::~bicubic_remap_2D_operator()
{
	delete [] lat_begindx_src;
	delete [] lat_begindx_dst;
	delete [] wgt_indx_src;
	delete [] wgt_values;
	delete [] coefficients;
	delete [] derivative_of_data_src;
}

double f_case2(double lon, double lat){
    	lon = lon*PI/(double)180.0;
	lat = lat*PI/(double)180.0;
	return ((double)2.0+cos(lat)*cos(lat)*cos(2*lon));
}

double f_case2_diff_lon(double lon, double lat){
    	lon = lon*PI/(double)180.0;
	lat = lat*PI/(double)180.0;
	return -((double)2.0*sin(2*lon)*cos(lat)*cos(lat));
}

double f_case2_diff_lat(double lon, double lat){
    	lon = lon*PI/(double)180.0;
	lat = lat*PI/(double)180.0;
    return -((double)2.0*cos(2*lon)*cos(lat)*sin(lat));
}

double f_case2_diff_lat_lon(double lon, double lat){
    	lon = lon*PI/(double)180.0;
	lat = lat*PI/(double)180.0;
    return ((double)4.0*sin(2*lon)*cos(lat)*sin(lat));
}

void bicubic_remap_2D_operator::wgt_bicubic( double lons_src[2],
											double lats_src[2],
											double lon_dst,
											double lat_dst,
											double *wgts) 
{
	double t;
	double u;
	
	if (!((lon_dst >= lons_src[0]) &&
		  (lons_src[1] >= lons_src[0]))) {
		if (lon_dst < lons_src[0])
			lon_dst += 360;
		if (lons_src[1] < lons_src[0])
			lons_src[1] += 360;
	}
	if (!((lon_dst >= lons_src[2]) &&
		  (lons_src[3] >= lons_src[2]))) {
		if (lon_dst < lons_src[2])
			lon_dst += 360;
		if (lons_src[3] < lons_src[2])
			lons_src[3] += 360;
	}
    
	t=(lon_dst-lons_src[0])/(lons_src[1]-lons_src[0]);
	u=(lat_dst-lats_src[0])/(lats_src[1]-lats_src[0]);
    	
	wgts[0] = (1-(u*u)*(3-2*u))*(1-t*t*(3-2*t));
	wgts[1] = (1-(u*u)*(3-2*u))*t*t*(3-2*t);
	wgts[2] = u*u*(3-2*u)*t*t*(3-2*t);
	wgts[3] = u*u*(3-2*u)*(1-t*t*(3-2*t));
	
	wgts[4] = (1-u*u*(3-2*u))*t*(t-1)*(t-1);	
	wgts[5] = (1-u*u*(3-2*u))*t*t*(t-1);
	wgts[6] = u*u*(3-2*u)*t*t*(t-1);
	wgts[7] = u*u*(3-2*u)*t*(t-1)*(t-1);
	
	wgts[8] = u*(u-1)*(u-1)*(1-t*t*(3-2*t));
	wgts[9] = u*(u-1)*(u-1)*t*t*(3-2*t);
	wgts[10] = u*u*(u-1)*t*t*(3-2*t);
	wgts[11] = u*u*(u-1)*(1-t*t*(3-2*t));
	
	wgts[12] = t*(t-1)*(t-1)*u*(u-1)*(u-1);
	wgts[13] = t*t*(t-1)*u*(u-1)*(u-1);
	wgts[14] = t*t*(t-1)*u*u*(u-1);
	wgts[15] = t*(t-1)*(t-1)*u*u*(u-1);
}			 

void bicubic_remap_2D_operator::coefs_bicubic(int point_dst, double *coefs) 
{
    int *local_point_indx_src = wgt_indx_src;
    coefs[0] = derivative_of_data_src[4*local_point_indx_src[4*point_dst]]; 
    coefs[1] = derivative_of_data_src[4*local_point_indx_src[4*point_dst+1]];
    coefs[2] = derivative_of_data_src[4*local_point_indx_src[4*point_dst+3]];
    coefs[3] = derivative_of_data_src[4*local_point_indx_src[4*point_dst+2]];
    
    coefs[4] = derivative_of_data_src[4*local_point_indx_src[4*point_dst]+1];
    coefs[5] = derivative_of_data_src[4*local_point_indx_src[4*point_dst+1]+1];
    coefs[6] = derivative_of_data_src[4*local_point_indx_src[4*point_dst+3]+1]; 
    coefs[7] = derivative_of_data_src[4*local_point_indx_src[4*point_dst+2]+1];
    
    coefs[8] = derivative_of_data_src[4*local_point_indx_src[4*point_dst]+2];
    coefs[9] = derivative_of_data_src[4*local_point_indx_src[4*point_dst+1]+2];
    coefs[10] = derivative_of_data_src[4*local_point_indx_src[4*point_dst+3]+2]; 
    coefs[11] = derivative_of_data_src[4*local_point_indx_src[4*point_dst+2]+2];
    
    coefs[12] = derivative_of_data_src[4*local_point_indx_src[4*point_dst]+3];
    coefs[13] = derivative_of_data_src[4*local_point_indx_src[4*point_dst+1]+3];
    coefs[14] = derivative_of_data_src[4*local_point_indx_src[4*point_dst+3]+3];
    coefs[15] = derivative_of_data_src[4*local_point_indx_src[4*point_dst+2]+3];

}

void bicubic_remap_2D_operator::init_remap(double *temp_data_src)
{
	int m,n;
    	int bug;
	int i, point_dst;
	int bound_lat_indx1, bound_lat_indx2;
	int bound_lon_indx1, bound_lon_indx2;
	int left_nbr, right_nbr;
	double left_lon, right_lon;
	double left_lat, right_lat;
	double bicubic_lons[2];
	double bicubic_lats[2];
	grid *grid_src, *grid_dst;
	bool *mask_src, *mask_dst;
	int type_grid_src, type_grid_dst;
	double *max_lat_each_row_src, *min_lat_each_row_src;
	double last_lat = -10000;
	int nlons_src;
	bool lat_in, lon_in;
	double lons_src_box[4], lats_src_box[4];
	int corner_indx[4];
	double max_lat, min_lat;
	double max_lon, min_lon, tmp_lon;
	double *bound_boxes_src;
	double *cornner_boxes_lat_src, *cornner_boxes_lon_src;
    
	grid_src = cpl_grids->get_grid(grid_name_src);
	grid_dst = cpl_grids->get_grid(grid_name_dst);
	mask_src = grid_src->get_mymask();
	mask_dst = grid_dst->get_mymask();
	type_grid_src = grid_src->get_type_grid();
	type_grid_dst = grid_dst->get_type_grid();
	max_lat_each_row_src = grid_src->get_max_lat_each_row();
	min_lat_each_row_src = grid_src->get_min_lat_each_row();
	nlons_src = grid_src->get_num_lons();
    //generate_derivatives(data_src);
    
    //init derivatives
    int piece=nlon_each_lat_src[1];
    
    /*
     2----------3
     |          |
     |          |
     |   cell   |
     |          |
     |          |
     0----------1
     */
	//four corner-cell points
    //left_bottom_point
    derivative_of_data_src[0*4] = temp_data_src[0];
    derivative_of_data_src[0*4+1] = (temp_data_src[1] - temp_data_src[0]) / (center_lons_src[1] - center_lons_src[0]);
    derivative_of_data_src[0*4+2] = (temp_data_src[1*piece] - temp_data_src[0]) / (center_lats_src[1*piece] - center_lats_src[0]);
    derivative_of_data_src[0*4+3] = (temp_data_src[1+1*piece] - temp_data_src[1] - temp_data_src[1*piece] + temp_data_src[0]) / 
    ((center_lons_src[1] - center_lons_src[0]) * (center_lats_src[1*piece] - center_lats_src[0]));
    
    //right_bottom_point
    derivative_of_data_src[(piece-1)*4] = temp_data_src[piece-1];
    derivative_of_data_src[(piece-1)*4+1] = (temp_data_src[piece-1] - temp_data_src[piece-2]) / (center_lons_src[piece-1] - center_lons_src[piece-2]);
    derivative_of_data_src[(piece-1)*4+2] = (temp_data_src[piece-1+1*piece] - temp_data_src[piece-1]) / (center_lats_src[piece-1+1*piece] - center_lats_src[piece-1]);
    derivative_of_data_src[(piece-1)*4+3] = (temp_data_src[piece-1+1*piece] - temp_data_src[piece-1] - temp_data_src[piece-2+1*piece] + temp_data_src[piece-2]) / 
    ((center_lons_src[piece-1] - center_lons_src[piece-2]) * (center_lats_src[piece-1+1*piece] - center_lats_src[piece-1]));
    
    //left_top_point
    derivative_of_data_src[(npts_src-piece)*4] = temp_data_src[npts_src-piece];
    derivative_of_data_src[(npts_src-piece)*4+1] = (temp_data_src[npts_src-piece+1] - temp_data_src[npts_src-piece]) / (center_lons_src[npts_src-piece+1] - center_lons_src[npts_src-piece]);
    derivative_of_data_src[(npts_src-piece)*4+2] = (temp_data_src[npts_src-piece] - temp_data_src[npts_src-piece-1*piece]) / 
    (center_lats_src[npts_src-piece] - center_lats_src[npts_src-piece-1*piece]);
    derivative_of_data_src[(npts_src-piece)*4+3] = (temp_data_src[npts_src-piece+1] - temp_data_src[npts_src-piece+1-1*piece] - temp_data_src[npts_src-piece] + temp_data_src[npts_src-piece-1*piece])/   
    ((center_lons_src[npts_src-piece+1] - center_lons_src[npts_src-piece]) * (center_lats_src[npts_src-piece] - center_lats_src[npts_src-piece-1*piece]));
    
    //right_top_point
    derivative_of_data_src[(npts_src-1)*4] = temp_data_src[npts_src-1];
    derivative_of_data_src[(npts_src-1)*4+1] = (temp_data_src[npts_src-1] - temp_data_src[npts_src-2]) / (center_lons_src[npts_src-1] - center_lons_src[npts_src-2]);
    derivative_of_data_src[(npts_src-1)*4+2] = (temp_data_src[npts_src-1] - temp_data_src[npts_src-1-1*piece]) / (center_lats_src[npts_src-1] - center_lats_src[npts_src-1-1*piece]);
    derivative_of_data_src[(npts_src-1)*4+3] = (temp_data_src[npts_src-1] - temp_data_src[npts_src-1-1*piece] - temp_data_src[npts_src-2] + temp_data_src[npts_src-2-1*piece]) / 
    ((center_lons_src  [npts_src-1] - center_lons_src[npts_src-2]) * (center_lats_src[npts_src-1] - center_lats_src[npts_src-1-1*piece]));
    
    //four corner-line
    //bottom line
    for (m = 1; m < piece - 1; m ++) 
    {
        derivative_of_data_src[m*4] = temp_data_src[m];
        derivative_of_data_src[m*4+1] = (temp_data_src[m+1] - temp_data_src[m-1]) / (center_lons_src[m+1] - center_lons_src[m-1]);
        derivative_of_data_src[m*4+2] = (temp_data_src[m+1*piece] - temp_data_src[m]) / (center_lats_src[m+1*piece] - center_lats_src[m]);
        derivative_of_data_src[m*4+3] = (temp_data_src[m+1+1*piece] - temp_data_src[m+1] - temp_data_src[m-1+1*piece] + temp_data_src[m-1]) / 
        ((center_lons_src[m+1] - center_lons_src[m-1]) * (center_lats_src[m+1*piece] - center_lats_src[m]));
    }
    
    //top line
    for (m = npts_src-piece + 1; m < npts_src - 1; m ++) 
    {
        derivative_of_data_src[m*4] = temp_data_src[m];
        derivative_of_data_src[m*4+1] = (temp_data_src[m+1] - temp_data_src[m-1]) / (center_lons_src[m+1] - center_lons_src[m-1]);
        derivative_of_data_src[m*4+2] = (temp_data_src[m] - temp_data_src[m-1*piece]) / (center_lats_src[m] - center_lats_src[m-1*piece]);
        derivative_of_data_src[m*4+3] = (temp_data_src[m+1] - temp_data_src[m+1-1*piece] - temp_data_src[m-1] + temp_data_src[m-1-1*piece]) / 
        ((center_lons_src[m+1] - center_lons_src[m-1]) * (center_lats_src[m] - center_lats_src[m-1*piece]));
    }
    
    //left & right(together with normal points)
    for (m = piece; m < npts_src - piece; m ++) {	
        
        //left line
        if ((m % piece) == 0) {
            derivative_of_data_src[m*4] = temp_data_src[m];
            derivative_of_data_src[m*4+1] = (temp_data_src[m+1] - temp_data_src[m]) / (center_lons_src[m+1] - center_lons_src[m]);
            derivative_of_data_src[m*4+2] = (temp_data_src[m+1*piece] - temp_data_src[m-1*piece]) / (center_lats_src[m+1*piece] - center_lats_src[m-1*piece]);
            derivative_of_data_src[m*4+3] = (temp_data_src[m+1+1*piece] - temp_data_src[m+1-1*piece] - temp_data_src[m+1*piece] + temp_data_src[m-1*piece]) / 
            ((center_lons_src[m+1] - center_lons_src[m]) * (center_lats_src[m+1*piece] - center_lats_src[m-1*piece]));
            continue;
        }
        
        //right line
        if (((m+1) % piece) == 0){
            derivative_of_data_src[m*4] = temp_data_src[m];
            derivative_of_data_src[m*4+1] = (temp_data_src[m] - temp_data_src[m-1]) / (center_lons_src[m] - center_lons_src[m-1]);
            derivative_of_data_src[m*4+2] = (temp_data_src[m+1*piece] - temp_data_src[m-1*piece]) / (center_lats_src[m+1*piece] - center_lats_src[m-1*piece]);
            derivative_of_data_src[m*4+3] = (temp_data_src[m+1*piece] - temp_data_src[m-1*piece] - temp_data_src[m-1+1*piece] + temp_data_src[m-1-1*piece]) / 
            ((center_lons_src[m] - center_lons_src[m-1]) * (center_lats_src[m+1*piece] - center_lats_src[m-1*piece]));
            continue;
        }
        
        //normal points
        derivative_of_data_src[m*4] = temp_data_src[m];
        derivative_of_data_src[m*4+1] = (temp_data_src[m+1] - temp_data_src[m-1]) / (center_lons_src[m+1] - center_lons_src[m-1]);
        derivative_of_data_src[m*4+2] = (temp_data_src[m+1*piece] - temp_data_src[m-1*piece]) / (center_lats_src[m+1*piece] - center_lats_src[m-1*piece]);
        derivative_of_data_src[m*4+3] = (temp_data_src[m+1+1*piece] - temp_data_src[m+1-1*piece] - temp_data_src[m-1+1*piece] + temp_data_src[m-1-1*piece]) / 
        ((center_lons_src[m+1] - center_lons_src[m-1]) * (center_lats_src[m+1*piece] - center_lats_src[m-1*piece]));
    }
    //generate derivatives finished
/*
    
    int count;

    for (count = 0; count < npts_src; count ++) {
            derivative_of_data_src[4*count] = f_case2(center_lons_src[count], center_lats_src[count]);
            derivative_of_data_src[4*count+1] = f_case2_diff_lon(center_lons_src[count], center_lats_src[count]);
            derivative_of_data_src[4*count+2] = f_case2_diff_lat(center_lons_src[count], center_lats_src[count]);
            derivative_of_data_src[4*count+3] = f_case2_diff_lat_lon(center_lons_src[count], center_lats_src[count]);
    }
 
    
    
	for(count = 0; count < npts_src*4; count ++){
		printf("the %d-th diff is %lf\n",count,derivative_of_data_src[count]);
	}
  */  
		for (point_dst = 0; point_dst < npts_dst; point_dst ++) {
			if (last_lat != center_lats_dst[point_dst]) {
				last_lat = center_lats_dst[point_dst];
				bound_lat_indx1 = -1;
				for(i = 0; i < num_lats_src; i ++)
                    if(last_lat - max_lat_each_row_src[i] < 0) {
						bound_lat_indx1 = i - 1;
						break;
					}
				bound_lat_indx2 = bound_lat_indx1 + 1;
			}
			if (!mask_dst[point_dst]) {
                
                wgt_indx_src[point_dst * 4 + 0] = 0;
				wgt_indx_src[point_dst * 4 + 1] = 0;
				wgt_indx_src[point_dst * 4 + 2] = 0;
				wgt_indx_src[point_dst * 4 + 3] = 0;
                
                for (bug = 0; bug < 16; bug ++) {
                    wgt_values[point_dst * 16 + bug] = 0;
                }
                
                for (bug = 0; bug < 16; bug ++) {
                    coefficients[point_dst * 16 + bug] = 0;
                }
                
                continue;
			}
            
			if (bound_lat_indx1 == -1) {
				get_nearest_nbrs(center_lats_dst[point_dst], 
                                 center_lons_dst[point_dst],
                                 grid_src, 
                                 120.0, 
                                 4, 
                                 wgt_indx_src + point_dst * 4,
                                 wgt_values + point_dst * 16, 
                                 1.0);
                for (bug = 4; bug < 16; bug ++) {
                    wgt_values[point_dst * 16 + bug] = 0;
                }
                
                for (bug = 0; bug < 16; bug ++) {
                    coefficients[point_dst * 16 + bug] = 0;
                }
				continue;
			}
            
			bound_lon_indx2 = 0;
			bound_lon_indx1 = 0;
			for(i = 0; i < nlons_src; i ++) {
				left_nbr = (i - 1 + nlons_src) % nlons_src;
				right_nbr = i + lat_begindx_src[bound_lat_indx1];
				left_nbr += lat_begindx_src[bound_lat_indx1];
				left_lon = center_lons_src[left_nbr];
				right_lon = center_lons_src[right_nbr];
				if (right_lon >= left_lon) {
					if (center_lons_dst[point_dst] - right_lon < 0 &&
					    center_lons_dst[point_dst] - left_lon >= 0) {
						bound_lon_indx1 = i - 1;
						break;
					}
				}
				else if (center_lons_dst[point_dst] - right_lon < 0) {
					bound_lon_indx1 = i - 1;
					break;
				}
			}
			bound_lon_indx2 = i;
			if (i == 0 || i == nlons_src) {
				bound_lon_indx1 = nlons_src - 1;
				bound_lon_indx2 = 0;
			} 
			
			if (!mask_src[bound_lon_indx1 + lat_begindx_src[bound_lat_indx1]] ||
			    !mask_src[bound_lon_indx2 + lat_begindx_src[bound_lat_indx1]] ||
			    !mask_src[bound_lon_indx1 + lat_begindx_src[bound_lat_indx2]] ||
			    !mask_src[bound_lon_indx2 + lat_begindx_src[bound_lat_indx2]]) {
				get_nearest_nbrs(center_lats_dst[point_dst], 
                                 center_lons_dst[point_dst],
                                 grid_src, 
                                 120.0, 
                                 4, 
                                 wgt_indx_src + point_dst * 4,
                                 wgt_values + point_dst * 16, 
                                 1.0);
                for (bug = 4; bug < 16; bug ++) {
                    wgt_values[point_dst * 16 + bug] = 0;
                }
                
                for (bug = 0; bug < 16; bug ++) {
                    coefficients[point_dst * 16 + bug] = 0;
                }
			}
            else {
				wgt_indx_src[point_dst * 4 + 0] = bound_lon_indx1 + lat_begindx_src[bound_lat_indx1];
				wgt_indx_src[point_dst * 4 + 1] = bound_lon_indx2 + lat_begindx_src[bound_lat_indx1];
				wgt_indx_src[point_dst * 4 + 2] = bound_lon_indx1 + lat_begindx_src[bound_lat_indx2];
				wgt_indx_src[point_dst * 4 + 3] = bound_lon_indx2 + lat_begindx_src[bound_lat_indx2];
                
                bicubic_lons[0]= center_lons_src[bound_lon_indx1 + lat_begindx_src[bound_lat_indx1]];
				bicubic_lons[1] = center_lons_src[bound_lon_indx2 + lat_begindx_src[bound_lat_indx1]];
				
                bicubic_lats[0] = center_lats_src[lat_begindx_src[bound_lat_indx1]]; 
				bicubic_lats[1] = center_lats_src[lat_begindx_src[bound_lat_indx2]];
                
                
                /*
                printf("num id %d\n",point_dst);
                printf("lon_dst is %lf\n",center_lons_dst[point_dst]);
                printf("lat_dst is %lf\n", center_lats_dst[point_dst]);
                
                printf("zuobiao_lon[0] is %lf\t zuobiao_lon[1] is %lf\n",bicubic_lons[0],bicubic_lons[1]);
                printf("zuobiao_lon[0]_indx is %d\t zuobiao_lon[1]_indx is %d\n",bound_lon_indx1 + lat_begindx_src[bound_lat_indx1],bound_lon_indx2 + lat_begindx_src[bound_lat_indx1]);
                
                printf("zuobiao_lat[0] is %lf\t zuobiao_lat[1] is %lf\n",bicubic_lats[0],bicubic_lats[1]);
                printf("zuobiao_lat[0]_indx is %d\t zuobiao_lat[1]_indx is %d\n",lat_begindx_src[bound_lat_indx1],lat_begindx_src[bound_lat_indx2]);
                */
                
				wgt_bicubic(bicubic_lons, 
					    	bicubic_lats, 
					    	center_lons_dst[point_dst], 
					    	center_lats_dst[point_dst], 
					    	wgt_values + 16 * point_dst);
                
				coefs_bicubic(point_dst, coefficients + 16 * point_dst);
			}
		}
}
void bicubic_remap_2D_operator::cal_remap(double *data_src, double *data_dst)
{
	int i, j = 0;
	double *local_wgts = wgt_values;
	double *local_coefs = coefficients;
	
	for (i = 0, j = 0; i < npts_dst; i ++) 
	{
		data_dst[i] = 	local_coefs[j] * local_wgts[j] + 
				local_coefs[j+1] * local_wgts[j+1] +
				local_coefs[j+2] * local_wgts[j+2] +
				local_coefs[j+3] * local_wgts[j+3] +
				local_coefs[j+4] * local_wgts[j+4] + 
				local_coefs[j+5] * local_wgts[j+5] + 
				local_coefs[j+6] * local_wgts[j+6] +
				local_coefs[j+7] * local_wgts[j+7] + 
				local_coefs[j+8] * local_wgts[j+8] + 
				local_coefs[j+9] * local_wgts[j+9] + 
				local_coefs[j+10] * local_wgts[j+10] +
				local_coefs[j+11] * local_wgts[j+11] +
				local_coefs[j+12] * local_wgts[j+12] + 
				local_coefs[j+13] * local_wgts[j+13] + 
				local_coefs[j+14] * local_wgts[j+14] + 
				local_coefs[j+15] * local_wgts[j+15];
		j += 16;
	}
}
