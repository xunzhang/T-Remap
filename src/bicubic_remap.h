/*
 *  bicubic_remap.h
 *  bicubic_remap
 *
 *  Created by wu hong on 10-10-10.
 *  Copyright 2010 Tsinghua university. All rights reserved.
 *
 */
#ifndef BICUBIC_REMAP_H
#define BICUBIC_REMAP_H

#include "remap_base.h"
#include "remap_func.h"
class bicubic_remap_2D_operator : public Common_remap
{   
private:
    int *nlon_each_lat_src;			
    int *nlon_each_lat_dst;			
    int *lat_begindx_src;				
    int *lat_begindx_dst;
	double *coefficients;
	double *derivative_of_data_src;
	
public:
	bicubic_remap_2D_operator(char*, char*, char*);
	~bicubic_remap_2D_operator();
    double f1_case2(double, double);
    double f2_case2_diff_lon(double, double);
    double f3_case2_diff_lat(double, double);
    double f4_case2_diff_lat_lon(double, double);
	void wgt_bicubic(double*, double*, double, double, double*);
	void coefs_bicubic(int, double*);
    void init_remap(double *);
	void cal_remap(double *, double *);
};
#endif
