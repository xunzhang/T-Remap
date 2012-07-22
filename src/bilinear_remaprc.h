#ifndef BILINEAR_REMAP_H
#define BILINEAR_REMAP_H

#include "remap_base.h"
#include "remap_func.h"
class bilinear_remap_2D_operator : public common_remap
{	
	private:
		double *center_lons_src;			//
		double *center_lats_src;			//
		double *center_lons_dst;			//
		double *center_lats_dst;			//		
		int *nlon_each_lat_src;			//
		int *nlon_each_lat_dst;			//
		int *lat_begindx_src;				//
		int *lat_begindx_dst;				//
		int num_lats_src;					//
		int num_lats_dst;				//
		int *point_indx_src;				// The index of each point in src grid for each operation
		double *weights;					// The weight for each operation

	public:
		bilinear_remap_2D_operator(char*, char*, char*);
		~bilinear_remap_2D_operator();
		void wgt_bilinear( double nlons_src[4],//x1a
						   double nlats_src[2],//x2a
						   double nlons_dst_coords,   //x1
						   double nlats_dst_coords,   //x2
						   double data_src[4]);
		void init_remap(); 
		void cal_remap(double *, double *);
};
/*sequence of grid points——
		
		(x3,y3)      (x4,y4)
			p3----------p4
			|			|
			|  .P(x,y)  |
			|			|
			|			|
			p1----------p2
		(x1,y1)     (x4,y2)
*/
#endif
