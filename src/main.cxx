/*
极地有问题
*/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "dist_remap.h"
#include "grid.h"
#include "grid_func.h"
#include "data_func.h"
#include "io.h"
#include "remap_func.h"
#include "remap_cfg.h"
#include <sys/time.h>
#include <string.h>
#include "models_cfg.h"

void wtime(double *t)
{
  static int sec = -1;
  struct timeval tv;
  gettimeofday(&tv, (__timezone_ptr_t)0);
  if (sec < 0) 
   sec = tv.tv_sec;
  *t = (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;
}

void Print_Dst_To_File(char* str, int npts_dst_para, double* lon_coords_dst_para, double* lat_coords_dst_para, double* data_dst_para)
{
	FILE* out;
	if ((out = fopen(str, "a+")) == NULL)
	{
		printf("Cannot open file!\n");
	}
	else
	{
		int count = 0;
		while (count < npts_dst_para)
		{
			char out1[20] , out2[20] , out3[20];
			const char* out4 = ",";
			const char* out5 = ":";
			const char* out6 = "\n";
			sprintf(out1 , "res-%d: %lf" , count, lon_coords_dst_para[count]);
			sprintf(out2 , "%lf" , lat_coords_dst_para[count]);
			sprintf(out3 , "%lf" , data_dst_para[count]);

			fputs((const char*)out1, out);
			fputs(out4 , out);

			fputs((const char*)out2 , out);	
			fputs(out5 , out); 

			fputs((const char*)out3 , out);
			fputs(out6 , out);

			count++;
		}
	}
	fclose(out);
}


int main (int argc, char **argv)
{
	int i, j;
	int read_size;
	double *data_src;				// the src grid_2D data to be remapped
	int size_input_data_file;
	double *data_dst;				// the dst grid_2D data after remapping				
	FILE *fp_data_src;
	FILE *fp_data_dst;
	grid *grid_src, *grid_dst;
	char *name_data_field = data_remaps_cfg[0].name_data_field;
	char *grid_name1 = data_remaps_cfg[0].name_grid_src;
	char *grid_name2 = data_remaps_cfg[0].name_grid_dst;
	double time1, time2;

	if (argc < 5) {
		printf("usage <grid1> <grid2> <input_data> <algorithm>\n");
		exit(0);
	}

	strcpy(atm_model_set.models_cfg[0].grid_2D.name_file, argv[1]);
	strcpy(ocn_model_set.models_cfg[0].grid_2D.name_file, argv[2]);
	strcpy(data_remaps_cfg[0].name_remap_alg, argv[4]);
	generate_all_cpl_grids();
	grid_src = cpl_grids->get_grid(grid_name1);
	grid_dst = cpl_grids->get_grid(grid_name2);

	if ((fp_data_src = fopen(argv[3], "rb")) == NULL) {
			printf("Cannot open input data file!");
			exit(0);
	}
	if ((fp_data_dst = fopen("our_output_data", "w")) == NULL) {
			printf("Cannot open!");
			exit(0);
	}

	fseek(fp_data_src, 0, SEEK_END);
	size_input_data_file = ftell(fp_data_src) / sizeof(double);
	fseek(fp_data_src, 0, SEEK_SET);
	data_src = new double [size_input_data_file];
	data_dst = new double [grid_dst->get_num_points()];	
	fread(data_src, sizeof(char), size_input_data_file*sizeof(double), fp_data_src);

#ifdef NEW_DATA
	FILE *fp_new_data = fopen("file_data_atm_0.5x0.5_sparse4_lat-45to45_lon90to270", "w+");
#endif
    
    wtime(&time1);	
	generate_remap_operators(data_src);
	wtime(&time2);
	printf("initialize_time %f\n", time2 - time1);
    
	for (i = 0; i < size_input_data_file; i += grid_src->get_num_points()) {
		wtime(&time1);
		do_remap(name_data_field, grid_name1, grid_name2, data_src + i, data_dst);
#ifdef NEW_DATA		
		fwrite(data_dst, sizeof(double), grid_dst->get_num_points(), fp_new_data);
#endif
		//do_remap(name_data_field, grid_name2, grid_name1, data_dst, data_src + i / sizeof(double));
		wtime(&time2);

		printf("remap_time %f\n", time2 - time1);
	
		//Print_Dst_To_File("data_dst1", grid_dst->get_num_points(), grid_dst->get_lon_coords(), grid_dst->get_lat_coords(), data_dst);
	}

	finalize_all_cpl_grids();
	finalize_remap_operators();
	fwrite(data_dst, sizeof(char), grid_dst->get_num_points() * sizeof(double), fp_data_dst);	
	
	return 0;
}
