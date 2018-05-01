#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "/panfs/panfs.maxeler/packages/gcc/fftw/3.3.4/include/fftw3.h"
#include <complex.h> //after fftw
#include <stdbool.h>
#include <time.h>
#include <float.h>

#include "/panfs/panfs.maxeler/home/HCEEC001/nnm06/nxj57-nnm06/wiener_workspace/maxpower/src/maxpower/kernel/random/runtime/random_mt.h"

#include "Maxfiles.h"                   // Includes .max files
#include <MaxSLiCInterface.h>   // Simple Live CPU interface


float two_over_randmax_single = 2.0/((float)RAND_MAX);
float two_over_randmax = 2.0/(RAND_MAX);

static void* malloc_wrapper(size_t size) {
	void * pointer = malloc(size);
	if (pointer == NULL) {
		fprintf(stderr, "Was not able to allocate enough memory");
		exit(-1);
	}
	return pointer;
}


// Using Box-Muller in the polar form, return a complex number, 
// where the real and imag are both drawn from N(0,1)
// this has been tested and is roughly 2x faster than c++'s std::normal_distribution
static
void normal_pointer (float * u1, float * u2) {
	float w = 2.0;
	while ( w >=1.0 ) {
		*u1 = rand() * two_over_randmax_single - 1.0;
		*u2 = rand() * two_over_randmax_single - 1.0;
		w = (*u1) * (*u1) + (*u2) * (*u2);
	}

	w = sqrt((-2.0 * log(w))/w);
	*u1 = (*u1) * w;
	*u2 = (*u2) * w;
} 

// alternative (slower) method:
float complex normal_complex () {

	float u1, u2;
	float w = 2.0;
	while ( w >=1.0 ) {
		
		u1 = rand() * two_over_randmax_single - 1.0;
		u2 = rand() * two_over_randmax_single - 1.0;
		w = u1 * u1 + u2 * u2;
	}

	w = sqrt((-2.0 * log(w))/w);
	return (u1 + I * u2) * w; 
}


static
bool crelequal(float  a, float  b, float tolerance)
{
    float diff = cabsf(a - b);
    if ((diff <= cabsf(a) * tolerance) & (diff <= cabsf(b)*tolerance))
        return true;
    return false;
}

// Get time in seconds
static
double timesec()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	return ((double) (t.tv_sec * 1e6L + t.tv_usec)) / 1.0e6;
}


static
void sampler_cpu(int iterationmax, int map_number, int gridwidth, int size, float gridsquare_inv, float *tauvalues, float *Tftinvvalues , float *data, float *Sft, float *Nbar, float complex *sandtCPU)
{


	double start_time = timesec();
	float complex_normal;

	fftwf_complex *mapCPU = malloc_wrapper(size*sizeof(fftwf_complex));
	fftwf_complex *map2CPU = malloc_wrapper(size*sizeof(fftwf_complex));

	float *gaussian1 = malloc_wrapper(size*sizeof(float));
	float *gaussian2 = malloc_wrapper(size*sizeof(float));

	int i, j, fullsizehere, fullsize, mapcount, back_one_iteration, streamindex;

	// For each iteration
	int iteration = 0;
	for(iteration=0 ; iteration < iterationmax ; iteration++ ) {
		fullsize = iteration*size*map_number; //total number of pixels including this iteration .i.e. final index
		back_one_iteration = (iteration-1)*size*map_number; //total number of pixels excluding this iteration i.e. preceding index

		// for each map
		mapcount = 0;
		for(mapcount = 0 ; mapcount < map_number ; mapcount++ ){
			fullsizehere = size*mapcount; 

			// Take current map values into mapCPU as the current map, with the first index as the map element,
			// and the second index for real or imaginary. 
			// This format is so that mapCPU can be taken in by the fftw function.
			j = 0;
			for(j=0;j< size ; j++) {
				mapCPU[j][0]=creal(sandtCPU[j+fullsizehere+back_one_iteration]);
				mapCPU[j][1]=cimag(sandtCPU[j+fullsizehere+back_one_iteration]); 

				// Generate pairs of gaussian random variates for the two dot product sections
				normal_pointer(&gaussian1[j], &gaussian2[j]);
				
			}

			// First dot product
			i = 0;

			for(i=0;i<(size); i++) {
				streamindex = i + fullsizehere;
				mapCPU[i][0] =  ( tauvalues[mapcount] * data[streamindex] 
						+ mapCPU[i][0] * Nbar[streamindex])    
						/   (tauvalues[mapcount] + Nbar[streamindex] ) //
						+  sqrt((tauvalues[mapcount] * Nbar[streamindex]) 
						/ (tauvalues[mapcount]+Nbar[streamindex] ) ) * gaussian1[i]; //real part

				mapCPU[i][1] = ( mapCPU[i][1] * Nbar[streamindex])  
						/  (tauvalues[mapcount]+Nbar[streamindex] ); //imag part

			}

			

			//FFT 1
			fftwf_plan plan = fftwf_plan_dft_2d(gridwidth, gridwidth  ,  mapCPU, map2CPU, FFTW_FORWARD,  FFTW_ESTIMATE);
			fftwf_execute(plan);

			// Second dot product
			i = 0;
			for(i=0;i<(size); i++) {
				streamindex = i + fullsizehere;
				//real part divided by normalisation:
				map2CPU[i][0] =  gridsquare_inv * (( Tftinvvalues[mapcount] * map2CPU[i][0] * Sft[streamindex])
						/ (1.0 + Tftinvvalues[mapcount] * Sft[streamindex] ) //  
				//map2CPU[i][0] =  (( Tftinvvalues[mapcount] * map2CPU[i][0] * Sft[streamindex])
					//	/ (1.0 + Tftinvvalues[mapcount] * Sft[streamindex] ) );//  
						+ sqrt(  Sft[streamindex] 
							/ (1.0 + Tftinvvalues[mapcount] * Sft[streamindex] ) ) * gaussian2[i] ); 
				//imag part conjugated (minus sign) divided by normalisation:
				
				map2CPU[i][1] = (-gridsquare_inv) * (Tftinvvalues[mapcount] * map2CPU[i][1] * Sft[streamindex])
							/ (1.0+Tftinvvalues[mapcount] * Sft[streamindex] ); 
				//map2CPU[i][1] =  (-Tftinvvalues[mapcount] * map2CPU[i][1] * Sft[streamindex])
				//			/ (1.0+Tftinvvalues[mapcount] * Sft[streamindex] ); 

			}



			//FFT 2
			fftwf_plan plan2 = fftwf_plan_dft_2d(gridwidth, gridwidth  ,  map2CPU, mapCPU, FFTW_FORWARD,  FFTW_ESTIMATE);
			fftwf_execute(plan2);

			// put the back in sandtCPU
			j = 0;
			for(j=0;j< size ; j++) {
				sandtCPU[fullsize+j+fullsizehere] = mapCPU[j][0]+I*mapCPU[j][1];
			} 
		}
	}

	double end_time = timesec();

	printf("%d, %.6f\n", iterationmax,  end_time - start_time);
	
}




static
void sampler_DFE(int iterationmax, int map_number, int size, float gridsquare_inv, float *tauvalues, float *Tftinvvalues , float *data, float *Sft, float *Nbar, float complex *sandtDFE, float complex *reconDFE){


	double start_time = timesec();
	const size_t transfer_size_real = sizeof(float) * size*map_number;
	const size_t transfer_size_complex = 2*transfer_size_real;

	
	max_file_t *maxfile = sampler_init(); //!
	max_engine_t *engine = max_load(maxfile, "*");
	max_actions_t* act = max_actions_init(maxfile, NULL);

	max_set_ticks(act, "samplerKernel", size/4*map_number*iterationmax); //size/4 as 4 pipes in DFE

	//max_set_double(act, "samplerKernel", "iteration_scalar", iterationmax);  // these are commented out if outputing all iterations

	max_set_double(act, "samplerKernel", "Tftinv0", Tftinvvalues[0]);
	//max_set_double(act, "samplerKernel","Tftinv1" , Tftinvvalues[1] );
	//max_set_double(act, "samplerKernel","Tftinv2" , Tftinvvalues[2]);
	//max_set_double(act, "samplerKernel","Tftinv3", Tftinvvalues[3]);
	//max_set_double(act, "samplerKernel","Tftinv4", Tftinvvalues[4] );
	max_set_double(act, "samplerKernel","gridsquare_inv", gridsquare_inv);
	max_set_double(act, "samplerKernel","tau0", tauvalues[0]);
	//max_set_double(act, "samplerKernel","tau1", tauvalues[1]);
	//max_set_double(act, "samplerKernel","tau2",tauvalues[2]);
	//max_set_double(act, "samplerKernel","tau3", tauvalues[3]);
	//max_set_double(act, "samplerKernel","tau4", tauvalues[4]);



	max_queue_input(act,  "Nbar",  Nbar,  transfer_size_real);
	max_queue_input(act,  "Sft",  Sft,  transfer_size_real);
	max_queue_input(act,  "data",  data,  transfer_size_real);
	max_queue_input(act,  "sandtDFE",  sandtDFE,  transfer_size_complex);

	max_queue_output(act, "reconDFE", reconDFE, iterationmax*transfer_size_complex); //CHANGE

	random_mt_init(maxfile, act, "samplerKernel", "seed1_mt", time(NULL));
	random_mt_init(maxfile, act, "samplerKernel", "seed2_mt", time(NULL)+42);

	max_run(engine, act);
	max_unload(engine);

	double end_time = timesec();

	printf("%d, %.6f\n", iterationmax,  end_time - start_time);


}

static
void reorder_signal_covariance(int gridwidth, int map_number, int size, float *Sft) {
	
	float *Sft2 = malloc(map_number * size * sizeof(float));
	int mapcount;
	// Reorder Sft in a way such that no Fourier k shifts are required:
	int k = 0;
	for ( k = 0; k < size/2; k++){
		mapcount =0;
		for(mapcount=0; mapcount < map_number; mapcount++){
			Sft2[k + size*mapcount] = Sft[k + size*mapcount + size/2];
			Sft2[k + size*mapcount + size/2] = Sft[k + size*mapcount];
		}
	}

	int flipcount =0;
	for ( flipcount = 0; flipcount < gridwidth; flipcount++){
		k =0;
		for ( k = 0; k < gridwidth/2; k++){
			mapcount =0;
			for(mapcount=0; mapcount < map_number; mapcount++){
				Sft[flipcount*gridwidth + k + size*mapcount] = Sft2[flipcount*gridwidth +k+gridwidth/2 + size*mapcount];
				Sft[flipcount*gridwidth +k+gridwidth/2 + size*mapcount] = Sft2[flipcount*gridwidth + k + size*mapcount];
			}
		}
	}

	free(Sft2);

}



int main(void)
{

	// -------------------------------------------------------------------------//
	// ---------------------------Define variables------------------------------//
	// -------------------------------------------------------------------------//
	const int repeat_number = 1;
	const int iterationmin = 3;
	const int iterationmax = 3; // at the moment this is maxed at 1^5 by the variable CountMax in the kernel
	const bool print_bool = true; // whether or not to save the output
	printf("sampler_iterationmax %d \n", iterationmax);

	//Define these variables in the parameters file:
	const int map_number = sampler_map_number;
	const int gridwidth = sampler_N;
	const int DFE_frequency = sampler_streamFrequency;
	const int size = gridwidth*gridwidth;
	int mapcount, i, j;

	printf("Filtering %d maps (%d by %d) with %d iterations \n", map_number, gridwidth, gridwidth, iterationmax);
	printf("DFE stream clock frequency =  %d\n", DFE_frequency);
	// need this for format to stop integer overflows for large iteration numbers:
	size_t size_complex_output, size_complex_input;
	size_complex_input = (map_number*size*sizeof(float complex));
	size_complex_output = size_complex_input*iterationmax;

	
	// Define the vector/array quantities:
	int sizeBytes = size * sizeof(float);
	float *tauvalues = malloc(map_number*sizeof(float));
	float *Tftinv = malloc(map_number*sizeBytes); // check times with Tftinv streaming in compared to individual scalars
	float *Tftinvvalues = malloc(map_number*sizeof(float));
	float *Nbar = malloc(map_number*sizeBytes);
	float *Sft = malloc(map_number*sizeBytes);
	float *data = malloc(map_number*sizeBytes);
	float complex *reconDFE = malloc(size_complex_output);
	float complex *sandtDFE = malloc(size_complex_input);
	// check this- putting the CPU allocation first leads to a seg fault for large numbers of iterations
	float complex *sandtCPU = malloc(size_complex_output);

	printf("Test line \n");
	// -------------------------------------------------------------------------//
	// --------------------------Read in file data------------------------------//
	// -------------------------------------------------------------------------//

	char line[32];
	char line2[32];
	char line3[32];
	char line4[32];
	char line5[32];

	// reading in Sft and data
	i=0;
	FILE *infile ;
	FILE *infile2 ;
	FILE *infile3 ;
	FILE *infile4 ;
	FILE *infile5 ;

	//infile = fopen("maps256/paper_Sft_diagonal_128.txt", "r");
	infile = fopen("maps256/Sft1", "r");
	infile2 = fopen("maps256/Sft2", "r");
	infile3 = fopen("maps256/Sft3", "r");
	infile4 = fopen("maps256/Sft4", "r");
	infile5 = fopen("maps256/Sft5", "r");

	while ((fgets(line,sizeof line, infile) != NULL) & (fgets(line2,sizeof line2, infile2) != NULL) & (fgets(line3,sizeof line3, infile3) != NULL )& (fgets(line4,sizeof line4, infile4) != NULL) & (fgets(line5,sizeof line5, infile5) != NULL)) {
		Sft[i]=(atof(line));
		Sft[i+size] = atof(line2);
		Sft[i+2*size] = atof(line3);
		Sft[i+3*size] = atof(line4);
		Sft[i+4*size] = atof(line5);
		/*Sft[i+size] = atof(line);
		Sft[i+2*size] = atof(line);
		Sft[i+3*size] = atof(line);
		Sft[i+4*size] = atof(line);*/
		i++;
	}
	fclose(infile);
	fclose(infile2);
	fclose(infile3);
	fclose(infile4);
	fclose(infile5);

	i=0;
	//infile = fopen("maps256/paper_data_diagonal_128.txt", "r");
	infile = fopen("maps256/data1", "r");
	infile2 = fopen("maps256/data2", "r");
	infile3 = fopen("maps256/data3", "r");
	infile4 = fopen("maps256/data4", "r");
	infile5 = fopen("maps256/data5", "r");

	while ((fgets(line,sizeof line, infile) != NULL) & (fgets(line2,sizeof line2, infile2) != NULL) & (fgets(line3,sizeof line3, infile3) != NULL )& (fgets(line4,sizeof line4, infile4) != NULL) & (fgets(line5,sizeof line5, infile5) != NULL)) {
		data[i]= (atof(line));
		data[i+size]= (atof(line2));
		data[i+size*2]= (atof(line3));
		data[i+size*3]= (atof(line4));
		data[i+size*4]= (atof(line5));
		/*data[i+size]= (atof(line));
		data[i+size*2]= (atof(line));
		data[i+size*3]= (atof(line));
		data[i+size*4]= (atof(line));*/
		
		i++;
	}
	fclose(infile);
	fclose(infile2);
	fclose(infile3);
	fclose(infile4);
	fclose(infile5);

	// Read in N and calculate the minimum of N ( = tau1):
	i = 0;
	//infile = fopen("maps256/paper_N_diagonal_128.txt", "r");
	infile = fopen("maps256/N1", "r");
	infile2 = fopen("maps256/N2", "r");
	infile3 = fopen("maps256/N3", "r");
	infile4 = fopen("maps256/N4", "r");
	infile5 = fopen("maps256/N5", "r");

	while ((fgets(line,sizeof line, infile) != NULL )& (fgets(line2,sizeof line2, infile2) != NULL) & (fgets(line3,sizeof line3, infile3) != NULL )& (fgets(line4,sizeof line4, infile4) != NULL) & ( fgets(line5,sizeof line5, infile5) != NULL)) {

		Nbar[i]= (atof(line));
		Nbar[i+size]= (atof(line2));
		Nbar[i+size*2]= (atof(line3));
		Nbar[i+size*3]= (atof(line4));
		Nbar[i+size*4]= (atof(line5));
		/*Nbar[i+size]= (atof(line));
		Nbar[i+size*2]= (atof(line));
		Nbar[i+size*3]= (atof(line));
		Nbar[i+size*4]= (atof(line));*/

		mapcount=0;
		for(mapcount=0 ; mapcount < map_number; mapcount++) {
			if (i == 0){
				tauvalues[0]=Nbar[0];
				tauvalues[1] = Nbar[size];
				tauvalues[2] = Nbar[size*2];
				tauvalues[3] = Nbar[size*3];
				tauvalues[4] = Nbar[size*4];
			}
			else {
				if(Nbar[i+size*mapcount] < tauvalues[mapcount]){
					tauvalues[mapcount] = Nbar[i+size*mapcount];
				}
			}
		}
		i++;
	}
	fclose(infile);
	fclose(infile2);
	fclose(infile3);
	fclose(infile4);
	fclose(infile5);


	// Create Nbar = N - min(N) , where min(N) = tau1:

	j =0;
	for ( j = 0; j < size; j++){
		mapcount= 0;
		for(mapcount=0; mapcount < map_number; mapcount++){
			Nbar[j+mapcount*size]= Nbar[j+mapcount*size]-tauvalues[mapcount]; // rescale Nbar as Nbar = N - tau
			Tftinv[j+size*mapcount] = ((float)size)/tauvalues[mapcount]; // T^{-1} \propto 1 /{\tau}
		}
	}


	// -------------------------------------------------------------------------//
	// ----------------------Reorder data for management------------------------//
	// -------------------------------------------------------------------------//

	reorder_signal_covariance(gridwidth, map_number, size, Sft);

	//  Initialise input data
	i=0;
	for( i = 0; i < size*map_number; ++i) {
		sandtDFE[i] =  0.0 + 0.0*I;
		sandtCPU[i] =  0.0 + 0.0*I;
	}

	
	/*//  Initialise input data randomly
	for( i = 0; i < size*map_number; ++i) {
		sandtDFE[i] =  200 * creal(normal()) + 0.0 * I;
		sandtCPU[i] =  100 * creal(normal()) + 0.0 * I;
	}*/
	



	// -------------------------------------------------------------------------//
	// -------------------------Define scalars----------------------------------//
	// -------------------------------------------------------------------------//

	printf("Define scalars for DFE (and CPU)... \n");
	float gridsquare_inv = 1.0/(gridwidth*gridwidth);

	i=0;
	for(i=0; i < map_number; i++) {
		Tftinvvalues[i]= ((float)size)/tauvalues[i];
	}



	// -------------------------------------------------------------------------//
	// ---------------------- DFE sampler Filter---------------------------------//
	// -------------------------------------------------------------------------//

	
	int current_iteration;

	printf("Running on DFE.\n");
	current_iteration = iterationmin;
	for( current_iteration = iterationmin; current_iteration < iterationmax+1; current_iteration = 10*current_iteration){
		j = 0;
		for( j = 0; j < repeat_number; j++) {
			sampler_DFE( current_iteration,
						map_number,
						size,
						gridsquare_inv,
						tauvalues,
						Tftinvvalues ,
						data,
						Sft,
						Nbar,
						sandtDFE,
						reconDFE);  
		}
	}

	// -------------------------------------------------------------------------//
	// --------------------------Write to File----------------------------------//
	// -------------------------------------------------------------------------//

	int whichmap = iterationmax;
	int step = (whichmap-1)*size*map_number;
	if( print_bool == true ) {
		printf("Writing (real) DFE results to file... \n");

		FILE *g1=fopen("DFEresultsreal1.txt", "w");
		FILE *g2=fopen("DFEresultsreal2.txt", "w");
		FILE *g3=fopen("DFEresultsreal3.txt", "w");
		FILE *g4=fopen("DFEresultsreal4.txt", "w");
		FILE *g5=fopen("DFEresultsreal5.txt", "w");
		mapcount=0;

		if ( (g1 == NULL) | (g2 == NULL) | (g3 == NULL) | (g4 == NULL) | (g5 == NULL))
		{
			printf("Error opening DFE results file! \n");
		}
		else
		{
			i=0;
			for(i = 0 ; i < size ; i++) {
				fprintf(g1, "%.12e \n ", creal(reconDFE[i+step] ) );
				fprintf(g2, "%.12e \n ", creal(reconDFE[i+size+step] ) );
				fprintf(g3, "%.12e \n ", creal(reconDFE[i+2*size+step] ) );
				fprintf(g4, "%.12e \n ", creal(reconDFE[i+3*size+step] ) );
				fprintf(g5, "%.12e \n ", creal(reconDFE[i+4*size+step] ) );
			}
		}

		fclose(g1);
		fclose(g2);
		fclose(g3);
		fclose(g4);
		fclose(g5);

		printf("DFE write complete. \n");
	}

	// -------------------------------------------------------------------------//
	// ----------------------CPU sampler filter----------------------------------//
	// -------------------------------------------------------------------------//
	
	printf("Running on CPU... \n");
	current_iteration = iterationmin;
	for( current_iteration = iterationmin; current_iteration < iterationmax+1; current_iteration = 10*current_iteration){
		j = 0;
		for( j = 0; j < repeat_number; j++) {
			sampler_cpu( current_iteration,
				map_number,
				gridwidth,
				size,
				gridsquare_inv,
				tauvalues,
				Tftinvvalues ,
				data,
				Sft,
				Nbar,
				sandtCPU);
		}
	}


	// -------------------------------------------------------------------------//
	// --------------------------Write to File----------------------------------//
	// -------------------------------------------------------------------------//

	if ( print_bool == true) {
		printf("Writing (real) CPU results to file... \n");

		FILE *f1=fopen("CPUresultsreal1.txt", "w");
		FILE *f2=fopen("CPUresultsreal2.txt", "w");
		FILE *f3=fopen("CPUresultsreal3.txt", "w");
		FILE *f4=fopen("CPUresultsreal4.txt", "w");
		FILE *f5=fopen("CPUresultsreal5.txt", "w");
		mapcount=0;

		if ( (f1 == NULL) | (f2 == NULL) | (f3 == NULL) | (f4 == NULL) | (f5 == NULL))
		{
			printf("Error opening CPU results file! \n");
		}
		else
		{
			i=0;
			for(i = 0 ; i < size ; i++) {
				fprintf(f1, "%.12e \n ", creal(sandtCPU[step+i]));
				fprintf(f2, "%.12e \n ", creal(sandtCPU[step+i+size]));
				fprintf(f3, "%.12e \n ", creal(sandtCPU[step+i+2*size]));
				fprintf(f4, "%.12e \n ", creal(sandtCPU[step+i+3*size]));
				fprintf(f5, "%.12e \n ", creal(sandtCPU[step+i+4*size]));
			}
		}

		fclose(f1);
		fclose(f2);
		fclose(f3);
		fclose(f4);
		fclose(f5);
		printf("CPU write complete. \n");
	}



	// -------------------------------------------------------------------------//
	// -----------------------Check Results-------------------------------------//
	// -------------------------------------------------------------------------//

	printf("Comparing CPU and DFE results \n");

	int checkerror1 = 0;
	int checkerror2 = 0;
	int checkerror3 = 0;
	int checkerror4 = 0;
	i=0;
	for( i = 0; i < size*map_number; ++i){

		if (!crelequal(creal(sandtCPU[step+i]), creal(reconDFE[i+step]), 1e-4)) {
			checkerror1++;
		}

		if (!crelequal(creal(sandtCPU[step+i]), creal(reconDFE[i+step]), 1e-3)) {
			checkerror2++;
		}
		if (!crelequal(creal(sandtCPU[step+i]), creal(reconDFE[i+step]), 1e-2)) {
			checkerror3++;
		}
		if (!crelequal(creal(sandtCPU[step+i]), creal(reconDFE[i+step]), 1e-1)) {
			checkerror4++;
		}
	}

	printf("%d data points found at > 0.01%% CPU and DFE mismatch \n", checkerror1);
	printf("%d data points found at > 0.1%% CPU and DFE mismatch \n", checkerror2);
	printf("%d data points found at > 1%% CPU and DFE mismatch \n", checkerror3);
	printf("%d data points found at > 10%% CPU and DFE mismatch \n", checkerror4);
	printf("Done.\n");
	
	/*
	FILE *out1=fopen("outputs.csv", "w");
	i=0;
	for(i=50; i < map_number*size*iterationmax; i=i+map_number*size) {
		fprintf(out1, "%e, %e, %e, %e, %e \n", 
			creal(sandtCPU[i]), 
			creal(sandtCPU[i+size]), 
			creal(sandtCPU[i+2*size]), 
			creal(sandtCPU[i+3*size]), 
			creal(sandtCPU[i+4*size]) );
		//printf("%d, %e, %e, %e\n", i/(map_number*size), creal(sandtCPU[i]), creal(sandtCPU[i+size]), creal(sandtCPU[i+2*size]) );
	}
	fclose(out1);*/

	return 0;
}
