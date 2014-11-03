
/********************************************************
 ********************* HEADER ***************************
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


// GSL - GNU Scientific Library  
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


#define FALSE 0
#define TRUE  1



// data set structure
typedef struct{
  int n; // # data points
  float x[13000]; // input data & error
  float x_err[13000];  
  float y[13000]; 
  float y_err[13000];  
  float x_synth[13000]; // synthed data OR bootstrapped/resampled data
  float y_synth[13000];
  float x_rank[13000]; // rank of data
  float y_rank[13000]; 
  char name[50]; // file name
} dataset_t;


typedef struct{
  int infile;
  char data_file[50];
  int iter;
  int outfile;
  int seed;
  int method;
} options_t;




int parse_options(int argc, char *argv[], options_t *opt);
void show_usage(void);
void credits(void);

int print_head(char *string, char symbol, int bl_start, int bl_end); 

void read_data(char *filename, dataset_t *data);
void synth(int flag, dataset_t *data);
void synth2(int flag, int boot, dataset_t *data);
void bootstrap(int flag, dataset_t *data);
void rank(dataset_t *data);
float spear(dataset_t *data);

void full_monte(int flag, int iter, int method);


float mean(float x[], int n); 
float std_dev(float x[], int n, float mean);

float fisher(float r);
float zscore(int n, float Fr);
float tscore(int n, float r);



/*************************************************** 
 *              Variable declarations 
 ***************************************************/

dataset_t data;

gsl_rng * ridum; // for use with GSL random num generator
gsl_rng * ridum2; // for use with GSL random num generator
