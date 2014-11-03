
/*************************************************************************
 *  
 *  For the calculation of Spearman Rank coefficient of data 
 *   and estimation of errors via Monte Carlo
 *
 *  Author: Peter Curran <peter.a.curran@gmail.com>
 *
 *  Version 0.0 (13/05/2011; CEA-Saclay)
 *  v0.2 Added bootstrap method (06/06/2011; CEA-Saclay)
 *  v0.3 NumRec replaced by  GSL - GNU Scientific Library (30/10/2014; ICRAR) 
 *       Added Student's t routine (not implemented in MC)
 *************************************************************************/

 /*************************************************************************
 * TO COMPILE:   
 *        gcc -Wall -o mcspearman mcspearman_v0.3.c -L/usr/lib -lgsl -lgslcblas -lm 
 * where -L/usr/lib is your GSL directory 
 **************************************************************************/


#define version 0.3 
#define date "20/10/2014"
#define email "peter.curran@curtin.edu.au"
#define institute "Curtin University/ICRAR"

#include "mcspearman.h"  // GLOBAL HEADER FILE



/***********************************************************************
 ************************ MAIN FUNCTION ********************************
 ***********************************************************************/

int main(int argc, char *argv[])
{

  int i; 
  float temp;

  options_t opt;
  parse_options(argc, argv, &opt);

  char txt[30] = " version ";
  char version_str[10];
  sprintf(version_str, "%.1f", version);
  strcat(txt, version_str);
  strcat(txt, " ");
 
  print_head(NULL, '*', 2, 0);
  print_head(" Monte Carlo Spearman Rank Test ", '*', 0, 0);
  print_head(txt, '*', 0, 0);
  print_head(NULL, '*', 0, 2);




  // GSL RANDOM NUMBER SET-UP
  /* create a generator chosen by the 
     environment variable GSL_RNG_TYPE */

  const gsl_rng_type * T;
  //gsl_rng * r; // GLOBAL VARIABLE

  gsl_rng_env_setup();
  T = gsl_rng_default;
  ridum = gsl_rng_alloc (T);
  ridum2 = gsl_rng_alloc (T);

  //gsl_rng_free (ridum);
  //gsl_rng_free (ridum2);


  // ********* DEPRECATED *********
  // negative RANDOM NUMBER GENERATOR SEED
  // idum is relegated with NumRec
  long idum; // for use with random num gen
  //long idum2 -1145; 
  if(opt.seed != FALSE){
    idum = -1*abs(opt.seed);
    fprintf(stdout, "\n\nRandom number generator seed = %ld [DEPRECATED]\n\n",idum);
  }
   

  // READ IN DATA
  if(opt.infile == FALSE){
    fprintf(stdout, "Enter data file (e.g. abc123.data): ");
    scanf("%s", data.name); 
  } else {
    strcpy (data.name, opt.data_file);
  }
  read_data(data.name, &data);



  synth2(FALSE, FALSE, &data);  // Fill data._synth enteries
  rank(&data);

  fprintf(stdout, "\nSpearman Rank coefficient, z-score & student's t\n");

  float spear_rank;
  spear_rank = spear(&data);
  fprintf(stdout, "SR =  %.4f \n", spear_rank);

  temp = fisher(spear_rank); 
  //fprintf(stdout, "Fr =  %.4f \n", temp);

  temp = zscore(data.n,temp);
  fprintf(stdout, "z  =  %.4f \n", temp);

  temp = zscore(data.n,spear_rank);
  fprintf(stdout, "t  =  %.4f \n", temp);
  

  // Monte Carlo Analysis
  if(opt.iter == FALSE){
    fprintf(stdout, "\nHow many MC iterations in error analysis? (<1 will exit)\n");
    scanf("%i", &i); } else{
    i=opt.iter; 
  }
  if (i>0){
    full_monte(opt.outfile, i, opt.method);
  } 
  


 /************************* OUTRO ******************************/
  fprintf(stdout, "\n Thank you, come again... \n\n");
  return 0;
}




/***************** full_monte ************************
 * Monte Carlo error analysis
 *************************************************/

void full_monte(int flag, int iter, int method)
{
  int i;
  char filename[50];
  FILE *logfile;  

  float SR[10000];
  float Fr[10000];
  float ZS[10000];

  float average; 
  float error;

  int synthflag = TRUE;
  int bootflag = FALSE;


  printf("\nFull Monte Carlo error analysis:\n");

  if (method == 2){
    printf("Method # %i (Pertubation method) \n", method);
    synthflag = TRUE;
    bootflag = FALSE;
  }else if (method == 1){
    printf("Method # %i (Bootstrap method) \n", method);
    synthflag = FALSE;
    bootflag = TRUE;
  } else if (method == 3){
    printf("Method # %i (Composite method) \n", method);
    synthflag = TRUE;
    bootflag = TRUE;
  } else {
    printf("Method # %i is not a valid method \n", method);    
    return;
  }



  if (flag == TRUE){
    // Opening of logfile:
    sprintf(filename, "full_monte_x%i_%s", iter, data.name);
    logfile = fopen (filename, "w");
    if ( logfile == NULL ){
      printf("Error: cannot open '%s' - sorry :-)\n", filename);
    }else{ 
      printf("File '%s' opened successfully... & fitting...\n", filename);
    }
    
    for(i=1; i <= iter; i++){
      synth2(synthflag, bootflag, &data);/////////synth(TRUE, &data);
      rank(&data);
      
      SR[i] = spear(&data); 
      Fr[i] = fisher(SR[i]); 
      ZS[i] = zscore(data.n,Fr[i]);
      
      fprintf(logfile, "%.4f\t%.4f\t \n", SR[i], ZS[i]);
      //fprintf(logfile, "%.4f\t%.4f\t%.4f\t \n", SR[i], Fr[i], ZS[i]);
      //fprintf(stdout, "%.4f\t%.4f\t%.4f\t \n", SR[i], Fr[i], ZS[i]);
    } 
    fclose(logfile);  
  }else {

    fprintf(stdout, "\nOutput flag is %i (FALSE): No log file \n", flag);
    for(i=1; i <= iter; i++){
     synth2(synthflag, bootflag, &data); /////////synth(TRUE, &data);
      rank(&data);
      
      SR[i] = spear(&data); 
      Fr[i] = fisher(SR[i]); 
      ZS[i] = zscore(data.n,Fr[i]);
    } 
  }
  


  // Calculate averages & deviations
  printf("\nResults for %i iterations: \n", iter);
  printf("(mean +/- std)\n");
  
  average = mean(SR, iter);  
  error = std_dev(SR, iter, average);
  fprintf(stdout, "SR = %.4f +/- %.4f (ZS = %.2f) \n", average, error, zscore(data.n, fisher(average)) );
 
  //average = mean(Fr, iter);  
  //error = std_dev(Fr, iter, average);
  //fprintf(stdout, "Fr = %.3f +/- %.3f \n", average, error);

  average = mean(ZS, iter);  
  error = std_dev(ZS, iter, average);
  fprintf(stdout, "ZS = %.3f +/- %.3f \n", average, error);

 
 

  return;
}



/*************** READ DATA  **********************
 ********* Reads in the data from file ***********
 *************************************************/

void read_data(char *filename, dataset_t *data)
{
  int i;
  FILE *datafile;

  datafile = fopen(filename, "r");
  
  if(datafile==NULL){
    fprintf(stdout, "ERROR: cannot open file: %s\n", filename);
    fprintf(stdout, "Please ensure that filename is correct and a data file exists \n\n\n");
    exit(0);
  }else{
    fprintf(stdout, "Reading Data from '%s'...\n", filename);
    i=1;

    while((!feof(datafile))   ){
      fscanf(datafile, "%f %f %f %f", &data->x[i], &data->x_err[i], &data->y[i], &data->y_err[i]);
      i++;
    }   
    data->n = i-1;
    fprintf(stdout, "done.\n");
  }
    
  // test if last line was blank
  if(data->x[data->n] == 0 && data->y[data->n] == 0 && data->x_err[data->n] == 0 && data->y_err[data->n] == 0){
    data->n = data->n - 1;
  }


  //fprintf(stdout, "X\tX err  \tY\tY err \n");
  //for(i=1; i<=data->n; i++) {
  //   fprintf(stdout, "%.3f\t%.3f\t%.4f\t%.4f \n", data->x[i], data->x_err[i], data->y[i], data->y_err[i]);
  //}
  
  fclose(datafile); 
  fprintf(stdout, "Number of data points =  %i \n", data->n);


  return;
}



/******************* SPEAR *************************
 * returns Spearman coeficient of ranked data
 ***************************************************/

float spear(dataset_t *data) 
{
  int i;
  float coef;
  float temp;
  float dsum = 0;

  for(i=1; i <= data->n; i++){
    temp = data->x_rank[i] - data->y_rank[i];
    dsum = dsum + temp*temp;
    // fprintf(stdout, "%.1f %.1f  \n", temp, data->d2[i]);
  }

  coef = 1 - (  (6 * dsum ) /  ( data->n*(data->n*data->n -1) ) )  ;
  // fprintf(stdout, "SR coeff is %.4f  \n", coef);

  return coef;
}


/******************* RANK *************************
 * ranks data set 
 ***************************************************/

void rank(dataset_t *data) 
{
  int i,j;
    
  for(i=1; i <=data->n; i++){
    data->x_rank[i] =  1;
    data->y_rank[i] =  1;


    for(j=1; j <=data->n; j++){
      if(data->x_synth[j] < data->x_synth[i]){
	data->x_rank[i] = data->x_rank[i] +1;	
      } else if (data->x_synth[j] == data->x_synth[i]  && i != j ){
	data->x_rank[i] = data->x_rank[i] + 0.5;	
      } 

      if(data->y_synth[j] < data->y_synth[i]){
	data->y_rank[i] = data->y_rank[i] +1;	
      } else if (data->y_synth[j] == data->y_synth[i]  && i != j ){
	data->y_rank[i] = data->y_rank[i] + 0.5;	
      }

    }

    //fprintf(stdout, "%.1f  %.1f  \n", data->x_rank[i], data->y_rank[i]);

  }

  return;
}



/******************* SYNTH *************************
 * synthesizes data set
 ***************************************************/

void synth(int flag, dataset_t *data) 
{
  int i;
  float gasdevran;
  
  //if(flag == TRUE){
  //  fprintf(stdout, "synth data\n");
  //} else if (flag == FALSE){
  //  fprintf(stdout, "non-synth data\n");
  //} else {
  //  fprintf(stdout, " synth-flag error!\n");
  //}
  
  for(i=1; i <=data->n; i++){
    // perturb the synthetic data
    //gasdevran = gasdev(&idum);
    gasdevran = gsl_ran_gaussian(ridum, 1.0);

    data->x_synth[i] = data->x[i] + flag*gasdevran*data->x_err[i];
    //gasdevran = gasdev(&idum); 
    gasdevran = gsl_ran_gaussian(ridum, 1.0);
    data->y_synth[i] = data->y[i] + flag*gasdevran*data->y_err[i];
    //fprintf(stdout, "%.3f  %.4f  \n", data->x_synth[i], data->y_synth[i]);
  }

  return;
}


void synth2(int flag, int boot, dataset_t *data) 
{
  int i;
  float random;
  int random_int;
  int n =   data->n;
  float gasdevran;

  
  if (boot == FALSE){ // the perturb only method
    for(i=1; i <=data->n; i++){
      // perturb the synthetic data
      //gasdevran = gasdev(&idum); 
      gasdevran = gsl_ran_gaussian(ridum, 1.0);

      data->x_synth[i] = data->x[i] + flag*gasdevran*data->x_err[i];
      //gasdevran = gasdev(&idum);
      gasdevran = gsl_ran_gaussian(ridum, 1.0);
      data->y_synth[i] = data->y[i] + flag*gasdevran*data->y_err[i];
      //fprintf(stdout, "%.3f  %.4f  \n", data->x_synth[i], data->y_synth[i]);
    }
  } else{ // boot=TRUE, the resampled bootstrap method
    for(i=1; i <=n; i++){
      
      //random = ran1(&idum2);
      random = gsl_rng_uniform (ridum2);
      random_int = (random*n)+1;
      
      //gasdevran = gasdev(&idum);
      gasdevran = gsl_ran_gaussian(ridum, 1.0);
      data->x_synth[i] = data->x[random_int]+ flag*gasdevran*data->x_err[random_int];
      
      //gasdevran = gasdev(&idum); 
      gasdevran = gsl_ran_gaussian(ridum, 1.0);
      data->y_synth[i] = data->y[random_int] + flag*gasdevran*data->y_err[random_int];
      
      //fprintf(stdout, "%i %.3f %.3f ; %i %.3f %.3f\n", i,  data->x[i], data->y[i], random_int, data->x_synth[i], data->y_synth[i] );
    }
    
  }
  
  
  return;
}


/******************* BOOTSTRAP *************************
 * Creates resampled [& perturbed] data set
 ***************************************************/

void bootstrap(int flag, dataset_t *data) 
{
  int i;
  float random;
  int random_int;
  int n =   data->n;
  float gasdevran;


  for(i=1; i <=n; i++){

    //random = ran1(&idum2); 
    random = gsl_rng_uniform (ridum2);
    random_int = (random*n)+1;

    //gasdevran = gasdev(&idum); 
    gasdevran = gsl_ran_gaussian(ridum, 1.0);
    data->x_synth[i] = data->x[random_int]+ flag*gasdevran*data->x_err[random_int];
    
    //gasdevran = gasdev(&idum); 
    gasdevran = gsl_ran_gaussian(ridum, 1.0);
    data->y_synth[i] = data->y[random_int] + flag*gasdevran*data->y_err[random_int];

    //fprintf(stdout, "%i %.3f %.3f ; %i %.3f %.3f\n", i,  data->x[i], data->y[i], random_int, data->x_synth[i], data->y_synth[i] );
  }

  return;
}




/***************** FISHER ************************
 * Fisher transformation
 *************************************************/

float fisher(float r)
{  
 return 0.5 * log10( (1+r)/(1-r) );
}


/***************** Z-SCORE ************************
 * Standard/Z Score
 *************************************************/

float zscore(int n, float Fr)
{  
 return sqrt( (n-3)/1.06 ) * Fr; 
}


/*************** STUDENTS T **********************
 * Student's t test
 *************************************************/

float tscore(int n, float r)
{  
  return r / sqrt((1-r*r)/(n-2)) ;
}





/***************** MEAN ************************
 * calculates mean value of an array
 *************************************************/

float mean(float x[], int n)
{
  float total = 0; 
  int i; 
  
  for(i=1; i<=n; i++){
    total = total + x[i];
  }
  
  return total/n;
}


/***************** STD_DEV ***********************
 * calculates standard deviation of array
 *************************************************/

float std_dev(float x[], int n, float mean)
{
  float total = 0;
  float sigma;
  int i; 

  for(i=1; i<=n; i++){
    total = total + x[i]*x[i];
  }  
  sigma = pow(total/n - mean*mean, 0.5);  
  

  return sigma;
}




/************************************************************
 * Prints standard header lines 
 ************************************************************/
int print_head(char *string, char symbol, int bl_start, int bl_end)
{
  const int max_len = 74;

  int nstar;
  int l_nstar, r_nstar;
  int i;

  nstar = max_len;
  if(string!=NULL){
    nstar = nstar - strlen(string);
  }
  l_nstar = 0.5*nstar;
  r_nstar = nstar - l_nstar;

  for(i=0; i<bl_start; i++){
    fprintf(stdout, "\n");
  }

  for(i=0; i<l_nstar; i++){
    fprintf(stdout, "%c", symbol);
  }

  if(string!=NULL){
    fprintf(stdout, "%s", string);
  }

  for(i=0; i<r_nstar; i++){
    fprintf(stdout, "%c", symbol);
  }
  fprintf(stdout, "\n");

  for(i=0; i<bl_end; i++){
    fprintf(stdout, "\n");
  }

  return 0;
}





/************************************************************
 * Reads in command line options
 ************************************************************/
int parse_options(int argc, char *argv[], options_t *opt)
{
  int i;

  //set the default values
  opt->infile = FALSE;
  opt->iter = FALSE; 
  opt->outfile = FALSE;
  opt->seed = FALSE;
  opt->method = 2;

  for(i=1; i<argc; i++){
    if(!(strcmp(argv[i], "-infile"))){
      opt->infile = TRUE;
      strcpy (opt->data_file, argv[i+1]);
      i++;
      
    }else if(!(strcmp(argv[i], "-output"))){
      opt->outfile = TRUE;
      i++;
   

    } else if(!(strcmp(argv[i], "-seed"))){
      opt->seed = atof(argv[i+1]);
      i++;

    } else if(!(strcmp(argv[i], "-method"))){
      opt->method = atof(argv[i+1]);
      i++; 
      
    }else if(!(strcmp(argv[i], "-i"))){
      opt->iter = atof(argv[i+1]);
      i++;

    }else if(!(strcmp(argv[i], "-h"))){
      credits();
      exit(0);
    }else if(!(strcmp(argv[i], "-o"))){
      show_usage();
      exit(0);
    }else{
      fprintf(stderr, "Unknown option %s\n", argv[i]);
      show_usage();
      exit(1);}
  }

  return 0;
}



/************************************************************
 * Lists different functions available 
 * (either if an incorrect one is tried, or option -h)
 ************************************************************/

void show_usage(void)
{
  fprintf(stdout, "USE: mcspearman [OPTIONS]\n");
  fprintf(stdout, "\nOPTIONS:\n");

  fprintf(stdout, "  -infile ####  - Input data file name\n");
  fprintf(stdout, "  -i ###        - Number of iterations for MC error analysis (<0 skips) \n");


  fprintf(stdout, "  -method #     - Method: bootstrap (1), perturb (2; default), composite (3)\n");
  fprintf(stdout, "  -output       - Outputs MC results to file \n");
  //fprintf(stdout, "  -seed ####    - Random number generator seed (integer) \n"); 
  // DEPRECATED option


  fprintf(stdout, "  -o            - options \n");
  fprintf(stdout, "  -h            - credits & help\n");

  fprintf(stdout, "\nInput file format: X X_err   Y Y_err \n\n");

  fprintf(stdout, "\nExamples:\n");
  fprintf(stdout, "  >mcspearman -infile test.data  -i 100 -method 2\n");
  fprintf(stdout, "  >nice mcspearman -infile test.data -i 100 >& outfile.txt & \n");




  
  fprintf(stdout, "\n");
}

void credits(void)
{
  // print_head(NULL, '*', 1, 0);
  //print_head(" GRBz ", '*', 0, 0);
  //print_head(NULL, '*', 0, 2);

  fprintf(stdout, " \n");
  fprintf(stdout, " MC Spearman \n");
  fprintf(stdout, " =========== \n");
  fprintf(stdout, " A program to compute the Spearman's rank correlation coefficient, \n");
  fprintf(stdout, " including a Monte Carlo error analysis. Produces a 'full_monte_***.data \n");
  fprintf(stdout, " output file for histograms. \n\n");
 
  fprintf(stdout, " Current Version: v%.1f (%s)\n", version, date);
  fprintf(stdout, " Author: Peter A. Curran, %s \n", institute);
  fprintf(stdout, " <%s> \n <peter.a.curran@gmail.com> \n", email);

  fprintf(stdout, "\n\n");
  
  show_usage();
  
}
