/// \file Main.cpp
/// \brief Main.
///
/// \mainpage 2D Amortized Noise Evaluator
///
/// This project measures the performance of amortized Perlin noise and old-school
/// Perlin noise. It will prompt the user via stdin/stdout
/// for the number of repeats, the smallest and largest texture sizes (which
/// must be a power of 2), and the size delta. The input values are checked for
/// errors. If they pass, the app will measure the average CPU
/// times for generating textures of size s, s+d, s+2d,..., S, where s is
/// the smallest size, S is the largest size, and d is the size delta.
///
/// <b>Experiment 1</b>: Time finite (classical, old school) Perlin noise
/// against finite amortized noise and infinite amortized noise. The results
/// are reported to the console and saved in TimeDataAmortized.txt as a tab-separated
/// text file for easy copy-and-paste into Microsoft(R) Excel(TM). The sizes are
/// listed in reverse order, which is convenient for me but inconvenient for you
/// unless you are familiar with the Sort button on the Data pane in Excel(TM).
///
/// <b>Experiment 2</b>: Time finite (classical, old school) Perlin noise against
/// infinite Perlin noise and infinite smooth Perlin noise. Be very careful not
/// to set the number of repeats too high for this experiment because
/// infinite and infinite smooth Perlin noise are both very, very slow. This is,
/// in fact, the point of the paper (infinite amortized noise is fast
/// when compared to this). The results are reported to the console and saved
/// in TimeDataPerlin.txt as a tab-separated text file for easy copy-and-paste
/// into Microsoft(R) Excel(TM). The sizes are once again listed in reverse order
/// (see above).
///
/// Increasing the number of repeats will increase the quality of the data but
/// also increase the running time. The number of repeats specified by the user
/// is used for the largest texture size only. Experiments on smaller texture 
/// sizes automatically scale the number of repeats so that the run-time is 
/// about the same for all texture sizes. This maximizes data quality without,
/// hopefully, annoying you too much. The actual number of repeats is reported
/// in the text files mentioned above.
///
/// For more details on amortized noise, see Ian Parberry, "Amortized Noise",
/// <em>Journal of Computer Graphics Techniques</em>, Vol. ?, No. ?, pp. ?, To Appear.

// Copyright Ian Parberry, January 2014.
//
// This file is made available under the GNU All-Permissive License.
//
// Copying and distribution of this file, with or without modification,
// are permitted in any medium without royalty provided the copyright
// notice and this notice are preserved.  This file is offered as-is,
// without any warranty.
//
// Created by Ian Parberry, December 2013.
// Last updated May 7, 2014.

#include "defines.h" //OS porting defines 

#include <stdlib.h> //for rand()
#include <stdio.h> //for printf()

#include "FiniteAmortizedNoise2D.h"
#include "InfiniteAmortizedNoise2D.h"
#include "perlin.h"
#include "infiniteperlin.h"
#include "infinitesmoothperlin.h"

#include "CPUtime.h"

const int MAXCELLSIZE2D = 8192; ///< Maximum cell size, MUST be a power of 2.

CFiniteAmortizedNoise2D FiniteAmortizedNoise(MAXCELLSIZE2D); ///< Amortized 2D Perlin noise generator.
CInfiniteAmortizedNoise2D InfiniteAmortizedNoise(MAXCELLSIZE2D, 0); ///< Amortized 2D noise generator.

/// Time a random cell of 2D Amortized noise.
/// \param cell Pointer to array in which to store the noise.
/// \param octave0 Smallest octave.
/// \param octave1 Largest octave.
/// \param size Texture size.
/// \param repeats Number of times to repeat the experiment.
/// \return CPU time in milliseconds.

unsigned int Time2DImprovedAmortizedNoise(float** cell, int octave0, int octave1, int size, int repeats){    
  unsigned int  t = CPUTimeInMilliseconds(); //start time
  for(int i=0; i<repeats; i++)
    InfiniteAmortizedNoise.generate(rand(), rand(), octave0, octave1, size, cell);
  return CPUTimeInMilliseconds() - t;
} //Time2DImprovedAmortizedNoise

/// Time a random cell of 2D Amortized Perlin noise.
/// \param cell Pointer to array in which to store the noise.
/// \param octave0 Smallest octave.
/// \param octave1 Largest octave.
/// \param size Texture size.
/// \param repeats Number of times to repeat the experiment.
/// \return CPU time in milliseconds.

unsigned int Time2DAmortizedNoise(float** cell, int octave0, int octave1, int size, int repeats){    
  unsigned int  t = CPUTimeInMilliseconds(); //start time
  for(int i=0; i<repeats; i++)
    FiniteAmortizedNoise.generate(rand(), rand(), octave0, octave1, size, cell);
  return CPUTimeInMilliseconds() - t;
} //Time2DAmortizedNoise

/// Time a random cell of 2D Perlin noise.
/// \param cell Pointer to array in which to store the noise.
/// \param octave0 Smallest octave.
/// \param octave1 Largest octave.
/// \param size Texture size.
/// \param repeats Number of times to repeat the experiment.
/// \param noise Noise function to be timed, will be either vanilla, infinite, or infinite smoothed Perlin noise.
/// \return CPU time in milliseconds.

unsigned int Time2DPerlinNoise(float** cell, int octave0, int octave1, int size, int repeats, float (*noise)(float, float, float, float, int)){
  initPerlin2D(); //initialize Perlin noise   

  float x0 = (float)rand(), z0 = (float)rand(); //random start location

  int n = 2*size;
  for(int i=0; i<octave0; i++) //skip over unused octaves
    n /= 2;
  const float delta = 1.0f/(n - 1); //distance between points

  unsigned int  t = CPUTimeInMilliseconds(); //start time
  
  for(int k=0; k<repeats; k++){
    float x = x0; //current x location
    for(int i=0; i<size; i++){
      float z = z0; //current z location
      for(int j=0; j<size; j++){
        cell[i][j] = (*noise)(x, z, 0.5f, 2.0f, octave1 - octave0 + 1);
        z += delta;
      } //for
      x += delta;
    } //for
  } //for
  
  return CPUTimeInMilliseconds() - t;
} //Time2DPerlinNoise

/// Time 2D Perlin noise against infinite and infinite smooth Perlin noise.

void TestPerlinNoiseVariants(){  
  printf("\nTiming Perlin noise against infinite and infinite smooth Perlin noise.\n");

  FILE* datafile = fopen("TimeDataPerlin.txt", "wt");

  if(datafile)
    fprintf(datafile, "n\tPerlin\tInf.\tRatio1\tInf. Sm.\tRatio2\tRepeats\n");

  int nRepeats=-1, nSmallSize=-1, nLargeSize=-1, nSizeDelta=-1;
  
  while(nRepeats <= 0){
    printf("Number of repeats:\n> ");
    scanf("%d", &nRepeats); 
    if(nRepeats <= 0)
      printf("Please enter a number greater than zero.");
  } //while
  
  while(nSmallSize <= 0){
    printf("Smallest size:\n> ");
    scanf("%d", &nSmallSize);
    if(nSmallSize <= 0)
      printf("Please enter a number greater than zero.");
  } //while
  
  while(nLargeSize < nSmallSize){
    printf("Largest size:\n> ");
    scanf("%d", &nLargeSize);
    if(nLargeSize < nSmallSize)
      printf("Please enter a number that is at least %d.", nSmallSize);
  } //while
  
  while(nSizeDelta <= 0){
    printf("Size delta:\n> ");
    scanf("%d", &nSizeDelta);   
    if(nSizeDelta <= 0)
      printf("Please enter a number greater than zero.");
  } //while
  
  printf("n\tFinite\tInfinite\tSmooth\tRatio1\tRatio2\tRepeats\n");

  unsigned int nVanillaTime, nInfiniteTime, nSmoothInfiniteTime;

  for(int size=nLargeSize; size>=nSmallSize; size-=nSizeDelta){
    printf("%d\t", size);
    if(datafile)fprintf(datafile, "%d\t", size);
    
    //allocate space to store noise cells
    float** cell = new float* [size];
    for(int i=0; i<size; i++)
      cell[i] = new float [size];

    //vanilla Perlin noise
    nVanillaTime = Time2DPerlinNoise(cell, 1, 1, size, nRepeats, PerlinNoise2D);
    printf("%0.4f\t", (float)nVanillaTime/nRepeats);
    if(datafile)fprintf(datafile, "%0.4f\t", (float)nVanillaTime/nRepeats);

    //infinite Perlin noise
    nInfiniteTime = Time2DPerlinNoise(cell, 1, 1, size, nRepeats, InfinitePerlinNoise2D);
    printf("%0.4f\t", (float)nInfiniteTime/nRepeats);
    if(datafile)fprintf(datafile, "%0.4f\t", (float)nInfiniteTime/nRepeats);

    //smooth infinite Perlin noise
    nSmoothInfiniteTime = Time2DPerlinNoise(cell, 1, 1, size, nRepeats, InfiniteSmoothPerlinNoise2D);
    printf("%0.4f\t", (float)nSmoothInfiniteTime/nRepeats);
    if(datafile)fprintf(datafile, "%0.4f\t", (float)nSmoothInfiniteTime/nRepeats);
    
    //report ratios
    printf("%0.1f\t%0.1f\t", (float)nInfiniteTime/(float)nVanillaTime, (float)nSmoothInfiniteTime/(float)nVanillaTime);
    if(datafile)fprintf(datafile, "%0.1f\t%0.1f\t", (float)nInfiniteTime/(float)nVanillaTime, (float)nSmoothInfiniteTime/(float)nVanillaTime);

    //report number of repeats
    printf("%d\n", nRepeats);
    if(datafile)fprintf(datafile, "%d\n", nRepeats);

    //compute next repeat count, which can and should be larger for smaller cells
    if(size != nSizeDelta)
      nRepeats = size*size*nRepeats/((size - nSizeDelta)*(size - nSizeDelta));
    else nRepeats = size*size*nRepeats;

    //deallocate noise cell space
    for(int i=0; i<size; i++)
      delete [] cell[i];
    delete [] cell;
  } //for

  fclose(datafile);
} //TestPerlinNoiseVariants


/// Time 2D Perlin noise against finite and infinite amortized noise.

void TestAmortizedNoiseVariants(){  
  printf("\nTiming Perlin noise against finite and infinite Amortized noise.\n");
  
  FILE* datafile = fopen("TimeDataAmortized.txt", "wt");

  if(datafile)fprintf(datafile, "n\tPerlin\tAmortized\tRatio\tInfinite Amortized\tRatio\tRepeats\n");

  int nRepeats=-1, nSmallSize=-1, nLargeSize=-1, nSizeDelta=-1;
  
  while(nRepeats <= 0){
    printf("Number of repeats:\n> ");
    scanf("%d", &nRepeats); 
    if(nRepeats <= 0)
      printf("Please enter a number greater than zero.");
  } //while
  
  while(nSmallSize <= 0){
    printf("Smallest size:\n> ");
    scanf("%d", &nSmallSize);
    if(nSmallSize <= 0)
      printf("Please enter a number greater than zero.");
  } //while
  
  while(nLargeSize < nSmallSize){
    printf("Largest size:\n> ");
    scanf("%d", &nLargeSize);
    if(nLargeSize < nSmallSize)
      printf("Please enter a number that is at least %d.", nSmallSize);
  } //while
  
  while(nSizeDelta <= 0){
    printf("Size delta:\n> ");
    scanf("%d", &nSizeDelta);   
    if(nSizeDelta <= 0)
      printf("Please enter a number greater than zero.");
  } //while
  
  printf("n\tPerlin\tAmort.\tRatio\tInf.\tRatio\tRepeats\n");

  unsigned int nVanillaTime, nAmortizedTime, nInfiniteTime;

  for(int size=nLargeSize; size>=nSmallSize; size-=nSizeDelta){
    printf("%d\t", size);
    if(datafile)fprintf(datafile, "%d\t", size);
    
    //allocate space to store noise cells
    float** cell = new float* [size];
    for(int i=0; i<size; i++)
      cell[i] = new float [size];

    //vanilla Perlin noise
    nVanillaTime = Time2DPerlinNoise(cell, 1, 1, size, nRepeats, PerlinNoise2D);
    printf("%0.4f\t", (float)nVanillaTime/nRepeats);
    if(datafile)fprintf(datafile, "%0.4f\t", (float)nVanillaTime/nRepeats);

    //finite amortized noise
    nAmortizedTime = Time2DAmortizedNoise(cell, 1, 1, size, nRepeats);
    printf("%0.4f\t", (float)nAmortizedTime/nRepeats);
    if(datafile)fprintf(datafile, "%0.4f\t", (float)nAmortizedTime/nRepeats);

     //report ratio
    printf("%0.1f\t", (float)nVanillaTime/(float)nAmortizedTime);
    if(datafile)fprintf(datafile, "%0.1f\t", (float)nVanillaTime/(float)nAmortizedTime);
  
    //infinite amortized noise
    nInfiniteTime = Time2DImprovedAmortizedNoise(cell, 1, 1, size, nRepeats);
    printf("%0.4f\t", (float)nInfiniteTime/nRepeats);
    if(datafile)fprintf(datafile, "%0.4f\t", (float)nInfiniteTime/nRepeats);

    //report ratio
    printf("%0.1f\t", (float)nVanillaTime/(float)nInfiniteTime);
    if(datafile)fprintf(datafile, "%0.1f\t", (float)nVanillaTime/(float)nInfiniteTime);

    //report number of repeats
    printf("%d\n", nRepeats);
    if(datafile)fprintf(datafile, "%d\n", nRepeats);

    //compute next repeat count, which can and should be larger for smaller cells
    if(size != nSizeDelta)
      nRepeats = size*size*nRepeats/((size - nSizeDelta)*(size - nSizeDelta));
    else nRepeats = size*size*nRepeats;

    //deallocate noise cell space
    for(int i=0; i<size; i++)
      delete [] cell[i];
    delete [] cell;
  } //for

  fclose(datafile);
} //TestAmortizedNoiseVariants

/// Main.
/// \param argc Argument count
/// \param argv Arguments.
/// \return 0 for success, 1 for failure.

int main(int argc, char *argv[]){ 
  //print banner to console
  printf("Amortized 2D Noise CPU Time Experiments, Ian Parberry, 2014.\n");
  printf("Timing the generation of 2D noise with persistence 0.5 and lacunarity 2.\n");
  printf("--------------------------------------------------------------\n");
  
  srand(0); //seed random number generator
  TestAmortizedNoiseVariants(); 
  TestPerlinNoiseVariants();  

#if defined(_MSC_VER) //Windows Visual Studio 
  //wait for user keystroke and exit
  printf("\nHit Almost Any Key to Exit...\n");
  _getch();
#endif

  return 0;
} //main
