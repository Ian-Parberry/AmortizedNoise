/// \file Main.cpp
/// \brief Main.
///
/// \mainpage 2D Amortized Noise Generator
///
/// This project demonstrates how the 2D amortized noise algorithm can be used to produce
/// an infinite greyscale tiled texture. It will prompt the user via stdin/stdout
/// for the texture size (which must be a power of 2), a seed value, the largest and smallest octaves,
/// and tile coordinates in row-column format. The input values are checked for errors. 
/// If they pass, the tile texture will be saved in a png file (windows) or tga file (other OS).
/// The amount of user CPU time used in generating the noise is reported to stdout.
/// Saved files will be named  i[F]s[N]o[LS]s[D]r[R]c[C].X (for example, i1s256o25s5432r54c33.png), where:
/// <center><table>
/// <tr><td>F<td>0 for finite or 1 for infinite
/// <tr><td>N<td>tile side
/// <tr><td>L<td>largest octave
/// <tr><td>S<td>smallest octave
/// <tr><td>D<td>seed
/// <tr><td>R<td>row
/// <tr><td>C<td>column
/// <tr><td>X<td>png or tga
/// </table></center>
///
/// For more details on amortized noise, see Ian Parberry, "Amortized Noise",
/// <em>Journal of Computer Graphics Techniques</em>, Vol. ?, No. ?, pp. ?, To Appear.

// Copyright Ian Parberry, September 2013.
//
// This file is made available under the GNU All-Permissive License.
// Copying and distribution of this file, with or without modification,
// are permitted in any medium without royalty provided the copyright
// notice and this notice are preserved.  This file is offered as-is,
// without any warranty.
//
// Created by Ian Parberry, September 2013.
// Last updated May 9, 2014.

#include "defines.h" //OS porting defines

#include <stdlib.h> //for rand()
#include <stdio.h> //for printf()

#include "InfiniteAmortizedNoise2D.h"
#include "FiniteAmortizedNoise2D.h"
#include "CPUtime.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION ///< For platform-neutral image file output.
#include "stb_image_write.h"

CFiniteAmortizedNoise2D* g_pFiniteAmortizedNoise = NULL; ///< Pointer to the finite amortized 2D noise generator. 
CInfiniteAmortizedNoise2D* g_pInfiniteAmortizedNoise = NULL; ///< Pointer to the infinite amortized 2D noise generator. 

bool g_bInfinite = true; ///< Whether to use infinite noise: true for infinite noise, false for finite noise.

/// Save a cell of 2D noise as a png file.
/// \param cell Pointer to array in which to store the noise.
/// \param scale Rescale by this to get in the range -1..1.
/// \param basefilename Root of the file name under which to store the noise.
/// \param n Cell size.

void Save2DNoise(float** cell, const int n, const float scale, const char* basefilename){ 
  char filename[MAX_PATH];
  sprintf(filename, "%s.png", basefilename); 
  printf("Saving to %dx%d png file %s\n\n", n, n, filename);

  unsigned char* pixelbuffer = new unsigned char[4*n*n]; 

  //convert from noise cell to pixel buffer
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j){
      unsigned int b = (unsigned int)(128.0f*(cell[i][j]*scale + 1.0f));
      ((unsigned int*)pixelbuffer)[i*n + j] = 0xFF000000U | b<<16 | b<<8 | b;       
    } //for 

  stbi_write_png(filename, n, n, 4, pixelbuffer, 4*n); //save byte buffer as a png file

  delete [] pixelbuffer;
} //Save2DNoise

/// Generate a random cell of 2D noise.
/// \param cell Pointer to array in which to store the noise.
/// \param x Tile column index.
/// \param y Tile row index.
/// \param m0 Largest octave.
/// \param m1 Smallest octave.
/// \param n Cell size.

float Generate2DNoise(float** cell, const int x, const int y, const int m0, const int m1, const int n){ 
  printf("Generating %d octaves of 2D noise...\n", m1 - m0 + 1);
  int t = CPUTimeInMilliseconds();
  float scale = 1.0f;
  if(g_bInfinite)
    scale = g_pInfiniteAmortizedNoise->generate(y, x, m0, m1, n, cell);
  else scale = g_pFiniteAmortizedNoise->generate(y, x, m0, m1, n, cell);

  t = CPUTimeInMilliseconds() - t;
  printf("Generated %d points in %0.2f seconds CPU time.\n", n*n, t/1000.0f);
  return scale;
} //Generate2DNoise

/// Generate a cell of 2D noise and save it as an image file.
/// \param nRow Tile row index.
/// \param nCol Tile column index.
/// \param m0 Largest octave.
/// \param m1 Smallest octave.
/// \param n Cell size.
/// \param seed Hash function seed.

void GenerateAndSave2DNoise(int nRow, int nCol, const int m0, const int m1, const int n, const int seed){ 
  g_pFiniteAmortizedNoise = new CFiniteAmortizedNoise2D(n);  
  g_pInfiniteAmortizedNoise = new CInfiniteAmortizedNoise2D(n, seed);

  //construct file name
  char filename[MAX_PATH];
  sprintf(filename, "i%1ds%do%1d%1ds%dr%dc%d", g_bInfinite, n, m0, m1, seed, nRow, nCol);

  //adjust for tile size of smallest octave
  for(int i=1; i<m0; i++){ 
    nCol *= 2; nRow *= 2;
  } //for

  //allocate space to store noise cells
  float** cell = new float* [n];
  for(int i=0; i<n; i++)
    cell[i] = new float [n]; 

  //generate and save noise
  float scale = Generate2DNoise(cell, nCol, nRow, m0, m1, n);
  Save2DNoise(cell, n, scale, filename);

  //deallocate noise cell space
  for(int i=0; i<n; i++)
    delete [] cell[i];
  delete [] cell;

  delete g_pFiniteAmortizedNoise;
  delete g_pInfiniteAmortizedNoise;
} //GenerateAndSave2DNoise

/// Main.
/// \param argc Argument count
/// \param argv Arguments.
/// \return 0 for success, 1 for failure.

int main(int argc, char *argv[]){ 
  printf("Amortized 2D Noise Generator, Ian Parberry, 2014\n");
  printf("--------------------------------------------------------------\n\n");

  int response = 1;
  do{
    printf("Enter 0 for finite noise and 1 for infinite noise:\n> ");
    scanf("%d", &response);
  }while(response != 0 && response != 1);
  g_bInfinite = response == 1;

  int n = 2048;
  printf("Texture size (must be a power of 2 and at least 2):\n> ");
  scanf("%d", &n);

  if((n >= 2) && !(n & (n - 1))){ //n is a power of 2 and at least 2
    unsigned int seed = 0;
    printf("Hash seed:\n> "); scanf("%d", &seed);
    srand(seed); //for finite noise

    int log2n = 0; //log base 2 of n
    for(int temp=n; temp>1; temp=temp>>1) //compute log base 2 of n
      log2n++;

    int m0 = 4, m1 = 7; //largest and smallest octaves
    printf("Largest octave (must be at least 1 and at most %d):\n> ", log2n);
    scanf("%d", &m0);

    if(m0 >= 1 && m0 <= log2n){ //largest octave is correct
      printf("Smallest octave (must be at least %d and at most %d):\n> ", m0, log2n);
      scanf("%d", &m1);

      if(m1 >= m0){ //both octaves are correct       
        //get tile coordinates
        int nCol = 0, nRow = 0;
        printf("Tile coordinates:\n");  
        printf("  Row: "); scanf("%d", &nRow); 
        printf("  Col: "); scanf("%d", &nCol);

        //generate the noise texture and save it
        GenerateAndSave2DNoise(nRow, nCol, m0, m1, n, seed);
      } //if
      else printf("That's not between %d and %d. Bailing out.\n", m0, log2n);
    } //if
    else printf("That's not between 1 and %d. Bailing out.\n", log2n);
  } //if
  else printf("That's not a power of 2 and at least 2. Bailing out.\n");
  
#if defined(_MSC_VER) //Windows Visual Studio 
  //wait for user keystroke and exit
  printf("\nHit Almost Any Key to Exit...\n");
  _getch();
#endif

  return 0;
} //main
