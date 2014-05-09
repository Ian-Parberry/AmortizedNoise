/// \file Main.cpp
/// \brief Main.
///
/// \mainpage 3D Amortized Noise Generator
///
/// This project demonstrates how the 3D amortized noise algorithm can be used to produce
/// an infinite greyscale tiled 3D texture. It will prompt the user via stdin/stdout
/// for the texture size (which must be a power of 2), a seed value, the largest and smallest octaves,
/// and tile coordinates in row-column-sheet format. The input values are checked for errors. 
/// If they pass, the tile texture will be saved in a png file (windows) or tga file (other OS).
/// The amount of user CPU time used in generating the noise is reported to stdout.
/// Saved files will be named  i[F]s[N]o[LS]s[D]r[R]c[C]-NNN.X (for example, i1s256o25s5432r54c33.png), where:
/// <center><table>
/// <tr><td>F<td>0 for finite or 1 for infinite
/// <tr><td>N<td>tile side
/// <tr><td>L<td>largest octave
/// <tr><td>S<td>smallest octave
/// <tr><td>D<td>seed
/// <tr><td>R<td>row
/// <tr><td>C<td>column
/// <tr><td>H<td>sheet
/// <tr><td>X<td>png or tga
/// </table></center>
///
/// For more details on amortized noise, see Ian Parberry, "Amortized Noise",
/// <em>Journal of Computer Graphics Techniques</em>, Vol. ?, No. ?, pp. ?, To Appear.

// Copyright Ian Parberry, September 2013.
//
// This file is made available under the GNU All-Permissive License.
//
// Copying and distribution of this file, with or without modification,
// are permitted in any medium without royalty provided the copyright
// notice and this notice are preserved.  This file is offered as-is,
// without any warranty.
//
// Created by Ian Parberry, September 2013.
// Last updated May 9, 2014.

#include <stdlib.h> //for rand()
#include <stdio.h> //for printf()

#include "defines.h"  //OS porting defines

#include "FiniteAmortizedNoise3D.h"
#include "InfiniteAmortizedNoise3D.h"
#include "CPUtime.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION ///< For platform-neutral image file output.
#include "stb_image_write.h"

CFiniteAmortizedNoise3D* g_pFiniteAmortizedNoise = NULL; ///< Pointer to the finite amortized 3D noise generator. 
CInfiniteAmortizedNoise3D* g_pInfiniteAmortizedNoise = NULL; ///< Pointer to the infinite amortized 3D noise generator. 

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

/// Generate a random cell of 3D noise.
/// \param cell Pointer to array in which to store the noise.
/// \param x Tile column index.
/// \param y Tile row index.
/// \param z Tile sheet index.
/// \param m0 Largest octave.
/// \param m1 Smallest octave.
/// \param n Cell size.

float Generate3DNoise(float*** cell, const int x, const int y, const int z, const int m0, const int m1, const int n){ 
  printf("Generating %d octaves of 3D noise with persistence 0.5 and lacunarity 2.0\n", m1 - m0 + 1);
  int t = CPUTimeInMilliseconds();
  float scale = 1.0f;
  if(g_bInfinite)
    scale = g_pInfiniteAmortizedNoise->generate(z, y, x, m0, m1, n, cell);
  else scale = g_pFiniteAmortizedNoise->generate(z, y, x, m0, m1, n, cell);

  t = CPUTimeInMilliseconds() - t;
  printf("Generated %d points in %0.2f seconds CPU time.\n", n*n, t/1000.0f);
  return scale;
} //Generate3DNoise

/// Save a cell of 3D noise as a set of png files.
/// \param cell Pointer to array in which to store the noise.
/// \param scale Rescale by this to get in the range -1..1.
/// \param n Cell size.
/// \param frames Number of files to save.
/// \param filename Root of the file name under which to store the noise.

void Save3DNoise(float*** cell, const int n, const int frames, float scale, char* filename){ 
#if defined(_MSC_VER) //Windows Visual Studio only 
  printf("Saving %d frames to %dx%d png files %s-NNN.png\n\n", frames, n, n, filename);
#else
  printf("Saving %d frames to %dx%d tga files %s-NNN.tga\n\n", frames, n, n, filename);
#endif

  int delta = n/frames;

  char buffer[MAX_PATH];
  for(int i=0; i<n; i+=delta){
    sprintf(buffer, "%s-%c%c%c", filename, '0' + (i%1000)/100, '0' + (i%100)/10, '0' + i%10);
    Save2DNoise(cell[i], n, 1.0f, buffer);
  } //for
} //Save3DNoise

/// Generate a cell of 3D noise and save it as an image file.
/// \param nRow Tile row index.
/// \param nCol Tile column index.
/// \param nSheet Tile sheet index.
/// \param m0 Largest octave.
/// \param m1 Smallest octave.
/// \param n Cell size.
/// \param seed Hash function seed.
/// \param frames Number of 2D slices across the 3D texture to actually save.

void GenerateAndSave3DNoise(int nRow, int nCol, int nSheet, const int m0, const int m1,
                            const int n, const int seed, const int frames){ 
  g_pFiniteAmortizedNoise = new CFiniteAmortizedNoise3D(n);  
  g_pInfiniteAmortizedNoise = new CInfiniteAmortizedNoise3D(n, seed);

  //construct file name
  char filename[MAX_PATH];
  sprintf(filename, "i%1ds%do%1d%1ds%dr%dc%dh%d", g_bInfinite, n, m0, m1, seed, nRow, nCol, nSheet);

  //adjust for tile size of smallest octave
  for(int i=1; i<m0; i++){ 
    nCol *= 2; nRow *= 2; nSheet *= 2;
  } //for

  //allocate space to store noise cells
  float*** cell = new float** [n];
  for(int i=0; i<n; i++){
    cell[i] = new float* [n];
    for(int j=0; j<n; j++)
      cell[i][j] = new float [n];
  } //for 

  //generate and save noise
  float scale = Generate3DNoise(cell, nCol, nRow, nSheet, m0, m1, n);
  Save3DNoise(cell, n, frames, scale, filename);

  //deallocate noise cell space
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++)
      delete [] cell[i][j];
    delete [] cell[i];
  } //for
  delete [] cell;

  delete g_pFiniteAmortizedNoise;
  delete g_pInfiniteAmortizedNoise;
} //GenerateAndSave3DNoise

/// Main.
/// \param argc Argument count
/// \param argv Arguments.
/// \return 0 for success, 1 for failure

int main(int argc, char *argv[]){ 
  printf("Amortized 3D Noise Generator, Ian Parberry, 2014\n");
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
        int nCol = 0, nRow = 0, nSheet = 0;
        printf("Tile coordinates:\n");  
        printf("  Row: "); scanf("%d", &nRow); 
        printf("  Col: "); scanf("%d", &nCol);
        printf("  Sheet: "); scanf("%d", &nSheet);

         int frames = 1;
         printf("Number of frames to save (must be a power of 2):\n> ");
         scanf("%d", &frames);

         if((frames >= 2) && !(frames & (frames - 1))) //n is a power of 2 and at least 2     
           GenerateAndSave3DNoise(nRow, nCol, nSheet, m0, m1, n, seed, frames); //generate the noise texture and save it
         else      
           printf("That's not a power of 2. Bailing out.\n");
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
