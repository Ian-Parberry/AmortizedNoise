/// \file FiniteAmortizedNoise2D.cpp
/// \brief Code file for the 2D amortized noise class CFiniteAmortizedNoise2D.

// Copyright Ian Parberry, April 2014.
//
// This file is made available under the GNU All-Permissive License.
//
// Copying and distribution of this file, with or without modification,
// are permitted in any medium without royalty provided the copyright
// notice and this notice are preserved.  This file is offered as-is,
// without any warranty.
//
// Created by Ian Parberry, April 2014.
// Last updated May 7, 2014.

#define _USE_MATH_DEFINES ///< Enable use of constant M_SQRT2 in math.h
#include <math.h> //for trig functions
#include <stdlib.h> //for rand()

#include "FiniteAmortizedNoise2D.h"

#define sqr(x) ((x)*(x)) ///< Square function.

/// Constructor.
/// \param n Cell size.

CFiniteAmortizedNoise2D::CFiniteAmortizedNoise2D(const unsigned int n){ 
  //Allocate space for amortized noise tables.
  uax = new float [n]; vax = new float [n]; //ax
  ubx = new float [n]; vbx = new float [n]; //bx
  uay = new float [n]; uby = new float [n]; //ay
  vay = new float [n]; vby = new float [n]; //by

  //Allocate space for spline table.
  spline = new float [n]; 

  //initialize the Perlin noise tables.
  initPerlinNoiseTables();
} //constructor

CFiniteAmortizedNoise2D::~CFiniteAmortizedNoise2D(){
  //Deallocate space for amortized noise tables.
  delete [] uax; delete [] vax; //ax
  delete [] ubx; delete [] vbx; //bx
  delete [] uay; delete [] uby; //ay
  delete [] vay; delete [] vby; //by
  
  //Deallocate space for spline table.
  delete [] spline; 
} //destructor

/// Fill amortized noise table bottom up.
/// \param t Amortized noise table.
/// \param s Initial value.
/// \param n Granularity.

void CFiniteAmortizedNoise2D::FillUp(float* t, const float s, const int n){
  const float d = s/n; //increment amount
  t[0] = 0.0f; //corner
  for(int i=1; i<n; i++) //edge values
    t[i] = t[i-1] + d;
} //FillUp

/// Fill amortized noise table top down.
/// \param t Amortized noise table.
/// \param s Initial value.
/// \param n Granularity.

void CFiniteAmortizedNoise2D::FillDn(float* t, const float s, const int n){
  const float d = -s/n;  //increment amount
  t[n-1] = d; //corner
  for(int i=n-2; i>=0; i--) //edge values
    t[i] = t[i+1] + d;
} //FillDn

/// Get a random unit floating point number.
/// \return Random floating point number >=-1 and <=1.

float CFiniteAmortizedNoise2D::getRandomUnitFloat(){
  return (float)((rand()%(B + B)) - B)/B;
} //getRandomUnitFloat

/// A 1D hash function.
/// Hash one dimension into a single unsigned int.
/// \param x  value to be hashed.
/// \return Hash of x.

unsigned int CFiniteAmortizedNoise2D::h(const unsigned int x){
  return p[x & BM];
} //h

/// A 2D hash function.
/// Hash two dimensions into a single unsigned int.
/// \param x X coordinate of value to be hashed.
/// \param y Y coordinate of value to be hashed.
/// \return Hash of (x, y).

unsigned int CFiniteAmortizedNoise2D::h(const unsigned int x, const unsigned int y){
  return h(h(x) + y); 
} //h

/// Initialize the Perlin noise tables.
/// Loads p with a random permutation and g2 with a randomly
/// chosen selection of unit gradients.

void CFiniteAmortizedNoise2D::initPerlinNoiseTables(){
  //random normalized gradient vectors
  for(int i=0; i<B; i++){
    g2[i][0] = getRandomUnitFloat();
    g2[i][1] = getRandomUnitFloat();
    float m = sqrt(sqr(g2[i][0]) + sqr(g2[i][1]));
    g2[i][0] /= m; g2[i][1] /= m;
  } //for
 
  //identity permutation
  for(int i=0; i<B; i++)
    p[i] = i; 

  //randomize permutation 
  for(int i=B-1; i>0; i--){
    int tmp = p[i];
    int j = rand()%(i + 1);
    p[i] = p[j];
    p[j] = tmp;
  } //for
} //initPerlinNoiseTables

/// Initialize the amortized noise tables.
/// \param x0 x coordinate of top left corner of cell.
/// \param y0 y coordinate of top left corner of cell.
/// \param n Granularity.

void CFiniteAmortizedNoise2D::initEdgeTables(const int x0, const int y0, const int n){
  //compute gradients at corner points
  unsigned int b00 = h(x0,   y0);
  unsigned int b01 = h(x0,   y0+1); 
  unsigned int b10 = h(x0+1, y0);
  unsigned int b11 = h(x0+1, y0+1);
  
  //fill inferred gradient tables from corner gradients
  FillUp(uax, g2[b00][0], n); FillDn(vax, g2[b01][0], n);
  FillUp(ubx, g2[b10][0], n); FillDn(vbx, g2[b11][0], n);
  FillUp(uay, g2[b00][1], n); FillUp(vay, g2[b01][1], n);
  FillDn(uby, g2[b10][1], n); FillDn(vby, g2[b11][1], n);
} //initEdgeTables

/// Initialize the spline table.
/// \param n Granularity.

void CFiniteAmortizedNoise2D::initSplineTable(const int n){
  for(int i=0; i<n; i++){ //for each table entry
    float t = (float)i/n; //offset between grid points
    spline[i] = s_curve(t); //cubic spline
  } //for
} //initSplineTable

/// Compute a single point of a single octave of Perlin noise. This is similar
/// to Perlin's noise2 function except that it substitutes table lookups for
/// floating point multiplication.
/// \param i x coordinate of point.
/// \param j y coordinate of point.
/// \return Noise value at (i, j).

float CFiniteAmortizedNoise2D::getNoise(const int i, const int j){  
  float u = uax[j] + uay[i];
  float v = vax[j] + vay[i];
  const float a = lerp(spline[j], u, v); 
  u = ubx[j] + uby[i];
  v = vbx[j] + vby[i];
  const float b = lerp(spline[j], u, v);   
  return lerp(spline[i], a, b);   
} //getNoise

/// Get a single octave of noise into a subcell.
/// This differs from CInfiniteAmortizedNoise2D::addNoise in that it copies the noise
/// to the cell instead of adding it in.
/// \param n Granularity.
/// \param i0 x offset of this subcell in cell.
/// \param j0 y offset of this subcell in cell.
/// \param cell Noise cell.

void CFiniteAmortizedNoise2D::getNoise(const int n, const int i0, const int j0, float** cell){  
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++) 
      cell[i0 + i][j0 + j] = getNoise(i, j); //the only line that differs from addNoise
} //getNoise

/// Add a single octave of noise into a subcell.
/// This differs from CInfiniteAmortizedNoise2D::getNoise in that it adds the noise
/// to the cell instead of copying it there.
/// \param n Granularity.
/// \param i0 x offset of this subcell in cell.
/// \param j0 y offset of this subcell in cell.
/// \param scale Noise is to be rescaled by this factor.
/// \param cell Noise cell.

void CFiniteAmortizedNoise2D::addNoise(const int n, const int i0, const int j0, const float scale, float** cell){  
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      cell[i0 + i][j0 + j] += 
        scale * getNoise(i, j); //the only line that differs from getNoise
} //addNoise

/// Generate a cell of 1/f amortized noise with persistence 0.5 and lacunarity 2.0.
/// \param x x coordinate of top left corner of cell.
/// \param y y coordinate of top left corner of cell.
/// \param m0 First octave.
/// \param m1 Last octave.
/// \param n Granularity.
/// \param cell Cell to put generated noise into.
/// \return Multiply noise by this to get into the range -1..1

float CFiniteAmortizedNoise2D::generate(int x, int y, const int m0, const int m1, int n, float** cell){
  int r = 1; //Side length of cell divided by side length of subcell.

  //Skip over unwanted octaves.
  for(int i=1; i<m0; i++){
    n /= 2; r += r;
  } //for

  if(n < 2)return 1.0f; //fail and bail - should not happen

  //Generate first octave directly into cell.
  //  We could add all octaves into the cell directly if we zero out the cell
  //  before we begin. However, that is a nontrivial overhead that Perlin noise
  //  does not have, and we can avoid it too by putting in the first octave and
  //  adding in the rest.

  initSplineTable(n); //initialize the spline table to cells of size n
  for(int i0=0; i0<r; i0++)
    for(int j0=0; j0<r; j0++){ //for each subcell
      initEdgeTables(x + i0, y + j0, n); //initialize the amortized noise tables
      getNoise(n, i0*n, j0*n, cell); //generate noise directly into cell
    } //for

  float scale = 1.0f; //scale factor

  //Generate the other octaves and add them into cell. See previous comment.
  for(int k=m0; k<m1 && n>=2; k++){ //for each octave after the first
    n /= 2; r += r;  x += x; y += y; scale *= 0.5f; //rescale for next octave
    initSplineTable(n); //initialize the spline table to cells of size n
    for(int i0=0; i0<r; i0++)
      for(int j0=0; j0<r; j0++){ //for each subcell
        initEdgeTables(x + i0, y + j0, n); //initialize the edge tables
        addNoise(n, i0*n, j0*n, scale, cell); //generate directly into cell
      } //for
  } //for each octave

  //Compute 1/magnitude and return it. 
  //  A single octave of Perlin noise returns a value of magnitude at most 
  //  1/sqrt(2). Adding magnitudes over all scaled octaves gives a total
  //  magnitude of (1 + 0.5 + 0.25 +...+ scale)/sqrt(2). This is equal to
  //  (2 - scale)/sqrt(2) (using the standard formula for the sum of a geometric
  //  progression). 1/magnitude is therefore sqrt(2)/(2-scale).

  return (float)M_SQRT2/(2.0f - scale); //multiply by this to bring noise to [-1,1]
} //generate
