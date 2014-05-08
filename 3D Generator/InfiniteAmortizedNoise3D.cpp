/// \file InfiniteAmortizedNoise3D.cpp 
/// \brief Code file for the 3D amortized noise class CInfiniteAmortizedNoise3D.

// Copyright Ian Parberry, January 2014.
//
// This file is made available under the GNU All-Permissive License.
//
// Copying and distribution of this file, with or without modification,
// are permitted in any medium without royalty provided the copyright
// notice and this notice are preserved.  This file is offered as-is,
// without any warranty.
//
// Created by Ian Parberry, January 2014.
// Last updated May 7, 2014.

#include <math.h> //for trig functions
#include <stdlib.h> //for rand()

#include "InfiniteAmortizedNoise3D.h"
#include "Common.h"
#include "MurmurHash3.h"

CInfiniteAmortizedNoise3D::CInfiniteAmortizedNoise3D(int n, int s): seed(s){ 
  //Allocate space for amortized noise tables.
  uax = new float [n+1]; vax = new float [n+1]; 
  ubx = new float [n+1]; vbx = new float [n+1]; 
  uay = new float [n+1]; vay = new float [n+1]; 
  uby = new float [n+1]; vby = new float [n+1];

  ucx = new float [n+1]; vcx = new float [n+1]; 
  udx = new float [n+1]; vdx = new float [n+1]; 
  ucy = new float [n+1]; vcy = new float [n+1]; 
  udy = new float [n+1]; vdy = new float [n+1];

  uaz = new float [n+1]; vaz = new float [n+1]; 
  ubz = new float [n+1]; vbz = new float [n+1]; 
  ucz = new float [n+1]; vcz = new float [n+1]; 
  udz = new float [n+1]; vdz = new float [n+1];

  //Allocate space for spline table.
  spline = new float [n+1]; 
} //constructor

CInfiniteAmortizedNoise3D::~CInfiniteAmortizedNoise3D(){
  //Deallocate space for amortized noise tables.
  delete [] uax; delete [] ubx; 
  delete [] vax; delete [] vbx; 
  delete [] uay; delete [] uby; 
  delete [] vay; delete [] vby;

  delete [] ucx; delete [] udx; 
  delete [] vcx; delete [] vdx; 
  delete [] ucy; delete [] udy; 
  delete [] vcy; delete [] vdy; 
  
  delete [] uaz; delete [] vaz; 
  delete [] ubz; delete [] vbz; 
  delete [] ucz; delete [] vcz; 
  delete [] udz; delete [] vdz; 

  //Deallocate space for spline table.
  delete [] spline; 
} //destructor

///< Fill amortized noise table bottom up.
///< \param t Amortized noise table.
///< \param s Initial value.
///< \param n Granularity.

void CInfiniteAmortizedNoise3D::FillUp(float* t, float s, int n){
  const float d = s/n;
  t[0] = 0.0f; t[1] = d;
  for(int i=2; i<=n; i++)
    t[i] = t[i-1] + d;
} //FillUp

/// Fill amortized noise table top down.
/// \param t Amortized noise table.
/// \param s Initial value.
/// \param n Granularity.

void CInfiniteAmortizedNoise3D::FillDn(float* t, float s, int n){
  const float d = -s/n;
  t[n] = 0.0f; t[n-1] = d;
  for(int i=n-2; i>=0; i--)
    t[i] = t[i+1] + d;
} //FillDn

/// Get a random unit floating point number.
/// \return Random floating point number >=-1 and <=1.

float CInfiniteAmortizedNoise3D::getRandomUnitFloat(){
  return (float)((rand()%(B + B)) - B)/B;
} //getRandomUnitFloat

/// Data structure to help construct a hash key for MurmurHash.

union KeyHelper{
  unsigned int integer[3]; ///< The key as an array of 3 ints.
  unsigned char byte[12]; ///< The key as an array of 12 bytes.
};

/// A 3D hash function.
/// Hash three dimensions into a single unsigned int.
/// \param x X coordinate of value to be hashed.
/// \param y Y coordinate of value to be hashed.
/// \param z Z coordinate of value to be hashed.
/// \return Hash of (x, y, z).

float CInfiniteAmortizedNoise3D::h1(const unsigned int x, const unsigned int y, const unsigned int z){ 
  KeyHelper key;
  key.integer[0] = x;
  key.integer[1] = y;
  key.integer[2] = z;
  unsigned int result;
  MurmurHash3_32(&key.byte, 12, seed, &result); //do the heavy lifting
  return (float)result;
} //h1

float CInfiniteAmortizedNoise3D::h2(const unsigned int x, const unsigned int y, const unsigned int z){ 
  KeyHelper key;
  key.integer[0] = x;
  key.integer[1] = y;
  key.integer[2] = z;
  unsigned int result;
  MurmurHash3_32(&key.byte, 12, seed*13, &result); //do the heavy lifting
  return (float)result;
} //h2

/// Initialize amortized noise tables.
/// \param x0 x coordinate of top left front corner of cell.
/// \param y0 y coordinate of top left front corner of cell.
/// \param z0 z coordinate of top left front corner of cell.
/// \param n Granularity.

void CInfiniteAmortizedNoise3D::initEdgeTables(const int x0, const int y0, const int z0, const int n){ 
  //compute gradients at corner points
  float b0  = h1(x0,   y0, z0),   b1  = h1(x0,   y0+1, z0  ),
        b2  = h1(x0+1, y0, z0),   b3  = h1(x0+1, y0+1, z0  ),
        b4  = h1(x0,   y0, z0+1), b5  = h1(x0,   y0+1, z0+1),
        b6  = h1(x0+1, y0, z0+1), b7  = h1(x0+1, y0+1, z0+1),

        b8  = h2(x0,   y0, z0),   b9  = h2(x0,   y0+1, z0  ),
        b10 = h2(x0+1, y0, z0),   b11 = h2(x0+1, y0+1, z0  ),
        b12 = h2(x0,   y0, z0+1), b13 = h2(x0,   y0+1, z0+1),
        b14 = h2(x0+1, y0, z0+1), b15 = h2(x0+1, y0+1, z0+1);
  
  //fill inferred gradient tables from corner gradients

  FillUp(uax,  sinf(b0)*sinf(b8),  n); FillDn(vax,  sinf(b4)*sinf(b12), n);
  FillUp(ubx,  sinf(b1)*sinf(b9),  n); FillDn(vbx,  sinf(b5)*sinf(b13), n);
  FillUp(uay, -sinf(b0)*cosf(b8),  n); FillUp(vay, -sinf(b4)*cosf(b12), n);
  FillDn(uby, -sinf(b1)*cosf(b9),  n); FillDn(vby, -sinf(b5)*cosf(b13), n);

  FillUp(uaz,      cosf(b0),      n); FillUp(vaz,      cosf(b4),      n);
  FillUp(ubz,      cosf(b1),      n); FillUp(vbz,      cosf(b5),      n);
  FillUp(ucx,  sinf(b2)*sinf(b10), n); FillDn(vcx,  sinf(b6)*sinf(b14), n);
  FillUp(udx,  sinf(b3)*sinf(b11), n); FillDn(vdx,  sinf(b7)*sinf(b15), n);

  FillUp(ucy, -sinf(b2)*cosf(b10), n); FillUp(vcy, -sinf(b6)*cosf(b14), n);
  FillDn(udy, -sinf(b3)*cosf(b11), n); FillDn(vdy, -sinf(b7)*cosf(b15), n);
  FillDn(ucz,      cosf(b2),      n); FillDn(vcz,      cosf(b6),      n);
  FillDn(udz,      cosf(b3),      n); FillDn(vdz,      cosf(b7),      n); 
} //initEdgeTables

/// Initialize the spline table.
/// \param n Granularity.

void CInfiniteAmortizedNoise3D::initSplineTable(const int n){
  for(int i=0; i<=n; i++){ //for each table entry
    float t = (float)i/n; //offset between grid points
    spline[i] = s_curve(t);  //cubic spline
  } //for
} //initSplineTable

/// Compute a single point of a single octave of Perlin noise. This is similar
/// to Perlin's noise3 function except that it substitutes table lookups for
/// floating point multiplication.
/// \param i x coordinate of point.
/// \param j y coordinate of point.
/// \param k z coordinate of point.
/// \return Noise value at (i, j, k).

float CInfiniteAmortizedNoise3D::getNoise(const int i, const int j, const int k){  
  float u = uax[k] + uay[j] + uaz[i];
  float v = vax[k] + vay[j] + vaz[i];
  float a = lerp(spline[k], u, v); 
  u = ubx[k] + uby[j] + ubz[i];
  v = vbx[k] + vby[j] + vbz[i];
  float b = lerp(spline[k], u, v);
  float e = lerp(spline[j], a, b);  
      
  u = ucx[k] + ucy[j] + ucz[i];
  v = vcx[k] + vcy[j] + vcz[i];
  float c = lerp(spline[k], u, v); 
  u = udx[k] + udy[j] + udz[i];
  v = vdx[k] + vdy[j] + vdz[i];
  float d = lerp(spline[k], u, v);
  float f = lerp(spline[j], c, d);   

  return lerp(spline[i], e, f); 
} //getNoise

/// Get a single octave of noise into a subcell.
/// \param n Granularity.
/// \param i0 X offset of this subcell in cell.
/// \param j0 Y offset of this subcell in cell.
/// \param k0 Z offset of this subcell in cell.
/// \param cell Noise cell.

void CInfiniteAmortizedNoise3D::getNoise(const int n, const int i0, const int j0, const int k0, float*** cell){  
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      for(int k=0; k<n; k++) 
        cell[i0 + i][j0 + j][k0 + k] = getNoise(i, j, k); //this is the only line that differs from addNoise
} //getNoise

/// Add a single octave of noise into a subcell.
/// \param n Granularity.
/// \param i0 X offset of this subcell in cell.
/// \param j0 Y offset of this subcell in cell.
/// \param k0 Z offset of this subcell in cell.
/// \param scale Noise is to be rescaled by this factor.
/// \param cell Noise cell.

void CInfiniteAmortizedNoise3D::addNoise(const int n, const int i0, const int j0, const int k0, float scale, float*** cell){  
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      for(int k=0; k<n; k++) 
        cell[i0 + i][j0 + j][k0 + k] += scale * getNoise(i, j, k); //this is the only line that differs from getNoise
} //addNoise

///< Get a cell of 1/f amortized noise. Assumes step 0.5 and lacunarity 2.0.
///< \param x x coordinate of top left corner of cell.
///< \param y y coordinate of top left corner of cell.
///< \param z z coordinate of top left corner of cell.
///< \param m0 First octave.
///< \param m1 Last octave.
///< \param n Granularity.
///< \param cell Generated noise.

float CInfiniteAmortizedNoise3D::generate(int x, int y, int z, const int m0, const int m1, int n, float*** cell){
  int r = 1; // Side length of cell divided by side length of subcell.

  //skip over unwanted octaves
  for(int i=1; i<m0; i++){
    n /= 2; r *= 2;
  } //for
  
  if(n < 2)return 1.0f; //fail and bail - should not happen
    
  //Generate first octave directly into cell.
  //  We could add all octaves into the cell directly if we zero out the cell
  //  before we begin. However, that is a nontrivial overhead that Perlin noise
  //  does not have, and we can avoid it too by putting in the first octave and
  //  adding in the rest.

  initSplineTable(n); //initialize the spline table to cells of size n
  for(int i0=0; i0<r; i0++)
    for(int j0=0; j0<r; j0++)
      for(int k0=0; k0<r; k0++){ //for each subcell
        initEdgeTables(x + i0, y + j0, z + k0, n); //initialize the amortized noise tables
        getNoise(n, i0*n, j0*n, k0*n, cell); //generate noise directly into cell
      } //for

  float scale = 1.0f; //scale factor
  
  //Generate the other octaves and add them into cell. See previous comment.
  for(int k=m0; k<m1 && n>=2; k++){ //for each octave after the first
    n /= 2; r += r;  x += x; y += y; z += z; scale *= 0.5f; //rescale for next octave
    initSplineTable(n); //initialize the spline table to cells of size n
    for(int i0=0; i0<r; i0++)
      for(int j0=0; j0<r; j0++)
        for(int k0=0; k0<r; k0++){ //for each subcell
          initEdgeTables(x + i0, y + j0, z + k0, n); //initialize the edge tables
          addNoise(n, i0*n, j0*n, k0*n, scale, cell); //generate directly into cell
        } //for
  } //for each octave

  //Compute 1/magnitude and return it. 
  //  A single octave of Perlin noise returns a value of magnitude at most 
  //  1/sqrt(3). Adding magnitudes over all scaled octaves gives a total
  //  magnitude of (1 + 0.5 + 0.25 +...+ scale)/sqrt(3). This is equal to
  //  (2 - scale)/sqrt(3) (using the standard formula for the sum of a geometric
  //  progression). 1/magnitude is therefore sqrt(3)/(2-scale).

  return 1.732f/(2.0f - scale); //scale factor
} //generate
