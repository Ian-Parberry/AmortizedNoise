/// \file FiniteAmortizedNoise3D.cpp 
/// \brief Code file for the 3D amortized noise class CFiniteAmortizedNoise3D.

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

#include "FiniteAmortizedNoise3D.h"
#include "Common.h"

#define sqr(x) ((x)*(x)) ///< Square function.

CFiniteAmortizedNoise3D::CFiniteAmortizedNoise3D(int n){ 
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
  
  //initialize the Perlin noise tables.
  initPerlinNoiseTables();
} //constructor

CFiniteAmortizedNoise3D::~CFiniteAmortizedNoise3D(){
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

void CFiniteAmortizedNoise3D::FillUp(float* t, float s, int n){
  const float d = s/n;
  t[0] = 0.0f; t[1] = d;
  for(int i=2; i<=n; i++)
    t[i] = t[i-1] + d;
} //FillUp

/// Fill amortized noise table top down.
/// \param t Amortized noise table.
/// \param s Initial value.
/// \param n Granularity.

void CFiniteAmortizedNoise3D::FillDn(float* t, float s, int n){
  const float d = -s/n;
  t[n] = 0.0f; t[n-1] = d;
  for(int i=n-2; i>=0; i--)
    t[i] = t[i+1] + d;
} //FillDn

/// Get a random unit floating point number.
/// \return Random floating point number >=-1 and <=1.

float CFiniteAmortizedNoise3D::getRandomUnitFloat(){
  return (float)((rand()%(B + B)) - B)/B;
} //getRandomUnitFloat


/// A 1D hash function.
/// Hash one dimension into a single unsigned int.
/// \param x  value to be hashed.
/// \return Hash of x.

unsigned int CFiniteAmortizedNoise3D::h(const unsigned int x){
  return p[x & BM];
} //h

/// A 3D hash function.
/// Hash three dimensions into a single unsigned int.
/// \param x X coordinate of value to be hashed.
/// \param y Y coordinate of value to be hashed.
/// \param z Z coordinate of value to be hashed.
/// \return Hash of (x, y, z).

unsigned int CFiniteAmortizedNoise3D::h(const unsigned int x, const unsigned int y, const unsigned int z){
  return h(h(h(x) + y) + z);
} //h

/// Initialize the Perlin noise tables.
/// Loads p with a random permutation and g3 with a randomly
/// chosen selection of unit gradients.

void CFiniteAmortizedNoise3D::initPerlinNoiseTables(){
  //random normalized gradient vectors
  for(int i=0; i<B; i++){
    g3[i][0] = getRandomUnitFloat();
    g3[i][1] = getRandomUnitFloat();
    g3[i][2] = getRandomUnitFloat();
    float m = sqrt(sqr(g3[i][0]) + sqr(g3[i][1]) + sqr(g3[i][2]));
    g3[i][0] /= m; g3[i][1] /= m; g3[i][2] /= m;
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

/// Initialize amortized noise tables.
/// \param x0 x coordinate of top left front corner of cell.
/// \param y0 y coordinate of top left front corner of cell.
/// \param z0 z coordinate of top left front corner of cell.
/// \param n Granularity.

void CFiniteAmortizedNoise3D::initEdgeTables(const int x0, const int y0, const int z0, const int n){ 
  //compute gradients at corner points
  unsigned int b000 = h(x0  , y0  , z0  );
  unsigned int b010 = h(x0  , y0+1, z0  );
  unsigned int b100 = h(x0+1, y0  , z0  );
  unsigned int b110 = h(x0+1, y0+1, z0  );
  unsigned int b001 = h(x0  , y0  , z0+1);
  unsigned int b011 = h(x0  , y0+1, z0+1);
  unsigned int b101 = h(x0+1, y0  , z0+1);
  unsigned int b111 = h(x0+1, y0+1, z0+1);
  
  //fill inferred gradient tables from corner gradients
  FillUp(uax, g3[b000][0], n); FillDn(vax, g3[b001][0], n);
  FillUp(ubx, g3[b010][0], n); FillDn(vbx, g3[b011][0], n);
  FillUp(uay, g3[b000][1], n); FillUp(vay, g3[b001][1], n);
  FillDn(uby, g3[b010][1], n); FillDn(vby, g3[b011][1], n);
  FillUp(uaz, g3[b000][2], n); FillUp(vaz, g3[b001][2], n);
  FillUp(ubz, g3[b010][2], n); FillUp(vbz, g3[b011][2], n);

  FillUp(ucx, g3[b100][0], n); FillDn(vcx, g3[b101][0], n);
  FillUp(udx, g3[b110][0], n); FillDn(vdx, g3[b111][0], n);
  FillUp(ucy, g3[b100][1], n); FillUp(vcy, g3[b101][1], n);
  FillDn(udy, g3[b110][1], n); FillDn(vdy, g3[b111][1], n);
  FillDn(ucz, g3[b100][2], n); FillDn(vcz, g3[b101][2], n);
  FillDn(udz, g3[b110][2], n); FillDn(vdz, g3[b111][2], n);
} //initEdgeTables

/// Initialize the spline table.
/// \param n Granularity.

void CFiniteAmortizedNoise3D::initSplineTable(const int n){
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

float CFiniteAmortizedNoise3D::getNoise(const int i, const int j, const int k){  
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

void CFiniteAmortizedNoise3D::getNoise(const int n, const int i0, const int j0, const int k0, float*** cell){  
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

void CFiniteAmortizedNoise3D::addNoise(const int n, const int i0, const int j0, const int k0, float scale, float*** cell){  
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

float CFiniteAmortizedNoise3D::generate(int x, int y, int z, const int m0, const int m1, int n, float*** cell){
  int r = 1; // Side length of cell divided by side length of subcell.

  //skip over unwanted octaves
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
