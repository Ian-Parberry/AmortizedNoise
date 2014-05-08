/// \file perlin.cpp
/// \brief Code file for 2D Perlin noise functions.

// Copyright Ian Parberry, January 2014.
//
// This file is made available under the GNU All-Permissive License.
//
// Copying and distribution of this file, with or without modification,
// are permitted in any medium without royalty provided the copyright
// notice and this notice are preserved.  This file is offered as-is,
// without any warranty.

// Created by Ian Parberry, January 2014.
// Last updated May 7, 2014.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Common.h"

static int p[B + B + 2]; ///< Ken Perlin's permutation table
static float g2[B + B + 2][2]; ///< Ken Perlin's 2D gradient table.

/// Compute a single octave at 2D noise at a single point.
/// \param vec Point at which to evaluate noise.
/// \return Noise value between -1.0 and 1.0.

float noise2(float vec[2]){
  int bx0, bx1, by0, by1, b00, b10, b01, b11;
  float rx0, rx1, ry0, ry1, *q, sx, sy, a, b, t, u, v;
  int i, j;

  setup(0, bx0,bx1, rx0,rx1);
  setup(1, by0,by1, ry0,ry1);

  i = p[bx0];
  j = p[bx1];

  b00 = p[i + by0];
  b10 = p[j + by0];
  b01 = p[i + by1];
  b11 = p[j + by1];

  sx = s_curve(rx0);
  sy = s_curve(ry0);

  q = g2[b00] ; u = at2(rx0, ry0);
  q = g2[b10] ; v = at2(rx1, ry0);
  a = lerp(sx, u, v);

  q = g2[b01] ; u = at2(rx0, ry1);
  q = g2[b11] ; v = at2(rx1, ry1);
  b = lerp(sx, u, v);

  return lerp(sy, a, b);
} //noise2

/// 2D vector normalize. Works by side-effect.
/// \param v 2D vector.

void normalize2(float v[2]){
  float s = sqrt(v[0]*v[0] + v[1]*v[1]);
  v[0] /= s; v[1] /= s;
} //normalize2

/// Swap two values.
/// \param x First value.
/// \param y Second value.

void swap(int& x, int& y){
  int k = x; x = y; y = k;
} //swap

/// Initialize Perlin's permutation and gradient tables.

void initPerlin2D(){
  for(int i=0; i<B; i++)
    p[i] = i;
  for(int i=B-1; i>0; i--)
   swap(p[i], p[rand()%(i+1)]);
    
  for(int i=0; i<B; i++){
    g2[i][0] = (float)((rand()%(B + B)) - B)/B;
    g2[i][1] = (float)((rand()%(B + B)) - B)/B;
    normalize2(g2[i]);
  } //for

  for(int i=0; i<B+2; i++){
    p[B + i] = p[i];
    for(int j=0; j<2; j++)
      g2[B + i][j] = g2[i][j];
  } //for
} //initPerlin2D

/// Compute turbulence, also known as 1/f noise.
/// \param x X coordinate.
/// \param y Y coordinate.
/// \param alpha Persistence.
/// \param beta Lacunarity.
/// \param n Side of square grid.
/// \return Scale factor required to bring it down to range [-1, 1].

float PerlinNoise2D(float x, float y, float alpha, float beta, int n){
  float sum=0.0f, p[2], scale=1.0f;
  p[0] = x; p[1] = y;

  for(int i=0; i<n; i++){
    sum += noise2(p)*scale;
    scale *= alpha;
    p[0] *= beta;	p[1] *= beta;
  } //for

  return sum/(2.0f - scale);
} //PerlinNoise2D
