/// \file FiniteAmortizedNoise2D.h
/// \brief Header file for the 2D amortized noise class CFiniteAmortizedNoise2D.
///
/// Copyright Ian Parberry, April 2014.
///
/// This file is made available under the GNU All-Permissive License.
///
/// Copying and distribution of this file, with or without modification,
/// are permitted in any medium without royalty provided the copyright
/// notice and this notice are preserved.  This file is offered as-is,
/// without any warranty.
///
/// Created by Ian Parberry, April 2014.
/// Last updated April 10, 2014.

#pragma once

#include "Common.h"

/// \brief The finite amortized 2D noise class.
///
/// The 2D amortized noise class implements the 2D amortized noise
/// algorithm using Perlin's finite hash table.

class CFiniteAmortizedNoise2D{  
  private: //Perlin noise tables
    float g2[B][2]; ///< Perlin gradient table.
    int p[B];  ///< Perlin permutation table.
    void initPerlinNoiseTables(); ///< Initialize Perlin noise tables.

  private: //Amortized noise stuff
    float *uax; ///< X coordinate of u used to compute a.
    float *vax; ///< X coordinate of v used to compute a.
    float *ubx; ///< X coordinate of u used to compute b.
    float *vbx; ///< X coordinate of v used to compute b.
    float *uay; ///< Y coordinate of u used to compute a.
    float *vay; ///< Y coordinate of v used to compute a.
    float *uby; ///< Y coordinate of u used to compute b.
    float *vby; ///< Y coordinate of v used to compute b.
    float* spline; ///< Spline table.

    float getRandomUnitFloat(); ///< Random float in range [-1, 1].

    unsigned int h(const unsigned int x); ///< 1D hash function.
    unsigned int h(const unsigned int x, const unsigned int y); ///< 2D hash function.

    void FillUp(float* t, const float s, const int n); ///< Fill amortized noise table bottom up.
    void FillDn(float* t, const float s, const int n); ///< Fill amortized noise table top down.

    void initEdgeTables(const int x, const int y, const int n); ///< Initialize the amortized noise tables.
    void initSplineTable(const int n); ///< Initialize the spline table.

    float getNoise(const int i, const int j);  ///< Get one point of amortized noise. 
    void getNoise(const int n, const int i0, const int j0, float** cell);  ///< Get 1 octave of amortized noise into cell.
    void addNoise(const int n, const int i0, const int j0, const float scale, float** cell);  ///< Add 1 octave of amortized noise into cell.

  public:
    CFiniteAmortizedNoise2D(const unsigned n); ///< Constructor.
    ~CFiniteAmortizedNoise2D(); ///< Destructor.
    float generate(int x, int y, const int m0, const int m1, const int n, float** cell); ///< Generate a cell of 2D amortized noise.
}; //CFiniteAmortizedNoise2D
