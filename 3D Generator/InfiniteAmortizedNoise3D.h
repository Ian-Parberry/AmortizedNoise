/// \file InfiniteAmortizedNoise3D.h
/// \brief Header file for the 3D amortized noise class CInfiniteAmortizedNoise3D.
///
/// Copyright Ian Parberry, January 2014.
///
/// This file is made available under the GNU All-Permissive License.
///
/// Copying and distribution of this file, with or without modification,
/// are permitted in any medium without royalty provided the copyright
/// notice and this notice are preserved.  This file is offered as-is,
/// without any warranty.
///
/// Created by Ian Parberry, January 2014.
/// Last updated May 7, 2014.

#pragma once

#include "Common.h"

/// \brief The infinite amortized 3D noise class.
///
/// The 3D amortized noise class implements the 3D amortized noise
/// algorithm using Perlin's finite hash table.

class CInfiniteAmortizedNoise3D{
  private: //Amortized noise stuff
    float *uax; ///< X coordinate of u used to compute a.
    float *vax; ///< X coordinate of v used to compute a.
    float *ubx; ///< X coordinate of u used to compute b.
    float *vbx; ///< X coordinate of v used to compute b.
    float *uay; ///< Y coordinate of u used to compute a.
    float *vay; ///< Y coordinate of v used to compute a.
    float *uby; ///< Y coordinate of u used to compute b.
    float *vby; ///< Y coordinate of v used to compute b.
    float *ucx; ///< X coordinate of u used to compute c.
    float *vcx; ///< X coordinate of v used to compute c.
    float *udx; ///< X coordinate of u used to compute d.
    float *vdx; ///< X coordinate of v used to compute d.
    float *ucy; ///< Y coordinate of u used to compute c.
    float *vcy; ///< Y coordinate of v used to compute c.
    float *udy; ///< Y coordinate of u used to compute d.
    float *vdy; ///< Y coordinate of v used to compute d.
    float *uaz; ///< Z coordinate of u used to compute a.
    float *vaz; ///< Z coordinate of v used to compute a.
    float *ubz; ///< Z coordinate of u used to compute b.
    float *vbz; ///< Z coordinate of v used to compute b.
    float *ucz; ///< Z coordinate of u used to compute c.
    float *vcz; ///< Z coordinate of v used to compute c.
    float *udz; ///< Z coordinate of u used to compute d.
    float *vdz; ///< Z coordinate of v used to compute d.
    float* spline; ///< Spline array.
    unsigned int seed; ///< Hash seed.

    float getRandomUnitFloat(); ///< Random float in range [-1, 1].
    
    float h1(const unsigned int x, const unsigned int y, const unsigned int z); ///< 3D hash function.  
    float h2(const unsigned int x, const unsigned int y, const unsigned int z); ///< 3D hash function.
    
    void FillUp(float* t, float s, int n); ///< Fill amortized noise table bottom up.
    void FillDn(float* t, float s, int n); ///< Fill amortized noise table top down.
    
    float getNoise(const int i, const int j, const int k);  ///< Get one point of amortized noise. 
    void getNoise(const int n, const int i0, const int j0,
      const int k0, float*** cell); /// Get 1 octave of amortized noise into cell.
    void addNoise(const int n, const int i0, const int j0,
      const int k0, float scale, float*** cell); /// Add 1 octave of amortized noise into cell.

    void initSplineTable(const int n); ///< Initialize the spline table.
    void initEdgeTables(const int x, const int y, const int z, const int n); ///< Initialize the amortized noise tables.

  public:
    CInfiniteAmortizedNoise3D(const int n, const int s); ///< Constructor.
    ~CInfiniteAmortizedNoise3D(); ///< Destructor.
    float generate(int x, int y, int z, const int m0, const int m1, int n, float*** cell); ///< Generate a cell of 3D amortized noise.
}; //CInfiniteAmortizedNoise3D
