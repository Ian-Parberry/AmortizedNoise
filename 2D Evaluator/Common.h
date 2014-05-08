/// \file Common.h

// Copyright Ian Parberry, September 2013.
//
// This file is made available under the GNU All-Permissive License.
//
// Copying and distribution of this file, with or without modification,
// are permitted in any medium without royalty provided the copyright
// notice and this notice are preserved.  This file is offered as-is,
// without any warranty.

// Created by Ian Parberry, September 2013.
// Last updated May 6, 2014.

#pragma once

#pragma once

#define MurmurHash3_32 MurmurHash3_x86_32 ///< MurmurHash3 for x86 architectures

#define B 0x100 ///< Perlin's B, a power of 2 usually equal to 256.
#define BM 0xff ///< A bit mask, one less than B.
#define N 0x1000 ///< Perlin's N.

#define lerp(t, a, b) (a + t*(b - a)) ///< Linear interpolation.
#define s_curve(t) (t * t * (3.0f - 2.0f*t)) ///< Cubic spline.
#define s_curve2(t) (t * t * t * (10.0f + 3.0f * t * (2.0f*t - 5.0f))) ///< Quintic spline.

/// Perlin's setup macro.
#define setup(i,b0,b1,r0,r1)\
        t = vec[i] + N;\
        b0 = ((int)t) & BM;\
        b1 = (b0+1) & BM;\
        r0 = t - (int)t;\
        r1 = r0 - 1.0f;

#define at2(rx, ry) (rx*q[0] + ry*q[1]) ///< Perlin's dot product macro.

unsigned int h2(unsigned int x, unsigned int y);
