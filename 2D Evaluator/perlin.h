/// \file perlin.h
/// \brief Header file for 2D Perlin noise functions.

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

#pragma once

void initPerlin2D(); ///< Initialize gradient and permutation tables.
float PerlinNoise2D(float, float, float, float, int); ///< Generate a cell of 2D Perlin noise.
