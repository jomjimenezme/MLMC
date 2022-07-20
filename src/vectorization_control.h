/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * 
 */

#ifndef VECTORIZATION_CONTROL_H
#define VECTORIZATION_CONTROL_H

#define SIMD_LENGTH_double 4

#ifdef SSE

// #define SIMD_LENGTH_double 4

#define OPTIMIZED_COARSE_NEIGHBOR_COUPLING_double
#define OPTIMIZED_COARSE_SELF_COUPLING_double
#define INTERPOLATION_OPERATOR_LAYOUT_OPTIMIZED_double
#define INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_double
#define OPTIMIZED_NEIGHBOR_COUPLING_double
#define OPTIMIZED_SELF_COUPLING_double
#define GRAM_SCHMIDT_VECTORIZED_double
#define OPTIMIZED_LINALG_double
#define OPTIMIZED_LINALG_double

#include "sse_complex_double_intrinsic.h"

#endif

#define OPERATOR_COMPONENT_OFFSET_double  (SIMD_LENGTH_double *((l->num_eig_vect+SIMD_LENGTH_double -1)/SIMD_LENGTH_double ))

#define OPERATOR_TYPE_double double

#endif // VECTORIZATION_CONTROL_H
