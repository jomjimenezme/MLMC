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

#ifndef LINALG_double_HEADER
  #define LINALG_double_HEADER
  
#ifdef _M10TV
  #define VECTOR_FOR( start, end, expression, update, l ) do{ \
    if ( l->depth == 0 ) { \
      for ( start; end; ) \
        FOR12( expression; update; ) \
    } else { \
      for ( start; end; ) \
        FOR20( expression; update; ) \
    } \
  } while(0)
#else
  #define VECTOR_FOR( start, end, expression, update, l ) do{ \
    if ( l->depth == 0 ) { \
      for ( start; end; ) \
        FOR12( expression; update; ) \
    } else { \
      for ( start; end; ) \
        FOR2( expression; update; ) \
    } \
  } while(0)
#endif


#ifdef _M10TV
  #define REAL_VECTOR_FOR( start, end, expression, update, l ) do{ \
    if ( l->depth == 0 ) { \
      for ( start; end; ) \
        FOR24( expression; update; ) \
    } else { \
      for ( start; end; ) \
        FOR40( expression; update; ) \
    } \
  } while(0)
#else
  #define REAL_VECTOR_FOR( start, end, expression, update, l ) do{ \
    if ( l->depth == 0 ) { \
      for ( start; end; ) \
        FOR24( expression; update; ) \
    } else { \
      for ( start; end; ) \
        FOR4( expression; update; ) \
    } \
  } while(0)
#endif


#ifdef _M10TV
  #define THREADED_VECTOR_FOR( i, start_index, end_index, expression, update, l, threading ) do{ \
    int thread_start, thread_end; \
    if ( l->depth == 0 ) { \
      compute_core_start_end_custom(start_index, end_index, &thread_start, &thread_end, l, threading, 12); \
      for ( i=thread_start; i<thread_end; ) \
        FOR12( expression; update; ) \
    } else { \
      compute_core_start_end_custom(start_index, end_index, &thread_start, &thread_end, l, threading, 20); \
      for ( i=thread_start; i<thread_end; ) \
        FOR20( expression; update; ) \
    } \
  } while(0)
#else
  #define THREADED_VECTOR_FOR( i, start_index, end_index, expression, update, l, threading ) do{ \
    int thread_start, thread_end; \
    if ( l->depth == 0 ) { \
      compute_core_start_end_custom(start_index, end_index, &thread_start, &thread_end, l, threading, 12); \
      for ( i=thread_start; i<thread_end; ) \
        FOR12( expression; update; ) \
    } else { \
      compute_core_start_end_custom(start_index, end_index, &thread_start, &thread_end, l, threading, 2); \
      for ( i=thread_start; i<thread_end; ) \
        FOR2( expression; update; ) \
    } \
  } while(0)
#endif
  

  struct Thread;

  complex_double global_inner_product_double( vector_double x, vector_double y, int start, int end, level_struct *l, struct Thread *threading );
  complex_double process_inner_product_double( vector_double phi, vector_double psi, int start, int end, level_struct *l, struct Thread *threading );

  void process_multi_inner_product_double( int count, complex_double *results, vector_double *phi, vector_double psi,
      int start, int end, level_struct *l, struct Thread *threading );

  double global_norm_double( vector_double phi, int start, int end, level_struct *l, struct Thread *threading );
  double process_norm_double( vector_double x, int start, int end, level_struct *l, struct Thread *threading );
  
  complex_double local_xy_over_xx_double( vector_double phi, vector_double psi, int start, int end, level_struct *l  );
  void vector_double_plus( vector_double z, vector_double x, vector_double y, int start, int end, level_struct *l ); // z := x + y
  void vector_double_minus( vector_double z, vector_double x, vector_double y, int start, int end, level_struct *l ); // z := x - y
  void vector_double_scale( vector_double z, vector_double x, complex_double alpha, int start, int end, level_struct *l ); // z := alpha*x
  void vector_double_real_scale( vector_double z, vector_double x, complex_double alpha,
                                    int start, int end, level_struct *l );
  void vector_double_saxpy( vector_double z, vector_double x, vector_double y, complex_double alpha, int start, int end, level_struct *l ); // z := x + alpha*y
  void vector_double_copy( vector_double z, vector_double x, int start, int end, level_struct *l ); // z := x
  void vector_double_projection( vector_double z, vector_double v, int k, vector_double *W, complex_double *diag, 
                                  int orthogonal, level_struct *l, Thread *threading );
  
  void gram_schmidt_on_aggregates_double( vector_double *V, const int num_vec, level_struct *l, struct Thread *threading );
  
  // Gram-Schmidt on a block of vectors, used by Block-Gram-Schmidt
  void aggregate_gram_schmidt_block_double( double *V,
      int num_vec, int leading_dimension, level_struct *l, struct Thread *threading );
  // used by Block-Gram-Schmidt
  void aggregate_orthogonalize_block_wrt_orthonormal_block_double( double *B, double *U,
      int num_vec, level_struct *l, struct Thread *threading );
  // used by Block-Gram-Schmidt
  void aggregate_block_dot_block_double( double *S, double *U, double *B,
      int num_vec, int leading_dimension, level_struct *l, struct Thread *threading );
  // used by Block-Gram-Schmidt
  void aggregate_block_minus_block_times_dot_double( double *B, double *U, double *S,
      int num_vec, int leading_dimension, level_struct *l, struct Thread *threading );

  void gram_schmidt_double( vector_double *V, complex_double *buffer, const int start, const int n, level_struct *l, struct Thread *threading );
  void setup_gram_schmidt_double( vector_double *V, vector_double g5v,
                                     complex_double *buffer, const int n, level_struct *l,
                                     struct Thread *threading );
  void spinwise_double_skalarmultiply( vector_double eta1, vector_double eta2,
                                          vector_double phi, complex_double alpha, int start, int end, level_struct *l );
  void set_boundary_double( vector_double phi, complex_double alpha, level_struct *l, struct Thread *threading );
  
#endif
