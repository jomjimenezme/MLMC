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

#include "main.h"

#ifdef SSE

#ifdef OPTIMIZED_LINALG_double
void vector_double_scale( vector_double z, vector_double x, complex_double alpha, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_double_START( _LA6 );
  
  __m128d alpha_re = _mm_set1_pd( creal_double(alpha) );
  __m128d alpha_im = _mm_set1_pd( cimag_double(alpha) );
  double *zd = (double*)(z+start);
  double *xd = (double*)(x+start);
  
  for( int i=start; i<end; ) {
    FOR6(
      {
        __m128d z_re; __m128d z_im;
        __m128d x_re; __m128d x_im;
        sse_complex_deinterleaved_load_pd( xd, &x_re, &x_im );
        
        cmul_pd( alpha_re, alpha_im, x_re, x_im, &z_re, &z_im );
        
        sse_complex_interleaved_store_pd( z_re, z_im, zd );
        zd += SIMD_LENGTH_double*2;
        xd += SIMD_LENGTH_double*2;
        i += SIMD_LENGTH_double;
      }
    )
  }
  
  if(thread == 0 && start != end)
  PROF_double_STOP( _LA6, (double)(end-start)/(double)l->inner_vector_size );
}
#endif

#ifdef OPTIMIZED_LINALG_double
void vector_double_scale( vector_double z, vector_double x, complex_double alpha, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_double_START( _LA6 );
  
  __m128 alpha_re = _mm_set1_ps( creal_double(alpha) );
  __m128 alpha_im = _mm_set1_ps( cimag_double(alpha) );
  double *zf = (double*)(z+start);
  double *xf = (double*)(x+start);
  
  if ( l->depth == 0 ) {
    for( int i=start; i<end; ) {
      FOR3(
        {
          __m128 z_re; __m128 z_im;
          __m128 x_re; __m128 x_im;
          sse_complex_deinterleaved_load( xf, &x_re, &x_im );
          
          cmul( alpha_re, alpha_im, x_re, x_im, &z_re, &z_im );
          
          sse_complex_interleaved_store( z_re, z_im, zf );
          zf += SIMD_LENGTH_double*2;
          xf += SIMD_LENGTH_double*2;
          i += SIMD_LENGTH_double;
        }
      )
    }
  } else {
    for( int i=start; i<end; ) {
      __m128 z_re; __m128 z_im;
      __m128 x_re; __m128 x_im;
      sse_complex_deinterleaved_load( xf, &x_re, &x_im );
      
      cmul( alpha_re, alpha_im, x_re, x_im, &z_re, &z_im );
      
      sse_complex_interleaved_store( z_re, z_im, zf );
      zf += SIMD_LENGTH_double*2;
      xf += SIMD_LENGTH_double*2;
      i += SIMD_LENGTH_double;
    }
  }
  
  if(thread == 0 && start != end)
  PROF_double_STOP( _LA6, (double)(end-start)/(double)l->inner_vector_size );
}
#endif

#ifdef OPTIMIZED_LINALG_double
void vector_double_saxpy( vector_double z, vector_double x, vector_double y, complex_double alpha, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if (thread == 0 && start != end )
  PROF_double_START( _LA8 );
  
  __m128 alpha_re = _mm_set1_ps( creal_double(alpha) );
  __m128 alpha_im = _mm_set1_ps( cimag_double(alpha) );
  
  if ( l->depth == 0 ) {
    for ( int i=start; i<end; ) {
      FOR3(
        {
          __m128 x_re; __m128 x_im; __m128 y_re; __m128 y_im;
          sse_complex_deinterleaved_load( (double*)(x+i), &x_re, &x_im );
          sse_complex_deinterleaved_load( (double*)(y+i), &y_re, &y_im );
          cfmadd(alpha_re, alpha_im, y_re, y_im, &x_re, &x_im);
          sse_complex_interleaved_store( x_re, x_im, (double*)(z+i) );
          i+=SIMD_LENGTH_double;
        }
      )
    }
  } else {
    for ( int i=start; i<end; ) {
      __m128 x_re; __m128 x_im; __m128 y_re; __m128 y_im;
      sse_complex_deinterleaved_load( (double*)(x+i), &x_re, &x_im );
      sse_complex_deinterleaved_load( (double*)(y+i), &y_re, &y_im );
      cfmadd(alpha_re, alpha_im, y_re, y_im, &x_re, &x_im);
      sse_complex_interleaved_store( x_re, x_im, (double*)(z+i) );
      i+=SIMD_LENGTH_double;
    }
  }
  
  if( thread == 0 && start != end )
  PROF_double_STOP( _LA8, (double)(end-start)/(double)l->inner_vector_size );
}
#endif

#ifdef OPTIMIZED_LINALG_double
void vector_double_saxpy( vector_double z, vector_double x, vector_double y, complex_double alpha, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if (thread == 0 && start != end )
  PROF_double_START( _LA8 );
  
  __m128d alpha_re = _mm_set1_pd( creal_double(alpha) );
  __m128d alpha_im = _mm_set1_pd( cimag_double(alpha) );
  
  for ( int i=start; i<end; ) {
    FOR6(
      {
        __m128d x_re; __m128d x_im; __m128d y_re; __m128d y_im;
        sse_complex_deinterleaved_load_pd( (double*)(x+i), &x_re, &x_im );
        sse_complex_deinterleaved_load_pd( (double*)(y+i), &y_re, &y_im );
        cfmadd_pd(alpha_re, alpha_im, y_re, y_im, &x_re, &x_im);
        sse_complex_interleaved_store_pd( x_re, x_im, (double*)(z+i) );
        i+=SIMD_LENGTH_double;
      }
    )
  }
  
  if( thread == 0 && start != end )
  PROF_double_STOP( _LA8, (double)(end-start)/(double)l->inner_vector_size );
}
#endif

#ifdef OPTIMIZED_LINALG_double
complex_double global_inner_product_double( vector_double phi, vector_double psi, int start, int end, level_struct *l, struct Thread *threading ) {
  
  PROF_double_START( _GIP, threading );
  complex_double local_alpha = 0, global_alpha = 0;

  int thread_start;
  int thread_end;
  compute_core_start_end(start, end, &thread_start, &thread_end, l, threading);
  
  SYNC_CORES(threading)
  
  __m128d alpha_re = _mm_setzero_pd();
  __m128d alpha_im = _mm_setzero_pd();
  
  if ( l->depth == 0 ) {
    for( int i=thread_start; i<thread_end; ) {
      FOR3(
        {
          __m128d phi_re; __m128d phi_im;
          __m128d psi_re; __m128d psi_im;
          sse_complex_deinterleaved_load_pd( (double*)(phi+i), &phi_re, &phi_im );
          sse_complex_deinterleaved_load_pd( (double*)(psi+i), &psi_re, &psi_im );
          cfmadd_conj_pd( phi_re, phi_im, psi_re, psi_im, &alpha_re, &alpha_im );
          i+=SIMD_LENGTH_double;
        }
      )
    }
  } else {
    for( int i=thread_start; i<thread_end; ) {
      __m128d phi_re; __m128d phi_im;
      __m128d psi_re; __m128d psi_im;
      sse_complex_deinterleaved_load_pd( (double*)(phi+i), &phi_re, &phi_im );
      sse_complex_deinterleaved_load_pd( (double*)(psi+i), &psi_re, &psi_im );
      cfmadd_conj_pd( phi_re, phi_im, psi_re, psi_im, &alpha_re, &alpha_im );
      i+=SIMD_LENGTH_double;
    }
  }
  
  local_alpha = sse_reduce_add_pd( alpha_re ) + I* sse_reduce_add_pd( alpha_im );

  // sum over cores
  START_NO_HYPERTHREADS(threading)
  ((complex_double *)threading->workspace)[threading->core] = local_alpha;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int i=1; i<threading->n_core; i++)
    ((complex_double *)threading->workspace)[0] += ((complex_double *)threading->workspace)[i];
  local_alpha = ((complex_double *)threading->workspace)[0];
  END_MASTER(threading)
  
  if ( g.num_processes > 1 ) {
    START_MASTER(threading)
    PROF_double_START( _ALLR );
    MPI_Allreduce( &local_alpha, &global_alpha, 1, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
    PROF_double_STOP( _ALLR, 1 );
    ((complex_double *)threading->workspace)[0] = global_alpha;
    END_MASTER(threading)
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    global_alpha = ((complex_double *)threading->workspace)[0];
    PROF_double_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return global_alpha;
  } else {
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    local_alpha = ((complex_double *)threading->workspace)[0];
    PROF_double_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return local_alpha;
  }
}
#endif

#ifdef OPTIMIZED_LINALG_double
complex_double global_inner_product_double( vector_double phi, vector_double psi, int start, int end, level_struct *l, struct Thread *threading ) {
  
  PROF_double_START( _GIP, threading );
  complex_double local_alpha = 0, global_alpha = 0;

  int thread_start;
  int thread_end;
  compute_core_start_end(start, end, &thread_start, &thread_end, l, threading);
  
  SYNC_CORES(threading)
  
  __m128 alpha_re = _mm_setzero_ps();
  __m128 alpha_im = _mm_setzero_ps();
  
  double *phif = (double*)(phi+thread_start);
  double *psif = (double*)(psi+thread_start);
  
  if ( l->depth == 0 ) {
    for( int i=thread_start; i<thread_end; ) {
      FOR3(
        {
          __m128 phi_re; __m128 phi_im;
          __m128 psi_re; __m128 psi_im;
          sse_complex_deinterleaved_load( phif, &phi_re, &phi_im );
          sse_complex_deinterleaved_load( psif, &psi_re, &psi_im );
          cfmadd_conj( phi_re, phi_im, psi_re, psi_im, &alpha_re, &alpha_im );
          phif += 8;
          psif += 8;
          i+=SIMD_LENGTH_double;
        }
      )
    }
  } else {
    for( int i=thread_start; i<thread_end; ) {
      __m128 phi_re; __m128 phi_im;
      __m128 psi_re; __m128 psi_im;
      sse_complex_deinterleaved_load( phif, &phi_re, &phi_im );
      sse_complex_deinterleaved_load( psif, &psi_re, &psi_im );
      cfmadd_conj( phi_re, phi_im, psi_re, psi_im, &alpha_re, &alpha_im );
      phif += 8;
      psif += 8;
      i+=SIMD_LENGTH_double;
    }
  }
  
  local_alpha = sse_reduce_add_ps( alpha_re ) + I* sse_reduce_add_ps( alpha_im );

  // sum over cores
  START_NO_HYPERTHREADS(threading)
  ((complex_double *)threading->workspace)[threading->core] = local_alpha;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int i=1; i<threading->n_core; i++)
    ((complex_double *)threading->workspace)[0] += ((complex_double *)threading->workspace)[i];
  local_alpha = ((complex_double *)threading->workspace)[0];
  END_MASTER(threading)
  
  if ( g.num_processes > 1 ) {
    START_MASTER(threading)
    PROF_double_START( _ALLR );
    MPI_Allreduce( &local_alpha, &global_alpha, 1, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
    PROF_double_STOP( _ALLR, 1 );
    ((complex_double *)threading->workspace)[0] = global_alpha;
    END_MASTER(threading)
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    global_alpha = ((complex_double *)threading->workspace)[0];
    PROF_double_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return global_alpha;
  } else {
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    local_alpha = ((complex_double *)threading->workspace)[0];
    PROF_double_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return local_alpha;
  }
}
#endif

#ifdef OPTIMIZED_LINALG_double
double global_norm_double( vector_double x, int start, int end, level_struct *l, struct Thread *threading ) {
  
  PROF_double_START( _GIP, threading );
  
  double local_alpha = 0, global_alpha = 0;

  int thread_start;
  int thread_end;
  compute_core_start_end(start, end, &thread_start, &thread_end, l, threading);
  
  SYNC_CORES(threading)
  
  VECTOR_FOR( int i=thread_start, i<thread_end, local_alpha += NORM_SQUARE_double(x[i]), i++, l );

  // sum over cores
  START_NO_HYPERTHREADS(threading)
  ((double *)threading->workspace)[threading->core] = local_alpha;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int i=1; i<threading->n_core; i++)
    ((double *)threading->workspace)[0] += ((double *)threading->workspace)[i];
  local_alpha = ((double *)threading->workspace)[0];
  END_MASTER(threading)

  if ( g.num_processes > 1 ) {
    START_MASTER(threading)
    PROF_double_START( _ALLR );
    MPI_Allreduce( &local_alpha, &global_alpha, 1, MPI_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
    PROF_double_STOP( _ALLR, 1 );
    ((double *)threading->workspace)[0] = global_alpha;
    END_MASTER(threading)
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    global_alpha = ((double *)threading->workspace)[0];
    PROF_double_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return (double)sqrt((double)global_alpha);
  } else {
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    local_alpha = ((double *)threading->workspace)[0];
    PROF_double_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return (double)sqrt((double)local_alpha);
  }
}
#endif

#ifdef OPTIMIZED_LINALG_double
double global_norm_double( vector_double x, int start, int end, level_struct *l, struct Thread *threading ) {
  
  PROF_double_START( _GIP, threading );
  
  double local_alpha = 0, global_alpha = 0;

  int thread_start;
  int thread_end;
  compute_core_start_end(start, end, &thread_start, &thread_end, l, threading);
  
  SYNC_CORES(threading)
  
  __m128 alpha = _mm_setzero_ps(); 
  
  if ( l->depth == 0 ) {       
    for( int i=thread_start; i<thread_end; ) {
      FOR6(
        {
          __m128 phi = _mm_loadu_ps((double*)(x+i));
          alpha = sse_fmadd( phi, phi, alpha );
          i += 2;
        }
      )
    }
  } else {
    for( int i=thread_start; i<thread_end; ) {
      __m128 phi = _mm_loadu_ps((double*)(x+i));
      alpha = sse_fmadd( phi, phi, alpha );
      i += 2;
    }
  }
  
  local_alpha = sse_reduce_add_ps( alpha );

  // sum over cores
  START_NO_HYPERTHREADS(threading)
  ((double *)threading->workspace)[threading->core] = local_alpha;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int i=1; i<threading->n_core; i++)
    ((double *)threading->workspace)[0] += ((double *)threading->workspace)[i];
  local_alpha = ((double *)threading->workspace)[0];
  END_MASTER(threading)

  if ( g.num_processes > 1 ) {
    START_MASTER(threading)
    PROF_double_START( _ALLR );
    MPI_Allreduce( &local_alpha, &global_alpha, 1, MPI_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
    PROF_double_STOP( _ALLR, 1 );
    ((double *)threading->workspace)[0] = global_alpha;
    END_MASTER(threading)
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    global_alpha = ((double *)threading->workspace)[0];
    PROF_double_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return (double)sqrt((double)global_alpha);
  } else {
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    local_alpha = ((double *)threading->workspace)[0];
    PROF_double_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return (double)sqrt((double)local_alpha);
  }
}
#endif

#ifdef OPTIMIZED_LINALG_double
void vector_double_multi_saxpy( vector_double z, vector_double *V, complex_double *alpha,
                                int sign, int count, int start, int end, level_struct *l ) {

  int thread = omp_get_thread_num();
  if (thread == 0 && start != end )
  PROF_double_START( _LA8 );
  
  int flag = 0;
  __m128d alpha_re[count]; __m128d alpha_im[count];
  for ( int c=0; c<count; c++ ) {
    alpha_re[c] = _mm_set1_pd( sign*creal_double(alpha[c]) );
    alpha_im[c] = _mm_set1_pd( sign*cimag_double(alpha[c]) );
    if ( cimag_double(alpha[c]) > EPS_double || -cimag_double(alpha[c]) > EPS_double )
      flag = 1;
  }
  
  if ( flag == 0 ) {
    for ( int c=0; c<count; c++ ) {
      for ( int i=start; i<end; ) {
        FOR12(
          {
            __m128d z_re = _mm_loadu_pd( (double*)(z+i) );
            __m128d V_re = _mm_loadu_pd( (double*)(V[c]+i) );
            z_re = sse_fmadd_pd( alpha_re[c], V_re, z_re );
            _mm_storeu_pd( (double*)(z+i), z_re );
            i++;
          }
        )
      }
    }
  } else {
    for ( int c=0; c<count; c++ ) {
      for ( int i=start; i<end; ) {
        FOR6(
          {
            __m128d z_re; __m128d z_im; __m128d V_re; __m128d V_im; 
            sse_complex_deinterleaved_load_pd( (double*)(z+i), &z_re, &z_im );
            sse_complex_deinterleaved_load_pd( (double*)(V[c]+i), &V_re, &V_im );
            cfmadd_pd(alpha_re[c], alpha_im[c], V_re, V_im, &z_re, &z_im);
            sse_complex_interleaved_store_pd( z_re, z_im, (double*)(z+i) );
            i += SIMD_LENGTH_double;
          }
        )
      }
    }
  }
  
  if( thread == 0 && start != end )
  PROF_double_STOP( _LA8, (double)(count) );
}
#endif

#ifdef OPTIMIZED_LINALG_double
void vector_double_multi_saxpy( vector_double z, vector_double *V, complex_double *alpha,
                               int sign, int count, int start, int end, level_struct *l ) {

  __m128 V_re; __m128 V_im;
  __m128 z_re; __m128 z_im;
  __m128 alpha_re[count]; __m128 alpha_im[count];
  int flag = 0;
  
  int thread = omp_get_thread_num();
  if (thread == 0 && start != end )
  PROF_double_START( _LA8 );
  
  for ( int c=0; c<count; c++ ) {
    alpha_re[c] = _mm_set1_ps( sign*creal_double(alpha[c]) );
    alpha_im[c] = _mm_set1_ps( sign*cimag_double(alpha[c]) );
    if ( cimag_double(alpha[c]) > EPS_double || -cimag_double(alpha[c]) > EPS_double )
      flag = 1;
  }
  
  if ( l->depth == 0 ) {
    if ( flag == 0 ) {
      for ( int c=0; c<count; c++ ) {
        for ( int i=start; i<end; ) {
          FOR6(
            {
              z_re = _mm_loadu_ps( (double*)(z+i) );
              V_re = _mm_loadu_ps( (double*)(V[c]+i) );
              z_re = sse_fmadd( alpha_re[c], V_re, z_re );
              _mm_storeu_ps( (double*)(z+i), z_re );
              i+=2;
            }
          )
        }
      }
    } else {
      for ( int c=0; c<count; c++ ) {
        for ( int i=start; i<end; ) {
          FOR3(
            {
              sse_complex_deinterleaved_load( (double*)(z+i), &z_re, &z_im );
              sse_complex_deinterleaved_load( (double*)(V[c]+i), &V_re, &V_im );
              cfmadd(alpha_re[c], alpha_im[c], V_re, V_im, &z_re, &z_im);
              sse_complex_interleaved_store( z_re, z_im, (double*)(z+i) );
              i+=SIMD_LENGTH_double;
            }
          )
        }
      }
    }
  } else {
    if ( flag == 0 ) {
      for ( int c=0; c<count; c++ ) {
        for ( int i=start; i<end; ) {
          z_re = _mm_loadu_ps( (double*)(z+i) );
          V_re = _mm_loadu_ps( (double*)(V[c]+i) );
          z_re = sse_fmadd( alpha_re[c], V_re, z_re );
          _mm_storeu_ps( (double*)(z+i), z_re );
          i+=2;
        }
      }
    } else {
      for ( int c=0; c<count; c++ ) {
        for ( int i=start; i<end; ) {
          sse_complex_deinterleaved_load( (double*)(z+i), &z_re, &z_im );
          sse_complex_deinterleaved_load( (double*)(V[c]+i), &V_re, &V_im );
          cfmadd(alpha_re[c], alpha_im[c], V_re, V_im, &z_re, &z_im);
          sse_complex_interleaved_store( z_re, z_im, (double*)(z+i) );
          i+=SIMD_LENGTH_double;
        }
      }
    }
  }
  
  if( thread == 0 && start != end )
  PROF_double_STOP( _LA8, (double)(count) );
}
#endif

#ifdef OPTIMIZED_LINALG_double
void process_multi_inner_product_MP( int count, complex_double *results, vector_double *phi,
                                     vector_double psi, int start, int end, level_struct *l,
                                     struct Thread *threading ) {

  PROF_double_START( _PIP, threading );
  int i;
  for(int c=0; c<count; c++)
    results[c] = 0.0;

  int thread_start;
  int thread_end;

  SYNC_CORES(threading)
  
  compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 12);
  
  __m128 result_re = _mm_setzero_ps();
  __m128 result_im = _mm_setzero_ps();
  
  for( int c=0; c<count; c++) {
    for ( i=thread_start; i<thread_end; i+=12 ) {
      __m128 psi_re; __m128 psi_im;
      __m128 phi_re; __m128 phi_im;
      // deinterleave complex numbers into 4 real parts and 4 imag parts
      sse_complex_deinterleaved_load( (double*)(psi+i), &psi_re, &psi_im );
      sse_complex_deinterleaved_load( (double*)(phi[c]+i), &phi_re, &phi_im );
      
      cmul_conj(phi_re, phi_im, psi_re, psi_im, &result_re, &result_im);
      
      sse_complex_deinterleaved_load( (double*)(psi+i+4), &psi_re, &psi_im );
      sse_complex_deinterleaved_load( (double*)(phi[c]+i+4), &phi_re, &phi_im );
      
      cfmadd_conj(phi_re, phi_im, psi_re, psi_im, &result_re, &result_im);
      
      sse_complex_deinterleaved_load( (double*)(psi+i+8), &psi_re, &psi_im );
      sse_complex_deinterleaved_load( (double*)(phi[c]+i+8), &phi_re, &phi_im );
      
      cfmadd_conj(phi_re, phi_im, psi_re, psi_im, &result_re, &result_im);
      
      results[c] += sse_reduce_add_ps(result_re) + I* sse_reduce_add_ps(result_im);
    }
  }

  START_NO_HYPERTHREADS(threading)
  ((complex_double **)threading->workspace)[threading->core] = results;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int c=0; c<count; c++)
    for(int i=1; i<threading->n_core; i++)
      ((complex_double **)threading->workspace)[0][c] += ((complex_double **)threading->workspace)[i][c];
  END_MASTER(threading)
  // all threads need the result of the norm
  SYNC_MASTER_TO_ALL(threading)
  for(int c=0; c<count; c++)
    results[c] = ((complex_double **)threading->workspace)[0][c];

  PROF_double_STOP( _PIP, (double)(end-start)/(double)l->inner_vector_size, threading );
}
#endif

#ifdef OPTIMIZED_LINALG_double
void process_multi_inner_product_double( int count, complex_double *results, vector_double *phi, vector_double psi,
    int start, int end, level_struct *l, struct Thread *threading ) {

  PROF_double_START( _PIP, threading );
  int i;
  for(int c=0; c<count; c++)
    results[c] = 0.0;

  int thread_start;
  int thread_end;

  SYNC_CORES(threading)
  
  if ( l->depth == 0 ) {
    compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 12);
    for(int c=0; c<count; c++) {
      __m128 result_re = _mm_setzero_ps();
      __m128 result_im = _mm_setzero_ps();
      for ( i=thread_start; i<thread_end; ) {
        FOR3(
          {
            __m128 phi_re; __m128 phi_im;
            __m128 psi_re; __m128 psi_im;
            
            // deinterleave complex numbers into 4 real parts and 4 imag parts        
            sse_complex_deinterleaved_load( (double*)(phi[c]+i), &phi_re, &phi_im );
            sse_complex_deinterleaved_load( (double*)(psi+i), &psi_re, &psi_im );

            cfmadd_conj(phi_re, phi_im, psi_re, psi_im, &result_re, &result_im);
            i+=SIMD_LENGTH_double;
          }
        )
      }
      results[c] += sse_reduce_add_ps(result_re) + I*sse_reduce_add_ps(result_im);
    }
  } else {
    compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 4);
    for(int c=0; c<count; c++) {
      __m128 result_re = _mm_setzero_ps();
      __m128 result_im = _mm_setzero_ps();
      for ( i=thread_start; i<thread_end; i+=SIMD_LENGTH_double ) {
        __m128 phi_re; __m128 phi_im;
        __m128 psi_re; __m128 psi_im;
        
        // deinterleave complex numbers into 4 real parts and 4 imag parts        
        sse_complex_deinterleaved_load( (double*)(phi[c]+i), &phi_re, &phi_im );
        sse_complex_deinterleaved_load( (double*)(psi+i), &psi_re, &psi_im );

        cfmadd_conj(phi_re, phi_im, psi_re, psi_im, &result_re, &result_im);
      }
      results[c] += sse_reduce_add_ps(result_re) + I*sse_reduce_add_ps(result_im);
    }
  }
  
  START_NO_HYPERTHREADS(threading)
  ((complex_double **)threading->workspace)[threading->core] = results;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int c=0; c<count; c++)
    for(int i=1; i<threading->n_core; i++)
      ((complex_double **)threading->workspace)[0][c] += ((complex_double **)threading->workspace)[i][c];
  END_MASTER(threading)
  // all threads need the result of the norm
  SYNC_MASTER_TO_ALL(threading)
  for(int c=0; c<count; c++)
    results[c] = ((complex_double **)threading->workspace)[0][c];

  PROF_double_STOP( _PIP, (double)(count*(end-start))/(double)l->inner_vector_size, threading );
}
#endif

#ifdef OPTIMIZED_LINALG_double
void process_multi_inner_product_double( int count, complex_double *results, vector_double *phi, vector_double psi,
    int start, int end, level_struct *l, struct Thread *threading ) {

  PROF_double_START( _PIP, threading );
  int i;
  for(int c=0; c<count; c++)
    results[c] = 0.0;

  int thread_start;
  int thread_end;
  
  SYNC_CORES(threading)

  if ( l->depth == 0 ) {
    compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 12);
    for(int c=0; c<count; c++) {
      __m128d result_re = _mm_setzero_pd();
      __m128d result_im = _mm_setzero_pd();
      for ( i=thread_start; i<thread_end; ) {
        FOR6(
          {
            __m128d phi_re; __m128d phi_im;
            __m128d pdi_re; __m128d pdi_im;
            
            // deinterleave complex numbers into 4 real parts and 4 imag parts        
            sse_complex_deinterleaved_load_pd( (double*)(phi[c]+i), &phi_re, &phi_im );
            sse_complex_deinterleaved_load_pd( (double*)(psi+i), &pdi_re, &pdi_im );

            cfmadd_conj_pd(phi_re, phi_im, pdi_re, pdi_im, &result_re, &result_im);
            i+=SIMD_LENGTH_double;
          }
        )
      }
      results[c] += sse_reduce_add_pd(result_re) + I*sse_reduce_add_pd(result_im);
    }
  } else {
    compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 2);
    for(int c=0; c<count; c++) {
      __m128d result_re = _mm_setzero_pd();
      __m128d result_im = _mm_setzero_pd();
      for ( i=thread_start; i<thread_end; i+=SIMD_LENGTH_double ) {
        __m128d phi_re; __m128d phi_im;
        __m128d pdi_re; __m128d pdi_im;
        
        // deinterleave complex numbers into 4 real parts and 4 imag parts        
        sse_complex_deinterleaved_load_pd( (double*)(phi[c]+i), &phi_re, &phi_im );
        sse_complex_deinterleaved_load_pd( (double*)(psi+i), &pdi_re, &pdi_im );

        cfmadd_conj_pd(phi_re, phi_im, pdi_re, pdi_im, &result_re, &result_im);
      }
      results[c] += sse_reduce_add_pd(result_re) + I*sse_reduce_add_pd(result_im);
    }
  }
  
  START_NO_HYPERTHREADS(threading)
  ((complex_double **)threading->workspace)[threading->core] = results;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int c=0; c<count; c++)
    for(int i=1; i<threading->n_core; i++)
      ((complex_double **)threading->workspace)[0][c] += ((complex_double **)threading->workspace)[i][c];
  END_MASTER(threading)
  // all threads need the result of the norm
  SYNC_MASTER_TO_ALL(threading)
  for(int c=0; c<count; c++)
    results[c] = ((complex_double **)threading->workspace)[0][c];

  PROF_double_STOP( _PIP, (double)(count*(end-start))/(double)l->inner_vector_size, threading );
}
#endif

#endif // SSE

