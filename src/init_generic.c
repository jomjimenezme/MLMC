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

void prof_double_init( level_struct *l ) {

/*********************************************************************************
* Initializes the profiling struct by specifying the name of every entry.
*********************************************************************************/
  
  if ( l != NULL ) {
    for ( int i=0; i<_NUM_PROF; i++ ) {
      l->prof_double.time[i] = 0.0;
      l->prof_double.count[i] = 0.0;
      l->prof_double.flop[i] = 0.0;
    }
    
    double level_ratio = 1;
    for ( int mu=0; mu<4; mu++ )
      level_ratio *= (double)g.global_lattice[l->depth][mu]/(double)g.global_lattice[0][mu];
    
    sprintf( l->prof_double.name[_GIP], "global inner product, double" );
    l->prof_double.flop[_GIP] = level_ratio*l->num_lattice_site_var*8.0;
    sprintf( l->prof_double.name[_PIP], "process inner product, double" );
    l->prof_double.flop[_PIP] = level_ratio*l->num_lattice_site_var*8.0;
    sprintf( l->prof_double.name[_LA2], "2 flop vector operations, double" );
    l->prof_double.flop[_LA2] = level_ratio*l->num_lattice_site_var*2.0;
    sprintf( l->prof_double.name[_LA6], "6 flop vector operations, double" );
    l->prof_double.flop[_LA6] = level_ratio*l->num_lattice_site_var*6.0;
    sprintf( l->prof_double.name[_LA8], "8 flop vector operations, double" );
    l->prof_double.flop[_LA8] = level_ratio*l->num_lattice_site_var*8.0;
    sprintf( l->prof_double.name[_LA], "other vector operations, double" );
    sprintf( l->prof_double.name[_GRAM_SCHMIDT], "Gram-Schmidt, double" );
    sprintf( l->prof_double.name[_GRAM_SCHMIDT_ON_AGGREGATES], "Gram-Schmidt on aggregates, double" );
    sprintf( l->prof_double.name[_CPY], "copy operations, double" );
    sprintf( l->prof_double.name[_SET], "set value operations, double" );
    sprintf( l->prof_double.name[_PR], "interpolation and restriction, double" );
    l->prof_double.flop[_PR] = level_ratio*l->num_lattice_site_var*8.0*(l->num_lattice_site_var/2);
    sprintf( l->prof_double.name[_SC], "self coupling, double" );
    l->prof_double.flop[_SC] = (l->depth==0)?552.0:level_ratio*SQUARE(l->num_lattice_site_var)*8.0;
    sprintf( l->prof_double.name[_NC], "neighbor coupling, double" );
    l->prof_double.flop[_NC] = (l->depth==0)?1368.0:level_ratio*8.0*SQUARE(l->num_lattice_site_var)*8.0;
    sprintf( l->prof_double.name[_SM], "smoother, double" );
    double ncflops = l->prof_double.flop[_SC];
    for ( int mu=0; mu<4; mu++ )
      ncflops += (l->prof_double.flop[_NC]/4.0)*((double)(l->block_lattice[mu]-1)/(double)l->block_lattice[mu]);
    l->prof_double.flop[_SM] = ncflops * (double)(g.odd_even?l->block_iter+1:l->block_iter);
    l->prof_double.flop[_SM] += (l->prof_double.flop[_NC] + l->prof_double.flop[_SC]);
    sprintf( l->prof_double.name[_OP_COMM], "operator comm init, double" );
    sprintf( l->prof_double.name[_OP_IDLE], "operator comm wait, double" );
    sprintf( l->prof_double.name[_ALLR], "allreduces, double" );
    sprintf( l->prof_double.name[_GD_COMM], "data re-distribution comm init, double" );
    sprintf( l->prof_double.name[_GD_IDLE], "data re-distribution comm wait, double" );
    sprintf( l->prof_double.name[_SM1], "smoother - pt 1, res no comm, double" );
    sprintf( l->prof_double.name[_SM2], "smoother - pt 2, solve no comm, double" );
    sprintf( l->prof_double.name[_SM3], "smoother - pt 3, res comm, double" );
    sprintf( l->prof_double.name[_SM4], "smoother - pt 4, solve comm, double" );
    
    sprintf( l->prof_double.name[_SMALL1], "Hessenberg: qr update double" );
    sprintf( l->prof_double.name[_SMALL2], "Hessenberg: bkwd subst double" );
  }
}


double prof_double_print( level_struct *l ) {
  double flop = 0;
  for ( int i=0; i<_NUM_PROF; i++ )
    if ( l->prof_double.count[i] > 0 ) {
      if ( l->prof_double.count[i] > 9999999 )
        printf0("| %37s: %8.2le(%7.1le) |\n", l->prof_double.name[i], l->prof_double.time[i], l->prof_double.count[i] );
      else
        printf0("| %37s: %8.2le(%7d) |\n", l->prof_double.name[i], l->prof_double.time[i], (int)l->prof_double.count[i] );
      flop += (double)l->prof_double.count[i] * l->prof_double.flop[i];
    }
  return flop;
}


void fine_level_double_alloc( level_struct *l ) {
  
  int n = 8;
#ifdef HAVE_TM1p1
  MALLOC( l->vbuf_double[0], complex_double, 2*n*l->vector_size );  
  for ( int i=1; i<n; i++ )
    l->vbuf_double[i] = l->vbuf_double[0] + 2*i*l->vector_size;
  MALLOC( l->p_double.b, complex_double, 2*2*l->inner_vector_size );
  l->p_double.x = l->p_double.b + 2*l->inner_vector_size;
#else
  MALLOC( l->vbuf_double[0], complex_double, n*l->vector_size );  
  for ( int i=1; i<n; i++ )
    l->vbuf_double[i] = l->vbuf_double[0] + i*l->vector_size;
  MALLOC( l->p_double.b, complex_double, 2*l->inner_vector_size );
  l->p_double.x = l->p_double.b + l->inner_vector_size;
#endif
}


void fine_level_double_free( level_struct *l ) {
  
  int n = 8;
  
#ifdef HAVE_TM1p1
  FREE( l->vbuf_double[0], complex_double, 2*n*l->vector_size );  
  for ( int i=1; i<n; i++ )
    l->vbuf_double[i] = NULL;
  FREE( l->p_double.b, complex_double, 2*2*l->inner_vector_size );
  l->p_double.x = NULL;
#else
  FREE( l->vbuf_double[0], complex_double, n*l->vector_size );  
  for ( int i=1; i<n; i++ )
    l->vbuf_double[i] = NULL;
  FREE( l->p_double.b, complex_double, 2*l->inner_vector_size );
  l->p_double.x = NULL;
#endif
}


void next_level_double_setup( level_struct *l ) {

  prof_double_init( l->next_level );
  prof_double_init( l->next_level );
  gathering_double_next_level_init( &(l->next_level->gs_double), l );  
  gathering_double_setup( &(l->next_level->gs_double), l->next_level );
  
  if ( !l->idle ) {
    coarsening_index_table_double_alloc( &(l->is_double), l );
    coarsening_index_table_double_define( &(l->is_double), &(l->s_double), l );

    if ( l->level == 1 && !l->next_level->idle ) {
      // coarsest-level solver

#if defined(GCRODR) && defined(POLYPREC) && defined(BLOCK_JACOBI) && defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
      flgcrodr_double_struct_alloc( g.coarse_iter, g.coarse_restart, l->next_level->vector_size, g.coarse_tol, 
                                       _COARSE_GMRES, _RIGHT, apply_polyprec_double,
                                       g.method==6?(g.odd_even?g5D_coarse_apply_schur_complement_double:g5D_apply_coarse_operator_double)
                                       :(g.odd_even?coarse_apply_schur_complement_double:apply_coarse_operator_double),
                                       &(l->next_level->p_double), l->next_level );
#elif defined(GCRODR) && defined(POLYPREC) && defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
      flgcrodr_double_struct_alloc( g.coarse_iter, g.coarse_restart, l->next_level->vector_size, g.coarse_tol, 
                                       _COARSE_GMRES, _RIGHT, apply_polyprec_double,
                                       g.method==6?(g.odd_even?g5D_coarse_apply_schur_complement_double:g5D_apply_coarse_operator_double)
                                       :(g.odd_even?coarse_apply_schur_complement_double:apply_coarse_operator_double),
                                       &(l->next_level->p_double), l->next_level );
#elif defined(GCRODR) && defined(MUMPS_ADDS) && defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
      flgcrodr_double_struct_alloc( g.coarse_iter, g.coarse_restart, l->next_level->vector_size, g.coarse_tol, 
                                       _COARSE_GMRES, _RIGHT, mumps_solve_double,
                                       g.method==6?(g.odd_even?g5D_coarse_apply_schur_complement_double:g5D_apply_coarse_operator_double)
                                       :(g.odd_even?coarse_apply_schur_complement_double:apply_coarse_operator_double),
                                       &(l->next_level->p_double), l->next_level );
#elif !defined(GCRODR) && defined(MUMPS_ADDS) && !defined(SINGLE_ALLREDUCE_ARNOLDI) && !defined(PIPELINED_ARNOLDI)
      fgmres_double_struct_alloc( g.coarse_iter, g.coarse_restart, l->next_level->vector_size, g.coarse_tol, 
                                     _COARSE_GMRES, _RIGHT, mumps_solve_double,
                                     g.method==6?(g.odd_even?g5D_coarse_apply_schur_complement_double:g5D_apply_coarse_operator_double)
                                     :(g.odd_even?coarse_apply_schur_complement_double:apply_coarse_operator_double),
                                     &(l->next_level->p_double), l->next_level );
#elif defined(GCRODR) && defined(MUMPS_ADDS) && !defined(SINGLE_ALLREDUCE_ARNOLDI) && !defined(PIPELINED_ARNOLDI)
      flgcrodr_double_struct_alloc( g.coarse_iter, g.coarse_restart, l->next_level->vector_size, g.coarse_tol,
                                       _COARSE_GMRES, _RIGHT, mumps_solve_double,
                                       g.method==6?(g.odd_even?g5D_coarse_apply_schur_complement_double:g5D_apply_coarse_operator_double)
                                       :(g.odd_even?coarse_apply_schur_complement_double:apply_coarse_operator_double),
                                       &(l->next_level->p_double), l->next_level );
#elif defined(GCRODR) && !defined(MUMPS_ADDS) && !defined(SINGLE_ALLREDUCE_ARNOLDI) && !defined(PIPELINED_ARNOLDI)
      flgcrodr_double_struct_alloc( g.coarse_iter, g.coarse_restart, l->next_level->vector_size, g.coarse_tol,
                                       _COARSE_GMRES, _NOTHING, NULL,
                                       g.method==6?(g.odd_even?g5D_coarse_apply_schur_complement_double:g5D_apply_coarse_operator_double)
                                       :(g.odd_even?coarse_apply_schur_complement_double:apply_coarse_operator_double),
                                       &(l->next_level->p_double), l->next_level );
#elif !defined(GCRODR) && !defined(POLYPREC) && !defined(BLOCK_JACOBI) && !defined(SINGLE_ALLREDUCE_ARNOLDI) && !defined(PIPELINED_ARNOLDI) && !defined(MUMPS_ADDS)
      fgmres_double_struct_alloc( g.coarse_iter, g.coarse_restart, l->next_level->vector_size, g.coarse_tol, 
                                     _COARSE_GMRES, _NOTHING, NULL,
                                     g.method==6?(g.odd_even?g5D_coarse_apply_schur_complement_double:g5D_apply_coarse_operator_double)
                                     :(g.odd_even?coarse_apply_schur_complement_double:apply_coarse_operator_double),
                                     &(l->next_level->p_double), l->next_level );
#else
      error0("Combination of coarsest-level algorithms not allowed!\n");
#endif

    } else {
      if ( g.kcycle ) {
        fgmres_double_struct_alloc( g.kcycle_restart, g.kcycle_max_restart, l->next_level->vector_size, g.kcycle_tol, 
                                       _K_CYCLE, _RIGHT, vcycle_double,
                                       g.method==6?g5D_apply_coarse_operator_double:apply_coarse_operator_double,
                                       &(l->next_level->p_double), l->next_level );
      } else {
#ifdef HAVE_TM1p1
        MALLOC( l->next_level->p_double.b, complex_double, 2*2*l->next_level->vector_size );
        l->next_level->p_double.x = l->next_level->p_double.b + 2*l->next_level->vector_size;
#else
        MALLOC( l->next_level->p_double.b, complex_double, 2*l->next_level->vector_size );
        l->next_level->p_double.x = l->next_level->p_double.b + l->next_level->vector_size;
#endif
        l->next_level->p_double.v_start = 0;
        l->next_level->p_double.v_end = l->next_level->inner_vector_size;
#ifdef BLOCK_JACOBI
        if ( l->next_level->level==0 ) {
          l->next_level->p_double.block_jacobi_double.local_p.v_start = 0;
          l->next_level->p_double.block_jacobi_double.local_p.v_end = l->next_level->inner_vector_size;
        }
#endif
      }
    }

    int i, n = (l->next_level->level>0)?6:4;
#ifdef HAVE_TM1p1
    MALLOC( l->next_level->vbuf_double[0], complex_double, 2*n*l->next_level->vector_size );
    for ( i=1; i<n; i++ )
      l->next_level->vbuf_double[i] = l->next_level->vbuf_double[0] + 2*i*l->next_level->vector_size;
#else
    MALLOC( l->next_level->vbuf_double[0], complex_double, n*l->next_level->vector_size );
    for ( i=1; i<n; i++ )
      l->next_level->vbuf_double[i] = l->next_level->vbuf_double[0] + i*l->next_level->vector_size;
#endif
  }
}


void next_level_double_free( level_struct *l ) {
  
  coarse_grid_correction_double_free( l );

  if ( !l->idle ) {
    if ( ( l->level == 1 && !l->next_level->idle ) || g.kcycle ) {
#ifdef GCRODR
      if ( l->level == 1 && !l->next_level->idle ) {
        flgcrodr_double_struct_free( &(l->next_level->p_double), l->next_level );
      } else {
        fgmres_double_struct_free( &(l->next_level->p_double), l->next_level );
      }
#else
      fgmres_double_struct_free( &(l->next_level->p_double), l->next_level );
#endif
    } else {
#ifdef HAVE_TM1p1
      FREE( l->next_level->p_double.b, complex_double, 2*2*l->next_level->vector_size );
#else
      FREE( l->next_level->p_double.b, complex_double, 2*l->next_level->vector_size );
#endif
    }
  
    int i, n = (l->next_level->level>0)?6:4;  
    for ( i=1; i<n; i++)
      l->next_level->vbuf_double[i] = NULL;
#ifdef HAVE_TM1p1
    FREE( l->next_level->vbuf_double[0], complex_double, 2*n*l->next_level->vector_size );
#else
    FREE( l->next_level->vbuf_double[0], complex_double, n*l->next_level->vector_size );
#endif
    coarsening_index_table_double_free( &(l->is_double), l );
  }

  gathering_double_free( &(l->next_level->gs_double), l->next_level );
}


void level_double_init( level_struct *l ) {

  for ( int i=0; i<9; i++ )
    l->vbuf_double[i] = NULL;
  
  operator_double_init( &(l->op_double) );
  operator_double_init( &(l->oe_op_double) );
  schwarz_double_init( &(l->s_double), l );
  interpolation_double_struct_init( &(l->is_double) );
#ifdef GCRODR
  if ( l->level==0 ) {
    flgcrodr_double_struct_init( &(l->p_double) );
    flgcrodr_double_struct_init( &(l->sp_double) );
  } else {
    fgmres_double_struct_init( &(l->p_double) );
    fgmres_double_struct_init( &(l->sp_double) );
  }
#else
  fgmres_double_struct_init( &(l->p_double) );
  fgmres_double_struct_init( &(l->sp_double) );
#endif
}


void vcycle_timing_double( int n, level_struct *l, struct Thread *threading ) {
  
  ASSERT( g.mixed_precision );
  vector_double v1 = NULL, v2 = NULL;
  double t0=0, t1=0;
  PUBLIC_MALLOC( v1, complex_double, l->inner_vector_size );
  PUBLIC_MALLOC( v2, complex_double, l->inner_vector_size );

  START_LOCKED_MASTER(threading)
  vector_double_define_random( v2, 0, l->inner_vector_size, l );
  END_LOCKED_MASTER(threading)
  
  START_MASTER(threading)
  t0 = MPI_Wtime();
  END_MASTER(threading)
  for ( int i=0; i<n; i++ ) {
    vcycle_double( v1, NULL, v2, _NO_RES, l, threading );
  }
  START_MASTER(threading)
  t1 = MPI_Wtime();
  printf0("100 v-cycles: %le seconds\n", t1-t0 );
  END_MASTER(threading)

  START_LOCKED_MASTER(threading)
  PUBLIC_FREE( v1, complex_double, l->inner_vector_size );
  PUBLIC_FREE( v2, complex_double, l->inner_vector_size );
  END_LOCKED_MASTER(threading)
}
