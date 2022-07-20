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

void inv_iter_2lvl_extension_setup_double( int setup_iter, level_struct *l, struct Thread *threading );
void inv_iter_inv_fcycle_double( int setup_iter, level_struct *l, struct Thread *threading );
void testvector_analysis_double( vector_double *test_vectors, level_struct *l, struct Thread *threading );
void read_tv_from_file_double( level_struct *l, struct Thread *threading );

void coarse_grid_correction_double_setup( level_struct *l, struct Thread *threading ) {
  
  if ( !l->idle ) {
    
    START_LOCKED_MASTER(threading)
    coarse_operator_double_alloc( l );
#ifndef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_double
    coarse_operator_double_setup( l->is_double.interpolation, l );
    END_LOCKED_MASTER(threading)
#else
    END_LOCKED_MASTER(threading)
    coarse_operator_double_setup_vectorized( l->is_double.operator, l, threading );
#endif
    
    START_LOCKED_MASTER(threading)
    if ( !l->next_level->idle ) {
      if ( l->next_level->level > 0 ) {
        schwarz_double_alloc( &(l->next_level->s_double), l->next_level );
        schwarz_layout_double_define( &(l->next_level->s_double), l->next_level );
      } else {
        operator_double_alloc( &(l->next_level->s_double.op), _ORDINARY, l->next_level );
        operator_double_define( &(l->next_level->s_double.op), l->next_level );
        interpolation_double_alloc( l->next_level );
      }
    } else {
      interpolation_double_dummy_alloc( l->next_level );
    }
    conf_double_gather( &(l->next_level->s_double.op), &(l->next_level->op_double), l->next_level );
    
    END_LOCKED_MASTER(threading)
    
    if ( !l->next_level->idle && l->next_level->level > 0 ) {
      START_LOCKED_MASTER(threading)
      schwarz_double_boundary_update( &(l->next_level->s_double), l->next_level );
      END_LOCKED_MASTER(threading)
      if ( g.method >= 4 && g.odd_even ) {
        START_LOCKED_MASTER(threading)
        coarse_oddeven_alloc_double( l->next_level );
        END_LOCKED_MASTER(threading)
        coarse_oddeven_setup_double( &(l->next_level->s_double.op), _REORDER, l->next_level, threading );
      }
      coarse_operator_double_set_couplings( &(l->next_level->s_double.op), l->next_level, threading );
      START_LOCKED_MASTER(threading)
      l->next_level->p_double.op = &(l->next_level->s_double.op);
#ifdef BLOCK_JACOBI
      if ( l->next_level->level==0 ) l->next_level->p_double.block_jacobi_double.local_p.op = &(l->next_level->s_double.op);
#endif
      END_LOCKED_MASTER(threading)
    }
    if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
      START_LOCKED_MASTER(threading)
      coarse_oddeven_alloc_double( l->next_level );
      END_LOCKED_MASTER(threading)
      coarse_oddeven_setup_double( &(l->next_level->s_double.op), _NO_REORDERING, l->next_level, threading );
    } else if ( !l->next_level->idle && l->next_level->level == 0 ) {
      coarse_operator_double_set_couplings( &(l->next_level->s_double.op), l->next_level, threading );
    }
  }

  if ( l->next_level->level > 0 ) {
    next_level_setup( NULL, l->next_level, threading );
    START_LOCKED_MASTER(threading)
    if ( !l->next_level->idle )
      interpolation_double_alloc( l->next_level );
    END_LOCKED_MASTER(threading)
    SYNC_HYPERTHREADS(threading)
    if ( !l->idle ) {
      for ( int i=0; i<MIN(l->next_level->num_eig_vect,l->num_eig_vect); i++ ) {
        restrict_double( l->next_level->is_double.test_vector[i], l->is_double.test_vector[i], l, threading );
      }
      START_LOCKED_MASTER(threading)
      for ( int i=MIN(l->next_level->num_eig_vect,l->num_eig_vect); i<l->next_level->num_eig_vect; i++ ) {
        if ( !l->next_level->idle )
          vector_double_define_random( l->next_level->is_double.test_vector[i], 0,
                                          l->next_level->inner_vector_size, l->next_level );
      }
      END_LOCKED_MASTER(threading)
    }
    if ( !l->next_level->idle )
      interpolation_double_define( NULL, l->next_level, threading );
    
    coarse_grid_correction_double_setup( l->next_level, threading );
  }
}


void iterative_double_setup( int setup_iter, level_struct *l, struct Thread *threading ) {
  if ( l->depth == 0 ) {
    switch ( g.interpolation ) {
      case 2: inv_iter_inv_fcycle_double( setup_iter, l, threading ); break;
      case 3: inv_iter_inv_fcycle_double( setup_iter, l, threading ); break;
      case 4: read_tv_from_file_double( l, threading ); break;
      default: inv_iter_2lvl_extension_setup_double( setup_iter, l, threading ); break;
    }
  }

  level_struct *lp = l;
  while( lp->level > 0 ) {
    testvector_analysis_double( lp->is_double.test_vector, lp, threading );
    lp = lp->next_level;
    if ( lp == NULL )
      break;
  }
}


void read_tv_from_file_double( level_struct *l, struct Thread *threading ) {
  
  if ( l->depth == 0 ) {
    if ( g.tv_io_single_file ) {
      START_LOCKED_MASTER(threading)
      vector_io_single_file( NULL, NULL, g.tv_io_file_name, _READ, l->num_eig_vect, "test vectors", l );
      END_LOCKED_MASTER(threading)
      re_setup_double( l, threading );
    } else {
      START_LOCKED_MASTER(threading)

      int n = l->num_eig_vect, i;
      char filename[STRINGLENGTH+1];
      vector_double tmp = NULL;
      
      MALLOC( tmp, complex_double, l->inner_vector_size );
      
      for ( i=0; i<n; i++ ) {
        sprintf( filename, "%s.%02d", g.tv_io_file_name, i );
        printf0("%s.%02d\n", g.tv_io_file_name, i );
        vector_io( (double*)tmp, filename, _READ, l );
        trans_double( l->is_double.test_vector[i], tmp, l->s_double.op.translation_table, l, no_threading );
      }
      
      FREE( tmp, complex_double, l->inner_vector_size );

      END_LOCKED_MASTER(threading)

      re_setup_double( l, threading );
    }
  }
}


void coarse_grid_correction_double_free( level_struct *l ) {
  
  next_level_free( l->next_level );
  
  if ( !l->idle ) {
    if ( !l->next_level->idle ) {
      if ( l->next_level->level > 0 ) {
        schwarz_double_free( &(l->next_level->s_double), l->next_level );
        if ( g.method >= 4 && g.odd_even ) {
          coarse_oddeven_free_double( l->next_level );
        }
      } else {
        operator_double_free( &(l->next_level->s_double.op), _ORDINARY, l->next_level );
        interpolation_double_free( l->next_level );
        if ( g.odd_even )
          coarse_oddeven_free_double( l->next_level );
      }
    } else {
      interpolation_double_dummy_free( l->next_level );
    }
    interpolation_double_free( l );
    coarse_operator_double_free( l );
  }  
}


void interpolation_double_define( vector_double *V, level_struct *l, struct Thread *threading ) {
  
  int k, i, n = l->num_eig_vect,
    pc = 0;
#ifdef DEBUG
  int pi = 1, pn = n*6;
#endif
  vector_double *buffer = NULL;
  int start = threading->start_index[l->depth];
  int end   = threading->end_index[l->depth];
    
  if ( V == NULL ) {
    
    PUBLIC_MALLOC( buffer, complex_double*, 3 );
    START_MASTER(threading)
    buffer[0] = NULL;
    END_MASTER(threading)
    PUBLIC_MALLOC( buffer[0], complex_double, l->vector_size*3 );
    
    START_MASTER(threading)
    for( i=1; i<3; i++)
      buffer[i] = buffer[0] + l->vector_size*i;
    if ( g.print > 0 ) printf0("initial definition --- depth: %d\n", l->depth );
#ifdef DEBUG
    if ( g.print > 0 ) { printf0("\033[0;42m\033[1;37m|"); fflush(0); }
#endif
    END_MASTER(threading)
    

    for ( k=0; k<n; k++ ) {
//       if ( l->depth == 0 ) {
        START_LOCKED_MASTER(threading)
        vector_double_define_random( l->is_double.test_vector[k], 0, l->inner_vector_size, l );
        END_LOCKED_MASTER(threading)
//       }
      
      smoother_double( buffer[0], NULL, l->is_double.test_vector[k], 1, _NO_RES, l, threading );
      vector_double_copy( l->is_double.test_vector[k], buffer[0], start, end, l );
      smoother_double( buffer[0], NULL, l->is_double.test_vector[k], g.method>=4?1:2, _NO_RES, l, threading );
      vector_double_copy( l->is_double.test_vector[k], buffer[0], start, end, l );
      smoother_double( buffer[0], NULL, l->is_double.test_vector[k], g.method>=4?1:3, _NO_RES, l, threading );
      vector_double_copy( l->is_double.test_vector[k], buffer[0], start, end, l );
        
      pc += 6;
#ifdef DEBUG
      START_MASTER(threading)
      if ( pc >= 0.2*pi*pn ) { if ( g.print > 0 ) printf0("%4d%% |", 20*pi); if ( g.my_rank == 0 ) fflush(0); pi++; }
      END_MASTER(threading)
#endif
    }
    
    PUBLIC_FREE( buffer[0], complex_double, l->vector_size*3 );
    PUBLIC_FREE( buffer, complex_double*, 3 );
    
    for ( k=0; k<n; k++ ) {
      vector_double_real_scale( l->is_double.test_vector[k], l->is_double.test_vector[k],
                                  1.0/global_norm_double( l->is_double.test_vector[k], 0, l->inner_vector_size, l, threading ),
                                  start, end, l );
    }
    
#ifdef DEBUG
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("\033[0m\n");
    END_MASTER(threading)
#endif
    
    } else {
    for ( i=0; i<n; i++ ) {
      trans_double( l->is_double.test_vector[i], V[i], l->s_double.op.translation_table, l, threading );
    }
  }

#ifndef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_double
  for ( k=0; k<n; k++ ) {
    vector_double_copy( l->is_double.interpolation[k], l->is_double.test_vector[k], start, end, l );
  }
#endif
  
  
    
  testvector_analysis_double( l->is_double.test_vector, l, threading );

#ifdef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_double
  define_interpolation_double_operator( l->is_double.test_vector, l, threading );  
  gram_schmidt_on_aggregates_double_vectorized( l->is_double.operator, n, l, threading );
#else
  gram_schmidt_on_aggregates_double( l->is_double.interpolation, n, l, threading );
  define_interpolation_double_operator( l->is_double.interpolation, l, threading );
#endif
  
}


void re_setup_double( level_struct *l, struct Thread *threading ) {
  
  if ( l->level > 0 ) {
    if ( !l->idle ) {
#ifdef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_double
      define_interpolation_double_operator( l->is_double.test_vector, l, threading );
      gram_schmidt_on_aggregates_double_vectorized( l->is_double.operator, l->num_eig_vect, l, threading );
      if ( l->depth > 0 )
        gram_schmidt_on_aggregates_double_vectorized( l->is_double.operator, l->num_eig_vect, l, threading );
      coarse_operator_double_setup_vectorized( l->is_double.operator, l, threading );
      START_LOCKED_MASTER(threading)
#else
      for ( int i=0; i<l->num_eig_vect; i++ ) {
        vector_double_copy( l->is_double.interpolation[i], l->is_double.test_vector[i],
                               threading->start_index[l->depth], threading->end_index[l->depth], l );
      }
      gram_schmidt_on_aggregates_double( l->is_double.interpolation, l->num_eig_vect, l, threading );
      if ( l->depth > 0 )
        gram_schmidt_on_aggregates_double( l->is_double.interpolation, l->num_eig_vect, l, threading );
      define_interpolation_double_operator( l->is_double.interpolation, l, threading );
      START_LOCKED_MASTER(threading)
      coarse_operator_double_setup( l->is_double.interpolation, l );
#endif
      conf_double_gather( &(l->next_level->s_double.op), &(l->next_level->op_double), l->next_level );
      END_LOCKED_MASTER(threading)
      if ( !l->next_level->idle && l->next_level->level > 0 ) {
        START_LOCKED_MASTER(threading)
        schwarz_double_boundary_update( &(l->next_level->s_double), l->next_level );
        END_LOCKED_MASTER(threading)
        if ( g.method >= 4 && g.odd_even ) {
          coarse_oddeven_setup_double( &(l->next_level->s_double.op), _REORDER, l->next_level, threading );
        } else {
          coarse_operator_double_set_couplings( &(l->next_level->s_double.op), l->next_level, threading );
        }
      }
      if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
        coarse_oddeven_setup_double( &(l->next_level->s_double.op), _NO_REORDERING, l->next_level, threading );
      } else if ( !l->next_level->idle && l->next_level->level == 0 ) {
        coarse_operator_double_set_couplings( &(l->next_level->s_double.op), l->next_level, threading );
      }
      re_setup_double( l->next_level, threading );
    }
  }

#if defined(POLYPREC) || defined(GCRODR) || defined(BLOCK_JACOBI) || defined(MUMPS_ADDS)

  else {
    // this runs on level 0 only
#ifdef POLYPREC
    l->p_double.polyprec_double.update_lejas = 1;
    l->p_double.polyprec_double.preconditioner = NULL;
#endif
#ifdef GCRODR
    l->p_double.gcrodr_double.update_CU = 1;
    l->p_double.gcrodr_double.upd_ctr = 0;
#endif
#ifdef BLOCK_JACOBI
    l->p_double.block_jacobi_double.local_p.polyprec_double.update_lejas = 1;
    l->p_double.block_jacobi_double.BJ_usable = 0;
#endif
#ifdef MUMPS_ADDS
    mumps_setup_double(l, threading);       //setup vals, Is, Js
    START_MASTER(threading)
    double t0,t1;
    t0 = MPI_Wtime();
    g.mumps_id.job = 2;        //factorize
    zmumps_c(&(g.mumps_id));
    t1 = MPI_Wtime();
    printf0("MUMPS factorize time (seconds) : %f\n",t1-t0);
    END_MASTER(threading)
    SYNC_CORES(threading)
#endif
  }
#endif
}


void inv_iter_2lvl_extension_setup_double( int setup_iter, level_struct *l, struct Thread *threading ) {

  printf0( "-----------------------------------------------WITHIN !!!\n" );

  if ( !l->idle ) {
    vector_double buf1 = NULL;
    gmres_double_struct gmres;
    
    // TODO: bugfix - threading, etc
    
    START_LOCKED_MASTER(threading)
    MALLOC( buf1, complex_double, l->vector_size );
    fgmres_double_struct_init( &gmres );
    fgmres_double_struct_alloc( g.coarse_iter, g.coarse_restart, l->next_level->vector_size, g.coarse_tol, 
                                   _COARSE_GMRES, _NOTHING, NULL, apply_coarse_operator_double, &gmres, l->next_level );
    
    if ( g.odd_even && l->next_level->level == 0 ){
      gmres.v_end = l->next_level->oe_op_double.num_even_sites*l->next_level->num_lattice_site_var;
      printf0("FIXME ASAP !\n");
    }
    END_LOCKED_MASTER(threading)
    
    for ( int k=0; k<setup_iter; k++ ) {
      int pc = 0;
#ifdef DEBUG
      int pi = 1, pn = l->num_eig_vect*l->post_smooth_iter;
#endif
      START_MASTER(threading)
      printf0("depth: %d, 2lvl correction step number %d...\n", l->depth, k+1 ); 
#ifdef DEBUG
      printf0("\033[0;42m\033[1;37m|"); fflush(0); 
#endif
      END_MASTER(threading)
      for ( int i=0; i<l->num_eig_vect; i++ ) {
        restrict_double( gmres.b, l->is_double.test_vector[i], l, threading );
        if ( !l->next_level->idle ) {
          if ( g.odd_even && l->next_level->level == 0 ) {
            coarse_solve_odd_even_double( &gmres, &(l->next_level->oe_op_double), l->next_level, threading );
          } else {
            fgmres_double( &gmres, l->next_level, threading );
          }
        }
        interpolate3_double( buf1, gmres.x, l, threading );
        smoother_double( buf1, NULL, l->is_double.test_vector[i], l->post_smooth_iter, _RES, l, threading );
        vector_double_real_scale( l->is_double.test_vector[i], buf1,
                                     1.0/global_norm_double( buf1, 0, l->inner_vector_size, l, threading ),
                                     threading->start_index[l->depth], threading->end_index[l->depth], l );
        pc += l->post_smooth_iter;
#ifdef DEBUG
        START_MASTER(threading)
        if ( pc >= 0.2*pi*pn ) { printf0("%4d%% |", 20*pi); fflush(0); pi++; }
        END_MASTER(threading)
#endif
      }
#ifdef DEBUG
      START_MASTER(threading)
      printf0("\033[0m\n");
      END_MASTER(threading)
#endif

#ifdef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_double
      define_interpolation_double_operator( l->is_double.test_vector, l, threading );
      gram_schmidt_on_aggregates_double_vectorized( l->is_double.operator, l->num_eig_vect, l, threading );
      if ( l->depth > 0 )
        gram_schmidt_on_aggregates_double_vectorized( l->is_double.operator, l->num_eig_vect, l, threading );
      coarse_operator_double_setup_vectorized( l->is_double.operator, l, threading );
      START_LOCKED_MASTER(threading)
#else
      for ( int i=0; i<l->num_eig_vect; i++ )
        vector_double_copy( l->is_double.interpolation[i], l->is_double.test_vector[i],
            threading->start_index[l->depth], threading->end_index[l->depth], l );
      gram_schmidt_on_aggregates_double( l->is_double.interpolation, l->num_eig_vect, l, threading );
      if ( l->depth > 0 )
        gram_schmidt_on_aggregates_double( l->is_double.interpolation, l->num_eig_vect, l, threading );
      define_interpolation_double_operator( l->is_double.interpolation, l, threading );
      START_LOCKED_MASTER(threading)
      coarse_operator_double_setup( l->is_double.interpolation, l );
#endif
      conf_double_gather( &(l->next_level->s_double.op), &(l->next_level->op_double), l->next_level );
      END_LOCKED_MASTER(threading)
      if ( !l->next_level->idle && l->next_level->level > 0 ) {
        START_LOCKED_MASTER(threading)
        schwarz_double_boundary_update( &(l->next_level->s_double), l->next_level );
        END_LOCKED_MASTER(threading)
        if ( g.method >= 4 && g.odd_even ) {
          coarse_oddeven_setup_double( &(l->next_level->s_double.op), _REORDER, l->next_level, threading );
        } else {
          coarse_operator_double_set_couplings( &(l->next_level->s_double.op), l->next_level, threading );
        }
      }
      if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
        coarse_oddeven_setup_double( &(l->next_level->s_double.op), _NO_REORDERING, l->next_level, threading );
      } else if ( !l->next_level->idle && l->next_level->level == 0 ) {
        coarse_operator_double_set_couplings( &(l->next_level->s_double.op), l->next_level, threading );
      }
    }
    
    if ( l->level > 1 )
      inv_iter_2lvl_extension_setup_double( setup_iter, l->next_level, threading );

    START_LOCKED_MASTER(threading)
    FREE( buf1, complex_double, l->vector_size );
    fgmres_double_struct_free( &gmres, l );
    END_LOCKED_MASTER(threading)
  }
}


void set_kcycle_tol_double( double tol, level_struct *l ) {
  
  if ( !l->idle ){
    l->p_double.tol = tol;
#ifdef BLOCK_JACOBI
    if ( l->level==0 ) l->p_double.block_jacobi_double.local_p.tol = tol;
#endif
  }
  
  if ( l->level > 1 )
    set_kcycle_tol_double( tol, l->next_level );
}


void test_vector_double_update( int i, level_struct *l, struct Thread *threading ) {
  
  if ( l->level > 1 )
    test_vector_double_update( i, l->next_level, threading );
  
  if ( !l->idle )
    vector_double_real_scale( l->is_double.test_vector[i], l->p_double.x,
                                 1.0/global_norm_double( l->p_double.x, 0, l->inner_vector_size, l, threading ),
                                 threading->start_index[l->depth], threading->end_index[l->depth], l );
}


void inv_iter_inv_fcycle_double( int setup_iter, level_struct *l, struct Thread *threading ) {
  
  vector_double v_buf = NULL;
  complex_double *buffer = NULL;
  
  PUBLIC_MALLOC( buffer, complex_double, 2*l->num_eig_vect );
      
  START_LOCKED_MASTER(threading)
  if ( l->depth == 0 )
    set_kcycle_tol_double( g.coarse_tol, l );
  END_LOCKED_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  
  PUBLIC_MALLOC( v_buf, complex_double, l->vector_size );
  
  if ( !l->idle ) {
    for ( int j=0; j<setup_iter; j++ ) {
      int pc = 0;
#ifdef DEBUG
      int pi = 1, pn = l->num_eig_vect*l->post_smooth_iter;
#endif
      
      START_LOCKED_MASTER(threading)
      if ( g.print > 0 ) printf0("depth: %d, bootstrap step number %d...\n", l->depth, j+1 );
#ifdef DEBUG
      if ( g.print > 0 ) { printf0("\033[0;42m\033[1;37m|"); if ( g.my_rank == 0 ) fflush(0); }
#endif
      END_LOCKED_MASTER(threading)
      
      gram_schmidt_double( l->is_double.test_vector, buffer, 0, l->num_eig_vect, l, threading );
      
      for ( int i=0; i<l->num_eig_vect; i++ ) {
        vcycle_double( l->p_double.x, NULL, l->is_double.test_vector[i], _NO_RES, l, threading );
        
        test_vector_double_update( i, l, threading );
        
        pc += l->post_smooth_iter;
#ifdef DEBUG
        START_MASTER(threading)
        if ( pc >= (int)((0.2*pi)*pn) ) { if ( g.print > 0 ) { printf0("%4d%% |", 20*pi); if ( g.my_rank == 0 ) fflush(0); } pi++; }
        END_MASTER(threading)
#endif
      }

#ifdef DEBUG
      START_MASTER(threading)
      if ( g.print > 0 ) printf0("\033[0m\n");
      END_MASTER(threading)
#endif      
      re_setup_double( l, threading );
      
      if ( l->depth == 0 && l->next_level->level > 0 ) {
        inv_iter_inv_fcycle_double( MAX(1,round( ((double)(j+1)*l->next_level->setup_iter)/
        ((double)setup_iter) )), l->next_level, threading );
      }
    }
    if ( l->depth > 0 && l->next_level->level > 0 ) {
      inv_iter_inv_fcycle_double( MAX(1, round((double)(l->next_level->setup_iter*setup_iter)/
      ((double)l->setup_iter))), l->next_level, threading );
    }
  }
  
  PUBLIC_FREE( v_buf, complex_double, l->vector_size );
  PUBLIC_FREE( buffer, complex_double, 2*l->num_eig_vect );
  
  if ( l->depth == 0 ) {
    START_LOCKED_MASTER(threading)
    set_kcycle_tol_double( g.kcycle_tol, l );
    END_LOCKED_MASTER(threading)
  }
}


void testvector_analysis_double( vector_double *test_vectors, level_struct *l, struct Thread *threading ) {
#ifdef TESTVECTOR_ANALYSIS
  START_UNTHREADED_FUNCTION(threading)
  if ( l->depth == 0 ) {
    
  complex_double lambda;
  double mu;
  printf0("--------------------------------------- depth: %d ----------------------------------------\n", l->depth );
  for ( int i=0; i<l->num_eig_vect; i++ ) {
    printf0("vector #%02d: ", i+1 );
    apply_operator_double( l->vbuf_double[3], test_vectors[i], &(l->p_double), l, no_threading );
    coarse_gamma5_double( l->vbuf_double[0], l->vbuf_double[3], 0, l->inner_vector_size, l );
    lambda = global_inner_product_double( test_vectors[i], l->vbuf_double[0], 0, l->inner_vector_size, l, no_threading );
    lambda /= global_inner_product_double( test_vectors[i], test_vectors[i], 0, l->inner_vector_size, l, no_threading );
    vector_double_saxpy( l->vbuf_double[1], l->vbuf_double[0], test_vectors[i], -lambda, 0, l->inner_vector_size, l );
    mu = global_norm_double( l->vbuf_double[1], 0, l->inner_vector_size, l, no_threading )/global_norm_double( test_vectors[i], 0, l->inner_vector_size, l, no_threading );
    printf0("singular value: %+lf%+lfi, singular vector precision: %le\n", (double)creal(lambda), (double)cimag(lambda), (double)mu );
  }
  printf0("--------------------------------------- depth: %d ----------------------------------------\n", l->depth );
  
  }
  END_UNTHREADED_FUNCTION(threading)
#endif
}
