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
#include "vcycle_PRECISION.h"

void smoother_PRECISION( vector_PRECISION phi, vector_PRECISION Dphi, vector_PRECISION eta,
                         int n, const int res, level_struct *l, struct Thread *threading ) {
  
  ASSERT( phi != eta );

  START_MASTER(threading);
  PROF_PRECISION_START( _SM );
  END_MASTER(threading);
  
  if ( g.method == 1 ) {
    additive_schwarz_PRECISION( phi, Dphi, eta, n, res, &(l->s_PRECISION), l, threading );
  } else if ( g.method == 2 ) {
    red_black_schwarz_PRECISION( phi, Dphi, eta, n, res, &(l->s_PRECISION), l, threading );
  } else if ( g.method == 3 ) {
    sixteen_color_schwarz_PRECISION( phi, Dphi, eta, n, res, &(l->s_PRECISION), l, threading );
  } else {
    int start = threading->start_index[l->depth];
    int end   = threading->end_index[l->depth];
    START_LOCKED_MASTER(threading)
    l->sp_PRECISION.initial_guess_zero = res;
    l->sp_PRECISION.num_restart = n;
    END_LOCKED_MASTER(threading)
    if ( g.method == 4 || g.method == 6 ) {
      if ( g.odd_even ) {
        if ( res == _RES ) {
          apply_operator_PRECISION( l->sp_PRECISION.x, phi, &(l->p_PRECISION), l, threading );
          vector_PRECISION_minus( l->sp_PRECISION.x, eta, l->sp_PRECISION.x, start, end, l );
        }
        block_to_oddeven_PRECISION( l->sp_PRECISION.b, res==_RES?l->sp_PRECISION.x:eta, l, threading );
        START_LOCKED_MASTER(threading)
        l->sp_PRECISION.initial_guess_zero = _NO_RES;
        END_LOCKED_MASTER(threading)
        if ( g.method == 6 ) {
          if ( l->depth == 0 ) g5D_solve_oddeven_PRECISION( &(l->sp_PRECISION), &(l->oe_op_PRECISION), l, threading );
          else g5D_coarse_solve_odd_even_PRECISION( &(l->sp_PRECISION), &(l->oe_op_PRECISION), l, threading );
        } else {
          if ( l->depth == 0 ) solve_oddeven_PRECISION( &(l->sp_PRECISION), &(l->oe_op_PRECISION), l, threading );
          else coarse_solve_odd_even_PRECISION( &(l->sp_PRECISION), &(l->oe_op_PRECISION), l, threading );
        }
        if ( res == _NO_RES ) {
          oddeven_to_block_PRECISION( phi, l->sp_PRECISION.x, l, threading );
        } else {
          oddeven_to_block_PRECISION( l->sp_PRECISION.b, l->sp_PRECISION.x, l, threading );
          vector_PRECISION_plus( phi, phi, l->sp_PRECISION.b, start, end, l );
        }
      } else {
        START_LOCKED_MASTER(threading)
        l->sp_PRECISION.x = phi; l->sp_PRECISION.b = eta;
        END_LOCKED_MASTER(threading)
        fgmres_PRECISION( &(l->sp_PRECISION), l, threading );
      }
    } else if ( g.method == 5 ) {
      vector_PRECISION_copy( l->sp_PRECISION.b, eta, start, end, l );
      bicgstab_PRECISION( &(l->sp_PRECISION), l, threading );
      vector_PRECISION_copy( phi, l->sp_PRECISION.x, start, end, l );
    }
    ASSERT( Dphi == NULL );
  }
  
  START_MASTER(threading);
  PROF_PRECISION_STOP( _SM, n );
  END_MASTER(threading);
}


void vcycle_PRECISION( vector_PRECISION phi, vector_PRECISION Dphi, vector_PRECISION eta,
                       int res, level_struct *l, struct Thread *threading ) {

  if ( g.interpolation && l->level>0 ) {
    for ( int i=0; i<l->n_cy; i++ ) {
      if ( i==0 && res == _NO_RES ) {
        restrict_PRECISION( l->next_level->p_PRECISION.b, eta, l, threading );
      } else {
        int start = threading->start_index[l->depth];
        int end   = threading->end_index[l->depth];
        apply_operator_PRECISION( l->vbuf_PRECISION[2], phi, &(l->p_PRECISION), l, threading );
        vector_PRECISION_minus( l->vbuf_PRECISION[3], eta, l->vbuf_PRECISION[2], start, end, l );
        restrict_PRECISION( l->next_level->p_PRECISION.b, l->vbuf_PRECISION[3], l, threading );
      }
      if ( !l->next_level->idle ) {
        START_MASTER(threading)
        if ( l->depth == 0 )
          g.coarse_time -= MPI_Wtime();
        END_MASTER(threading)
        if ( l->level > 1 ) {
          if ( g.kcycle )
            fgmres_PRECISION( &(l->next_level->p_PRECISION), l->next_level, threading );
          else
            vcycle_PRECISION( l->next_level->p_PRECISION.x, NULL, l->next_level->p_PRECISION.b, _NO_RES, l->next_level, threading );
        } else {
#ifndef MUMPS_ADDS
          if ( g.odd_even ) {
            if ( g.method == 6 ) {
              g5D_coarse_solve_odd_even_PRECISION( &(l->next_level->p_PRECISION), &(l->next_level->oe_op_PRECISION), l->next_level, threading );
            } else {

              START_MASTER(threading)
              g.coarsest_time -= MPI_Wtime();
              END_MASTER(threading)

              coarse_solve_odd_even_PRECISION( &(l->next_level->p_PRECISION), &(l->next_level->oe_op_PRECISION), l->next_level, threading );

              START_MASTER(threading)
              g.coarsest_time += MPI_Wtime();
              END_MASTER(threading)

            }
          } else {
#endif
            gmres_PRECISION_struct* px = &(l->next_level->p_PRECISION);
            level_struct* lx = l->next_level;

            int start, end;
            compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);

#ifdef BLOCK_JACOBI
            if ( lx->level==0 && lx->p_PRECISION.block_jacobi_PRECISION.local_p.polyprec_PRECISION.update_lejas == 1 ) {
              // re-construct Lejas
              local_re_construct_lejas_PRECISION( lx, threading );
            }
#endif

#ifdef POLYPREC
            px->preconditioner = px->polyprec_PRECISION.preconditioner;
#endif

#ifdef BLOCK_JACOBI
            // if Block Jacobi is enabled, solve the problem : M^{-1}Ax = M^{-1}b
              if ( px->block_jacobi_PRECISION.BJ_usable == 1 ) {
                // create a backup of b
                vector_PRECISION_copy( px->block_jacobi_PRECISION.b_backup, px->b, start, end, l );
                block_jacobi_apply_PRECISION( px->b, px->block_jacobi_PRECISION.b_backup, px, lx, threading );
              }
#endif

            int fgmres_iters;

            START_MASTER(threading)
            g.coarsest_time -= MPI_Wtime();
            END_MASTER(threading)

#ifdef GCRODR
            fgmres_iters = flgcrodr_PRECISION( px, lx, threading );
#else
            fgmres_iters = fgmres_PRECISION( px, lx, threading );
#endif

            START_MASTER(threading)
            g.coarsest_time += MPI_Wtime();
            END_MASTER(threading)

            START_MASTER(threading)
            printf0("gmres iters = %d\n", fgmres_iters);
            END_MASTER(threading)

#ifdef BLOCK_JACOBI
            // restore the rhs
            if ( px->block_jacobi_PRECISION.BJ_usable == 1 ) {
              vector_PRECISION_copy( px->b, px->block_jacobi_PRECISION.b_backup, start, end, l );
            }
#endif

#ifdef POLYPREC
            if ( lx->level==0 && lx->p_PRECISION.polyprec_PRECISION.update_lejas == 1 ) {
              if ( fgmres_iters >= 1.5*px->polyprec_PRECISION.d_poly ) {
                // re-construct Lejas
                re_construct_lejas_PRECISION( lx, threading );
              }
            }
#endif

#ifndef MUMPS_ADDS
          }
#endif
        }
        START_MASTER(threading)
        if ( l->depth == 0 )
          g.coarse_time += MPI_Wtime();
        END_MASTER(threading)
      }
      if( i == 0 && res == _NO_RES )
        interpolate3_PRECISION( phi, l->next_level->p_PRECISION.x, l, threading );
      else
        interpolate_PRECISION( phi, l->next_level->p_PRECISION.x, l, threading );
      smoother_PRECISION( phi, Dphi, eta, l->post_smooth_iter, _RES, l, threading );
      res = _RES;
    }
  } else {
    smoother_PRECISION( phi, Dphi, eta, (l->depth==0)?l->n_cy:l->post_smooth_iter, res, l, threading );
  }
}
