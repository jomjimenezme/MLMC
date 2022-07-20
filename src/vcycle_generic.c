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
#include "vcycle_double.h"

void smoother_double( vector_double phi, vector_double Dphi, vector_double eta,
                         int n, const int res, level_struct *l, struct Thread *threading ) {
  
  ASSERT( phi != eta );

  START_MASTER(threading);
  PROF_double_START( _SM );
  END_MASTER(threading);
  
  if ( g.method == 1 ) {
    additive_schwarz_double( phi, Dphi, eta, n, res, &(l->s_double), l, threading );
  } else if ( g.method == 2 ) {
    red_black_schwarz_double( phi, Dphi, eta, n, res, &(l->s_double), l, threading );
  } else if ( g.method == 3 ) {
    sixteen_color_schwarz_double( phi, Dphi, eta, n, res, &(l->s_double), l, threading );
  } else {
    int start = threading->start_index[l->depth];
    int end   = threading->end_index[l->depth];
    START_LOCKED_MASTER(threading)
    l->sp_double.initial_guess_zero = res;
    l->sp_double.num_restart = n;
    END_LOCKED_MASTER(threading)
    if ( g.method == 4 || g.method == 6 ) {
      if ( g.odd_even ) {
        if ( res == _RES ) {
          apply_operator_double( l->sp_double.x, phi, &(l->p_double), l, threading );
          vector_double_minus( l->sp_double.x, eta, l->sp_double.x, start, end, l );
        }
        block_to_oddeven_double( l->sp_double.b, res==_RES?l->sp_double.x:eta, l, threading );
        START_LOCKED_MASTER(threading)
        l->sp_double.initial_guess_zero = _NO_RES;
        END_LOCKED_MASTER(threading)
        if ( g.method == 6 ) {
          if ( l->depth == 0 ) g5D_solve_oddeven_double( &(l->sp_double), &(l->oe_op_double), l, threading );
          else g5D_coarse_solve_odd_even_double( &(l->sp_double), &(l->oe_op_double), l, threading );
        } else {
          if ( l->depth == 0 ) solve_oddeven_double( &(l->sp_double), &(l->oe_op_double), l, threading );
          else coarse_solve_odd_even_double( &(l->sp_double), &(l->oe_op_double), l, threading );
        }
        if ( res == _NO_RES ) {
          oddeven_to_block_double( phi, l->sp_double.x, l, threading );
        } else {
          oddeven_to_block_double( l->sp_double.b, l->sp_double.x, l, threading );
          vector_double_plus( phi, phi, l->sp_double.b, start, end, l );
        }
      } else {
        START_LOCKED_MASTER(threading)
        l->sp_double.x = phi; l->sp_double.b = eta;
        END_LOCKED_MASTER(threading)
        fgmres_double( &(l->sp_double), l, threading );
      }
    } else if ( g.method == 5 ) {
      vector_double_copy( l->sp_double.b, eta, start, end, l );
      bicgstab_double( &(l->sp_double), l, threading );
      vector_double_copy( phi, l->sp_double.x, start, end, l );
    }
    ASSERT( Dphi == NULL );
  }
  
  START_MASTER(threading);
  PROF_double_STOP( _SM, n );
  END_MASTER(threading);
}


void vcycle_double( vector_double phi, vector_double Dphi, vector_double eta,
                       int res, level_struct *l, struct Thread *threading ) {

  if ( g.interpolation && l->level>0 ) {
    for ( int i=0; i<l->n_cy; i++ ) {
      if ( i==0 && res == _NO_RES ) {
        restrict_double( l->next_level->p_double.b, eta, l, threading );
      } else {
        int start = threading->start_index[l->depth];
        int end   = threading->end_index[l->depth];
        apply_operator_double( l->vbuf_double[2], phi, &(l->p_double), l, threading );
        vector_double_minus( l->vbuf_double[3], eta, l->vbuf_double[2], start, end, l );
        restrict_double( l->next_level->p_double.b, l->vbuf_double[3], l, threading );
      }
      if ( !l->next_level->idle ) {
        START_MASTER(threading)
        if ( l->depth == 0 )
          g.coarse_time -= MPI_Wtime();
        END_MASTER(threading)
        if ( l->level > 1 ) {
          if ( g.kcycle )
            fgmres_double( &(l->next_level->p_double), l->next_level, threading );
          else
            vcycle_double( l->next_level->p_double.x, NULL, l->next_level->p_double.b, _NO_RES, l->next_level, threading );
        } else {
          if ( g.odd_even ) {
            if ( g.method == 6 ) {
              g5D_coarse_solve_odd_even_double( &(l->next_level->p_double), &(l->next_level->oe_op_double), l->next_level, threading );
            } else {

              START_MASTER(threading)
              g.coarsest_time -= MPI_Wtime();
              END_MASTER(threading)

              coarse_solve_odd_even_double( &(l->next_level->p_double), &(l->next_level->oe_op_double), l->next_level, threading );

              START_MASTER(threading)
              g.coarsest_time += MPI_Wtime();
              END_MASTER(threading)

            }
          } else {
            gmres_double_struct *px = &(l->next_level->p_double);
            level_struct *lx = l->next_level;
            int start, end, coarsest_iters;
            compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);
#ifdef BLOCK_JACOBI
            if ( lx->level==0 && lx->p_double.block_jacobi_double.local_p.polyprec_double.update_lejas == 1 ) {
              // re-construct Lejas
              //local_re_construct_lejas_double( lx, threading );
            }
#endif
#ifdef POLYPREC
            px->preconditioner = px->polyprec_double.preconditioner;
#endif
#ifdef BLOCK_JACOBI
            // if Block Jacobi is enabled, solve the problem : M^{-1}Ax = M^{-1}b
            if ( px->block_jacobi_double.BJ_usable == 1 ) {
              // create a backup of b
              vector_double_copy( px->block_jacobi_double.b_backup, px->b, start, end, lx );
              block_jacobi_apply_double( px->b, px->block_jacobi_double.b_backup, px, lx, threading );
            }
#endif
#ifdef GCRODR
            coarsest_iters = flgcrodr_double( px, lx, threading );
#else
            coarsest_iters = fgmres_double( px, lx, threading );
#endif
            START_MASTER(threading)
            g.avg_b1 += coarsest_iters;
            g.avg_b2 += 1;
            g.avg_crst = g.avg_b1/g.avg_b2;
            END_MASTER(threading)
            SYNC_MASTER_TO_ALL(threading)
            SYNC_CORES(threading)
#ifdef BLOCK_JACOBI
            // restore the rhs
            if ( px->block_jacobi_double.BJ_usable == 1 ) {
              vector_double_copy( px->b, px->block_jacobi_double.b_backup, start, end, lx );
            }
#endif
            SYNC_MASTER_TO_ALL(threading)
            SYNC_CORES(threading)

            //printf0("RIGHT BEFORE LEJAS CONSTRUCTION ... !!\n");
            //MPI_Finalize();
            //exit(0);

#ifdef POLYPREC
            if ( lx->level==0 && lx->p_double.polyprec_double.update_lejas == 1 ) {
              if ( coarsest_iters >= 1.5*px->polyprec_double.d_poly ) {
                // re-construct Lejas
                re_construct_lejas_double( lx, threading );
              }
            }
#endif
            START_MASTER(threading)
            printf0( "coarsest iters = %d\n",coarsest_iters );
            END_MASTER(threading)
          }
        }
        START_MASTER(threading)
        if ( l->depth == 0 )
          g.coarse_time += MPI_Wtime();
        END_MASTER(threading)
      }
      if( i == 0 && res == _NO_RES )
        interpolate3_double( phi, l->next_level->p_double.x, l, threading );
      else
        interpolate_double( phi, l->next_level->p_double.x, l, threading );
      smoother_double( phi, Dphi, eta, l->post_smooth_iter, _RES, l, threading );
      res = _RES;
    }
  } else {
    smoother_double( phi, Dphi, eta, (l->depth==0)?l->n_cy:l->post_smooth_iter, res, l, threading );
  }
}
