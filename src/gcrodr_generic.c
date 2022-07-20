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

#ifdef GCRODR


// declarations of aux functions
int  fgmresx_double( gmres_double_struct*, level_struct*, struct Thread* );
void gev_buildAB_double(complex_double**, complex_double**, complex_double**, vector_double*, vector_double*,
                           int, gmres_double_struct*, level_struct*, struct Thread*);
void order_pairs_double( vector_double, complex_double*, int*, int );
void build_CU_double( complex_double**, vector_double*, vector_double*,
                         gmres_double_struct*, level_struct*, struct Thread*, int );
int  arnoldix_step_double( vector_double *V, vector_double *Z, vector_double w,
                              complex_double **H, complex_double* buffer, int j, void (*prec)(),
                              gmres_double_struct *p, level_struct *l, struct Thread *threading );
void re_scale_Uk_double( gmres_double_struct *p, level_struct *l, struct Thread *threading );


// -------------------------------------------------------------



void flgcrodr_double_struct_init( gmres_double_struct *p ) {
  fgmres_double_struct_init( p );

  p->gcrodr_double.eigslvr.ordr_idxs = NULL;
  p->gcrodr_double.eigslvr.ordr_keyscpy = NULL;

  p->gcrodr_double.C = NULL;
  p->gcrodr_double.Cc = NULL;
  p->gcrodr_double.U = NULL;
  p->gcrodr_double.gev_A = NULL;
  p->gcrodr_double.gev_B = NULL;
  p->gcrodr_double.Bbuff = NULL;
  p->gcrodr_double.eigslvr.w = NULL;
  p->gcrodr_double.eigslvr.beta = NULL;
  p->gcrodr_double.eigslvr.vl = NULL;
  p->gcrodr_double.eigslvr.vr = NULL;
  p->gcrodr_double.eigslvr.qr_tau = NULL;

  p->gcrodr_double.Yk = NULL;
  p->gcrodr_double.Pk = NULL;

  p->gcrodr_double.QR = NULL;
  p->gcrodr_double.Q = NULL;
  p->gcrodr_double.R = NULL;
  p->gcrodr_double.Rinv = NULL;

  p->gcrodr_double.ort_B = NULL;
  p->gcrodr_double.G = NULL;
  p->gcrodr_double.Gc = NULL;
  p->gcrodr_double.hatZ = NULL;
  p->gcrodr_double.hatW = NULL;
#ifdef BLOCK_JACOBI
  p->gcrodr_double.r_aux = NULL;
#endif

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  p->gcrodr_double.PC = NULL;
  p->gcrodr_double.DPC = NULL;
#endif
}


void flgcrodr_double_struct_alloc( int m, int n, long int vl, double tol, const int type, const int prec_kind,
                                      void (*precond)(), void (*eval_op)(), gmres_double_struct *p, level_struct *l ) {

  fgmres_double_struct_alloc( m, n, vl, tol, type, prec_kind, precond, eval_op, p, l );

#ifdef HAVE_TM1p1
  vl*=2;
#endif

  if ( l->level==0 ) {
    if ( g.gcrodr_k >= p->restart_length ) {
      error0("The value of k in GCRO-DR needs to be smaller than the restart length m\n");
    }
    p->gcrodr_double.k = g.gcrodr_k;

    // g_ln is the length m+k of subspaces used in FL-GCRO-DR
    int g_ln = p->restart_length + p->gcrodr_double.k;
    int i;

    MALLOC( p->gcrodr_double.Bbuff, complex_double*, g_ln );

    p->gcrodr_double.Bbuff[0] = NULL;
    MALLOC( p->gcrodr_double.Bbuff[0], complex_double, g_ln*(g_ln+1)*2 );
    for ( i=1; i<g_ln; i++ ) {
      p->gcrodr_double.Bbuff[i] = p->gcrodr_double.Bbuff[0] + i*(g_ln+1);
    }

    MALLOC( p->gcrodr_double.eigslvr.ordr_idxs, int, g_ln );
    MALLOC( p->gcrodr_double.eigslvr.ordr_keyscpy, complex_double, g_ln );

    p->gcrodr_double.CU_usable = 0;

    // allocating C and U, which will contain the info associated to the recycling subspace
    MALLOC( p->gcrodr_double.C, vector_double, p->gcrodr_double.k );
    p->gcrodr_double.C[0] = NULL;
    MALLOC( p->gcrodr_double.C[0], complex_double, vl * p->gcrodr_double.k );
    for ( i=1; i<p->gcrodr_double.k; i++ ) {
      p->gcrodr_double.C[i] = p->gcrodr_double.C[0] + i*vl;
    }
    MALLOC( p->gcrodr_double.Cc, vector_double, p->gcrodr_double.k );
    p->gcrodr_double.Cc[0] = NULL;
    MALLOC( p->gcrodr_double.Cc[0], complex_double, vl * p->gcrodr_double.k );
    for ( i=1; i<p->gcrodr_double.k; i++ ) {
      p->gcrodr_double.Cc[i] = p->gcrodr_double.Cc[0] + i*vl;
    }
    MALLOC( p->gcrodr_double.U, vector_double, p->gcrodr_double.k );
    p->gcrodr_double.U[0] = NULL;
    MALLOC( p->gcrodr_double.U[0], complex_double, vl * p->gcrodr_double.k );
    for ( i=1; i<p->gcrodr_double.k; i++ ) {
      p->gcrodr_double.U[i] = p->gcrodr_double.U[0] + i*vl;
    }

#ifdef BLOCK_JACOBI
    MALLOC( p->gcrodr_double.r_aux, complex_double, vl );
#endif

    MALLOC( p->gcrodr_double.gev_A, complex_double*, g_ln );
    MALLOC( p->gcrodr_double.gev_B, complex_double*, g_ln );
    p->gcrodr_double.gev_A[0] = NULL;
    p->gcrodr_double.gev_B[0] = NULL;
    MALLOC( p->gcrodr_double.gev_A[0], complex_double, g_ln*g_ln );
    MALLOC( p->gcrodr_double.gev_B[0], complex_double, g_ln*g_ln );
    for ( i=1; i<g_ln; i++ ) {
      p->gcrodr_double.gev_A[i] = p->gcrodr_double.gev_A[0] + i*g_ln;
      p->gcrodr_double.gev_B[i] = p->gcrodr_double.gev_B[0] + i*g_ln;
    }

    MALLOC( p->gcrodr_double.eigslvr.w, complex_double, g_ln );
    MALLOC( p->gcrodr_double.eigslvr.beta, complex_double, g_ln );

    MALLOC( p->gcrodr_double.eigslvr.vr, complex_double, g_ln*g_ln );

    // setting values for using the Generalized Eigenvalue Solver from LAPACK
    p->gcrodr_double.eigslvr.jobvl = 'N';
    p->gcrodr_double.eigslvr.jobvr = 'V';
    p->gcrodr_double.eigslvr.A = p->gcrodr_double.gev_A[0];
    p->gcrodr_double.eigslvr.lda = p->restart_length + p->gcrodr_double.k;
    p->gcrodr_double.eigslvr.B = p->gcrodr_double.gev_B[0];
    p->gcrodr_double.eigslvr.ldb = p->restart_length + p->gcrodr_double.k;
    p->gcrodr_double.eigslvr.ldvl = p->restart_length + p->gcrodr_double.k;
    p->gcrodr_double.eigslvr.ldvr = p->restart_length + p->gcrodr_double.k;
  
    // matrix Y containing Yk = Zm * Pk
    MALLOC( p->gcrodr_double.Yk, complex_double*, p->gcrodr_double.k );
    p->gcrodr_double.Yk[0] = NULL;
    MALLOC( p->gcrodr_double.Yk[0], complex_double, vl * p->gcrodr_double.k );
    for ( i=1; i<p->gcrodr_double.k; i++ ) {
      p->gcrodr_double.Yk[i] = p->gcrodr_double.Yk[0] + i*vl;
    }

    MALLOC( p->gcrodr_double.Pk, complex_double*, p->gcrodr_double.k );

    p->gcrodr_double.syst_size = vl;

    MALLOC( p->gcrodr_double.QR, complex_double*, p->gcrodr_double.k );
    p->gcrodr_double.QR[0] = NULL;
    MALLOC( p->gcrodr_double.QR[0], complex_double, (g_ln+1) * p->gcrodr_double.k );
    for ( i=1; i<p->gcrodr_double.k; i++ ) {
      p->gcrodr_double.QR[i] = p->gcrodr_double.QR[0] + i*(g_ln+1);
    }

    p->gcrodr_double.eigslvr.qr_QR = p->gcrodr_double.QR;
    p->gcrodr_double.eigslvr.qr_lda = g_ln+1;

    p->gcrodr_double.Q = p->gcrodr_double.QR;
    p->gcrodr_double.eigslvr.qr_Q = p->gcrodr_double.Q;

    MALLOC( p->gcrodr_double.R, complex_double*, p->gcrodr_double.k );
    p->gcrodr_double.R[0] = NULL;
    MALLOC( p->gcrodr_double.R[0], complex_double, p->gcrodr_double.k * p->gcrodr_double.k );
    for ( i=1; i<p->gcrodr_double.k; i++ ) {
      p->gcrodr_double.R[i] = p->gcrodr_double.R[0] + i * p->gcrodr_double.k;
    }

    p->gcrodr_double.Rinv = p->gcrodr_double.R;
    p->gcrodr_double.eigslvr.qr_R = p->gcrodr_double.R;
    p->gcrodr_double.eigslvr.qr_Rinv = p->gcrodr_double.R;

    MALLOC( p->gcrodr_double.eigslvr.qr_tau, complex_double, p->gcrodr_double.k );
  
    MALLOC( p->gcrodr_double.ort_B, complex_double*, g_ln );
    p->gcrodr_double.ort_B[0] = NULL;
    MALLOC( p->gcrodr_double.ort_B[0], complex_double, p->gcrodr_double.k * g_ln );
    for ( i=1; i<g_ln; i++ ) {
      p->gcrodr_double.ort_B[i] = p->gcrodr_double.ort_B[0] + i * p->gcrodr_double.k;
    }

    MALLOC( p->gcrodr_double.G, complex_double*, g_ln );
    p->gcrodr_double.G[0] = NULL;
    MALLOC( p->gcrodr_double.G[0], complex_double, (g_ln+1) * g_ln );
    for ( i=1; i<g_ln; i++ ) {
      p->gcrodr_double.G[i] = p->gcrodr_double.G[0] + i * (g_ln+1);
    }
    memset( p->gcrodr_double.G[0], 0.0, sizeof(complex_double)*g_ln*(g_ln+1) );
    MALLOC( p->gcrodr_double.Gc, complex_double*, g_ln );
    p->gcrodr_double.Gc[0] = NULL;
    MALLOC( p->gcrodr_double.Gc[0], complex_double, (g_ln+1) * g_ln );
    for ( i=1; i<g_ln; i++ ) {
      p->gcrodr_double.Gc[i] = p->gcrodr_double.Gc[0] + i * (g_ln+1);
    }
    memset( p->gcrodr_double.Gc[0], 0.0, sizeof(complex_double)*g_ln*(g_ln+1) );

    MALLOC( p->gcrodr_double.hatZ, complex_double*, g_ln );
    MALLOC( p->gcrodr_double.hatW, complex_double*, g_ln+1 );

    for ( i=0; i<p->gcrodr_double.k; i++ ) {
      p->gcrodr_double.hatZ[i] = p->gcrodr_double.U[i];
      p->gcrodr_double.hatW[i] = p->gcrodr_double.C[i];
    }
    for ( i=p->gcrodr_double.k; i<g_ln; i++ ) {
      if ( precond==NULL ) {
        p->gcrodr_double.hatZ[i] = p->V[i - p->gcrodr_double.k];
      } else {
        p->gcrodr_double.hatZ[i] = p->Z[i - p->gcrodr_double.k];
      }
      p->gcrodr_double.hatW[i] = p->V[i - p->gcrodr_double.k];
    }
    p->gcrodr_double.hatW[g_ln] = p->V[p->restart_length];
  
    p->gcrodr_double.orth_against_Ck = 0;
    p->gcrodr_double.finish = 0;

    p->gcrodr_double.update_CU = 1;

    p->gcrodr_double.upd_ctr = 0;

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    MALLOC( p->gcrodr_double.PC, vector_double, p->gcrodr_double.k );
    p->gcrodr_double.PC[0] = NULL;
    MALLOC( p->gcrodr_double.PC[0], complex_double, vl * p->gcrodr_double.k );

    MALLOC( p->gcrodr_double.DPC, vector_double, p->gcrodr_double.k );
    p->gcrodr_double.DPC[0] = NULL;
    MALLOC( p->gcrodr_double.DPC[0], complex_double, vl * p->gcrodr_double.k );

    for ( i=1; i<p->gcrodr_double.k; i++ )
    {
      p->gcrodr_double.PC[i] = p->gcrodr_double.PC[0] + i*vl;
      p->gcrodr_double.DPC[i] = p->gcrodr_double.DPC[0] + i*vl;
    }    

    p->gcrodr_double.recompute_DPCk_poly = 0;
    p->gcrodr_double.recompute_DPCk_plain = 0;
#endif
  }

}


void flgcrodr_double_struct_free( gmres_double_struct *p, level_struct *l ) {

  fgmres_double_struct_free( p, l );

  if ( l->level==0 ) {
    // g_ln is the length m+k of subspaces used in FL-GCRO-DR
    int g_ln = p->restart_length + p->gcrodr_double.k;

    // GCRO-DR specific
    FREE( p->gcrodr_double.Yk[0], complex_double, p->gcrodr_double.syst_size * p->gcrodr_double.k );
    FREE( p->gcrodr_double.Yk, complex_double*, p->gcrodr_double.k );
    FREE( p->gcrodr_double.C[0], complex_double, p->gcrodr_double.syst_size * p->gcrodr_double.k );
    FREE( p->gcrodr_double.C, vector_double, p->gcrodr_double.k );
    FREE( p->gcrodr_double.Cc[0], complex_double, p->gcrodr_double.syst_size * p->gcrodr_double.k );
    FREE( p->gcrodr_double.Cc, vector_double, p->gcrodr_double.k );
    FREE( p->gcrodr_double.U[0], complex_double, p->gcrodr_double.syst_size * p->gcrodr_double.k );
    FREE( p->gcrodr_double.U, vector_double, p->gcrodr_double.k );
    FREE( p->gcrodr_double.Pk, complex_double*, p->gcrodr_double.k );
    FREE( p->gcrodr_double.ort_B[0], complex_double, p->gcrodr_double.k * g_ln );
    FREE( p->gcrodr_double.ort_B, complex_double*, g_ln );
    FREE( p->gcrodr_double.G[0], complex_double, (g_ln+1) * g_ln );
    FREE( p->gcrodr_double.G, complex_double*, g_ln );
    FREE( p->gcrodr_double.Gc[0], complex_double, (g_ln+1) * g_ln );
    FREE( p->gcrodr_double.Gc, complex_double*, g_ln );
    FREE( p->gcrodr_double.hatZ, complex_double*, g_ln );
    FREE( p->gcrodr_double.hatW, complex_double*, g_ln+1 );
#ifdef BLOCK_JACOBI
    FREE( p->gcrodr_double.r_aux, complex_double, p->gcrodr_double.syst_size );
#endif

    // for eigensolving
    FREE( p->gcrodr_double.Bbuff[0], complex_double, g_ln*(g_ln+1)*2 );
    FREE( p->gcrodr_double.Bbuff, complex_double*, g_ln );
    FREE( p->gcrodr_double.gev_A[0], complex_double, g_ln*g_ln );
    FREE( p->gcrodr_double.gev_B[0], complex_double, g_ln*g_ln );
    FREE( p->gcrodr_double.gev_A, complex_double*, g_ln );
    FREE( p->gcrodr_double.gev_B, complex_double*, g_ln );

    // for QR
    FREE( p->gcrodr_double.QR[0], complex_double, (g_ln+1) * p->gcrodr_double.k );
    FREE( p->gcrodr_double.QR, complex_double*, p->gcrodr_double.k );
    FREE( p->gcrodr_double.R[0], complex_double, p->gcrodr_double.k * p->gcrodr_double.k );
    FREE( p->gcrodr_double.R, complex_double*, p->gcrodr_double.k );
    FREE( p->gcrodr_double.eigslvr.qr_tau, complex_double, p->gcrodr_double.k );

    // vectors - QR
    FREE( p->gcrodr_double.eigslvr.vr, complex_double, g_ln*g_ln );
    FREE( p->gcrodr_double.eigslvr.w, complex_double, g_ln );
    FREE( p->gcrodr_double.eigslvr.beta, complex_double, g_ln );

    // ints - ordering
    FREE( p->gcrodr_double.eigslvr.ordr_idxs, int, g_ln );
    FREE( p->gcrodr_double.eigslvr.ordr_keyscpy, complex_double, g_ln );

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    FREE( p->gcrodr_double.PC[0], complex_double, p->gcrodr_double.syst_size * p->gcrodr_double.k );
    FREE( p->gcrodr_double.PC, vector_double, p->gcrodr_double.k );

    FREE( p->gcrodr_double.DPC[0], complex_double, p->gcrodr_double.syst_size * p->gcrodr_double.k );
    FREE( p->gcrodr_double.DPC, vector_double, p->gcrodr_double.k );
#endif    
  }
}


// ASSUMING : right preconditioner (or no preconditioner)
int flgcrodr_double( gmres_double_struct *p, level_struct *l, struct Thread *threading ){

  // TODO : double-check if the following assignments of hatZ and hatW are general enough ...
  int i, g_ln = p->restart_length + p->gcrodr_double.k;
  START_MASTER(threading)
  for ( i=0; i<p->gcrodr_double.k; i++ ) {
    p->gcrodr_double.hatZ[i] = p->gcrodr_double.U[i];
    p->gcrodr_double.hatW[i] = p->gcrodr_double.C[i];
  }
  for ( i=p->gcrodr_double.k; i<g_ln; i++ ) {
    if ( p->preconditioner==NULL ) {
      p->gcrodr_double.hatZ[i] = p->V[i - p->gcrodr_double.k];
    } else {
      p->gcrodr_double.hatZ[i] = p->Z[i - p->gcrodr_double.k];
    }
    p->gcrodr_double.hatW[i] = p->V[i - p->gcrodr_double.k];
  }
  p->gcrodr_double.hatW[g_ln] = p->V[p->restart_length];
  END_MASTER(threading)

  // TODO : add <extra> profiling (this extra profiling is to be added to FGMRES) ?

  // start and end indices for vector functions depending on thread
  // NOTE : in this context, <m> changes! It is not (necessarily) p->restart_length
  int start, end, fgmresx_iter=0, m, j, ol, k=p->gcrodr_double.k;

  double beta=0, norm_r0=1;

  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  SYNC_MASTER_TO_ALL(threading);
  SYNC_CORES(threading)

  // compute initial residual
  if ( p->initial_guess_zero == 0 ) {
    apply_operator_double( p->w, p->x, p, l, threading ); // compute w = D*x
    vector_double_minus( p->r, p->b, p->w, start, end, l ); // compute r = b - w
  } else {
    vector_double_copy( p->r, p->b, start, end, l ); // compute r = b - w
  }

  SYNC_MASTER_TO_ALL(threading);
  SYNC_CORES(threading)

  beta = global_norm_double( p->r, p->v_start, p->v_end, l, threading ); // gamma_0 = norm(r)
  norm_r0 = beta;

  START_MASTER(threading);
  p->gcrodr_double.norm_r0 = norm_r0;
  END_MASTER(threading);

  START_MASTER(threading)
  // setting the following line for the upcoming call to fgmresx_double(...)
  p->gamma[0] = beta;
  END_MASTER(threading);

  if (!p->initial_guess_zero) {
    p->gcrodr_double.b_norm = global_norm_double( p->b, p->v_start, p->v_end, l, threading );
  } else {
    START_MASTER(threading)
    p->gcrodr_double.b_norm = beta;
    END_MASTER(threading);
    SYNC_MASTER_TO_ALL(threading);
  }

  START_MASTER(threading)
  p->gcrodr_double.finish = 0;
  END_MASTER(threading);

  SYNC_MASTER_TO_ALL(threading);
  SYNC_CORES(threading)

  START_MASTER(threading)
  //printf0("beginning of gcrodr ...\n");
  END_MASTER(threading)


  if ( p->gcrodr_double.CU_usable==1 ) {
    vector_double *Uk;

    if ( p->initial_guess_zero == 1 )
      vector_double_define( p->x, 0, start, end, l );

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)

    if ( p->gcrodr_double.update_CU == 1 ) {

      // Yk = copy(Uk)
      START_MASTER(threading);
      complex_double **tmp_ptr = p->gcrodr_double.U;
      p->gcrodr_double.U = p->gcrodr_double.Yk;
      p->gcrodr_double.Yk = tmp_ptr;
      for ( i=0; i<k; i++ ) {
        p->gcrodr_double.hatZ[i] = p->gcrodr_double.U[i];
      }
      END_MASTER(threading);
      SYNC_MASTER_TO_ALL(threading);

      // ------ QR = A*Yk

      START_MASTER(threading)
      //printf0("before ...\n");
      END_MASTER(threading)

      // QR = A*Yk
      for ( i=0; i<k; i++ ) {
        apply_operator_double( p->gcrodr_double.C[i], p->gcrodr_double.Yk[i], p, l, threading );
      }

      START_MASTER(threading)
      //printf0("after ...\n");
      END_MASTER(threading)

      int i_length = p->v_end - p->v_start;
      pqr_double( i_length, k, p->gcrodr_double.C, p->gcrodr_double.R, p, l, threading );
      // re-orthogonalize ?
      //pqr_double( i_length, k, p->gcrodr_double.C, p->gcrodr_double.R, p, l, threading );

      SYNC_MASTER_TO_ALL(threading);
      SYNC_CORES(threading)

      START_MASTER(threading);
      inv_tri_double( &(p->gcrodr_double.eigslvr) );
      END_MASTER(threading);
      SYNC_MASTER_TO_ALL(threading);
      complex_double **Rinv = p->gcrodr_double.Rinv;

      // by this point : Ck = Q

      // Uk = Yk * R^{-1}
      Uk = p->gcrodr_double.U;
      vector_double *Yk = p->gcrodr_double.Yk;
      for ( i=0; i<k; i++ ) {
        // set all vectors in Yk to zero, to accumulate
        vector_double_define( Uk[i], 0, start, end, l );
        // and then, multi saxpy to obtain Yk
        // (the <i+1> in the 5th parameter is due to the triangular nature of Rinv)
        vector_double_multi_saxpy( Uk[i], Yk, Rinv[i], 1, i+1, start, end, l );
      }

      SYNC_MASTER_TO_ALL(threading);
      SYNC_CORES(threading)

      p->gcrodr_double.update_CU = 0;

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
      p->gcrodr_double.recompute_DPCk_poly = 1;
      p->gcrodr_double.recompute_DPCk_plain = 1;
#endif
    }

    // x  +=  Uk * Ck^{H} * r
    // r  -=  Ck * Ck^{H} * r

    complex_double *bf = p->gcrodr_double.Bbuff[0];
    complex_double *buffer = p->gcrodr_double.Bbuff[0] + k;
    vector_double *Ck = p->gcrodr_double.C;

    Uk = p->gcrodr_double.U;

    complex_double tmpx[k];
    process_multi_inner_product_double( k, tmpx, Ck, p->r, p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    // buffer is of length m, and k<m
    for ( i=0; i<k; i++ )
      buffer[i] = tmpx[i];
    if ( g.num_processes > 1 ) {
      PROF_double_START( _ALLR );
      MPI_Allreduce( buffer, bf, k, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
      PROF_double_STOP( _ALLR, 1 );
    } else {
      for( i=0; i<k; i++ )
        bf[i] = buffer[i];
    }
    END_MASTER(threading)

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)

    for( i=0; i<k; i++ )
      vector_double_saxpy( p->x, p->x, Uk[i], bf[i], start, end, l );

    apply_operator_double( p->w, p->x, p, l, threading ); // compute w = D*x
    vector_double_minus( p->r, p->b, p->w, start, end, l ); // compute r = b - w

    // TODO : restore (?) this portion of code, but with a flag that indicates that the first k diagonal entries
    //        of G need to be recomputed
    /*
    START_MASTER(threading)
    {
      int g_ln = k + p->restart_length;
      memset( p->gcrodr_double.G[0],  0.0, sizeof(complex_double)*( g_ln*(g_ln+1) ) );
      memset( p->gcrodr_double.Gc[0], 0.0, sizeof(complex_double)*( g_ln*(g_ln+1) ) );
    }
    END_MASTER(threading)

    // compute \tilde{ Uk }
    for ( i=0; i<k; i++ ) {
      complex_double diag_term = global_norm_double( Uk[i], p->v_start, p->v_end, l, threading );
      START_MASTER(threading)
      diag_term = 1.0 / diag_term;
      p->gcrodr_double.G[i][i]  = diag_term;
      p->gcrodr_double.Gc[i][i] = diag_term;
      END_MASTER(threading)
    }
    */

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)
  } else if ( p->gcrodr_double.CU_usable==0 ) {
    // call one cycle of FGMRES

    if ( p->initial_guess_zero == 1 )
      vector_double_define( p->x, 0, start, end, l );

    START_MASTER(threading)
    l->dup_H = 1;
    END_MASTER(threading);
    SYNC_MASTER_TO_ALL(threading);

    m = fgmresx_double(p, l, threading);
    fgmresx_iter += m;

    START_MASTER(threading)
    l->dup_H = 0;
    END_MASTER(threading);

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)

    // update the solution p->x (this, from the inside, applies back-substitution on p->y)
    compute_solution_double( p->x, (p->preconditioner&&p->kind==_RIGHT)?p->Z:p->V,
                                p->y, p->gamma, p->H, m-1, 1, p, l, threading );

#ifdef BLOCK_JACOBI

    // computing the actual residual in case of Block Jacobi
    /*
    {
      p->eval_operator( p->w, p->x, p->op, l, threading ); // compute w = D*x
      vector_double_minus( p->gcrodr_double.r_aux, p->block_jacobi_double.b_backup, p->w, start, end, l ); // compute r = b - w

      double norm_r0xx = global_norm_double( p->block_jacobi_double.b_backup, start, end, l, threading );
      double betaxx = global_norm_double( p->gcrodr_double.r_aux, start, end, l, threading );

      printf0("(proc=%d) 'real' rel residual = %f\n", g.my_rank, betaxx/norm_r0xx);
      printf0("(proc=%d) norm_r0xx = %f\n", g.my_rank, norm_r0xx);

      if ( betaxx/norm_r0xx > p->tol ) {
        p->gcrodr_double.finish = 0;
      }
    }

    if ( l->level==0 && p->block_jacobi_double.BJ_usable==1 ) {
      block_jacobi_apply_double( l->p_double.block_jacobi_double.xtmp, p->w, p, l, threading );
      vector_double_minus( p->r, p->b, l->p_double.block_jacobi_double.xtmp, start, end, l ); // compute r = b - w
    } else {
      apply_operator_double( p->w, p->x, p, l, threading ); // compute w = D*x
      vector_double_minus( p->r, p->b, p->w, start, end, l ); // compute r = b - w
    }
    */

    apply_operator_double( p->w, p->x, p, l, threading ); // compute w = D*x
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)
    vector_double_minus( p->r, p->b, p->w, start, end, l ); // compute r = b - w
    double norm_r0xx = global_norm_double( p->b, start, end, l, threading );
    double betaxx = global_norm_double( p->r, start, end, l, threading );
    if ( betaxx/norm_r0xx > p->tol ) {
      p->gcrodr_double.finish = 0;
    }

#else

    apply_operator_double( p->w, p->x, p, l, threading ); // compute w = D*x
    vector_double_minus( p->r, p->b, p->w, start, end, l ); // compute r = b - w

#endif

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)

    // if m<k, there's not enough information to build the recycling subspace
    if ( m<k ) {
      return m;
    }

    if ( p->preconditioner==NULL ) {
      // build the matrices A and B used for generalized-eigensolving
      gev_buildAB_double( p->gcrodr_double.gev_A, p->gcrodr_double.gev_B, p->gcrodr_double.eigslvr.Hc,
                             p->V, p->V, m, p, l, threading );

      SYNC_MASTER_TO_ALL(threading);
      SYNC_CORES(threading)

      // build C and U
      build_CU_double( p->gcrodr_double.eigslvr.Hc, p->V, p->V, p, l, threading, m );
    } else {
      // build the matrices A and B used for generalized-eigensolving
      gev_buildAB_double( p->gcrodr_double.gev_A, p->gcrodr_double.gev_B, p->gcrodr_double.eigslvr.Hc,
                             p->V, p->Z, m, p, l, threading );

      SYNC_MASTER_TO_ALL(threading);
      SYNC_CORES(threading)

      // build C and U
      build_CU_double( p->gcrodr_double.eigslvr.Hc, p->V, p->Z, p, l, threading, m );
    }

    // FIXME : issue when disabling this ...
    START_MASTER(threading)
    p->gcrodr_double.CU_usable=1;
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading);

    //printf0( "Updated CU!\n" );

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
      p->gcrodr_double.recompute_DPCk_poly = 1;
      p->gcrodr_double.recompute_DPCk_plain = 1;
#endif

    // check if this first call to fgmresx_double was enough
    if ( p->gcrodr_double.finish==1 ) {
      return m;
    }

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)

  } else{ error0("Invalid value for p->gcrodr_double.CU_usable \n"); }

#ifdef BLOCK_JACOBI
  double norm_r0xx = global_norm_double( p->block_jacobi_double.b_backup, start, end, l, threading );
#endif

  for ( ol=0; ol < p->num_restart; ol++ )  {

    beta = global_norm_double( p->r, p->v_start, p->v_end, l, threading ); // gamma_0 = norm(r)

    START_MASTER(threading)
    // setting the following line for the upcoming call to fgmresx_double(...)
    p->gamma[0] = beta;
    END_MASTER(threading);
    SYNC_MASTER_TO_ALL(threading);
    
    START_MASTER(threading)
    l->dup_H = 1;
    p->gcrodr_double.orth_against_Ck = 1;
    END_MASTER(threading);

    // NOTE the value of m from here onwards
    m = fgmresx_double(p, l, threading);
    fgmresx_iter += m;

    START_MASTER(threading)
    l->dup_H = 0;
    p->gcrodr_double.orth_against_Ck = 0;
    END_MASTER(threading);

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)

    // After calling fgmresx_double(...), G is of
    // size ( k+(m+1) ) x ( k+m )
    START_MASTER(threading)
    // set B within G
    for ( j=0; j<m; j++ ) {
      memcpy( p->gcrodr_double.G[k+j] , p->gcrodr_double.ort_B[j], sizeof(complex_double)*k );
      memcpy( p->gcrodr_double.Gc[k+j], p->gcrodr_double.ort_B[j], sizeof(complex_double)*k );
    }
    // set H within G
    for ( j=0; j<m; j++ ) {
      memcpy( p->gcrodr_double.G[k+j]+k , p->gcrodr_double.eigslvr.Hc[j], sizeof(complex_double)*(j+2) );
      memcpy( p->gcrodr_double.Gc[k+j]+k, p->gcrodr_double.eigslvr.Hc[j], sizeof(complex_double)*(j+2) );
    }
    END_MASTER(threading);

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)

    // and, the last ingredient for the least squares problem is : <bf> as the rhs, <G> as the matrix
    
    complex_double *bf = p->gcrodr_double.Bbuff[0];
    START_MASTER(threading)
    memset(bf, 0.0, sizeof(complex_double)*(k+m+1));
    bf[k] = beta;
    END_MASTER(threading)

    START_MASTER(threading)
    gels_double( LAPACK_COL_MAJOR, 'N', k+m+1, k+m, 1, p->gcrodr_double.G[0], k+p->restart_length+1, bf, k+p->restart_length+1);
    END_MASTER(threading)

    START_MASTER(threading)
    // restoring G from Gc
    {
      int g_ln = k + p->restart_length;
      memset( p->gcrodr_double.G[0], 0.0, sizeof(complex_double)*( g_ln*(g_ln+1) ) );
    }
    for ( i=0; i<k; i++ ) { p->gcrodr_double.G[i][i] = p->gcrodr_double.Gc[i][i]; }
    END_MASTER(threading)

    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    // update the solution
    for( i=0; i<k; i++ ) {
      vector_double_saxpy( p->x, p->x, p->gcrodr_double.hatZ[i], bf[i]*p->gcrodr_double.Gc[i][i], start, end, l );
    }
    for( i=k; i<(k+m); i++ )
      vector_double_saxpy( p->x, p->x, p->gcrodr_double.hatZ[i], bf[i], start, end, l );

    // updating p->r

#ifdef BLOCK_JACOBI

    // computing the actual residual in case of Block Jacobi
    /*
    {
      p->eval_operator( p->w, p->x, p->op, l, threading ); // compute w = D*x
      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)
      vector_double_minus( p->gcrodr_double.r_aux, p->block_jacobi_double.b_backup, p->w, start, end, l ); // compute r = b - w

      //double norm_r0xx = global_norm_double( p->block_jacobi_double.b_backup, start, end, l, threading );
      double betaxx = global_norm_double( p->gcrodr_double.r_aux, start, end, l, threading );

      //START_MASTER(threading)
      //printf0("(proc=%d) 'real' rel residual = %f\n", g.my_rank, betaxx/norm_r0xx);
      //END_MASTER(threading)

      if ( betaxx/norm_r0xx > p->tol ) {
        p->gcrodr_double.finish = 0;
      }
    }

    if ( l->level==0 && p->block_jacobi_double.BJ_usable==1 ) {
      block_jacobi_apply_double( l->p_double.block_jacobi_double.xtmp, p->w, p, l, threading );
      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)
      vector_double_minus( p->r, p->b, l->p_double.block_jacobi_double.xtmp, start, end, l ); // compute r = b - w
    } else {
      apply_operator_double( p->w, p->x, p, l, threading ); // compute w = D*x
      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)
      vector_double_minus( p->r, p->b, p->w, start, end, l ); // compute r = b - w
    }
    */

    apply_operator_double( p->w, p->x, p, l, threading ); // compute w = D*x
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)
    vector_double_minus( p->r, p->b, p->w, start, end, l ); // compute r = b - w
    double norm_r0xx = global_norm_double( p->b, start, end, l, threading );
    double betaxx = global_norm_double( p->r, start, end, l, threading );
    if ( betaxx/norm_r0xx > p->tol ) {
      p->gcrodr_double.finish = 0;
    }

#else

    apply_operator_double( p->w, p->x, p, l, threading ); // compute w = D*x
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)
    vector_double_minus( p->r, p->b, p->w, start, end, l ); // compute r = b - w

#endif

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)

    if ( p->gcrodr_double.upd_ctr < g.gcrodr_upd_itrs ) {
      // build the matrices A and B used for generalized-eigensolving
      gev_buildAB_double( p->gcrodr_double.gev_A, p->gcrodr_double.gev_B, p->gcrodr_double.Gc,
                             p->gcrodr_double.hatW, p->gcrodr_double.hatZ, k+m, p, l, threading );

      SYNC_MASTER_TO_ALL(threading);
      SYNC_CORES(threading)

      // build C and U
      build_CU_double( p->gcrodr_double.Gc, p->gcrodr_double.hatW, p->gcrodr_double.hatZ, p, l, threading, k+m );

      START_MASTER(threading)
      p->gcrodr_double.upd_ctr++;
      END_MASTER(threading)
    }

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)

    // check if tolerance has been reached
    if ( p->gcrodr_double.finish==1 ) {
      return fgmresx_iter;
    }

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)

    //START_MASTER(threading)
    //if( ol==3 ){
    //  //printf0("\n\nexiting (%d)\n\n\n", iter);
    //  MPI_Barrier( MPI_COMM_WORLD );
    //  exit(0);
    //}
    //END_MASTER(threading)

  }

  SYNC_MASTER_TO_ALL(threading);
  SYNC_CORES(threading)

  return fgmresx_iter;
}



// ---------- AUXILIARY FUNCTIONS



// build A and B for the generalized eigenvalue problem within FL-GCRO-DR
void gev_buildAB_double( complex_double **A, complex_double **B, complex_double **G,
                            vector_double *W, vector_double *Z, int mk, gmres_double_struct *p,
                            level_struct *l, struct Thread *threading ){

  int start, end, i, j, k;
  int g_ln = p->restart_length + p->gcrodr_double.k;

  complex_double **Bbuff = p->gcrodr_double.Bbuff;
  complex_double *Bbuffx = Bbuff[0] + g_ln*(g_ln+1);
  // tmpy is stored in column-major
  complex_double tmpy[(mk+1)*(mk)];

  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  // -------- building B

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  for ( j=0; j<mk; j++ ) {
    process_multi_inner_product_double( mk+1, tmpy+j*(mk+1), W, Z[j], p->v_start, p->v_end, l, threading );
  }

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  START_MASTER(threading)
  if ( g.num_processes > 1 ) {
    PROF_double_START( _ALLR );
    MPI_Allreduce( tmpy, Bbuffx, (mk+1)*mk, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
    PROF_double_STOP( _ALLR, 1 );
  } else {
    for( i=0; i<(mk+1)*mk; i++ )
      Bbuffx[i] = tmpy[i];
  }

  // from Bbuffx to Bbuff
  for ( j=0; j<mk; j++ ) {
    for ( i=0; i<(mk+1); i++ ) {
      Bbuff[j][i] = (Bbuffx+j*(mk+1))[i];
    }
  }

  for ( j=0; j<mk; j++ ) {
    if (p->gcrodr_double.CU_usable == 1) {
      if (j<p->gcrodr_double.k) for (i=0; i<(mk+1); i++) Bbuff[j][i] *= G[j][j];
    }
  }
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  // FIXME : improve the following matrix-matrix multiplication by using more threads than <master>
  START_MASTER(threading)
  for ( j=0; j<mk; j++ ) {
    for ( i=0; i<mk; i++ ) {
      B[j][i] = 0.0;
      for ( k=0; k<mk+1; k++ ) {
        B[j][i] += conj_double(G[i][k])*(Bbuff[j][k]);
      }
    }
  }
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading);
  SYNC_CORES(threading)

  // -------- building A

  // FIXME : improve the following matrix-matrix multiplication by using more threads than <master>
  START_MASTER(threading)
  for ( j=0; j<mk; j++ ) {
    for ( i=0; i<mk; i++ ) {
      A[j][i] = 0.0;
      for ( k=0; k<mk+1; k++ ) {
        A[j][i] += conj_double(G[i][k])*(G[j][k]);
      }
    }
  }
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading);
  SYNC_CORES(threading)
}


int fgmresx_double( gmres_double_struct *p, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* This is a trimmed version of FGMRES, as a helper to flgcrodr_double(...)
* This function assumes:
*	-- p->r already contains the computed residual
*********************************************************************************/  

  SYNC_MASTER_TO_ALL(threading);
  SYNC_CORES(threading)

  // start and end indices for vector functions depending on thread
  int start, end, j=-1, finish=0, iter=0, il;

  double norm_r0=1, gamma_jp1=1;

  START_LOCKED_MASTER(threading)
  if ( l->level==0 && g.num_levels > 1 && g.interpolation ) p->tol = g.coarse_tol;
  END_LOCKED_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);
  
  norm_r0 = p->gcrodr_double.norm_r0;

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  vector_double_real_scale( p->V[0], p->r, 1/p->gamma[0], start, end, l ); // v_0 = r / gamma_0

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  if ( l->level == 0 && l->depth > 0 ) {
    arnoldi_step_double( p->V, p->Z, p->w, p->H, p->y, 0, p->preconditioner, p, l, threading );
  }
#endif

#ifdef BLOCK_JACOBI
  int rewinder = p->restart_length/10;
  if (rewinder<6) rewinder = 6;
  if (rewinder>10) rewinder = 10;
  // two more iterations might be enough to cover the difference between 'real' and 'fake' residuals
  for( il=0; (il < rewinder) || (il<p->restart_length && finish==0); il++) {
#else
  for( il=0; il<p->restart_length && finish==0; il++) {
#endif
    j = il; iter++;

    //START_MASTER(threading)
    //printf0( "IL : j = %d\n",j );
    //END_MASTER(threading)

    // one step of Arnoldi
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    if ( l->level == 0 && l->depth > 0 ) {
      if ( !arnoldi_step_double( p->V, p->Z, p->w, p->H, p->y, j+1, p->preconditioner, p, l, threading ) ) {
        printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+2, j+1 );
        break;
      }
    } else {
      if ( !arnoldi_step_double( p->V, p->Z, p->w, p->H, p->y, j, p->preconditioner, p, l, threading ) ) {
        printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+1, j );
        break;
      }
    }
#else
    if ( !arnoldi_step_double( p->V, p->Z, p->w, p->H, p->y, j, p->preconditioner, p, l, threading ) ) {
      printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+1, j );
      break;
    }
#endif

    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)
      
    if ( cabs( p->H[j][j+1] ) > p->tol/10 ) {
      qr_update_double( p->H, p->s, p->c, p->gamma, j, l, threading );
      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)
      gamma_jp1 = cabs( p->gamma[j+1] );

      START_MASTER(threading)
      //printf0("g (proc=%d,j=%d) rel residual (gcro-dr) = %f\n\n", g.my_rank, j, gamma_jp1/norm_r0);
      END_MASTER(threading)

      if( gamma_jp1/norm_r0 < p->tol || gamma_jp1/norm_r0 > 1E+5 ) { // if satisfied ... stop

        START_MASTER(threading)
        p->gcrodr_double.finish = 1;
        END_MASTER(threading)
        SYNC_MASTER_TO_ALL(threading)
        
        finish = 1;
        
        START_MASTER(threading)
        if ( gamma_jp1/norm_r0 > 1E+5 ) printf0("Divergence of fgmresx_double, iter = %d, level=%d\n", iter, l->level );
        END_MASTER(threading)

      }
    } else {
      START_MASTER(threading)
      printf0("from gcrodr : depth: %d, iter: %d, p->H(%d,%d) = %+lf+%lfi\n", l->depth, iter, j+1, j, CSPLIT( p->H[j][j+1] ) );
      END_MASTER(threading)
      finish = 1;
      break;
    }
  } // end of the (only and) single restart
  
  if ( l->level == 0 ) {
    START_LOCKED_MASTER(threading)
    g.coarse_iter_count += iter;
    END_LOCKED_MASTER(threading)
  }

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  return iter;
}


// order the eigenpairs returned by the generalized eigenvalue solver from LAPACK
void order_pairs_double( vector_double keys, complex_double *keys_cpy, int *output, int n ){
  int i, j, buff1;

  complex_double buff2;

  memcpy(keys_cpy, keys, sizeof(complex_double)*n);

  for ( i=0; i<n; i++ ) {
    output[i] = i;
  }

  for ( i=0; i<n; i++ ) {
    for ( j=i+1; j<n; j++ ) {
      if ( cabs(keys_cpy[i]) > cabs(keys_cpy[j]) ) {
        buff2 =  keys_cpy[i];
        keys_cpy[i] = keys_cpy[j];
        keys_cpy[j] = buff2;

        buff1 =  output[i];
        output[i] = output[j];
        output[j] = buff1;
      }
    }
  }
}


void build_CU_double( complex_double **G, vector_double *W, vector_double *Z,
                         gmres_double_struct *p, level_struct *l, struct Thread *threading, int m ){

  // --------- eigensolving first

  // FIXME : improve the following eigensolve by using more threads than <master>
  START_MASTER(threading)
  // calling LAPACK's generalized eigenvalue solver through LAPACKE
  p->gcrodr_double.eigslvr.N = m;
  gen_eigslvr_double( &(p->gcrodr_double.eigslvr) );
  // the actual eigenvalues coming of LAPACK's gen-eigensolver are w/beta
  for (int i=0; i<m; i++) {
    p->gcrodr_double.eigslvr.w[i] /= p->gcrodr_double.eigslvr.beta[i];
  }
  // p->gcrodr_double.eigslvr.ordr_idxs contains the indices to access w and vr in ascending magnitude
  order_pairs_double( p->gcrodr_double.eigslvr.w, p->gcrodr_double.eigslvr.ordr_keyscpy,
                         p->gcrodr_double.eigslvr.ordr_idxs, m );
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  // ---------------- then, computing C and U

  int i, j, kl, start, end, g_ln=p->restart_length + p->gcrodr_double.k, k=p->gcrodr_double.k;

  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  complex_double **QR = p->gcrodr_double.QR;
  complex_double **Q = p->gcrodr_double.Q;
  complex_double **R = p->gcrodr_double.R;
  complex_double **Rinv = p->gcrodr_double.eigslvr.qr_Rinv;

  vector_double *Yk = p->gcrodr_double.Yk;
  vector_double *Ck2 = p->gcrodr_double.Cc;
  vector_double *Uk = p->gcrodr_double.U;
  vector_double vr  = p->gcrodr_double.eigslvr.vr;

  int *idxs = p->gcrodr_double.eigslvr.ordr_idxs;

  // for each new eigensolution, we have a new mapping for Pk
  vector_double *Pk = p->gcrodr_double.Pk;
  START_MASTER(threading)
  for ( i=0; i<k; i++ ) {
    Pk[i] = vr + idxs[i] * g_ln;
  }
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  // compute Yk
  for ( i=0; i<k; i++ ) {
    // set all vectors in Yk to zero, to accumulate
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)
    vector_double_define( Yk[i], 0, start, end, l );
    // and then, multi saxpy to obtain Yk
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    if (p->gcrodr_double.CU_usable == 1) {
      for (j=0; j<k; j++) vector_double_saxpy( Yk[i], Yk[i], Z[j], Pk[i][j] * G[j][j], start, end, l );
      SYNC_MASTER_TO_ALL(threading);
      SYNC_CORES(threading)
      vector_double_multi_saxpy( Yk[i], Z+k, Pk[i]+k, 1, m-k, start, end, l );
    } else {
      vector_double_multi_saxpy( Yk[i], Z, Pk[i], 1, m, start, end, l );
    }
  }

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  // build the matrix for QR
  // FIXME : improve the following matrix-matrix multiplication by using more threads than <master>
  START_MASTER(threading)
  for ( j=0; j<k; j++ ) {
    // set all column j of QR to zero
    memset( QR[j], 0.0, sizeof(complex_double)*(m+1) );
    // and then do accumulations over that column
    for ( i=0; i<m; i++ ) {
      for ( kl=0; kl<(m+1); kl++ ) {
        QR[j][kl] += G[i][kl] * Pk[j][i];
      }
    }
  }
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  // QR decomposition of G*Pk
  // FIXME : improve the following matrix-matrix multiplication by using more threads than <master>
  START_MASTER(threading)
  // QR
  p->gcrodr_double.eigslvr.qr_m = m+1;
  p->gcrodr_double.eigslvr.qr_n = k;
  qr_double( &(p->gcrodr_double.eigslvr) );
  p->gcrodr_double.eigslvr.qr_k = k;
  memset( R[0], 0.0, sizeof(complex_double)*k*k );
  // compute R^{-1}
  for ( j=0; j<k; j++ ) {
    for ( i=0; i<k; i++ ) {
      R[j][i] = QR[j][i];
    }
  }
  inv_tri_double( &(p->gcrodr_double.eigslvr) );
  // compute Q
  q_from_qr_double( &(p->gcrodr_double.eigslvr) );
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  // compute Ck
  for ( i=0; i<k; i++ ) {
    // set all vectors in Yk to zero, to accumulate
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)
    vector_double_define( Ck2[i], 0, start, end, l );
    // and then, multi saxpy to obtain Yk
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)
    vector_double_multi_saxpy( Ck2[i], W, Q[i], 1, m+1, start, end, l );
  }

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  // and then swap pointers
  START_MASTER(threading)
  complex_double **tmp_ptr;
  tmp_ptr = p->gcrodr_double.C;
  p->gcrodr_double.C = p->gcrodr_double.Cc;
  p->gcrodr_double.Cc = tmp_ptr;
  for ( i=0; i<k; i++ ) {
    p->gcrodr_double.hatW[i] = p->gcrodr_double.C[i];
  }
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  // compute Uk
  for ( i=0; i<k; i++ ) {
    // set all vectors in Yk to zero, to accumulate
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)
    vector_double_define( Uk[i], 0, start, end, l );
    // and then, multi saxpy to obtain Yk
    // (the <i+1> in the 5th parameter is due to the triangular nature of Rinv)
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)
    vector_double_multi_saxpy( Uk[i], Yk, Rinv[i], 1, i+1, start, end, l );
  }

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  // compute \tilde{ Uk }
  complex_double tmpy[k+1];
  complex_double **Bbuff = p->gcrodr_double.Bbuff;

  for ( i=0; i<k; i++ ) {
    process_multi_inner_product_double( 1, tmpy+i, &(Uk[i]), Uk[i], p->v_start, p->v_end, l, threading );
  }

  START_MASTER(threading)
  if ( g.num_processes > 1 ) {
    PROF_double_START( _ALLR );
    MPI_Allreduce( tmpy, Bbuff[0], k, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
    PROF_double_STOP( _ALLR, 1 );
  } else {
    for( i=0; i<k; i++ )
      Bbuff[0][i] = tmpy[i];
  }

  for ( i=0; i<k; i++ ) {
    complex_double diag_term = sqrt( creal(Bbuff[0][i]) );
    START_MASTER(threading)
    diag_term = 1.0 / diag_term;
    p->gcrodr_double.G[i][i]  = diag_term;
    p->gcrodr_double.Gc[i][i] = diag_term;
    END_MASTER(threading)
  }
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  p->gcrodr_double.recompute_DPCk_poly = 1;
  p->gcrodr_double.recompute_DPCk_plain = 1;
#endif

  printf0( "Updated CU!\n" );

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)
}


int arnoldix_step_double( vector_double *V, vector_double *Z, vector_double w,
                             complex_double **H, complex_double* buffer, int j, void (*prec)(),
                             gmres_double_struct *p, level_struct *l, struct Thread *threading ) {

  int return_val;
  return_val = arnoldi_step_double( V, Z, w, H, buffer, j, prec, p, l, threading );
  return return_val;

}

void re_scale_Uk_double( gmres_double_struct *p, level_struct *l, struct Thread *threading ){
  int i;
  int k = p->gcrodr_double.k;
  int start, end;
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  vector_double *Uk = p->gcrodr_double.U;

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  for ( i=0; i<k; i++ ) {
    complex_double diag_term = 1.0 / p->gcrodr_double.Gc[i][i];
    vector_double_scale( Uk[i], Uk[i], diag_term, start, end, l );
  }

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

}


#endif
