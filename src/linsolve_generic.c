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

void print_Matrix_double(complex_double** A, int mv, int mh )
{
  int i,j;

  // printf("\n\n");
  // for (i=0; i < mv; i++)
  // {
  //     for(j=0; j < mh; j++)
  //     {
  //             fprintf(stdout, "%6.6f +i%6.6f\t", creal(A[i*mh + j]), cimag(A[i*mh+j]));
  //     }
  //     fprintf(stdout, "\n");
  // }
  // printf("--\n");

  printf("\n\n");
  for (i=0; i < mh; i++)
  {
    for(j=0; j < mv; j++)
    {
      // fprintf(stdout, "%6.6f +i%6.6f\t", creal(A[j*mh + i]), cimag(A[j*mh+i]));
      fprintf(stdout, "%6.6f +i%6.6f\t", creal(A[j][i]), cimag(A[j][i]));
    }
    fprintf(stdout, "\n");
  }
  printf("--\n");
  printf("\n\n");
}



void print_double_eigenvalues( char* desc, int n, double complex* w) {
    int j;
    printf( "\n %s\n", desc );

    for( j = 0; j < n; j++ )
    {
        printf( " (%6.2f,%6.2f)", creal(w[j]), cimag(w[j]) );
    }
    printf( "\n" );
}


void fgmres_double_struct_init( gmres_double_struct *p ) {

/*********************************************************************************
* Initializes all declared pointers with NULL.                              
*********************************************************************************/

  p->Z = NULL;
  p->V = NULL;
  p->H = NULL;
  p->x = NULL;
  p->b = NULL;
  p->r = NULL;
  p->w = NULL;
  p->y = NULL;
  p->gamma = NULL;
  p->c = NULL;
  p->s = NULL;
  p->preconditioner = NULL;
  p->eval_operator = NULL;

  // copy of Hesselnberg matrix
#if defined(GCRODR) && defined(POLYPREC)
  p->gcrodr_double.eigslvr.Hc = NULL;
  p->polyprec_double.eigslvr.Hc = NULL;
#elif defined(GCRODR)
  p->gcrodr_double.eigslvr.Hc = NULL;
#elif defined(POLYPREC)
  p->polyprec_double.eigslvr.Hc = NULL;
#endif

#ifdef POLYPREC
  p->polyprec_double.Hcc = NULL; 
  p->polyprec_double.L = NULL;
  p->polyprec_double.col_prods = NULL;
  p->polyprec_double.accum_prod = NULL;
  p->polyprec_double.product = NULL;
  p->polyprec_double.temp = NULL;
  p->polyprec_double.h_ritz = NULL;
  p->polyprec_double.lejas = NULL;
  p->polyprec_double.random_rhs = NULL;
  p->polyprec_double.xtmp = NULL;

  p->polyprec_double.eigslvr.vl = NULL;
  p->polyprec_double.eigslvr.vr = NULL;
  p->polyprec_double.dirctslvr.ipiv = NULL;
  p->polyprec_double.dirctslvr.x = NULL;
  p->polyprec_double.dirctslvr.b = NULL;
#endif

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  p->Va = NULL;
  p->Za = NULL;
#endif

#ifdef BLOCK_JACOBI
  p->block_jacobi_double.b_backup = NULL;
  local_fgmres_double_struct_init( &(p->block_jacobi_double.local_p) );
#endif

#ifdef MUMPS_ADDS
    p->mumps_vals = NULL;
    p->mumps_Is = NULL;
    p->mumps_Js = NULL;

    p->mumps_rhs_loc = NULL;
    p->mumps_irhs_loc = NULL;
    p->mumps_SOL = NULL;
#endif
}


void fgmres_double_struct_alloc( int m, int n, long int vl, double tol, const int type, const int prec_kind,
                                    void (*precond)(), void (*eval_op)(), gmres_double_struct *p, level_struct *l ) {

/*********************************************************************************
* Allocates memory for the fgmres struct and sets its values.                  
* int m: Restart length                                                            
* int n: Number of restarts                                                        
* int vl: System size                                                              
* double tol: Tolerance for relative residual                                         
* const int type: Specifies the problem for which fgmres should be applied               
*                 (_GLOBAL_FGMRES, _K_CYCLE, _COARSE_GMRES)                              
* const int prec_kind: type of preconditioning: _RIGHT (flexible preconditioner),
*                                               _LEFT (stationary preconditioner)
*                                               or _NOTHING                                     
* void (*precond): Function pointer to the preconditioner                           
*********************************************************************************/  
  
  long int total=0; 
  int i, k=0;
  
  p->restart_length = m;
  p->num_restart = n;

  p->preconditioner = precond;

  p->eval_operator = eval_op; 
  p->tol = tol;
  p->kind = prec_kind;

#ifdef HAVE_TM1p1
  vl*=2;
#endif
  
  if(m > 0) {
  total += (m+1)*m; // Hessenberg matrix
  MALLOC( p->H, complex_double*, m );
  
  total += (5+m)*vl; // x, r, b, w, V
  MALLOC( p->V, complex_double*, m+1 );

  if ( precond != NULL ) {
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    if ( l->level==0 && l->depth>0 ) {
      total += (m+2)*vl;
      k = m+2;
      MALLOC( p->Z, complex_double*, k );
    } else {
#endif
      if ( prec_kind == _RIGHT ) {
        total += (m+1)*vl; // Z
        k = m+1;
      } else {
        total += vl;
        k = 1;
      }
      MALLOC( p->Z, complex_double*, k );
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    }
#endif
  } else {
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    if ( l->level == 0 && l->depth > 0 ) {
      total += (m+2)*vl;
      k = m+2;
      MALLOC( p->Z, complex_double*, k );
    }
#else
    k = 0;
#endif
  }

#ifdef PERS_COMMS
  g.pers_comms_nrZs = k;
#endif

  total += 4*(m+1); // y, gamma, c, s
  
  p->H[0] = NULL; // allocate connected memory
  MALLOC( p->H[0], complex_double, total );
  
  p->total_storage = total;
  total = 0;
  
  // ordering: H, y, gamma, c, s, w, V, Z, x, r, b
  // H
  for ( i=1; i<m; i++ )
    p->H[i] = p->H[0] + i*(m+1);
  total += m*(m+1);
  
  // y
  p->y = p->H[0] + total; total += m+1;
  // gamma
  p->gamma = p->H[0] + total; total += m+1;
  // c
  p->c = p->H[0] + total; total += m+1;
  // s
  p->s = p->H[0] + total; total += m+1;
  // w
  p->w = p->H[0] + total; total += vl;
  // V
  for ( i=0; i<m+1; i++ ) {
    p->V[i] = p->H[0] + total; total += vl;
  }
  // Z
  for ( i=0; i<k; i++ ) {
    p->Z[i] = p->H[0] + total; total += vl;
  }

  // x
  p->x = p->H[0] + total; total += vl;
  // r
  p->r = p->H[0] + total; total += vl;
  // b
  p->b = p->H[0] + total; total += vl;
  
  ASSERT( p->total_storage == total );
  }
  
  if ( type == _GLOBAL_FGMRES ) {    
    p->timing = 1;
    p->print = g.vt.evaluation?0:1;
    p->initial_guess_zero = 1;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    p->op = &(g.op_double);
  } else if ( type == _K_CYCLE ) {
    // these settings also work for GMRES as a smoother
    p->timing = 0;
    p->print = 0;
    p->initial_guess_zero = 1;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    p->op = &(l->s_double.op);
  } else if ( type == _COARSE_GMRES ) {
    p->timing = 0;
    p->print = 0;
    p->initial_guess_zero = 1;
    p->layout = -1;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    if ( g.odd_even )
      p->op = &(l->oe_op_double);
    else  
      p->op = &(l->s_double.op);
  } else {
    ASSERT( type < 3 );
  }

#if defined(GCRODR) || defined(POLYPREC)
  if (l->level==0) {
#endif

  // FIXME : is this function-pointer-assignment really necessary ?
#if defined(GCRODR) || defined(POLYPREC)
  //p->polyprec_double.eigslvr.eigslvr_double = eigslvr_double;
#endif

  // copy of Hesselnberg matrix
#if defined(GCRODR) && defined(POLYPREC)
  MALLOC(p->gcrodr_double.eigslvr.Hc, complex_double*, m);
  p->polyprec_double.eigslvr.Hc = p->gcrodr_double.eigslvr.Hc;
  p->gcrodr_double.eigslvr.Hc[0] = NULL; // allocate connected memory
  MALLOC( p->gcrodr_double.eigslvr.Hc[0], complex_double, m*(m+1) );
  for ( i=1; i<m; i++ )
    p->gcrodr_double.eigslvr.Hc[i] = p->gcrodr_double.eigslvr.Hc[0] + i*(m+1);
  p->polyprec_double.eigslvr.Hc[0] = p->gcrodr_double.eigslvr.Hc[0];
#elif defined(GCRODR)
  MALLOC(p->gcrodr_double.eigslvr.Hc, complex_double*, m);
  p->gcrodr_double.eigslvr.Hc[0] = NULL; // allocate connected memory
  MALLOC( p->gcrodr_double.eigslvr.Hc[0], complex_double, m*(m+1) );
  for ( i=1; i<m; i++ )
    p->gcrodr_double.eigslvr.Hc[i] = p->gcrodr_double.eigslvr.Hc[0] + i*(m+1);
#elif defined(POLYPREC)
  MALLOC(p->polyprec_double.eigslvr.Hc, complex_double*, m);
  p->polyprec_double.eigslvr.Hc[0] = NULL; // allocate connected memory
  MALLOC( p->polyprec_double.eigslvr.Hc[0], complex_double, m*(m+1) );
  for ( i=1; i<m; i++ )
    p->polyprec_double.eigslvr.Hc[i] = p->polyprec_double.eigslvr.Hc[0] + i*(m+1);
#endif

#ifdef POLYPREC
  p->polyprec_double.d_poly = g.polyprec_d;
  int d_poly=p->polyprec_double.d_poly;

  MALLOC( p->polyprec_double.col_prods, complex_double, d_poly);
  MALLOC( p->polyprec_double.h_ritz, complex_double, d_poly);
  MALLOC( p->polyprec_double.lejas, complex_double, d_poly);
  MALLOC( p->polyprec_double.random_rhs, complex_double, vl );
  MALLOC( p->polyprec_double.accum_prod, complex_double, vl );
  MALLOC( p->polyprec_double.product, complex_double, vl );
  MALLOC( p->polyprec_double.temp, complex_double, vl );

  MALLOC( p->polyprec_double.xtmp, complex_double, vl );

  MALLOC( p->polyprec_double.Hcc, complex_double, d_poly*d_poly );
  MALLOC( p->polyprec_double.L, complex_double*, d_poly+ 1);

  p->polyprec_double.L[0] = NULL;

  MALLOC( p->polyprec_double.L[0], complex_double, (d_poly+1)*d_poly );

  for (i=1; i<d_poly+1; i++)
  {
    p->polyprec_double.L[i] = p->polyprec_double.L[0] + i*d_poly;
  }

  MALLOC( p->polyprec_double.dirctslvr.ipiv, int, d_poly);
  MALLOC( p->polyprec_double.dirctslvr.x, complex_double, d_poly);
  MALLOC( p->polyprec_double.dirctslvr.b, complex_double, d_poly);

  p->polyprec_double.dirctslvr.N = d_poly;
  p->polyprec_double.dirctslvr.lda = d_poly; // m here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  p->polyprec_double.dirctslvr.ldb = d_poly;
  p->polyprec_double.dirctslvr.nrhs = 1;
  p->polyprec_double.dirctslvr.Hcc = p->polyprec_double.Hcc;
  p->polyprec_double.dirctslvr.dirctslvr_double = dirctslvr_double;

  MALLOC( p->polyprec_double.eigslvr.vl, complex_double, d_poly*d_poly );
  MALLOC( p->polyprec_double.eigslvr.vr, complex_double, d_poly*d_poly );

  p->polyprec_double.eigslvr.jobvl = 'N';
  p->polyprec_double.eigslvr.jobvr = 'N';

  p->polyprec_double.eigslvr.N = d_poly;
  p->polyprec_double.eigslvr.lda = p->restart_length + 1;
  p->polyprec_double.eigslvr.ldvl = d_poly;
  p->polyprec_double.eigslvr.ldvr = d_poly;
  p->polyprec_double.eigslvr.w = p->polyprec_double.h_ritz;
  p->polyprec_double.Hc = p->polyprec_double.eigslvr.Hc;
  p->polyprec_double.eigslvr.eigslvr_double = eigslvr_double;    

  p->polyprec_double.update_lejas = 1;
  p->polyprec_double.preconditioner = NULL;
  p->polyprec_double.preconditioner_bare = p->preconditioner;
  p->polyprec_double.syst_size = vl;

  p->polyprec_double.eigslvr.A = p->polyprec_double.Hc[0];
#endif

#if defined(GCRODR) || defined(POLYPREC)
  }
#endif

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  p->syst_size = vl;
  MALLOC( p->Va, complex_double*, m+2 );
  MALLOC( p->Za, complex_double*, m+2 );
  p->Va[0] = NULL;
  p->Za[0] = NULL;
  MALLOC( p->Va[0], complex_double, (m+2)*vl );
  MALLOC( p->Za[0], complex_double, (m+2)*vl );

  for ( i=1; i<m+2; i++ )
  {
    p->Va[i] = p->Va[0] + i*vl;
    p->Za[i] = p->Za[0] + i*vl;
  }

#ifdef PERS_COMMS
  g.pers_comms_nrZas = m+2;
#endif
#endif

#ifdef BLOCK_JACOBI
  p->block_jacobi_double.syst_size = vl;

  if (l->level==0) {
    // these two always go together
    p->block_jacobi_double.BJ_usable = 0;
    p->block_jacobi_double.local_p.polyprec_double.update_lejas = 1;

    MALLOC( p->block_jacobi_double.b_backup, complex_double, vl );
    MALLOC( p->block_jacobi_double.xtmp, complex_double, vl );

    p->block_jacobi_double.local_p.polyprec_double.d_poly = g.local_polyprec_d;

     local_fgmres_double_struct_alloc( g.local_polyprec_d, 1, vl, g.coarse_tol, 
                                          _COARSE_GMRES, _NOTHING, NULL,
                                         coarse_local_apply_schur_complement_double,
#ifdef MUMPS_ADDS
                                         local_apply_coarse_operator_double,
#else
                                         g.odd_even?coarse_local_apply_schur_complement_double:local_apply_coarse_operator_double,
#endif
                                          &(p->block_jacobi_double.local_p), l );
   }
 #endif

#ifdef MUMPS_ADDS
  if (l->level==0) {

    int site_var = l->num_lattice_site_var;
    int nr_nodes = l->num_inner_lattice_sites;
    MALLOC( p->mumps_vals,complex_double,SQUARE(site_var)*nr_nodes * 9);
    MALLOC( p->mumps_Is,int,SQUARE(site_var)*nr_nodes * 9); // nr. of el per node * nr. of nodes * 9   //9 = self + T+ + T- + Z+ + Z- + Y+ ...
    MALLOC(p->mumps_Js,int,SQUARE(site_var)*nr_nodes * 9);
    memset(l->p_double.mumps_Is, 0, SQUARE(site_var)*nr_nodes * 9 * sizeof(int));
    memset(l->p_double.mumps_Js, 0, SQUARE(site_var)*nr_nodes * 9 * sizeof(int));
    memset(l->p_double.mumps_vals, 0, SQUARE(site_var)*nr_nodes * 9 * sizeof(complex_double));

    int mumps_n = site_var * nr_nodes * g.num_processes;       //order of Matrix
    int nnz = SQUARE(site_var) * nr_nodes *9 * g.num_processes;        //number of nonzero elements
    int nnz_loc = SQUARE(site_var) * nr_nodes *9;
//SOLUTION
    if (g.my_rank == 0){
      MALLOC(l->p_double.mumps_SOL, complex_double, mumps_n);
      memset(l->p_double.mumps_SOL, 0, mumps_n * sizeof(complex_double));
    }

//######### SET UP RHS #############
    int rhs_len = l->p_double.v_end-l->p_double.v_start;  //entire vector eta
    MALLOC(l->p_double.mumps_irhs_loc, int, rhs_len);
    MALLOC(l->p_double.mumps_rhs_loc, complex_double, rhs_len);
    memset(l->p_double.mumps_rhs_loc, 0, rhs_len * sizeof(complex_double));
    memset(l->p_double.mumps_irhs_loc, 0, rhs_len * sizeof(int));

       //confige MUMPS_struct
    g.mumps_id.job = JOB_INIT;
    g.mumps_id.par = 1;
    g.mumps_id.sym = 0;
    g.mumps_id.comm_fortran = USE_COMM_WORLD;
    zmumps_c(&(g.mumps_id));

    g.mumps_id.ICNTL(5) = 0;   //assembled matrix
    g.mumps_id.ICNTL(18) = 3;  //distributed local triplets for analysis and factorization
    g.mumps_id.ICNTL(20) = 10; //distributed RHS. compare to inctl(20) = 11
    g.mumps_id.ICNTL(35) = 2;  //BLR feature is activated during factorization and solution phase
//          mumps_id.ICNTL(35) = 3;   //BLR feature is activablrted during factorization, not used in solve
    g.mumps_id.cntl[6] = g.mumps_drop_tol;    //dropping parameter ε   (absolute error)        //original 7 but in c 6

//LHS
    g.mumps_id.n = mumps_n;    //needed at least on P0
    g.mumps_id.nnz_loc = nnz_loc;
    g.mumps_id.irn_loc = p->mumps_Is;
    g.mumps_id.jcn_loc = p->mumps_Js;
    g.mumps_id.a_loc = p->mumps_vals;

//RHS
    g.mumps_id.nloc_rhs = rhs_len;
    g.mumps_id.rhs_loc = p->mumps_rhs_loc;
    g.mumps_id.irhs_loc = p->mumps_irhs_loc;
    g.mumps_id.lrhs_loc = rhs_len; //leading dimension
    if (g.my_rank == 0){
      g.mumps_id.rhs = p->mumps_SOL;
    }

//outputs
    g.mumps_id.ICNTL(1) = 0;//6;       //error messages
    g.mumps_id.ICNTL(2) = 0;//1;       //diagnostic printing and statistics local to each MPI process
    g.mumps_id.ICNTL(3) = 0;//6;       //global information, collected on host (default 6)
    g.mumps_id.ICNTL(4) = 0;//2;       //level of printing for error, warning, and diagnostic messages (default 2)
  }
#endif

}


void fgmres_double_struct_free( gmres_double_struct *p, level_struct *l ) {

/*********************************************************************************
* Frees the allocated space for the gmres struct p.                            
*********************************************************************************/ 
  
  int k=0;

  int m = p->restart_length;
  if ( p->preconditioner != NULL ) {
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    if ( l->level==0 && l->depth>0 ) {
      k = m+2;
    } else {
#endif
      if ( p->kind == _RIGHT ) {
        k = m+1;
      } else {
        k = 1;
      }
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    }
#endif
  } else {
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    if ( l->level == 0 && l->depth > 0 ) {
      k = m+2;
    }
#else
    k = 0;
#endif
  }

  if(p->restart_length > 0) {
  FREE( p->H[0], complex_double, p->total_storage );
  FREE( p->H, complex_double*, p->restart_length );
  FREE( p->V, complex_double*, p->restart_length+1 );
  
  if ( p->Z != NULL )
    FREE( p->Z, complex_double*, k );
  }
  
  p->D = NULL;
  p->clover = NULL;

  // --- COARSEST-LEVEL IMPROVEMENTS

#if defined(GCRODR) || defined(POLYPREC)
  if (l->level==0) {
#endif

  // copy of Hesselnberg matrix
#if defined(GCRODR) && defined(POLYPREC)
  int m = p->restart_length;
  FREE( p->gcrodr_double.eigslvr.Hc[0], complex_double, m*(m+1) );
  FREE(p->gcrodr_double.eigslvr.Hc, complex_double*, m);
#elif defined(GCRODR)
  int m = p->restart_length;
  FREE( p->gcrodr_double.eigslvr.Hc[0], complex_double, m*(m+1) );
  FREE(p->gcrodr_double.eigslvr.Hc, complex_double*, m);
#elif defined(POLYPREC)
  int m = p->restart_length;
  FREE( p->polyprec_double.eigslvr.Hc[0], complex_double, m*(m+1) );
  FREE(p->polyprec_double.eigslvr.Hc, complex_double*, m);
#endif

#ifdef POLYPREC
  int d_poly = 10;
  int vl = p->polyprec_double.syst_size;
  FREE( p->polyprec_double.Hcc, complex_double, d_poly*d_poly );
  FREE( p->polyprec_double.L[0], complex_double, (d_poly+1)*d_poly );
  FREE( p->polyprec_double.L, complex_double*, d_poly+1 );
  FREE( p->polyprec_double.h_ritz,complex_double, d_poly );
  FREE( p->polyprec_double.lejas,complex_double, d_poly );
  FREE( p->polyprec_double.accum_prod, complex_double, vl );
  FREE( p->polyprec_double.product, complex_double, vl );    
  FREE( p->polyprec_double.temp, complex_double, vl );    
  FREE( p->polyprec_double.xtmp, complex_double, vl );
  FREE( p->polyprec_double.random_rhs, complex_double, vl );
  FREE( p->polyprec_double.col_prods, complex_double, d_poly );

  FREE( p->polyprec_double.eigslvr.vl,complex_double, d_poly*d_poly );
  FREE( p->polyprec_double.eigslvr.vr,complex_double, d_poly*d_poly );  

  FREE( p->polyprec_double.dirctslvr.ipiv, int, d_poly );
  FREE( p->polyprec_double.dirctslvr.x, complex_double, d_poly );
  FREE( p->polyprec_double.dirctslvr.b, complex_double, d_poly );
#endif

#if defined(GCRODR) || defined(POLYPREC)
  }
#endif

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    FREE( p->Va[0], complex_double*, (p->restart_length+2)*p->syst_size );
    FREE( p->Za[0], complex_double*, (p->restart_length+2)*p->syst_size );
    FREE( p->Va, complex_double, p->restart_length+2 );
    FREE( p->Za, complex_double, p->restart_length+2 );
#endif

#ifdef BLOCK_JACOBI
  if (l->level==0) {
    FREE( p->block_jacobi_double.b_backup, complex_double, p->block_jacobi_double.syst_size );
    FREE( p->block_jacobi_double.xtmp, complex_double, p->block_jacobi_double.syst_size );

    local_fgmres_double_struct_free( &(p->block_jacobi_double.local_p), l );
  }
#endif

#ifdef MUMPS_ADDS
       //release MUMPS_struct
  g.mumps_id.job = JOB_END;
  zmumps_c(&(g.mumps_id));
  int site_var = l->num_lattice_site_var;
  int nr_nodes = l->num_inner_lattice_sites * g.num_processes;
  FREE( p->mumps_vals,complex_double,SQUARE(site_var)*nr_nodes );
  FREE( p->mumps_Is,int,SQUARE(site_var)*nr_nodes );
  FREE( p->mumps_Js,int,SQUARE(site_var)*nr_nodes );
  FREE( p->mumps_irhs_loc, int, l->p_double.v_end-l->p_double.v_start);
  FREE( p->mumps_rhs_loc, complex_double, l->p_double.v_end-l->p_double.v_start);
  FREE( p->mumps_SOL, complex_double, site_var * nr_nodes * g.num_processes);       //order of Matrix
#endif

}


int fgmres_double( gmres_double_struct *p, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Uses FGMRES to solve the system D x = b, where b is taken from p->b and x is 
* stored in p->x.                                                              
*********************************************************************************/  

  // start and end indices for vector functions depending on thread
  int start;
  int end;

  int j=-1, finish=0, iter=0, il, ol, res;
  complex_double gamma0 = 0;

  complex_double beta = 0;

  double norm_r0=1, gamma_jp1=1, t0=0, t1=0;
  START_LOCKED_MASTER(threading)

  if ( l->depth==0 && ( p->timing || p->print ) ) prof_init( l );

  if ( l->level==0 && g.num_levels > 1 && g.interpolation ) p->tol = g.coarse_tol;
  if ( l->depth > 0 ) p->timing = 1;
  if ( l->depth == 0 ) t0 = MPI_Wtime();
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
  if ( p->print && g.print > 0 ) printf0("+----------------------------------------------------------+\n");
#endif
  END_LOCKED_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  START_MASTER(threading)
//   if (l->level==0) printf0( "g ---> start=%d, end=%d, diff=%d, m0=%f, op=%p\n", p->v_start, p->v_end, p->v_end-p->v_start, p->op->m0, p->op );
  END_MASTER(threading)

  for( ol=0; ol<p->num_restart && finish==0; ol++ )  {
  
    if( ol == 0 && p->initial_guess_zero ) {
      res = _NO_RES;
      vector_double_copy( p->r, p->b, start, end, l );
    } else {
      res = _RES;
      if ( p->kind == _LEFT && p->preconditioner ) {
        apply_operator_double( p->Z[0], p->x, p, l, threading );
        if ( g.method == 5 ) {
          START_LOCKED_MASTER(threading)
          g.bicgstab_tol = (!g.mixed_precision)?p->tol:MAX( 1E-3, (p->tol/(gamma_jp1/norm_r0))*5E-1 );
          END_LOCKED_MASTER(threading)
        }
        p->preconditioner( p->w, NULL, p->Z[0], _NO_RES, l, threading );
      } else {
        apply_operator_double( p->w, p->x, p, l, threading ); // compute w = D*x
      }
      vector_double_minus( p->r, p->b, p->w, start, end, l ); // compute r = b - w
    }
    gamma0 = global_norm_double( p->r, p->v_start, p->v_end, l, threading ); // gamma_0 = norm(r)
    START_MASTER(threading)
    p->gamma[0] = gamma0;
    END_MASTER(threading);
    SYNC_MASTER_TO_ALL(threading);

    if( cabs(gamma0)/norm_r0 < p->tol ){ break; }

    if ( ol == 0 ) {
     if (l->depth == 0 && !p->initial_guess_zero) {
       norm_r0 = global_norm_double( p->b, p->v_start, p->v_end, l, threading );
       printf0("| initial guess relative residual:            %le |\n", creal(gamma0)/norm_r0);
     } else {
       norm_r0 = creal(p->gamma[0]);
     }
    }

    //complex_double actual_norm = p->gamma[0]/norm_r0;
    //printf0( "actual residual = %.15f \n",CSPLIT(actual_norm) );

    vector_double_real_scale( p->V[0], p->r, 1/p->gamma[0], start, end, l ); // v_0 = r / gamma_0
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    if ( l->level == 0 && l->depth > 0 ) {
      arnoldi_step_double( p->V, p->Z, p->w, p->H, p->y, 0, p->preconditioner, p, l, threading );
    }
#endif   
    
    for( il=0; il<p->restart_length && finish==0; il++) {

      j = il; iter++;
      if ( g.method == 5 ) {
        START_LOCKED_MASTER(threading)
        g.bicgstab_tol = (!g.mixed_precision)?p->tol:MAX( 1E-3, (p->tol/(gamma_jp1/norm_r0))*5E-1 );
        END_LOCKED_MASTER(threading)
      }

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

      if ( cabs( p->H[j][j+1] ) > p->tol/10 ) {
        qr_update_double( p->H, p->s, p->c, p->gamma, j, l, threading );
        gamma_jp1 = cabs( p->gamma[j+1] );
        
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
        if ( iter%10 == 0 || p->preconditioner != NULL || l->depth > 0 ) {
          START_MASTER(threading)
          if ( p->print && g.print > 0 )
            printf0("| approx. rel. res. after  %-6d iterations: %e |\n", iter, gamma_jp1/norm_r0 );
          END_MASTER(threading)
        }
#endif

        //if ( l->level==0 ) printf0("g rel residual (gmres) = %.16f\n", gamma_jp1/norm_r0);

        if( gamma_jp1/norm_r0 < p->tol || gamma_jp1/norm_r0 > 1E+5 ) { // if satisfied ... stop

          finish = 1;
          START_MASTER(threading)
          if ( gamma_jp1/norm_r0 > 1E+5 ) printf0("Divergence of fgmres_double, iter = %d, level=%d\n", iter, l->level );
          END_MASTER(threading)

/*
#ifdef BLOCK_JACOBI
          if ( l->level==0 )
          {
            // backup of p->x, just in case tol hasn't been reached we need to restore ...
            vector_double_copy( p->block_jacobi_double.xtmp, p->x, start, end, l );

            compute_solution_double( p->x, (p->preconditioner&&p->kind==_RIGHT)?p->Z:p->V,
                                        p->y, p->gamma, p->H, j, (res==_NO_RES)?ol:1, p, l, threading );

            p->eval_operator( p->w, p->x, p->op, l, threading );
            vector_double_minus( p->r, p->block_jacobi_double.b_backup, p->w, start, end, l ); // compute r = b - w
            double norm_r0xx = global_norm_double( p->block_jacobi_double.b_backup, start, end, l, threading );
            double betaxx = global_norm_double( p->r, start, end, l, threading );
            if ( betaxx/norm_r0xx < p->tol ) {
              finish = 1;
            } else {
              // restore p->x
              vector_double_copy( p->x, p->block_jacobi_double.xtmp, start, end, l );
            }
            START_MASTER(threading)
            if ( betaxx/norm_r0xx > 1E+5 ) printf0("Divergence of fgmres_double, iter = %d, level=%d\n", iter, l->level );
            END_MASTER(threading)
          } else {
            finish = 1;
            START_MASTER(threading)
            if ( gamma_jp1/norm_r0 > 1E+5 ) printf0("Divergence of fgmres_double, iter = %d, level=%d\n", iter, l->level );
            END_MASTER(threading)
          }
#else
          finish = 1;
          START_MASTER(threading)
          if ( gamma_jp1/norm_r0 > 1E+5 ) printf0("Divergence of fgmres_double, iter = %d, level=%d\n", iter, l->level );
          END_MASTER(threading)
#endif
*/
        }
      } else {
        START_MASTER(threading)
//         printf0("from fgmres : depth: %d, iter: %d, p->H(%d,%d) = %+lf+%lfi\n", l->depth, iter, j+1, j, CSPLIT( p->H[j][j+1] ) );
        END_MASTER(threading)
        finish = 1;
        break;
      }
    } // end of a single restart
#ifdef BLOCK_JACOBI
    if ( l->level==0 ) {
      if ( finish==0 ) {
        compute_solution_double( p->x, (p->preconditioner&&p->kind==_RIGHT)?p->Z:p->V,
                                    p->y, p->gamma, p->H, j, (res==_NO_RES)?ol:1, p, l, threading );
      }
    } else {
      compute_solution_double( p->x, (p->preconditioner&&p->kind==_RIGHT)?p->Z:p->V,
                                  p->y, p->gamma, p->H, j, (res==_NO_RES)?ol:1, p, l, threading );
    }
#else
    compute_solution_double( p->x, (p->preconditioner&&p->kind==_RIGHT)?p->Z:p->V,
                                p->y, p->gamma, p->H, j, (res==_NO_RES)?ol:1, p, l, threading );
#endif

  } // end of fgmres

  START_LOCKED_MASTER(threading)
  if ( l->depth == 0 ) { t1 = MPI_Wtime(); g.total_time = t1-t0; g.iter_count = iter; g.norm_res = gamma_jp1/norm_r0; }
  END_LOCKED_MASTER(threading)
  
  if ( p->print ) {
#ifdef FGMRES_RESTEST
    apply_operator_double( p->w, p->x, p, l, threading );
    vector_double_minus( p->r, p->b, p->w, start, end, l );
    beta = global_norm_double( p->r, p->v_start, p->v_end, l, threading );
#else
    beta = gamma_jp1;
#endif
    START_MASTER(threading)
    g.norm_res = creal(beta)/norm_r0;
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    if ( g.print > 0 ) printf0("+----------------------------------------------------------+\n\n");
#endif
    printf0("+----------------------------------------------------------+\n");
    printf0("|       FGMRES iterations: %-6d coarse average: %-6.2lf   |\n", iter,
            ((double)g.coarse_iter_count)/((double)iter) );
    printf0("| exact relative residual: ||r||/||b|| = %e      |\n", creal(beta)/norm_r0 );
    printf0("| elapsed wall clock time: %-8.4lf seconds                |\n", t1-t0 );
    if ( g.coarse_time > 0 ) 
      printf0("|        coarse grid time: %-8.4lf seconds (%04.1lf%%)        |\n",
              g.coarse_time, 100*(g.coarse_time/(t1-t0)) );
    printf0("|        coarsest grid time: %-8.4lf seconds (%04.1lf%%)        |\n",
              g.coarsest_time, 100*(g.coarsest_time/(t1-t0)) );
    printf0("|  consumed core minutes*: %-8.2le (solve only)           |\n", ((t1-t0)*g.num_processes*MAX(1,threading->n_core))/60.0 );
    printf0("|    max used mem/MPIproc: %-8.2le GB                     |\n", g.max_storage/1024.0 );
    printf0("+----------------------------------------------------------+\n");
    printf0("*: only correct if #MPIprocs*#threads == #CPUs\n\n");
    END_MASTER(threading)
  }
  
#ifdef COARSE_RES
  if ( l->depth > 0 ) {
    START_MASTER(threading)
    char number[3]; sprintf( number, "%2d", 31+l->depth ); printf0("\033[1;%2sm|", number );
    printf0(" - depth: %d, gmres iter: %2d, approx rel res: %le |", l->depth, iter, gamma_jp1/norm_r0 );
    printf0("\033[0m\n"); fflush(0);
    END_MASTER(threading)
  }
#endif

  if ( l->level == 0 ) {
    START_LOCKED_MASTER(threading)
    g.coarse_iter_count += iter;
    END_LOCKED_MASTER(threading)
  }
    
  if ( l->depth == 0 && g.vt.p_end != NULL  ) {
    if ( g.vt.p_end != NULL ) {
      START_LOCKED_MASTER(threading)
      printf0("solve iter: %d\n", iter );
      printf0("solve time: %le seconds\n", t1-t0 );
      g.vt.p_end->values[_SLV_TIME] += (t1-t0)/((double)g.vt.average_over);
      g.vt.p_end->values[_SLV_ITER] += iter/((double)g.vt.average_over);
      g.vt.p_end->values[_CRS_ITER] += (((double)g.coarse_iter_count)/((double)iter))/((double)g.vt.average_over);
      g.vt.p_end->values[_CRS_TIME] += g.coarse_time/((double)g.vt.average_over);
    END_LOCKED_MASTER(threading)
    }
  }
  if ( l->depth == 0 && ( p->timing || p->print ) && !(g.vt.p_end != NULL )  ) {
    START_MASTER(threading)
    if ( g.method != 6 ) prof_print( l );
    END_MASTER(threading)
  }
  
  return iter;
}


void bicgstab_double( gmres_double_struct *ps, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Uses BiCGstab to solve the system D x = b, where b is taken from ps->b and x is 
* stored in ps->x.                                                              
*********************************************************************************/
  
  vector_double x, b, r, r_tilde, p, pp, v, s, t; // Krylov subspace size: 5
  complex_double alpha=1, beta=1, rho=1, rho_old=1, omega=1;
  int iter=0, maxiter;
  double tol, b_norm, r_norm, s_norm;
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(ps->v_start, ps->v_end, &start, &end, l, threading);
  
  tol = (l->level==0 && g.num_levels > 1 && g.interpolation )?g.coarse_tol:g.bicgstab_tol;
  maxiter = 1000000; r = ps->r; b = ps->b; x = ps->x; p = ps->w;
  pp = ps->V[0]; r_tilde = ps->V[1]; v = ps->V[2]; s = ps->V[3]; t = ps->V[4];
  
  vector_double_copy( r, b, start, end, l );
  vector_double_copy( r_tilde, b, start, end, l );
  vector_double_define( x, 0, start, end, l );
  vector_double_define( v, 0, start, end, l );
  vector_double_define( s, 0, start, end, l );
  vector_double_define( t, 0, start, end, l );
  b_norm = global_norm_double( b, ps->v_start, ps->v_end, l, threading );

  r_norm = b_norm;
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)  
  START_MASTER(threading)
  printf0("+----------------------------------------------------------+\n");
  END_MASTER(threading)
#endif
  while ( r_norm/b_norm > tol && iter < maxiter ) {
    iter++;
    
    rho_old = rho;
    rho = global_inner_product_double( r_tilde, r, ps->v_start, ps->v_end, l, threading );

    if ( rho == 0 ) {
      START_MASTER(threading)
      printf0("rho = 0: BiCGstab did not converge.\n");
      END_MASTER(threading)
      break;
    }
    
    if ( iter == 1 ) {
      vector_double_copy( p, r, start, end, l );
    } else {
      beta = (rho/rho_old)*(alpha/omega);
      vector_double_saxpy( pp, p,  v, -omega, start, end, l );
      vector_double_saxpy( p,  r, pp,   beta, start, end, l );
    }    
    apply_operator_double( v, p, ps, l, threading );
    alpha = rho / global_inner_product_double( r_tilde, v, ps->v_start, ps->v_end, l, threading );
    vector_double_saxpy( s, r, v, -alpha, start, end, l );
    s_norm = global_norm_double( s, ps->v_start, ps->v_end, l, threading );

    if ( s_norm/b_norm < tol ) {
      vector_double_saxpy( x, x, p, alpha, start, end, l );
      break;
    }
    
    apply_operator_double( t, s, ps, l, threading );
    omega = global_inner_product_double( t, s, ps->v_start, ps->v_end, l, threading )
          / global_inner_product_double( t, t, ps->v_start, ps->v_end, l, threading );
    
    vector_double_saxpy( x, x, p,  alpha, start, end, l );
    vector_double_saxpy( x, x, s,  omega, start, end, l );
    vector_double_saxpy( r, s, t, -omega, start, end, l );

    r_norm = global_norm_double( r, ps->v_start, ps->v_end, l, threading );

#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    START_MASTER(threading)
    if ( iter % 100 == 0 ) printf0("| biCGstab relres: %12.6le,  iterations: %-8d     |\n", r_norm/b_norm, iter );
    END_MASTER(threading)
#endif
  }
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
  START_MASTER(threading)
  printf0("| biCGstab relres: %12.6le,  iterations: %-8d     |\n", r_norm/b_norm, iter );
  printf0("+----------------------------------------------------------+\n");
  END_MASTER(threading)
#endif
}


void cgn_double( gmres_double_struct *ps, level_struct *l, struct Thread *threading ) {
  
/*********************************************************************************
* Uses CGN to solve the system D x = b, where b is taken from ps->b and x is 
* stored in ps->x.                                                              
*********************************************************************************/

  vector_double r_old, r_new, r_true, p, pp, Dp, x, b;
  complex_double alpha, beta=0, gamma;
  int maxiter, iter=0;
  double tol, r0_norm, r_norm, prod_rr_old, t0=0, t1=0;
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  
  b = ps->b; x = ps->x;
  r_old = ps->V[2]; r_new = ps->V[3]; r_true = ps->r;
  p = ps->w; pp = ps->V[0]; Dp = ps->V[1];
  tol = (l->level==0 && g.num_levels > 1 && g.interpolation )?g.coarse_tol:ps->tol;
  maxiter = ps->num_restart;
  
  START_MASTER(threading)
  if ( ps->timing || ps->print ) t0 = MPI_Wtime();
  END_MASTER(threading)

  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(ps->v_start, ps->v_end, &start, &end, l, threading);

  vector_double_define( x, 0, start, end, l );
  apply_operator_double( Dp, x, ps, l, threading );
  vector_double_minus( pp, b, Dp, start, end, l );
  apply_operator_dagger_double( r_old, pp, ps, l, threading );
  
  vector_double_copy( p, r_old, start, end, l );
  r0_norm = global_norm_double( r_old, ps->v_start, ps->v_end, l, threading );
  //  prod_rr_old = global_inner_product_double( r_old, r_old, ps->v_start, ps->v_end, l, threading );
  prod_rr_old = r0_norm*r0_norm;

#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
  if ( ps->print ) {
    START_MASTER(threading)
    printf0("\n+----------------------------------------------------------+\n");
    END_MASTER(threading)
  }  
#endif
  while ( sqrt(prod_rr_old) / r0_norm > tol && iter < maxiter ) {
    iter++;
    
    apply_operator_double( pp, p, ps, l, threading );
    apply_operator_dagger_double( Dp, pp, ps, l, threading );
    
    gamma = global_inner_product_double( p, Dp, ps->v_start, ps->v_end, l, threading );
    alpha = prod_rr_old / gamma;
    vector_double_saxpy( x, x, p, alpha, start, end, l );
    vector_double_saxpy( r_new, r_old, Dp, -alpha, start, end, l );
    
    gamma = global_inner_product_double( r_new, r_new, ps->v_start, ps->v_end, l, threading );
    beta = gamma / prod_rr_old;
    
    vector_double_saxpy( p, r_new, p, beta, start, end, l );
    vector_double_copy( r_old, r_new, start, end, l );
    prod_rr_old = gamma;
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)    
    if ( iter%100 == 0 && ps->print >=1 ) {
      START_MASTER(threading)
      printf0("|      NE rel. res. after  %-6d iterations: %e |\n", iter, sqrt(prod_rr_old)/r0_norm );
      END_MASTER(threading)
    }
#endif
  }
  
  r0_norm = global_norm_double( b, ps->v_start, ps->v_end, l, threading );
  apply_operator_double( Dp, x, ps, l, threading );
  vector_double_minus( r_true, b, Dp, start, end, l );
  r_norm = global_norm_double( r_true, ps->v_start, ps->v_end, l, threading );

#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)  
  if ( ps->print ) {
    START_MASTER(threading)
    printf0("+----------------------------------------------------------+\n");
    printf0("| switching to CGNR, iter  %-6d true r res: %e |\n", iter, r_norm/r0_norm );
    printf0("+----------------------------------------------------------+\n");
    END_MASTER(threading)
  }
#endif
  
  while ( r_norm / r0_norm > tol && iter < maxiter ) {
    iter++;
    
    apply_operator_double( pp, p, ps, l, threading );
    apply_operator_dagger_double( Dp, pp, ps, l, threading );
    
    gamma = global_inner_product_double( p, Dp, ps->v_start, ps->v_end, l, threading );
    alpha = prod_rr_old / gamma;
    vector_double_saxpy( x, x, p, alpha, start, end, l );
    vector_double_saxpy( r_new, r_old, Dp, -alpha, start, end, l );
    
    // residual update
    vector_double_saxpy( r_true, r_true, pp, -alpha, start, end, l );
    r_norm = global_norm_double( r_true, ps->v_start, ps->v_end, l, threading );
    gamma = global_inner_product_double( r_new, r_new, ps->v_start, ps->v_end, l, threading );
    beta = gamma / prod_rr_old;
    
    vector_double_saxpy( p, r_new, p, beta, start, end, l );
    vector_double_copy( r_old, r_new, start, end, l );
    prod_rr_old = gamma;
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)    
    if ( iter%100 ==  0 && ps->print >=1 ) {
      START_MASTER(threading)
      printf0("|         rel. res. after  %-6d iterations: %e |\n", iter, r_norm/r0_norm );
      END_MASTER(threading)
    }
#endif
  }
  
  if ( ps->timing || ps->print ) t1 = MPI_Wtime();
  if ( ps->print ) {
    START_MASTER(threading)
    printf0("+----------------------------------------------------------+\n");
    printf0("|          CGN iterations: %-6d                          |\n", iter );
    END_MASTER(threading)
    apply_operator_double( Dp, x, ps, l, threading );
    vector_double_minus( pp, b, Dp, start, end, l );

    beta = global_norm_double( pp, ps->v_start, ps->v_end, l, threading );
    START_MASTER(threading)
    if ( ps->timing ) printf0("| exact relative residual: ||r||/||b|| = %e      |\n", creal(beta/r0_norm) );
    printf0("| elapsed wall clock time: %-12g seconds            |\n", t1-t0 );
    printf0("|  consumed core minutes*: %-8.2le (solve only)           |\n", ((t1-t0)*g.num_processes*MAX(1,threading->n_core))/60.0 );
    printf0("|    max used mem/MPIproc: %-8.2le GB                     |\n", g.max_storage/1024.0 );
    printf0("+----------------------------------------------------------+\n");
    printf0("*: only correct if #MPIprocs*#threads == #CPUs\n\n");
    END_MASTER(threading)
  }
  
  START_LOCKED_MASTER(threading)
  if ( l->level == 0 )
    g.coarse_iter_count += iter;

  if ( l->depth == 0 && g.vt.p_end != NULL  ) {
    if ( g.vt.p_end != NULL ) {
      g.vt.p_end->values[_SLV_TIME] += (t1-t0)/((double)g.vt.average_over);
      g.vt.p_end->values[_SLV_ITER] += (iter)/((double)g.vt.average_over);
    }
  }
  END_LOCKED_MASTER(threading)
}


int arnoldi_step_double( vector_double *V, vector_double *Z, vector_double w,
                            complex_double **H, complex_double* buffer, int j, void (*prec)(),
                            gmres_double_struct *p, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Extends the Arnoldi basis by one vector.
* - vector_double *V: Contains the Arnoldi basis vectors.
* - vector_double *Z: If a right precond. P is used, contains P*V[j] for all j.
* - vector_double w: Will be appended to existing Arnoldi basis at 
*   position j+1.
* - complex_double **H: Contains full Hessenberg matrix from the Arnoldi 
*   decomposition (columnmajor!)
* - complex_double* buffer: Buffer for local inner products.
* - int j: index of the new Arnoldi vector to be orthonormalized
*   against all previous ones.
* - void (*prec)(): Function pointer to preconditioner (can be NULL if no 
*   preconditioning is used).
*********************************************************************************/

// TODO :: add restriction to have always both SINGLE_ALLREDUCE_ARNOLDI and PIPELINED_ARNOLDI together
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  if ( l->level == 0 && l->depth > 0 ) {


    //if ( l->level==0 ) printf0("AHA1!\n");


    if (prec == NULL)
    {
#ifdef GCRODR
      int k = p->gcrodr_double.k;
      vector_double *Ck = p->gcrodr_double.C;
      complex_double **B = p->gcrodr_double.ort_B;
      //complex_double *bf = p->gcrodr_double.Bbuff[0];
      vector_double *DPCk = p->gcrodr_double.DPC;
#endif

      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)
      MPI_Request req;
      MPI_Status stat;
      int start, end, i;
      //const complex_double sigma = 0;
      compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

      if ( j == 0 ){
#ifdef PERS_COMMS
        g.pers_comms_id2 = 0;
        g.use_pers_comms1 = 1;
#endif
        apply_operator_double( Z[0], V[0], p, l, threading );
#ifdef PERS_COMMS
        g.pers_comms_id2 = -1;
        g.use_pers_comms1 = 0;
#endif

#ifdef GCRODR
        if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 ) {
          if (p->gcrodr_double.recompute_DPCk_plain == 1) {
            for( i=0; i<p->gcrodr_double.k; i++ ) {

#ifdef PERS_COMMS
              // g.pers_comms_nrZxs should be zero here anyways
              //g.pers_comms_id2 = p->restart_length + g.pers_comms_nrZxs + i;
              //g.use_pers_comms1 = 1;
#endif
              apply_operator_double( DPCk[i], Ck[i], p, l, threading );
#ifdef PERS_COMMS
              //g.pers_comms_id2 = -1;
              //g.use_pers_comms1 = 0;
#endif
            }
            p->gcrodr_double.recompute_DPCk_plain = 0;
          }
        }
#endif

        // TODO : re-enable !
        //if ( sigma ) vector_double_saxpy( Z[j+1], Z[j+1], Z[j], -sigma, start, end, l );
        return 1;
      } else {
        vector_double_copy( V[j], Z[j-1], start, end, l );
      }

#ifdef GCRODR
      // orthogonalize against Ck whenever necessary
      if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 ) {
        SYNC_MASTER_TO_ALL(threading)
        SYNC_CORES(threading)

        // buffer
        complex_double *bf1 = p->gcrodr_double.Bbuff[0];
        complex_double *bf2 = p->gcrodr_double.Bbuff[0]+(k+j+1);

        // space to orthogonalize against
        vector_double ort_sp[k+j+1];
        for (i=0; i<k; i++) ort_sp[i] = Ck[i];
        for (i=0; i<(j+1); i++) ort_sp[k+i] = V[i];

        complex_double tmpx[k+j+1];
        process_multi_inner_product_double( k+j+1, tmpx, ort_sp, V[j], p->v_start, p->v_end, l, threading );
        START_MASTER(threading)
        // buffer is of length m, and k<m
        for ( i=0; i<(k+j+1); i++ )
          bf2[i] = tmpx[i];
        if ( g.num_processes > 1 ) {
          PROF_double_START( _ALLR );
          MPI_Iallreduce( bf2, bf1, k+j+1, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm, &req );
          PROF_double_STOP( _ALLR, 1 );
        } else {
          for( i=0; i<(k+j+1); i++ )
            bf1[i] = bf2[i];
        }

        //memcpy( B[j-1], bf1, sizeof(complex_double)*k );
        //memcpy( H[j-1], bf1+k, sizeof(complex_double)*(j+1) );
        END_MASTER(threading)
      } else {
        complex_double tmp[j+1];
        process_multi_inner_product_double( j+1, tmp, V, V[j], p->v_start, p->v_end, l, threading );
        START_MASTER(threading)
        PROF_double_START( _ALLR );
        for( i=0; i<=j; i++ )
          buffer[i] = tmp[i];
        if ( g.num_processes > 1 ) {
          MPI_Iallreduce( buffer, H[j-1], j+1, MPI_COMPLEX_double, MPI_SUM,
                          (l->depth==0)?g.comm_cart:l->gs_double.level_comm, &req );
        } else {
          for( i=0; i<=j; i++ )
            H[j-1][i] = buffer[i];
        }
        PROF_double_STOP( _ALLR, 1 );
        END_MASTER(threading)
      }
#else
      complex_double tmp[j+1];
      process_multi_inner_product_double( j+1, tmp, V, V[j], p->v_start, p->v_end, l, threading );
      START_MASTER(threading)
      PROF_double_START( _ALLR );
      for( i=0; i<=j; i++ )
        buffer[i] = tmp[i];
      if ( g.num_processes > 1 ) {
        MPI_Iallreduce( buffer, H[j-1], j+1, MPI_COMPLEX_double, MPI_SUM,
                        (l->depth==0)?g.comm_cart:l->gs_double.level_comm, &req );
      } else {
        for( i=0; i<=j; i++ )
          H[j-1][i] = buffer[i];
      }
      PROF_double_STOP( _ALLR, 1 );
      END_MASTER(threading)
#endif

#ifdef PERS_COMMS
      g.pers_comms_id2 = j;
      g.use_pers_comms1 = 1;
#endif
      apply_operator_double( Z[j], V[j], p, l, threading );
#ifdef PERS_COMMS
      g.pers_comms_id2 = -1;
      g.use_pers_comms1 = 0;
#endif

      START_MASTER(threading)
      PROF_double_START( _ALLR );
      if ( g.num_processes > 1 ) {
        MPI_Wait( &req, &stat );
      }
      PROF_double_STOP( _ALLR, 0 );
      END_MASTER(threading)

#ifdef GCRODR
      if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 ) {
        SYNC_MASTER_TO_ALL(threading)
        SYNC_CORES(threading)

        complex_double *bf1 = p->gcrodr_double.Bbuff[0];
        START_MASTER(threading)
        // copy the B and H coefficients to the corresponding matrix
        memcpy( B[j-1], bf1, sizeof(complex_double)*k );
        memcpy( H[j-1], bf1+k, sizeof(complex_double)*(j+1) );
        END_MASTER(threading)
      }
#endif

      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)

#ifdef GCRODR
      if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 ) {
        for( i=0; i<k; i++ )
          vector_double_saxpy( V[j], V[j], Ck[i], -B[j-1][i], start, end, l );
      }
#endif

      for( i=0; i<j; i++ )
        vector_double_saxpy( V[j], V[j], V[i], -H[j-1][i], start, end, l );


//=============================================================
       START_MASTER(threading)
       complex_double tmp = H[j-1][j];
       for ( i=0; i<=j-1; i++ )
         tmp -= conj_double( H[j-1][i] )*H[j-1][i];
 #ifdef GCRODR
       if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 )
       {      
         for( i=0; i<k; i++ )
           tmp -= conj_double( B[j-1][i] )*B[j-1][i];
       }
 #endif
//       tmp = sqrt( creal_double( tmp ) );
      
      H[j-1][j] = tmp;
      ((complex_double*)threading->workspace)[0] = creal_double(tmp);
      END_MASTER(threading)
      
      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)
      if ( (creal_double(((complex_double*)threading->workspace)[0]) < 0.) || (sqrt(creal_double(H[j-1][j])) <= p->tol/10) )
      //if ( (creal_double(((complex_double*)threading->workspace)[0]) < 0.) )
      {
          H[j-1][j] = global_norm_double( V[j], p->v_start, p->v_end, l, threading );
          START_MASTER(threading)
          printf0( "FULL NORM !!\n" );
          END_MASTER(threading)
      }
      else
      {
          START_MASTER(threading)
          H[j-1][j] = sqrt(creal_double(H[j-1][j]));
          END_MASTER(threading)
          SYNC_MASTER_TO_ALL(threading)
          SYNC_CORES(threading)
      }
//=============================================================
      //H[j-1][j] = global_norm_double( V[j], p->v_start, p->v_end, l, threading );

//printf("H[j-1][j] = %f + i%f, tmp = %f + i%f\n", creal(H[j-1][j]), cimag(H[j-1][j]), creal(tmp), cimag(tmp));

      vector_double_real_scale( V[j], V[j], 1/H[j-1][j], start, end, l );

      //START_MASTER(threading)
      //if ( j > 0 ) {
      //  H[j-1][j-1] += sigma;
      //}
      //END_MASTER(threading)

      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)

      if ( j == 0 ) {
        //if ( sigma ) vector_double_saxpy( Z[j+1], Z[j+1], Z[j], -sigma, start, end, l );
      } else {
#ifdef GCRODR
        if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 ) {
          for( i=0; i<k; i++ )
            vector_double_saxpy( Z[j], Z[j], DPCk[i], -B[j-1][i], start, end, l );
        }
#endif
        for( i=0; i<j; i++ )
          vector_double_saxpy( Z[j], Z[j], Z[i], -H[j-1][i], start, end, l );
      }
      vector_double_real_scale( Z[j], Z[j], 1/H[j-1][j], start, end, l );
    }
    else
    {

      // ------------------------------------------------------------------------------------
      vector_double *Va = p->Va;
      vector_double *Za = p->Za;

      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)
      MPI_Request req;
      MPI_Status stat;
      int start, end, i;

      const complex_double sigma = 0;
      compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

#ifdef GCRODR
      int k = p->gcrodr_double.k;
      vector_double *Ck = p->gcrodr_double.C;
      complex_double **B = p->gcrodr_double.ort_B;
      //complex_double *bf = p->gcrodr_double.Bbuff[0];
      vector_double *PCk = p->gcrodr_double.PC;
      vector_double *DPCk = p->gcrodr_double.DPC;
#endif

      if ( j == 0 ) {
        if (prec == NULL) vector_double_copy( Z[0], V[0], start, end, l );
        else prec( Z[0], NULL, V[0], _NO_RES, l, threading );
#ifdef PERS_COMMS
        //g.pers_comms_id2 = p->restart_length + 0;
        //g.use_pers_comms1 = 1;
#endif
        apply_operator_double( Va[0], Z[0], p, l, threading );
#ifdef PERS_COMMS
        //g.pers_comms_id2 = -1;
        //g.use_pers_comms1 = 0;
#endif
        if ( sigma ) vector_double_saxpy( Va[j], Va[j], Va[j], -sigma, start, end, l );

#ifdef GCRODR
        if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 ) {
          if (p->gcrodr_double.recompute_DPCk_poly == 1) {
            for( i=0; i<p->gcrodr_double.k; i++ ) {
              if (prec == NULL) vector_double_copy( PCk[i], Ck[i], start, end, l );
              else prec( PCk[i], NULL, Ck[i], _NO_RES, l, threading );
#ifdef PERS_COMMS
              //g.pers_comms_id2 = p->restart_length + g.pers_comms_nrZxs + g.gcrodr_k + i;
              //g.use_pers_comms1 = 1;
#endif
              apply_operator_double( DPCk[i], PCk[i], p, l, threading );
#ifdef PERS_COMMS
              //g.pers_comms_id2 = -1;
              //g.use_pers_comms1 = 0;
#endif
            }
            p->gcrodr_double.recompute_DPCk_poly = 0;
          }
        }
#endif        
        return 1;
      }
      else
        vector_double_copy( V[j], Va[j-1], start, end, l );

#ifdef GCRODR
      // orthogonalize against Ck whenever necessary
      if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 ) {

        // buffer
        complex_double *bf1 = p->gcrodr_double.Bbuff[0];
        complex_double *bf2 = p->gcrodr_double.Bbuff[0]+(k+j+1);

        // space to orthogonalize against
        vector_double ort_sp[k+j+1];
        for (i=0; i<k; i++) ort_sp[i] = Ck[i];
        for (i=0; i<(j+1); i++) ort_sp[k+i] = V[i];

        complex_double tmpx[k+j+1];
        process_multi_inner_product_double( k+j+1, tmpx, ort_sp, V[j], p->v_start, p->v_end, l, threading );
        START_MASTER(threading)
        // buffer is of length m, and k<m
        for ( i=0; i<(k+j+1); i++ )
          bf2[i] = tmpx[i];
        if ( g.num_processes > 1 ) {
          PROF_double_START( _ALLR );
          MPI_Iallreduce( bf2, bf1, k+j+1, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm, &req );
          PROF_double_STOP( _ALLR, 1 );
        } else {
          for( i=0; i<(k+j+1); i++ )
            bf1[i] = bf2[i];
        }

        //memcpy( B[j-1], bf1, sizeof(complex_double)*k );
        //memcpy( H[j-1], bf1+k, sizeof(complex_double)*(j+1) );
        END_MASTER(threading)

        SYNC_MASTER_TO_ALL(threading)
        SYNC_CORES(threading)
      } else {
        complex_double tmp[j+1];
        process_multi_inner_product_double( j+1, tmp, V, V[j], p->v_start, p->v_end, l, threading );

        START_MASTER(threading)
        PROF_double_START( _ALLR );
        for( i=0; i<=j; i++ )
          buffer[i] = tmp[i];
        if ( g.num_processes > 1 ) {
          MPI_Iallreduce( buffer, H[j-1], j+1, MPI_COMPLEX_double, MPI_SUM,
                          (l->depth==0)?g.comm_cart:l->gs_double.level_comm, &req );
        } else {
            for( i=0; i<=j; i++ )
                H[j-1][i] = buffer[i];
        }
        PROF_double_STOP( _ALLR, 1 );
        END_MASTER(threading)
      }
#else
      complex_double tmp[j+1];
      process_multi_inner_product_double( j+1, tmp, V, V[j], p->v_start, p->v_end, l, threading );

      START_MASTER(threading)
      PROF_double_START( _ALLR );
      for( i=0; i<=j; i++ )
        buffer[i] = tmp[i];
      if ( g.num_processes > 1 ) {
        MPI_Iallreduce( buffer, H[j-1], j+1, MPI_COMPLEX_double, MPI_SUM,
                        (l->depth==0)?g.comm_cart:l->gs_double.level_comm, &req );
      } else {
          for( i=0; i<=j; i++ )
              H[j-1][i] = buffer[i];
      }
      PROF_double_STOP( _ALLR, 1 );
      END_MASTER(threading)
#endif

      if (prec == NULL) vector_double_copy( Za[j-1], Va[j-1], start, end, l );
      else prec( Za[j-1], NULL, Va[j-1], _NO_RES, l, threading );

#ifdef PERS_COMMS
      g.pers_comms_id2 = p->restart_length + j-1;
      g.use_pers_comms1 = 1;
#endif
      apply_operator_double( Va[j], Za[j-1], p, l, threading );
#ifdef PERS_COMMS
      g.pers_comms_id2 = -1;
      g.use_pers_comms1 = 0;
#endif

      START_MASTER(threading)
      PROF_double_START( _ALLR );
      if ( g.num_processes > 1 ) {
        MPI_Wait( &req, &stat );
      }
      PROF_double_STOP( _ALLR, 0 );
      END_MASTER(threading)

#ifdef GCRODR
      if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 ) {
        SYNC_MASTER_TO_ALL(threading)
        SYNC_CORES(threading)

        complex_double *bf1 = p->gcrodr_double.Bbuff[0];
        START_MASTER(threading)
        // copy the B and H coefficients to the corresponding matrix
        memcpy( B[j-1], bf1, sizeof(complex_double)*k );
        memcpy( H[j-1], bf1+k, sizeof(complex_double)*(j+1) );
        END_MASTER(threading)
      }
#endif

      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)

#ifdef GCRODR
      if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 ) {
        for( i=0; i<k; i++ )
          vector_double_saxpy( V[j], V[j], Ck[i], -B[j-1][i], start, end, l );
      }
#endif
      for( i=0; i<=j-1; i++ )
        vector_double_saxpy( V[j], V[j], V[i], -H[j-1][i], start, end, l );
      
//=============================================================
       START_MASTER(threading)
       complex_double tmp = H[j-1][j];
       for ( i=0; i<=j-1; i++ )
         tmp -= conj_double( H[j-1][i] )*H[j-1][i];
 #ifdef GCRODR
       if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 )
       {      
         for( i=0; i<k; i++ )
           tmp -= conj_double( B[j-1][i] )*B[j-1][i];
       }
 #endif
//       tmp = sqrt( creal_double( tmp ) );
      
      H[j-1][j] = tmp;
      ((complex_double*)threading->workspace)[0] = creal_double(tmp);
      END_MASTER(threading)
      
      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)
      if ( (creal_double(((complex_double*)threading->workspace)[0]) < 0.) || (sqrt(creal_double(H[j-1][j])) <= p->tol/10))
      //if ( (creal_double(((complex_double*)threading->workspace)[0]) < 0.) )
      {
          H[j-1][j] = global_norm_double( V[j], p->v_start, p->v_end, l, threading );
          START_MASTER(threading)
          printf0( "FULL NORM !!\n" );
          END_MASTER(threading)
      }
      else
      {
          START_MASTER(threading)
          H[j-1][j] = sqrt(creal_double(H[j-1][j]));
          END_MASTER(threading)
          SYNC_MASTER_TO_ALL(threading)
          SYNC_CORES(threading)
      }
//=============================================================
//      double tmp2 = global_norm_double( V[j], p->v_start, p->v_end, l, threading );
//      START_MASTER(threading)
//      H[j-1][j] = tmp2;
//      END_MASTER(threading)

//      SYNC_MASTER_TO_ALL(threading)
//      SYNC_CORES(threading)

      vector_double_real_scale( V[j], V[j], 1/H[j-1][j], start, end, l );

      vector_double_copy( Z[j], Za[j-1], start, end, l );
#ifdef GCRODR
      if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 )
      {      
        for( i=0; i<k; i++ )
          vector_double_saxpy( Z[j], Z[j], PCk[i], -B[j-1][i], start, end, l );
      }
#endif

      for( i=0; i<=j-1; i++ )
        vector_double_saxpy( Z[j], Z[j], Z[i], -H[j-1][i], start, end, l );
      vector_double_real_scale( Z[j], Z[j], 1/H[j-1][j], start, end, l );

      START_MASTER(threading)
        H[j-1][j-1] += sigma;
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)

#ifdef GCRODR
      if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 )
      {      
        for( i=0; i<k; i++ )
          vector_double_saxpy( Va[j], Va[j], DPCk[i], -B[j-1][i], start, end, l );
      }
#endif
      for( i=0; i<=j-1; i++ )
        vector_double_saxpy( Va[j], Va[j], Va[i], -H[j-1][i], start, end, l );
      //if ( cabs_double( H[j-1][j] ) > 1e-15 )
      vector_double_real_scale( Va[j], Va[j], 1/H[j-1][j], start, end, l );

    }

  } else {
#endif

// --------------------------------------------------------------------------------



  //if ( l->level==0 ) printf0("AHA2!\n");



  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)
  int i;
  // start and end indices for vector functions depending on thread
  int start, end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  if ( prec != NULL ) {
    if ( p->kind == _LEFT ) {
#ifdef PERS_COMMS
      //g.pers_comms_id2 = j;
      //g.use_pers_comms1 = 1;
#endif
      apply_operator_double( Z[0], V[j], p, l, threading );
#ifdef PERS_COMMS
      //g.pers_comms_id2 = -1;
      //g.use_pers_comms1 = 0;
#endif

      //MPI_Barrier(MPI_COMM_WORLD);
      //SYNC_MASTER_TO_ALL(threading);
      //SYNC_CORES(threading)

      //MPI_Barrier(MPI_COMM_WORLD);
      //SYNC_MASTER_TO_ALL(threading);
      //SYNC_CORES(threading)

      prec( w, NULL, Z[0], _NO_RES, l, threading );
    } else {
      if ( l->level == 0 ) {
        prec( Z[j], NULL, V[j], _NO_RES, l, threading );
#ifdef PERS_COMMS
        //g.pers_comms_id2 = p->restart_length + j;
        //g.use_pers_comms1 = 1;
#endif
        apply_operator_double( w, Z[j], p, l, threading );
#ifdef PERS_COMMS
        //g.pers_comms_id2 = -1;
        //g.use_pers_comms1 = 0;
#endif
      } else {
        if ( g.mixed_precision == 2 && (g.method >= 1 && g.method <= 2 ) ) {
          prec( Z[j], w, V[j], _NO_RES, l, threading );
          // obtains w = D * Z[j] from Schwarz
        } else {
          prec( Z[j], NULL, V[j], _NO_RES, l, threading );
#ifdef PERS_COMMS
          //g.pers_comms_id2 = p->restart_length + j;
          //g.use_pers_comms1 = 1;
#endif
          apply_operator_double( w, Z[j], p, l, threading ); // w = D*Z[j]
#ifdef PERS_COMMS
          //g.pers_comms_id2 = -1;
          //g.use_pers_comms1 = 0;
#endif
        }
      }
    }
  } else {
#ifdef PERS_COMMS
    //g.pers_comms_id2 = j;
    //g.use_pers_comms1 = 1;
#endif
    apply_operator_double( w, V[j], p, l, threading ); // w = D*V[j]
#ifdef PERS_COMMS
    //g.pers_comms_id2 = -1;
    //g.use_pers_comms1 = 0;
#endif
  }

#ifdef GCRODR
  // orthogonalize against Ck whenever necessary
  if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 ) {
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    int k = p->gcrodr_double.k;
    vector_double *Ck = p->gcrodr_double.C;
    complex_double **B = p->gcrodr_double.ort_B;
    // buffer
    complex_double *bf = p->gcrodr_double.Bbuff[0];

    complex_double tmpx[k+j+2];

    // merging Ck and V into a single vector of pointers
    complex_double* VCk[k+j+2];
    for( i=0;i<k;i++ ){ VCk[i] = Ck[i]; }
    for( i=0;i<(j+1);i++ ){ VCk[k+i] = V[i]; }
    VCk[k+j+1] = w; 

    process_multi_inner_product_double( k+j+2, tmpx, VCk, w, p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    // buffer is of length m, and k<m
    for ( i=0; i<(k+j+2); i++ )
      buffer[i] = tmpx[i];
    if ( g.num_processes > 1 ) {
      //if ( l->level==0 ) printf0("CALLING MPI_Allreduce(...) !!!\n");
      PROF_double_START( _ALLR );
      MPI_Allreduce( buffer, bf, k+j+2, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
      PROF_double_STOP( _ALLR, 1 );
    } else {
      for( i=0; i<(k+j+2); i++ )
        bf[i] = buffer[i];
    }

    // copy the B coefficients to the corresponding matrix
    memcpy( B[j], bf, sizeof(complex_double)*k );
    END_MASTER(threading)

    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    for( i=0; i<k; i++ )
      vector_double_saxpy( w, w, Ck[i], -B[j][i], start, end, l );

    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    START_MASTER(threading)
    // copy the H coefficients to the corresponding matrix
    memcpy( H[j], bf+k, sizeof(complex_double)*(j+2) );
    END_MASTER(threading)

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)

    for( i=0; i<=j; i++ )
      vector_double_saxpy( w, w, V[i], -H[j][i], start, end, l );
#ifdef REORTH
    // re-orthogonalization
    process_multi_inner_product_double( j+1, tmp, V, w, p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    for( i=0; i<=j; i++ )
      buffer[i] = tmp[i];
    if ( g.num_processes > 1 ) {
      //if ( l->level==0 ) printf0("CALLING MPI_Allreduce(...) !!!\n");
      PROF_double_START( _ALLR );
      MPI_Allreduce( buffer, tmp, j+1, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
      PROF_double_STOP( _ALLR, 1 );
    }

    for( i=0; i<=j; i++ )
      H[j][i] += tmp[i];

    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    for( i=0; i<=j; i++ )
      vector_double_saxpy( w, w, V[i], -tmp[i], start, end, l );
#endif

  } else {

    // orthogonalization
    complex_double tmp[j+2];
    complex_double *V_buff[j+2];

    if ( l->level==0 ) {
      for (i=0; i < j+1; i++) V_buff[i] = V[i];
      V_buff[j+1] = w;
    }

    if ( l->level==0 ) process_multi_inner_product_double( j+2, tmp, V_buff, w, p->v_start, p->v_end, l, threading );
    else process_multi_inner_product_double( j+1, tmp, V, w, p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    if ( l->level==0 ) {
      for( i=0; i<=j+1; i++ )
        buffer[i] = tmp[i];
    } else {
      for( i=0; i<=j; i++ )
        buffer[i] = tmp[i];
    }
    if ( g.num_processes > 1 ) {
      //if ( l->level==0 ) printf0("CALLING MPI_Allreduce(...) !!!\n");
      PROF_double_START( _ALLR );
      if ( l->level==0 ) MPI_Allreduce( buffer, H[j], j+2, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
      else MPI_Allreduce( buffer, H[j], j+1, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
      PROF_double_STOP( _ALLR, 1 );
    } else {
      if ( l->level==0 ) {
        for( i=0; i<=j+1; i++ )
          H[j][i] = buffer[i];
      } else {
        for( i=0; i<=j; i++ )
          H[j][i] = buffer[i];
      }
    }
    END_MASTER(threading)

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading)

    for( i=0; i<=j; i++ )
      vector_double_saxpy( w, w, V[i], -H[j][i], start, end, l );
#ifdef REORTH
    // re-orthogonalization
    process_multi_inner_product_double( j+1, tmp, V, w, p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    for( i=0; i<=j; i++ )
      buffer[i] = tmp[i];
    if ( g.num_processes > 1 ) {
      //if ( l->level==0 ) printf0("CALLING MPI_Allreduce(...) !!!\n");
      PROF_double_START( _ALLR );
      MPI_Allreduce( buffer, tmp, j+1, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
      PROF_double_STOP( _ALLR, 1 );
    }

    for( i=0; i<=j; i++ )
      H[j][i] += tmp[i];

    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    for( i=0; i<=j; i++ )
      vector_double_saxpy( w, w, V[i], -tmp[i], start, end, l );
#endif

  }

#else

  // orthogonalization
  complex_double tmp[j+2];
  complex_double *V_buff[j+2];

  if ( l->level==0 ) {
    for (i=0; i < j+1; i++) V_buff[i] = V[i];
    V_buff[j+1] = w;
  }

  if ( l->level==0 ) process_multi_inner_product_double( j+2, tmp, V_buff, w, p->v_start, p->v_end, l, threading );
  else process_multi_inner_product_double( j+1, tmp, V, w, p->v_start, p->v_end, l, threading );
  START_MASTER(threading)
  if ( l->level==0 ) {
    for( i=0; i<=j+1; i++ )
      buffer[i] = tmp[i];
  } else {
    for( i=0; i<=j; i++ )
      buffer[i] = tmp[i];
  }
  if ( g.num_processes > 1 ) {
    //if ( l->level==0 ) printf0("CALLING MPI_Allreduce(...) !!!\n");
    PROF_double_START( _ALLR );
    if ( l->level==0 ) MPI_Allreduce( buffer, H[j], j+2, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
    else MPI_Allreduce( buffer, H[j], j+1, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
    PROF_double_STOP( _ALLR, 1 );
  } else {
    if ( l->level==0 ) {
      for( i=0; i<=j+1; i++ )
        H[j][i] = buffer[i];
    } else {
      for( i=0; i<=j; i++ )
        H[j][i] = buffer[i];
    }
  }
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading);
  SYNC_CORES(threading)

  for( i=0; i<=j; i++ )
    vector_double_saxpy( w, w, V[i], -H[j][i], start, end, l );
#ifdef REORTH
  // re-orthogonalization
  process_multi_inner_product_double( j+1, tmp, V, w, p->v_start, p->v_end, l, threading );
  START_MASTER(threading)
  for( i=0; i<=j; i++ )
    buffer[i] = tmp[i];
  if ( g.num_processes > 1 ) {
    //if ( l->level==0 ) printf0("CALLING MPI_Allreduce(...) !!!\n");
    PROF_double_START( _ALLR );
    MPI_Allreduce( buffer, tmp, j+1, MPI_COMPLEX_double, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_double.level_comm );
    PROF_double_STOP( _ALLR, 1 );
  }
  
  for( i=0; i<=j; i++ )
    H[j][i] += tmp[i];

  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  for( i=0; i<=j; i++ )
    vector_double_saxpy( w, w, V[i], -tmp[i], start, end, l );
#endif

#endif // from GCRO-DR


//=============================================================
  if ( l->level==0 ) {
    START_MASTER(threading)
    complex_double tmp = H[j][j+1];
    for ( i=0; i<=j; i++ )
      tmp -= conj_double( H[j][i] )*H[j][i];
 #ifdef GCRODR
    if ( l->level==0 && p->gcrodr_double.orth_against_Ck == 1 )
    {      
      int k = p->gcrodr_double.k;
      complex_double **B = p->gcrodr_double.ort_B;

      for( i=0; i<k; i++ )
        tmp -= conj_double( B[j][i] )*B[j][i];
    }
 #endif
      
    H[j][j+1] = tmp;
    ((complex_double*)threading->workspace)[0] = creal_double(tmp);
    END_MASTER(threading)
      
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)
    if ( creal_double(((complex_double*)threading->workspace)[0]) < 0.)
    {
      H[j][j+1] = global_norm_double( w, p->v_start, p->v_end, l, threading );
    }
    else
    {
      START_MASTER(threading)
      H[j][j+1] = sqrt(creal_double(H[j][j+1]));
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)
    }

  } else {
    double tmp2 = global_norm_double( w, p->v_start, p->v_end, l, threading );
    START_MASTER(threading)
    H[j][j+1] = tmp2;
    END_MASTER(threading)
  }
//=============================================================

  SYNC_MASTER_TO_ALL(threading);
  SYNC_CORES(threading)
  
  // V_j+1 = w / H_j+1,j
  if ( cabs_double( H[j][j+1] ) > 1e-15 )
    vector_double_real_scale( V[j+1], w, 1/H[j][j+1], start, end, l );

#ifdef PIPELINED_ARNOLDI
  }
#endif

#if defined(GCRODR) || defined(MUMPS_ADDS)
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  //int jx = j-1;
  int jx;
  if ( j==0 ) {
    jx = j;
  } else {
    jx = j-1;
  }
#else
  int jx = j;
#endif
#endif

  // copy of Hesselnberg matrix (only level=0 currently)
#if defined(GCRODR) && defined(MUMPS_ADDS)
  if (l->dup_H==1 && l->level==0)
  {
    memcpy( p->gcrodr_double.eigslvr.Hc[jx], H[jx], sizeof(complex_double)*(jx+2) );
    memset( p->gcrodr_double.eigslvr.Hc[jx]+jx+2, 0.0, sizeof(complex_double)*(p->restart_length + 1 - (jx+2)) );
  }
#elif defined(GCRODR)
  if (l->dup_H==1 && l->level==0)
  {
    memcpy( p->gcrodr_double.eigslvr.Hc[jx], H[jx], sizeof(complex_double)*(jx+2) );
    memset( p->gcrodr_double.eigslvr.Hc[jx]+jx+2, 0.0, sizeof(complex_double)*(p->restart_length + 1 - (jx+2)) );
  }
#elif defined(POLYPREC)
  if (l->dup_H==1 && l->level==0)
  {
    memcpy( p->polyprec_double.eigslvr.Hc[jx], H[jx], sizeof(complex_double)*(jx+2) );
    memset( p->polyprec_double.eigslvr.Hc[jx]+jx+2, 0.0, sizeof(complex_double)*(p->restart_length + 1 - (jx+2)) );
  }
#endif

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  return 1;
}


void qr_update_double( complex_double **H, complex_double *s,
                          complex_double *c, complex_double *gamma, int j,
                          level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Applies one Givens rotation to the Hessenberg matrix H in order to solve the 
* least squares problem in (F)GMRES for computing the solution.
* - complex_double **H: Hessenberg matrix from Arnoldi decomposition
* - complex_double *s: sin values from givens rotations
* - complex_double *c: cos valies from givens rotations
* - complex_double *gamma: Approximation to residual from every step
* - int j: Denotes current iteration.
*********************************************************************************/  

  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)
  START_MASTER(threading)
  
  PROF_double_START( _SMALL1 );
  
  int i;
  complex_double beta;
  
  // update QR factorization
  // apply previous Givens rotation
  for ( i=0; i<j; i++ ) {
    beta = (-s[i])*H[j][i] + (c[i])*H[j][i+1];
    H[j][i] = conj_double(c[i])*H[j][i] + conj_double(s[i])*H[j][i+1];
    H[j][i+1] = beta;
  }
  // compute current Givens rotation
  beta = (complex_double) sqrt( NORM_SQUARE_double(H[j][j]) + NORM_SQUARE_double(H[j][j+1]) );
  s[j] = H[j][j+1]/beta; c[j] = H[j][j]/beta;
  // update right column
  gamma[j+1] = (-s[j])*gamma[j]; gamma[j] = conj_double(c[j])*gamma[j];
  // apply current Givens rotation
  H[j][j] = beta; H[j][j+1] = 0;
  
  PROF_double_STOP( _SMALL1, 6*j+6 );
  
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading);
  SYNC_CORES(threading)
}


void compute_solution_double( vector_double x, vector_double *V, complex_double *y,
                                 complex_double *gamma, complex_double **H, int j, int ol,
                                 gmres_double_struct *p, level_struct *l, struct Thread *threading ) {
  
  int i, k;
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  START_MASTER(threading)
  
  PROF_double_START( _SMALL2 );
  
  // backward substitution
  for ( i=j; i>=0; i-- ) {
    y[i] = gamma[i];
    for ( k=i+1; k<=j; k++ ) {
      y[i] -= H[k][i]*y[k];
    }
    y[i] /= H[i][i];
  }
  
  PROF_double_STOP( _SMALL2, ((j+1)*(j+2))/2 + j+1 );
  
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  
  // x = x + V*y
  if ( ol ) {
    for ( i=0; i<=j; i++ ) {
      vector_double_saxpy( x, x, V[i], y[i], start, end, l );
    }
  } else {
    vector_double_scale( x, V[0], y[0], start, end, l );
    for ( i=1; i<=j; i++ ) {
      vector_double_saxpy( x, x, V[i], y[i], start, end, l );
    }
  }
}


void local_minres_double( vector_double phi, vector_double eta, vector_double latest_iter,
                             int start, schwarz_double_struct *s, level_struct *l, struct Thread *threading ) {
  
/*********************************************************************************
* Minimal Residual iteration solver used to solve the block systems
*     blockD phi = eta
* within the Schwarz method, phi contains an initial guess and its updated version
* is returned after the block solve has been performed.
* eta is overwritten by the block residual r.
* To calculate the missing contributions to r on the current Schwarz block
* coming from outside of the block, an update "phi_new - phi_old" is returned in
* latest_iter -> cheaper residual update in the Schwarz method
*********************************************************************************/
  
  START_UNTHREADED_FUNCTION(threading)

  int i, nv = l->num_lattice_site_var, n = l->block_iter,
    end = (g.odd_even&&l->depth==0)?(start+nv*s->num_block_even_sites):(start+s->block_vector_size);
  vector_double Dr = s->local_minres_buffer[0];
  vector_double r = s->local_minres_buffer[1];
  vector_double lphi = s->local_minres_buffer[2];
  complex_double alpha;
  void (*block_op)() = (l->depth==0)?(g.odd_even?apply_block_schur_complement_double:block_d_plus_clover_double)
                                    :coarse_block_operator_double;

  vector_double_copy( r, eta, start, end, l );
  vector_double_define( lphi, 0, start, end, l );
  
  for ( i=0; i<n; i++ ) {
    // Dr = blockD*r
    block_op( Dr, r, start, s, l, no_threading );
    // alpha = <Dr,r>/<Dr,Dr>
    alpha = local_xy_over_xx_double( Dr, r, start, end, l );
    // phi += alpha * r
    vector_double_saxpy( lphi, lphi, r, alpha, start, end, l );
    // r -= alpha * Dr
    vector_double_saxpy( r, r, Dr, -alpha, start, end, l );
  }
  
  if ( latest_iter != NULL ) vector_double_copy( latest_iter, lphi, start, end, l );
  if ( phi != NULL ) vector_double_plus( phi, phi, lphi, start, end, l );
  vector_double_copy( eta, r, start, end, l );

  END_UNTHREADED_FUNCTION(threading)
}


void fgcr_double( gmres_double_struct *p, level_struct *l ) { 

/*********************************************************************************
* Uses FGCR to solve the system D x = b, where b is taken from p->b and x is 
* stored in p->x.                                                              
*********************************************************************************/

  int i, j=-1, finish=0, iter=0, il, ol;
  complex_double beta = 0, alpha;
  double r0_norm=0, t0=0, t1=0;
  
  if ( p->timing || p->print ) t0 = MPI_Wtime();
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)  
  if ( p->print ) printf0("+----------------------------------------------------------+\n");
#endif
  for( ol=0; ol<p->num_restart && finish==0; ol++ )  {
  
    if( ol == 0 && p->initial_guess_zero ) {
      vector_double_copy( p->r, p->b, p->v_start, p->v_end, l );

    } else {
      apply_operator_double( p->w, p->x, p, l, no_threading ); // compute w = D*x
      vector_double_minus( p->r, p->b, p->w, p->v_start, p->v_end, l ); // compute r = b - w
    }
    
    if( ol == 0) {
      r0_norm = global_norm_double( p->r, p->v_start, p->v_end, l, no_threading );
    }
    
    for( il=0; il<p->restart_length && finish==0; il++ ) {
      
      j = il; iter++;
      
      p->preconditioner( p->V[j], p->r, _NO_RES, l, no_threading );
      apply_operator_double( p->Z[j], p->V[j], p, l, no_threading );
      
      for( i=0; i<j; i++ ) {
        beta = global_inner_product_double( p->Z[i], p->Z[j], p->v_start, p->v_end, l, no_threading ) / p->gamma[i];
        vector_double_saxpy( p->V[j], p->V[j], p->V[i], -beta, p->v_start, p->v_end, l );
        vector_double_saxpy( p->Z[j], p->Z[j], p->Z[i], -beta, p->v_start, p->v_end, l );
      }
      
      p->gamma[j] = global_inner_product_double( p->Z[j], p->Z[j], p->v_start, p->v_end, l, no_threading );
      alpha = global_inner_product_double( p->Z[j], p->r, p->v_start, p->v_end, l, no_threading ) / p->gamma[j];
      vector_double_saxpy( p->x, p->x, p->V[j], alpha, p->v_start, p->v_end, l );
      vector_double_saxpy( p->r, p->r, p->Z[j], -alpha, p->v_start, p->v_end, l );
      
      alpha = global_norm_double( p->r, p->v_start, p->v_end, l, no_threading ) / r0_norm;
      if ( creal(alpha) < p->tol ) {
        finish = 1;
        break;
      } else {
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
        if ( iter%10 == 0 && p->print  ) printf0("| approx. rel. res. after  %-6d iterations: %e |\n", iter, alpha );
#endif
      }
    } // end of restart
  } // end of fgcr
  
  if ( p->timing || p->print ) t1 = MPI_Wtime();
  if ( p->print ) {
    apply_operator_double( p->w, p->x, p, l, no_threading );
    vector_double_minus( p->r, p->b, p->w, p->v_start, p->v_end, l );
    beta = global_norm_double( p->r, p->v_start, p->v_end, l, no_threading );
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
    printf0("+----------------------------------------------------------+\n");
    printf0("\n");
#endif
    printf0("+----------------------------------------------------------+\n");
    printf0("|         FGCR iterations: %-6d                          |\n", iter );
    printf0("| exact relative residual: ||r||/||b|| = %e      |\n", creal(beta)/r0_norm );
    printf0("| elapsed wall clock time: %-7lf seconds                |\n", t1-t0 );
    if ( g.coarse_time > 0 ) 
      printf0("|        coarse grid time: %-7lf seconds (%04.1lf%%)        |\n",
              g.coarse_time, 100*(g.coarse_time/(t1-t0)) );
    printf0("+----------------------------------------------------------+\n\n");
  }
}
