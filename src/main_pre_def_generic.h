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

#ifndef MAIN_PRE_DEF_double_HEADER
  #define MAIN_PRE_DEF_double_HEADER
  
  typedef double complex complex_double;
  typedef double complex *config_double;
  typedef double complex *vector_double;

  typedef struct {
    int length[8], *boundary_table[8], max_length[4],
        comm_start[8], in_use[8], offset, comm,
        num_even_boundary_sites[8], num_odd_boundary_sites[8],
        num_boundary_sites[8];
    vector_double buffer[8];
    MPI_Request sreqs[8], rreqs[8];
  } comm_double_struct;
  
  typedef struct {
    int ilde, dist_local_lattice[4], dist_inner_lattice_sites,
        *permutation, *gather_list, gather_list_length;
    vector_double buffer, transfer_buffer;
    MPI_Request *reqs;
    MPI_Group level_comm_group;
    MPI_Comm level_comm;
  } gathering_double_struct;
  
  typedef struct {
    double m0;
    config_double D, clover, clover_oo_inv;
    config_double odd_proj; //identity on the odd sites
    int oe_offset, self_coupling, num_even_sites, num_odd_sites,
        *index_table, *neighbor_table, *translation_table, table_dim[4],
        *backward_neighbor_table,
        table_mod_dim[4], *config_boundary_table[4];
    vector_double *buffer, prnT, prnZ, prnY, prnX, prpT, prpZ, prpY, prpX;
    comm_double_struct c;
    OPERATOR_TYPE_double *D_vectorized;
    OPERATOR_TYPE_double *D_transformed_vectorized;
    OPERATOR_TYPE_double *clover_vectorized;
    OPERATOR_TYPE_double *clover_oo_inv_vectorized;
#ifdef HAVE_TM
    double mu, mu_odd_shift, mu_even_shift;
    config_double tm_term;
#endif
#ifdef HAVE_TM1p1
    double epsbar, epsbar_ig5_odd_shift, epsbar_ig5_even_shift;
    config_double epsbar_term, clover_doublet_oo_inv;
    OPERATOR_TYPE_double *clover_doublet_vectorized;
    OPERATOR_TYPE_double *clover_doublet_oo_inv_vectorized;
#endif
  } operator_double_struct;

#if defined(POLYPREC) || defined(GCRODR)
  typedef struct
  {
    int N, nrhs, lda, ldb, info;

    int *ipiv;
    vector_double x, b;
    complex_double *Hcc;  

    void (*dirctslvr_double)();

  } dirctslvr_double_struct;
#endif

#if defined(GCRODR) || defined(POLYPREC)
  // this is both eigensolver and generalized eigensolver
  typedef struct {
    char jobvl, jobvr;

    int N, lda, ldb, ldvl, ldvr, info, qr_m, qr_n, qr_lda, qr_k;

    int *ordr_idxs;

    complex_double *ordr_keyscpy, *qr_tau;
    vector_double vl, vr, w, beta, A, B;

    complex_double **qr_QR, **qr_Q, **qr_R, **qr_Rinv;
    complex_double **Hc;

    void (*eigslvr_double)();
    void (*gen_eigslvr_double)();
  } eigslvr_double_struct;
#endif

#ifdef GCRODR
  typedef struct {
    int i, k, CU_usable, syst_size, finish, orth_against_Ck, update_CU, recompute_DPCk_poly, recompute_DPCk_plain, upd_ctr;

    double b_norm, norm_r0;

    vector_double *Pk, *C, *Cc, *U, *Yk, *hatZ, *hatW;
#ifdef BLOCK_JACOBI
    vector_double r_aux;
#endif
    // Gc is used to copy G
    complex_double **gev_A, **gev_B, **Bbuff, **QR, **Q, **R, **Rinv, **ort_B, **G, **Gc;

    eigslvr_double_struct eigslvr;

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    vector_double *PC, *DPC;
#endif
  } gcrodr_double_struct;
#endif

#ifdef POLYPREC
  typedef struct
  {
    int update_lejas;
    int d_poly;
    int syst_size;
      
    complex_double **Hc;
    complex_double *Hcc;
    complex_double **L;
    complex_double *col_prods;
    vector_double h_ritz;
    vector_double lejas;
    vector_double random_rhs;
    vector_double accum_prod, product, temp, xtmp;

    void (*preconditioner)();
    void (*preconditioner_bare)();

    eigslvr_double_struct eigslvr;
    dirctslvr_double_struct dirctslvr;
  } polyprec_double_struct;
#endif

#ifdef BLOCK_JACOBI
  typedef struct {
    vector_double x, b, r, w, *V, *Z;
    complex_double **H, *y, *gamma, *c, *s;
    config_double *D, *clover;
    operator_double_struct *op;
    double tol;
    int num_restart, restart_length, timing, print, kind,
      initial_guess_zero, layout, v_start, v_end;
    long int total_storage;
    void (*eval_operator)();

    polyprec_double_struct polyprec_double;
  } local_gmres_double_struct;

  typedef struct {
    int BJ_usable, syst_size;
    vector_double b_backup;
    vector_double xtmp;
    local_gmres_double_struct local_p;

    // for direct solves
    OPERATOR_TYPE_double* bj_op_inv_vectorized;
    OPERATOR_TYPE_double* bj_op_vectorized;
    OPERATOR_TYPE_double* bj_doublet_op_inv_vectorized;
    OPERATOR_TYPE_double* bj_doublet_op_vectorized;

    //vector_double xxxtmp[4];

  } block_jacobi_double_struct;
#endif

  typedef struct {
    vector_double x, b, r, w, *V, *Z;
    complex_double **H, *y, *gamma, *c, *s;
    config_double *D, *clover;
    operator_double_struct *op;
    double tol;
    int num_restart, restart_length, timing, print, kind,
      initial_guess_zero, layout, v_start, v_end;
    long int total_storage;
    void (*preconditioner)();
    void (*eval_operator)();
    
#ifdef GCRODR
    gcrodr_double_struct gcrodr_double;
#endif
#ifdef POLYPREC
    polyprec_double_struct polyprec_double;
#endif
#ifdef BLOCK_JACOBI
    block_jacobi_double_struct block_jacobi_double;
#endif
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    int syst_size;
    vector_double *Va, *Za;
#endif

#ifdef PERS_COMMS
    vector_double* pers_comms_ins;
    vector_double* pers_comms_outs;
#endif

#ifdef MUMPS_ADDS
    vector_double mumps_vals;
    int *mumps_Is, *mumps_Js;
    vector_double mumps_rhs_loc;
    int *mumps_irhs_loc;
    vector_double mumps_SOL;
#endif

  } gmres_double_struct;

  typedef struct {
    vector_double buffer1; //solution
    vector_double buffer2;  //source
    vector_double block_buffer;
    vector_double* X;
    vector_double sample;
       
    complex_double rough_trace;
    double total_variance;
    
    vector_double mlmc_b1;
    vector_double mlmc_b1_double;
    vector_double* X_double;
    
    complex_double rt;
    complex_double trace;
       
    int nr_rough_ests, block_size;
    int max_iters, min_iters;
    
    double trace_tol;
    
  } hutchinson_double_struct;

  typedef struct {
    operator_double_struct op;
    vector_double buf1, buf2, buf3, buf4, buf5;
    vector_double oe_buf[4];
    vector_double local_minres_buffer[3];
    int block_oe_offset, *index[4], dir_length[4], num_blocks, num_colors,
        dir_length_even[4], dir_length_odd[4], *oe_index[4],
        num_block_even_sites, num_block_odd_sites, num_aggregates,
        block_vector_size, num_block_sites, block_boundary_length[9],
        **block_list, *block_list_length;
    block_struct *block;
  } schwarz_double_struct;
  
  typedef struct {
    int num_agg, *agg_index[4], agg_length[4], *agg_boundary_index[4],
        *agg_boundary_neighbor[4], agg_boundary_length[4], num_bootstrap_vect;
    vector_double *test_vector, *interpolation, *bootstrap_vector, tmp;
    complex_double *operator, *eigenvalues, *bootstrap_eigenvalues;
  } interpolation_double_struct;
  
  typedef struct {
    double time[_NUM_PROF];
    double flop[_NUM_PROF];
    double count[_NUM_PROF];
    char name[_NUM_PROF][50];
  } profiling_double_struct;
  
  #ifdef PROFILING
    #define PROF_double_START_UNTHREADED( TYPE ) do{ l->prof_double.time[TYPE] -= MPI_Wtime(); }while(0)
    #define PROF_double_START_THREADED( TYPE, threading ) do{ if(threading->core + threading->thread == 0) l->prof_double.time[TYPE] -= MPI_Wtime(); }while(0)
  #else
    #define PROF_double_START_UNTHREADED( TYPE )
    #define PROF_double_START_THREADED( TYPE, threading )
  #endif
  
  #ifdef PROFILING
    #define PROF_double_STOP_UNTHREADED( TYPE, COUNT ) do{ l->prof_double.time[TYPE] += MPI_Wtime(); \
    l->prof_double.count[TYPE] += COUNT; }while(0)
    #define PROF_double_STOP_THREADED( TYPE, COUNT, threading ) do{ if(threading->core + threading->thread == 0) { l->prof_double.time[TYPE] += MPI_Wtime(); \
    l->prof_double.count[TYPE] += COUNT; } }while(0)
  #else
    #define PROF_double_STOP_UNTHREADED( TYPE, COUNT )
    #define PROF_double_STOP_THREADED( TYPE, COUNT, threading )
  #endif

  #define GET_MACRO2(_1,_2,NAME,...) NAME
  #define GET_MACRO3(_1,_2,_3,NAME,...) NAME
  #define PROF_double_START(...) GET_MACRO2(__VA_ARGS__, PROF_double_START_THREADED, PROF_double_START_UNTHREADED, padding)(__VA_ARGS__)
  #define PROF_double_STOP(...) GET_MACRO3(__VA_ARGS__, PROF_double_STOP_THREADED, PROF_double_STOP_UNTHREADED, padding)(__VA_ARGS__)
  
#endif
