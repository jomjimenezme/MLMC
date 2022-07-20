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

#ifdef BLOCK_JACOBI


  // declarations of aux functions
  // TODO

  // -------------------------------------------------------------

  // main functions here -- IMPORTANT : all Block Jacobi functions are, for now,
  //			                exclusively per-process operations


  void bj_direct_op_apply_double( vector_double out, vector_double in, level_struct *l, struct Thread *threading ){

    int start, end;
    compute_core_start_end_custom( 0, l->p_double.op->num_even_sites, &start, &end, l, threading, 1 );

    int site_size = l->num_lattice_site_var;
    int lda = SIMD_LENGTH_double*((site_size+SIMD_LENGTH_double-1)/SIMD_LENGTH_double);
#ifdef HAVE_TM1p1
    OPERATOR_TYPE_double *clover = 
                            (g.n_flavours == 2) ? l->p_double.block_jacobi_double.bj_doublet_op_inv_vectorized:l->p_double.block_jacobi_double.bj_op_inv_vectorized;
#else
    OPERATOR_TYPE_double *clover = l->p_double.block_jacobi_double.bj_op_inv_vectorized;
#endif
    for(int i=start; i<end; i++) {
      for(int j=0; j<site_size; j++)
        out[i*site_size+j] = 0.0;
      cgemv(site_size, clover+i*2*site_size*lda, lda, (double *)(in+i*site_size), (double *)(out+i*site_size));
    }
  }


  void block_jacobi_apply_double( vector_double out, vector_double in, gmres_double_struct *p, level_struct *l, struct Thread *threading ){
    if ( out==in ) { return; }

    int start,end;
    compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);
    vector_double_copy( out, in, start, end, l );

    //----------------------------------------------------
    /*
    //START_MASTER(threading)

    printf("(proc=%d) PRINTING MATRICES! (level=%d, nr_procs=%d) \n", g.my_rank, l->level, g.num_processes);

    {
      int start = p->v_start;
      int end = p->v_end;
      int i;

      printf("(proc=%d) GLOBAL ... \n", g.my_rank);

      int exp_fctr = 10;

      vector_double right  = (vector_double) malloc( exp_fctr*(end-start)*sizeof(complex_double) );
      vector_double left  = (vector_double) malloc( exp_fctr*(end-start)*sizeof(complex_double) );

      vector_double gather_buff;
      if ( g.my_rank==0 ) gather_buff  = (vector_double) malloc( g.num_processes*(end-start)*sizeof(complex_double) );

      local_gmres_double_struct *loc_p = &(p->block_jacobi_double.local_p);

      FILE *fp;

      if ( g.my_rank==0 ) fp = fopen("global_coarsest_mat.txt", "w");

      //int nr_elems = (end-start)*g.num_processes;
      for (i=0;i<g.num_processes;i++) {

        // i
        for (int j=0;j<(end-start);j++) {
          for ( int k=0; k<(end-start); k++ ) { right[k] = 0; }
          if ( g.my_rank==i ) right[j] = 1.0;
          //printf("(proc=%d) %d\n", g.my_rank, j);
          MPI_Barrier( MPI_COMM_WORLD );

          //p->eval_operator( out, in, p->op, l, threading );
          coarse_apply_schur_complement_double( left,right,p->op,l,threading );

          MPI_Gather( left, end-start, MPI_COMPLEX_double, gather_buff, end-start, MPI_COMPLEX_double, 0, g.comm_cart);

          if ( g.my_rank==0 ) {
            for ( int k=0; k<(end-start)*g.num_processes; k++ ) {
              fprintf(fp, "%f+%fj -- ", creal(gather_buff[k]), cimag(gather_buff[k]));
            }
            fprintf(fp, "\n");
          }
        }

      }

      free(right);
      free(left);
      if ( g.my_rank==0 ) free(gather_buff);
      if ( g.my_rank==0 ) fclose(fp);

      printf("(proc=%d) ... done\n", g.my_rank);
    }

    {
      int start = p->v_start;
      int end = p->v_end;
      int i;

      printf("(proc=%d) LOCAL ... \n", g.my_rank);

      int exp_fctr = 10;

      vector_double right  = (vector_double) malloc( exp_fctr*(end-start)*sizeof(complex_double) );
      vector_double left  = (vector_double) malloc( exp_fctr*(end-start)*sizeof(complex_double) );

      vector_double gather_buff;
      if ( g.my_rank==0 ) gather_buff  = (vector_double) malloc( g.num_processes*(end-start)*sizeof(complex_double) );

      local_gmres_double_struct *loc_p = &(p->block_jacobi_double.local_p);

      FILE *fp;

      if ( g.my_rank==0 ) fp = fopen("local_coarsest_mat.txt", "w");

      //int nr_elems = (end-start)*g.num_processes;
      for (i=0;i<g.num_processes;i++) {

        // i
        for (int j=0;j<(end-start);j++) {
          for ( int k=0; k<(end-start); k++ ) { right[k] = 0; }
          if ( g.my_rank==i ) right[j] = 1.0;
          //printf("(proc=%d) %d\n", g.my_rank, j);
          MPI_Barrier( MPI_COMM_WORLD );

          //p->eval_operator( out, in, p->op, l, threading );
          //coarse_apply_schur_complement_double( left,right,p->op,l,threading );
          coarse_local_apply_schur_complement_double( left,right,p->op,l,threading );

          MPI_Gather( left, end-start, MPI_COMPLEX_double, gather_buff, end-start, MPI_COMPLEX_double, 0, g.comm_cart);

          if ( g.my_rank==0 ) {
            for ( int k=0; k<(end-start)*g.num_processes; k++ ) {
              fprintf(fp, "%f+%fj -- ", creal(gather_buff[k]), cimag(gather_buff[k]));
            }
            fprintf(fp, "\n");
          }
        }

      }

      free(right);
      free(left);
      if ( g.my_rank==0 ) free(gather_buff);
      if ( g.my_rank==0 ) fclose(fp);

      printf("(proc=%d) ... done\n", g.my_rank);
    }

    printf("(proc=%d) EXITING! \n", g.my_rank);

    MPI_Finalize();

    //END_MASTER(threading)

    exit(0);
    */
    //----------------------------------------------------

    //SYNC_MASTER_TO_ALL(threading)
    //SYNC_CORES(threading)

#ifndef OPTIMIZED_COARSE_SELF_COUPLING_double

    //START_MASTER(threading)
    //if ( p->block_jacobi_double.BJ_usable==1 ) {
      //local_apply_polyprec_double( out, NULL, in, 0, l, threading );
      /*
      {
        double tmpx1, tmpx2;
        int start = p->v_start;
        int end = p->v_end;
        local_gmres_double_struct *loc_p = &(p->block_jacobi_double.local_p);

        int size = end-start;

        vector_double solution = (vector_double) malloc( size*size*sizeof(complex_double) );
        vector_double rhs = (vector_double) malloc( size*size*sizeof(complex_double) );
        vector_double x = (vector_double) malloc( size*size*sizeof(complex_double) );

        vector_double_define_random( solution, start, end, l );

        loc_p->eval_operator(rhs, solution, loc_p->op, l, threading);

        // x ~ w
        local_apply_polyprec_double( x, NULL, rhs, 0, l, threading );

        vector_double diff_sol = (vector_double) malloc( size*size*sizeof(complex_double) );

        vector_double_minus( diff_sol, x, solution, start, end, l );

        tmpx1 = 0.0;
        tmpx2 = 0.0;
        VECTOR_FOR( int i=start, i<end, tmpx1 += NORM_SQUARE_double(diff_sol[i]), i++, l );
        tmpx1= (double)sqrt((double)tmpx1);
        VECTOR_FOR( int i=start, i<end, tmpx2 += NORM_SQUARE_double(solution[i]), i++, l );
        tmpx2= (double)sqrt((double)tmpx2);
        printf0("l (proc=%d) ---> approx rel error local POLYPREC = %f\n", g.my_rank, tmpx1/tmpx2);

        free(solution);
        free(rhs);
        free(x);
        free(diff_sol);
      }
      */
    //}
    //END_MASTER(threading)

#else

#ifdef BJ_DIR_SOLVS
    //bj_direct_op_apply_double( out, in, l, threading );
#else
    //if ( p->block_jacobi_double.BJ_usable==1 ) {
    //  local_apply_polyprec_double( out, NULL, in, 0, l, threading );
    //}
#endif

#endif

    //SYNC_MASTER_TO_ALL(threading)
    //SYNC_CORES(threading)

    //START_MASTER(threading)
    //MPI_Barrier( MPI_COMM_WORLD );
    //MPI_Finzalize();
    //END_MASTER(threading)

    //exit(0);
  }

  /*
  void block_jacobi_restore_from_buffer_double( vector_double out, gmres_double_struct *p, level_struct *l, struct Thread *threading ) {
    START_MASTER(threading)
    //vector_double out_x=0;
    // TODO :
    // 		1. assign out_x to the BJ buffer (the same buffer used by block_jacobi_apply_double(...))
    //		2. copy <out_x> into <out>
    END_MASTER(threading)
  }
  */

#endif
