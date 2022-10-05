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

global_struct g;
#ifdef HAVE_HDF5
Hdf5_fileinfo h5info;
#endif
struct common_thread_data *commonthreaddata;
struct Thread *no_threading;

int main( int argc, char **argv ) {
    
#ifdef HAVE_HDF5
  h5info.filename=NULL;
  h5info.file_id=-1; 
  h5info.rootgroup_id=-1; 
  h5info.configgroup_id=-1;
  h5info.eigenmodegroup_id=-1;
  h5info.thiseigenmodegroup_id=-1;
  h5info.isOpen=0;
  h5info.mode=-1;
#endif
  level_struct l;
  config_double hopp = NULL;
  
  MPI_Init( &argc, &argv );
  
  predefine_rank( MPI_COMM_WORLD );
  if ( g.my_rank == 0 ) {
    printf("\n\n+----------------------------------------------------------+\n");
    printf("| The DDalphaAMG solver library.                           |\n");
    printf("| Copyright (C) 2016, Matthias Rottmann, Artur Strebel,    |\n");
    printf("|       Simon Heybrock, Simone Bacchio, Bjoern Leder.      |\n");
    printf("|                                                          |\n");
    printf("| This program comes with ABSOLUTELY NO WARRANTY.          |\n");
    printf("+----------------------------------------------------------+\n\n");
  }
  
  method_init( &argc, &argv, &l );
  
  no_threading = (struct Thread *)malloc(sizeof(struct Thread));
  setup_no_threading(no_threading, &l);
  
  MALLOC( hopp, complex_double, 3*l.inner_vector_size );

  if(g.in_format == _LIME)
    lime_read_conf( (double*)(hopp), g.in, &(g.plaq_hopp) );
  else 
    read_conf( (double*)(hopp), g.in, &(g.plaq_hopp), &l );

  // store configuration, compute clover term
  dirac_setup( hopp, &l );
  FREE( hopp, complex_double, 3*l.inner_vector_size );

  commonthreaddata = (struct common_thread_data *)malloc(sizeof(struct common_thread_data));
  init_common_thread_data(commonthreaddata);
  
  THREADED(g.num_openmp_processes)
  {
    g.if_rademacher=0;

    struct Thread threading;
    setup_threading(&threading, commonthreaddata, &l);
    setup_no_threading(no_threading, &l);

    // setup up initial MG hierarchy
    method_setup( NULL, &l, &threading );
    // iterative phase
    method_update( l.setup_iter, &l, &threading );
    
    /*
    //init hutchinson
    l.h_double.block_size =12;   //set Before allocating BLOCK stuff!  
    l.h_double.max_iters = 500000; //set Before allocating BLOCK stuff!      
    hutchinson_diver_double_init( &l, &threading );  
    hutchinson_diver_double_alloc( &l, &threading );
    complex_double trace, rtrace;
    struct Thread *threadingx = &threading;  

    // ------------ROUGH with Plain-------------------------------------------
    START_MASTER(threadingx)
    if(g.my_rank==0) printf("Computing ROUGH trace through PLAIN Hutchinson ...\n");
    END_MASTER(threadingx)
    SYNC_MASTER_TO_ALL(threadingx)
    l.h_double.max_iters = 5;
    l.h_double.min_iters = 5;
    l.h_double.trace_tol=1e-4;

    rtrace = hutchinson_driver_double( &l, &threading );
            
    START_MASTER(threadingx)
    if(g.my_rank==0) printf("\n... done\n\n");
    END_MASTER(threadingx)
    // -------------------------------------------------------------------



    // ------------ Trace with Plain-------------------------------------------
    START_MASTER(threadingx)
    if(g.my_rank==0) printf("Computing trace through PLAIN Hutchinson ...\n");
    END_MASTER(threadingx)

    l.h_double.rt = rtrace;
    l.h_double.max_iters = 500000;
    l.h_double.min_iters = 5;
    l.h_double.trace_tol=1e-3;
    
    //  trace = hutchinson_driver_double( &l, &threading );
 	
    START_MASTER(threadingx)
    printf0("\n... done\n\n");
    END_MASTER(threadingx)
    SYNC_MASTER_TO_ALL(threadingx)
    // ---------------------------------------------------------------

    // ------------TRACE with BLOCK-------------------------------------------
    
    l.h_double.min_iters = 10;
    l.h_double.trace_tol = 1e-5;

   // block_hutchinson_driver_double( &l, &threading );
    
    // ---------------------------------------------------------------

    // ------------ Trace with MLMC------------------------------------
    
    START_MASTER(threadingx)
    printf0("Computing trace through MLMC Hutchinson ...\n");
    END_MASTER(threadingx)

    // get actual trace
    l.h_double.rt = rtrace;
    l.h_double.max_iters = 500000;
    l.h_double.min_iters = 5;
    l.h_double.trace_tol = 1.0e-3;
    

  //       trace = mlmc_hutchinson_diver_double( &l, &threading );
    for(int i=1; i<=100; i++)
     trace = mlmc_hutchinson_diver_double( &l, &threading );

    START_MASTER(threadingx)
    if(g.my_rank==0) printf("\n... done\n\n");
    END_MASTER(threadingx)
    
    // -------------------------------------------------------  

   // ------------ Trace with BLOCK MLMC------------------------------------
    
    START_MASTER(threadingx)
    if(g.my_rank==0) printf("Computing trace through BLOCK MLMC Hutchinson ...\n");
    END_MASTER(threadingx)
    
    l.h_double.rt = rtrace;
    l.h_double.max_iters = 1000;
    l.h_double.min_iters = 5;
    l.h_double.trace_tol = 1.0e-3;
    
    
      //trace = mlmc_block_hutchinson_diver_double( &l, &threading );
        
    START_MASTER(threadingx)
    if(g.my_rank==0) printf("\n... done\n\n");
    END_MASTER(threadingx)
   
   
   
   
    hutchinson_diver_double_free( &l, &threading );
    */


    // calling block power iteration on some operators

    // bare levels i.e. non-difference levels
    //char op_name[50] = "non-difference";
    // difference levels
    char op_name[50] = "difference";

    // depth at which to apply power iteration
    int depth_bp_op = 1;
    // number of block power iteration vectors
    int nr_bp_vecs = 24;
    // relative-residual tolerance for the operators in power iteration
    double bp_tol = 1.0e-5;
    // number of power iteration cycles
    int nr_bpi_cycles = 5;

    // IMPORTANT : always call this operation with the finest-level l
    if( strcmp(op_name,"difference") && depth_bp_op>g.num_levels ){
      error0("The depth cannot be larger than the total number of levels\n");
    }
    block_powerit_double_init_and_alloc( op_name, depth_bp_op, nr_bp_vecs, nr_bpi_cycles, bp_tol, &l, &threading );
    block_powerit_double( op_name, depth_bp_op, &l, &threading );
    block_powerit_double_free( op_name, depth_bp_op, &l, &threading );

  }

  finalize_common_thread_data(commonthreaddata);
  finalize_no_threading(no_threading);
  free(commonthreaddata);
  free(no_threading);

  method_free( &l );
  method_finalize( &l );
  
  MPI_Finalize();
  
  return 0;
}
