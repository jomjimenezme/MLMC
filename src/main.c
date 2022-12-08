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

    //for( int i=0;i<g.num_levels;i++ ){
    //  printf0( "g.trace_deflation_type[i] = %d\n",g.trace_deflation_type[i] );
    //  printf0( "g.trace_deflation_nr_vectors[i] = %d\n",g.trace_deflation_nr_vectors[i] );
    //  printf0( "g.trace_powerit_solver_tol[i] = %f\n",g.trace_powerit_solver_tol[i] );
    //  printf0( "g.trace_powerit_cycles[i] = %d\n",g.trace_powerit_cycles[i] );
    //  printf0( "g.trace_powerit_spectrum_type[i] = %d\n",g.trace_powerit_spectrum_type[i] );
    //}
    //char op_name[50];
    //strcpy( op_name,"difference" );
    //printf0( "comparison = %d\n",strcmp( op_name,"non-difference" ) );
    //exit(0);

    struct Thread threading;
    setup_threading(&threading, commonthreaddata, &l);
    setup_no_threading(no_threading, &l);

    // setup up initial MG hierarchy
    method_setup( NULL, &l, &threading );
    // iterative phase
    method_update( l.setup_iter, &l, &threading );
    
    //block_powerit_driver_double( &l, &threading );
    
    hutchinson_diver_double_init( &l, &threading );  
    hutchinson_diver_double_alloc( &l, &threading );
    complex_double trace, rtrace;
    struct Thread *threadingx = &threading;  

    
    START_MASTER(threadingx)
    if(g.my_rank==0) printf("Computing ROUGH trace through PLAIN Hutchinson ...\n");
    END_MASTER(threadingx)
    SYNC_MASTER_TO_ALL(threadingx)
    //rtrace = hutchinson_driver_double( &l, &threading );
    rtrace=778.0;      
    START_MASTER(threadingx)
    if(g.my_rank==0) printf("\n... done\n\n");
    END_MASTER(threadingx)
    SYNC_MASTER_TO_ALL(threadingx)


    
    // ------------ Trace with MLMC------------------------------------
    
    START_MASTER(threadingx)
    printf("Computing trace through MLMC Hutchinson ...\n");
    END_MASTER(threadingx)

    // get actual trace
    l.h_double.rough_trace = rtrace;
    l.h_double.rt= rtrace;
    l.h_double.max_iters = 10;
    l.h_double.min_iters = 5;
    l.h_double.trace_tol = 1.0e-3;
    SYNC_MASTER_TO_ALL(threadingx)
    
    //for(int i=1; i<=100; i++)
    //trace = mlmc_hutchinson_diver_double( &l, &threading );
    trace = split_mlmc_hutchinson_diver_double( &l, &threading );

    START_MASTER(threadingx)
    if(g.my_rank==0) printf("\n... done\n\n");
    printf("Resulting trace = %f+i%f\n", CSPLIT(trace));
    END_MASTER(threadingx)
    // -------------------------------------------------------  

   
   
   
   
    hutchinson_diver_double_free( &l, &threading );
    
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
