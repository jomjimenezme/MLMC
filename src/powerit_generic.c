
#include "main.h"



void bp_op_double_apply( char* op_name, level_struct* lx, struct Thread* threading );
void bp_qr_double( level_struct* lx, struct Thread* threading );
void test_powerit_quality( char* op_name, level_struct* lx, struct Thread* threading );



void block_powerit_double_init_and_alloc( int spec_type, int op_id, int depth_bp_op, int nr_vecs, int nr_bpi_cycles, double bp_tol, level_struct* l, struct Thread* threading ){

  START_MASTER(threading)

  int i;

  // access l at the right level
  level_struct* lx = l;
  for( i=0;i<g.num_levels;i++ ){
    if( i==depth_bp_op ){ break; }
    lx = lx->next_level;
  }

  lx->powerit.nr_vecs = nr_vecs;
  lx->powerit.bp_tol = bp_tol;
  lx->powerit.nr_cycles = nr_bpi_cycles;

  lx->powerit.spec_type = spec_type;

  lx->powerit.gs_buffer = NULL;
  MALLOC( lx->powerit.gs_buffer, complex_double, 2*lx->powerit.nr_vecs );

  lx->powerit.vecs = NULL;
  MALLOC( lx->powerit.vecs, complex_double*, lx->powerit.nr_vecs );
  lx->powerit.vecs[0] = NULL;
  MALLOC( lx->powerit.vecs[0], complex_double, lx->powerit.nr_vecs*lx->vector_size );
  for( i=1;i<lx->powerit.nr_vecs;i++ ){
    lx->powerit.vecs[i] = lx->powerit.vecs[0] + i*lx->vector_size;
  }

  lx->powerit.vecs_buff1 = NULL;
  MALLOC( lx->powerit.vecs_buff1, complex_double, lx->vector_size );

  lx->powerit.vecs_buff2 = NULL;
  MALLOC( lx->powerit.vecs_buff2, complex_double, lx->vector_size );

  END_MASTER(threading)
}


void block_powerit_double_free( level_struct* l, struct Thread* threading ){

  START_MASTER(threading)

  int i,j;

  for( j=0;j<g.num_levels;j++ ){
    // in case no deflation is requested
    if( g.trace_deflation_type[j]==3 ){ continue; }

    // access l at the right level
    level_struct* lx = l;
    for( i=0;i<g.num_levels;i++ ){
      if( i==j ){ break; }
      lx = lx->next_level;
    }

    FREE( lx->powerit.vecs[0], complex_double, lx->powerit.nr_vecs*lx->vector_size );
    FREE( lx->powerit.vecs, complex_double*, lx->powerit.nr_vecs );
    FREE( lx->powerit.vecs_buff1, complex_double, lx->vector_size );
    FREE( lx->powerit.vecs_buff2, complex_double, lx->vector_size );

    FREE( lx->powerit.gs_buffer, complex_double, 2*lx->powerit.nr_vecs );
  }

  END_MASTER(threading)
}


void block_powerit_driver_double( level_struct* l, struct Thread* threading ){

  int i,op_id,spec_type;

  // specify the following in the .ini input file, at different levels
  // dx trace deflation type: 0   // 0 is difference, 1 is non-difference, 2 is split orthogonal, 3 is no deflation
  // dx trace deflation nr vectors: 10

  for( i=0;i<g.num_levels;i++ ){

    // in case no deflation is requested
    if( g.trace_deflation_type[i]==3 ){ continue; }

    switch(g.trace_deflation_type[i]){
      case 0:
        op_id = _DIFF_OP;
        break;
      case 1:
        op_id = _NON_DIFF_OP;
        break;
      case 2:
        op_id = _SPLIT_OP;
      default:
        error0("Uknown type for operator in block power iteration\n");
    }

    int depth_bp_op = i;
    int nr_bp_vecs = g.trace_deflation_nr_vectors[i];
    double bp_tol = g.trace_powerit_solver_tol[i];
    int nr_bpi_cycles = g.trace_powerit_cycles[i];

    switch(g.trace_powerit_spectrum_type[i]){
      case 0:
        spec_type = _EVs;
        break;
      case 1:
        spec_type = _SVs;
        break;
      default:
        error0("Uknown type of spectrum to be extracted\n");
    }

    // IMPORTANT :
    //		   -- always call this operation with the finest-level l
    //		   -- after calling power iteration, the result is in lx->powerit.vecs, with lx
    //		      the level struct of the chosen level
    if( depth_bp_op==(g.num_levels-1) && op_id==_DIFF_OP ){
      error0("There is no difference level operator at the coarsest level\n");
    }

    block_powerit_double_init_and_alloc( spec_type, op_id, depth_bp_op, nr_bp_vecs, nr_bpi_cycles, bp_tol, l, threading );
    //block_powerit_double( op_id, depth_bp_op, l, threading );
  }
}


void block_powerit_double( int op_id, int depth_bp_op, level_struct *l, struct Thread *threading ){

  // access l at the right level
  level_struct* lx = l;
  int i;
  for( i=0;i<g.num_levels;i++ ){
    if( i==depth_bp_op ){ break; }
    lx = lx->next_level;
  }

  // set the power iteration vectors to random
  START_LOCKED_MASTER(threading)
  vector_double_define_random( lx->powerit.vecs[0], 0, lx->powerit.nr_vecs*lx->vector_size, lx );
  END_LOCKED_MASTER(threading)
  SYNC_CORES(threading)

  for( i=0;i<lx->powerit.nr_cycles;i++ ){
    // apply the operator on the vectors ...
    //bp_op_double_apply( op_name, lx, threading );
    // ... and the resulting vectors are in lx->powerit.vecs

    bp_qr_double( lx, threading );
  }

  // in the SVs case, this tests the eigenvectors coming out of the Hermitian problem
  //test_powerit_quality( op_name, lx, threading );

  // apply gamma5 to the final result, if singular vectors are wanted
  if( strcmp(lx->powerit.spec_type,"SVs")==0 ){
    for( i=0;i<lx->powerit.nr_vecs;i++ ){
      if( lx->depth==0 ){
        gamma5_double( lx->powerit.vecs[i], lx->powerit.vecs[i], lx, threading );
      }
      else{
        int startg5, endg5;
        compute_core_start_end_custom(0, lx->inner_vector_size, &startg5, &endg5, lx, threading, lx->num_lattice_site_var );
        coarse_gamma5_double( lx->powerit.vecs[i], lx->powerit.vecs[i], startg5, endg5, lx );
      }
    }
    SYNC_CORES(threading)
  }
}



// auxiliary functions

void bp_op_double_apply( char* op_name, level_struct* lx, struct Thread* threading ){

  int i;

  // apply gamma5 before either operator
  if( strcmp(lx->powerit.spec_type,"SVs")==0 ){
    for( i=0;i<lx->powerit.nr_vecs;i++ ){
      if( lx->depth==0 ){
        gamma5_double( lx->powerit.vecs[i], lx->powerit.vecs[i], lx, threading );
      }
      else{
        int startg5, endg5;
        compute_core_start_end_custom(0, lx->inner_vector_size, &startg5, &endg5, lx, threading, lx->num_lattice_site_var );
        coarse_gamma5_double( lx->powerit.vecs[i], lx->powerit.vecs[i], startg5, endg5, lx );
      }
    }
    SYNC_CORES(threading)
  }

  // FIXME : all these if statements should be changed, and we should use pointers to functions here ...
  //         or not? Think about this a bit more

  if( strcmp(op_name,"non-difference")==0 ){
    // TODO : include threading for setting some values, in this non-difference case
    
    int start, end;
    double buff_tol;
    complex_double* buff_b;
    complex_double* buff_x;
    gmres_double_struct* px;
    if( lx->depth==0 ){ px = &(g.p); } else { px = &(lx->p_double); }

    compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);

    buff_tol = px->tol;
    buff_b = px->b;
    buff_x = px->x;

    START_MASTER(threading)
    px->tol = lx->powerit.bp_tol;
    END_MASTER(threading)

    for( i=0;i<lx->powerit.nr_vecs;i++ ){
      START_MASTER(threading)
      px->b = lx->powerit.vecs[i];
      px->x = lx->powerit.vecs_buff1;
      END_MASTER(threading)

      fgmres_double( px, lx, threading );
      
      vector_double_copy( lx->powerit.vecs[i], px->x, start, end, lx );

      START_MASTER(threading)
      if (g.my_rank==0) printf(".");
      END_MASTER(threading)
    }
    START_MASTER(threading)
    if (g.my_rank==0) printf("\n");
    END_MASTER(threading)

    // restore values
    START_MASTER(threading)
    px->tol = buff_tol;
    px->b = buff_b;
    px->x = buff_x;
    END_MASTER(threading)
  } else if( strcmp(op_name,"difference")==0 ) {
    double buff_tol;
    int buff_print1, buff_print2;
    int start, end;

    complex_double* buff_b;
    complex_double* buff_x;
    gmres_double_struct* px;
    if( lx->depth==0 ){ px = &(g.p); } else { px = &(lx->p_double); }
    level_struct* lxc = lx->next_level;
    gmres_double_struct* pxc = &(lxc->p_double);

    compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);

    // fine
    buff_tol = px->tol;
    buff_b = px->b;
    buff_x = px->x;
    buff_print1 = px->print;
    buff_print2 = g.print;
    START_MASTER(threading)
    px->tol = lx->powerit.bp_tol;
    px->print = 0;
    g.print = 0;
    END_MASTER(threading)

    for( i=0;i<lx->powerit.nr_vecs;i++ ){
      START_MASTER(threading)
      px->b = lx->powerit.vecs[i];
      END_MASTER(threading)
      SYNC_CORES(threading)

      // coarse
      if( lx->depth==0 ){
        trans_double( lx->sbuf_double[0], px->b, lx->s_double.op.translation_table, lx, threading );     
        restrict_double( pxc->b, lx->sbuf_double[0], lx, threading );
      } else {
        restrict_double( pxc->b, px->b, lx, threading );
      }
      double buff_coarsest_tol=0, buff_coarse_tol;
      if( lx->next_level->level==0 ){
        buff_coarsest_tol = g.coarse_tol;
        START_MASTER(threading)
        g.coarse_tol = lx->powerit.bp_tol;
        END_MASTER(threading)
      }
      buff_coarse_tol = pxc->tol;
      START_MASTER(threading)
      pxc->tol = lx->powerit.bp_tol;
      END_MASTER(threading)
      SYNC_CORES(threading)
      fgmres_double( pxc, lxc, threading );
      if( lx->next_level->level==0 ){
        START_MASTER(threading)
        g.coarse_tol = buff_coarsest_tol;
        END_MASTER(threading)
      }
      START_MASTER(threading)
      pxc->tol = buff_coarse_tol;
      END_MASTER(threading)
      if( lx->depth==0 ){
        interpolate3_double( lx->sbuf_double[1], pxc->x, lx, threading );
        trans_back_double( lx->powerit.vecs_buff2, lx->sbuf_double[1], lx->s_double.op.translation_table, lx, threading );
      } else {
        interpolate3_double( lx->powerit.vecs_buff2, pxc->x, lx, threading );
      }

      // fine
      START_MASTER(threading)
      px->x = lx->powerit.vecs_buff1;
      END_MASTER(threading)
      SYNC_CORES(threading)
      fgmres_double( px, lx, threading );

      vector_double_minus( lx->powerit.vecs[i], px->x, lx->powerit.vecs_buff2, start, end, lx );

      START_MASTER(threading)
      if (g.my_rank==0) printf(".");
      END_MASTER(threading)
    }
    START_MASTER(threading)
    if (g.my_rank==0) printf("\n");
    END_MASTER(threading)

    // restore values
    px->tol = buff_tol;
    px->b = buff_b;
    px->x = buff_x;
    px->print = buff_print1;
    g.print = buff_print2;
  } else if( strcmp(op_name,"split-orthogonal")==0 ) {
    // TODO : this needs to be implemented!
    error0("Block power iteration for the split orthogonal operator still needs to be implemented\n");
  } else {
    error0("Unrecognized operator to apply block power iteration\n");
  }
}


void bp_qr_double( level_struct* lx, struct Thread* threading ){

  gram_schmidt_double( lx->powerit.vecs, lx->powerit.gs_buffer, 0, lx->powerit.nr_vecs, lx, threading );
  gram_schmidt_double( lx->powerit.vecs, lx->powerit.gs_buffer, 0, lx->powerit.nr_vecs, lx, threading );
}


void test_powerit_quality( char* op_name, level_struct* lx, struct Thread* threading ){

  int i, start, end;
  complex_double** vecs_buff1;
  complex_double** vecs_buff2;

  gmres_double_struct* px;
  if( lx->depth==0 ){ px = &(g.p); } else { px = &(lx->p_double); }
  compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);

  vecs_buff1 = NULL;
  START_MASTER(threading)
  MALLOC( vecs_buff1, complex_double*, lx->powerit.nr_vecs );
  vecs_buff1[0] = NULL;
  MALLOC( vecs_buff1[0], complex_double, lx->powerit.nr_vecs*lx->vector_size );
  for( i=1;i<lx->powerit.nr_vecs;i++ ){
    vecs_buff1[i] = vecs_buff1[0] + i*lx->vector_size;
  }
  ((vector_double *)threading->workspace)[0] = vecs_buff1;
  END_MASTER(threading)
  SYNC_CORES(threading)
  vecs_buff1 = ((vector_double *)threading->workspace)[0];

  vecs_buff2 = NULL;
  START_MASTER(threading)
  MALLOC( vecs_buff2, complex_double*, lx->powerit.nr_vecs );
  vecs_buff2[0] = NULL;
  MALLOC( vecs_buff2[0], complex_double, lx->powerit.nr_vecs*lx->vector_size );
  for( i=1;i<lx->powerit.nr_vecs;i++ ){
    vecs_buff2[i] = vecs_buff2[0] + i*lx->vector_size;
  }
  ((vector_double *)threading->workspace)[0] = vecs_buff1;
  END_MASTER(threading)
  SYNC_CORES(threading)
  vecs_buff1 = ((vector_double *)threading->workspace)[0];

  // backup of lx->powerit.vecs
  for( i=0;i<lx->powerit.nr_vecs;i++ ){
    vector_double_copy( vecs_buff1[i], lx->powerit.vecs[i], start, end, lx );
  }

  // apply the operator
  bp_op_double_apply( op_name, lx, threading );

  // swap pointers
  START_MASTER(threading)
  {
    complex_double** buff_ptr = lx->powerit.vecs;
    lx->powerit.vecs = vecs_buff1;
    vecs_buff1 = buff_ptr;
  }
  END_MASTER(threading)

  for( i=0;i<lx->powerit.nr_vecs;i++ ){
    // compute the Rayleigh quotient
    complex_double rq;
    rq = global_inner_product_double( lx->powerit.vecs[i], vecs_buff1[i], px->v_start, px->v_end, lx, threading );
    double norm = global_norm_double( lx->powerit.vecs[i], 0, lx->inner_vector_size, lx, threading );
    rq /= norm;

    // print the Rayleigh quotient
    START_MASTER(threading)
    if (g.my_rank==0) printf( "Rayleigh quotient = %.16f+i%.16f \t",CSPLIT(rq) );
    END_MASTER(threading)

    // compute the eigenvalue residual
    vector_double_scale( vecs_buff2[i], lx->powerit.vecs[i], rq, start, end, lx );
    vector_double_minus( vecs_buff1[i], vecs_buff1[i], vecs_buff2[i], start, end, lx );
    double resx = global_norm_double( vecs_buff1[i], 0, lx->inner_vector_size, lx, threading );
    
    // print the residuals
    START_MASTER(threading)
    if (g.my_rank==0) printf( "Eigenvalue residual = %.16f\n",resx );
    END_MASTER(threading)
  }

  START_MASTER(threading)
  FREE( vecs_buff1[0], complex_double, lx->powerit.nr_vecs*lx->vector_size );
  FREE( vecs_buff2[0], complex_double, lx->powerit.nr_vecs*lx->vector_size );
  END_MASTER(threading)
}
