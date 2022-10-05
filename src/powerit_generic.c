
#include "main.h"



void bp_op_double_apply( char* op_name, level_struct* lx, struct Thread* threading );
void bp_qr_double( level_struct* lx, struct Thread* threading );


void block_powerit_double_init_and_alloc( char* op_name, int depth_bp_op, int nr_vecs, int nr_bpi_cycles, double bp_tol, level_struct* l, struct Thread* threading ){

  START_MASTER(threading)

  int i;

  // access l at the right level
  level_struct* lx = l;
  for( i=0;i<g.num_levels;i++ ){
    if( i==depth_bp_op ){ break; }
    lx = l->next_level;
  }

  lx->powerit.nr_vecs = nr_vecs;
  lx->powerit.bp_tol = bp_tol;
  lx->powerit.nr_cycles = nr_bpi_cycles;

  lx->powerit.vecs = NULL;
  MALLOC(lx->powerit.vecs, complex_double*, lx->powerit.nr_vecs );
  lx->powerit.vecs[0] = NULL;
  MALLOC(lx->powerit.vecs[0], complex_double, lx->powerit.nr_vecs*lx->vector_size );
  for( i=1;i<lx->powerit.nr_vecs;i++ ){
    lx->powerit.vecs[i] = lx->powerit.vecs[0] + i*lx->vector_size;
  }

  lx->powerit.vecs_buff1 = NULL;
  MALLOC(lx->powerit.vecs_buff1, complex_double*, lx->powerit.nr_vecs );
  lx->powerit.vecs_buff1[0] = NULL;
  MALLOC(lx->powerit.vecs_buff1[0], complex_double, lx->powerit.nr_vecs*lx->vector_size );
  for( i=1;i<lx->powerit.nr_vecs;i++ ){
    lx->powerit.vecs_buff1[i] = lx->powerit.vecs_buff1[0] + i*lx->vector_size;
  }

  lx->powerit.vecs_buff2 = NULL;
  MALLOC(lx->powerit.vecs_buff2, complex_double*, lx->powerit.nr_vecs );
  lx->powerit.vecs_buff2[0] = NULL;
  MALLOC(lx->powerit.vecs_buff2[0], complex_double, lx->powerit.nr_vecs*lx->vector_size );
  for( i=1;i<lx->powerit.nr_vecs;i++ ){
    lx->powerit.vecs_buff2[i] = lx->powerit.vecs_buff2[0] + i*lx->vector_size;
  }

  END_MASTER(threading)
}


void block_powerit_double_free( char* op_name, int depth_bp_op, level_struct* l, struct Thread* threading ){

  START_MASTER(threading)

  // access l at the right level
  level_struct* lx = l;
  int i;
  for( i=0;i<g.num_levels;i++ ){
    if( i==depth_bp_op ){ break; }
    lx = l->next_level;
  }

  FREE(lx->powerit.vecs[0], complex_double, lx->powerit.nr_vecs*lx->vector_size );
  FREE(lx->powerit.vecs_buff1[0], complex_double, lx->powerit.nr_vecs*lx->vector_size );
  FREE(lx->powerit.vecs_buff2[0], complex_double, lx->powerit.nr_vecs*lx->vector_size );

  END_MASTER(threading)
}


void block_powerit_double( char* op_name, int depth_bp_op, level_struct *l, struct Thread *threading ){

  // access l at the right level
  level_struct* lx = l;
  int i;
  for( i=0;i<g.num_levels;i++ ){
    if( i==depth_bp_op ){ break; }
    lx = l->next_level;
  }

  // set the power iteration vectors to random
  START_LOCKED_MASTER(threading)
  vector_double_define_random( lx->powerit.vecs[0], 0, lx->powerit.nr_vecs*lx->vector_size, lx );
  END_LOCKED_MASTER(threading)

  for( i=0;i<lx->powerit.nr_cycles;i++ ){
    // apply the operator on the vectors ...
    bp_op_double_apply( op_name, lx, threading );
    // ... and the resulting vectors are in lx->powerit.vecs

    bp_qr_double( lx, threading );
  }
}


// auxiliary functions

void bp_op_double_apply( char* op_name, level_struct* lx, struct Thread* threading ){

  int i;

  if( strcmp(op_name,"non-difference")==0 ){
    // TODO : include threading for setting some values
  
    double buff_tol;
    complex_double* buff_b;
    complex_double* buff_x;
    gmres_double_struct* px;
    if( lx->depth==0 ){ px = &(g.p); } else { px = &(lx->p_double); }

    buff_tol = px->tol;
    px->tol = lx->powerit.bp_tol;
    buff_b = px->b;
    buff_x = px->x;

    for( i=0;i<lx->powerit.nr_vecs;i++ ){
      px->b = lx->powerit.vecs[i];
      px->x = lx->powerit.vecs_buff1[i];

      fgmres_double( px, lx, threading );
    }

    // restore values
    px->tol = buff_tol;
    px->b = buff_b;
    px->x = buff_x;
  } else if( strcmp(op_name,"difference")==0 ) {
    double buff_tol;
    int buff_print1, buff_print2;
    complex_double* buff_b;
    complex_double* buff_x;
    gmres_double_struct* px;
    if( lx->depth==0 ){ px = &(g.p); } else { px = &(lx->p_double); }
    level_struct* lxc = lx->next_level;
    gmres_double_struct* pxc = &(lxc->p_double);

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
      double buff_coarsest_tol, buff_coarse_tol;
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
        trans_back_double( lx->powerit.vecs_buff2[i], lx->sbuf_double[1], lx->s_double.op.translation_table, lx, threading );
      } else {
        interpolate3_double( lx->powerit.vecs_buff2[i], pxc->x, lx, threading );
      }

      // fine
      START_MASTER(threading)
      px->x = lx->powerit.vecs_buff1[i];
      END_MASTER(threading)
      SYNC_CORES(threading)
      fgmres_double( px, lx, threading );

      int start, end;
      compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);
      vector_double_minus( px->x, px->x, lx->powerit.vecs_buff2[i], start, end, lx );

      START_MASTER(threading)
      printf0(".");
      END_MASTER(threading)
    }
    START_MASTER(threading)
    printf0("\n");
    END_MASTER(threading)

    // restore values
    px->tol = buff_tol;
    px->b = buff_b;
    px->x = buff_x;
    px->print = buff_print1;
    g.print = buff_print2;
  } else {
    error0("Unrecognized operator to apply block power iteration\n");
  }

  // swap pointers
  START_MASTER(threading)
  complex_double* buff = lx->powerit.vecs_buff1[0];
  lx->powerit.vecs_buff1[0] = lx->powerit.vecs[0];
  lx->powerit.vecs[0] = buff;
  END_MASTER(threading)
}

void bp_qr_double( level_struct* lx, struct Thread* threading ){

  printf0("QR is under construction ...\n");
}
