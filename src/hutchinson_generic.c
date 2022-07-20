
#include "main.h"



void rhs_define_double( vector_double rhs, level_struct *l, struct Thread *threading ) {
  
  // no hyperthreading here
  if(threading->thread != 0)
    return;

  int start = threading->start_index[l->depth];
  int end = threading->end_index[l->depth];

  if ( g.rhs == 0 ) {
    vector_double_define( rhs, 1, start, end, l );
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = ones\n");
    END_MASTER(threading)
  } else if ( g.rhs == 1 )  {
    vector_double_define( rhs, 0, start, end, l );
    if ( g.my_rank == 0 ) {
      START_LOCKED_MASTER(threading)
      rhs[0] = 1.0;
      END_LOCKED_MASTER(threading)
    }
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = first unit vector\n");
    END_MASTER(threading)
  } else if ( g.rhs == 2 ) {
    // this would yield different results if we threaded it, so we don't
    START_LOCKED_MASTER(threading)
    if ( g.if_rademacher==1 )
      vector_double_define_random_rademacher( rhs, 0, l->inner_vector_size, l );
    else
      vector_double_define_random( rhs, 0, l->inner_vector_size, l );
    END_LOCKED_MASTER(threading)
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = random\n");
    END_MASTER(threading)
  //} else if ( g.rhs == 4 ) {
  //  // this would yield different results if we threaded it, so we don't
  //  START_LOCKED_MASTER(threading)
  //  vector_double_define_random_rademacher( rhs, 0, l->inner_vector_size, l );
  //  END_LOCKED_MASTER(threading)
  //  START_MASTER(threading)
  //  if ( g.print > 0 ) printf0("rhs = random\n");
  //  END_MASTER(threading)
  } else if ( g.rhs == 3 ) {
    vector_double_define( rhs, 0, start, end, l );
  } else {
    ASSERT( g.rhs >= 0 && g.rhs <= 4 );
  }
}

void solve_double( vector_double solution, vector_double source, level_struct *l, struct Thread *threading ){

  int iter = 0, start = threading->start_index[l->depth], end = threading->end_index[l->depth];

  vector_double rhs = g.mixed_precision==2?g.p_MP.dp.b:g.p.b;
  vector_double sol = g.mixed_precision==2?g.p_MP.dp.x:g.p.x;

  vector_double_copy( rhs, source, start, end, l );  
  iter = fgmres_double( &(g.p), l, threading );
  vector_double_copy( solution, sol, start, end, l );
}


complex_double hutchinson_driver_double( level_struct *l, struct Thread *threading ) {
  
#ifdef POLYPREC
  {
    // setting flag to re-update lejas
    level_struct *lx = l;
    while (1) {
      if ( lx->level==0 ) {
        if ( g.mixed_precision==0 ) {
          lx->p_double.polyprec_double.update_lejas = 1;
          lx->p_double.polyprec_double.preconditioner = NULL;
        }
        else {
          lx->p_double.polyprec_double.update_lejas = 1;
          lx->p_double.polyprec_double.preconditioner = NULL;
        }
        break;
      }
      else { lx = lx->next_level; }
    }
  }
#endif

#ifdef GCRODR
  {
    // setting flag to re-update recycling subspace
    level_struct *lx = l;
    while (1) {
      if ( lx->level==0 ) {
        if ( g.mixed_precision==0 ) {
          //lx->p_double.gcrodr_double.CU_usable = 0;
          lx->p_double.gcrodr_double.update_CU = 1;
          lx->p_double.gcrodr_double.upd_ctr = 0;
        }
        else {
          //lx->p_double.gcrodr_double.CU_usable = 0;
          lx->p_double.gcrodr_double.update_CU = 1;
          lx->p_double.gcrodr_double.upd_ctr = 0;
        }
        break;
      }
      else { lx = lx->next_level; }
    }
  }
#endif

  hutchinson_double_struct* h = &(l->h_double);
  vector_double solution = h->buffer1, source = h->buffer2;

  int i,j, nr_ests = h->max_iters;
  
  double trace_tol = l->h_double.trace_tol;
  complex_double trace=0.0;
  complex_double rough_trace=0.0;
  complex_double aux=0.0;
  complex_double sample[nr_ests];
  double variance = 0.0;
  double RMSD = 0.0;
  rough_trace = l->h_double.rt;
  
  gmres_double_struct* p = &(g.p); //g accesesible from any func.  
  double tt0 = MPI_Wtime();
  for( i=0; i<nr_ests  ; i++ ) {
    double t0 = MPI_Wtime();
    g.if_rademacher = 1; //Compute random vector
    rhs_define_double( source, l, threading );
    g.if_rademacher = 0;

    solve_double( solution, source, l, threading ); //get A⁻¹x

    sample[i] = global_inner_product_double( source, solution, p->v_start, p->v_end, l, threading ); //compute x'A⁻¹x

    aux += sample[i];    trace = aux/ (i+1); //average the samples

    variance = 0.0;
    for (j=0; j<i; j++) {
      variance += ( sample[j]-trace )*conj( sample[j]-trace );
    }
    variance /= (j+1);

    RMSD = sqrt(variance/(j+1));

    START_MASTER(threading)
    //printf( "%d \t var %.15f \t RMSD %.15f < %.15f, \t Trace: %.15f + i%.15f \n ", i, variance, RMSD, cabs(rough_trace)*trace_tol, CSPLIT(trace)  );
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)

    double t1 = MPI_Wtime();         
    printf0("iters = %d (time = %f)\n", i, t1-t0);
    if( i!=0 && RMSD < cabs(rough_trace)*trace_tol && i>=l->h_double.min_iters-1) break;
  }
  
  START_MASTER(threading)
  if(g.my_rank==0 && rough_trace!=0) printf( "\n %.15f+i%.15f  \t %d  \n ",CSPLIT(trace), i );
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  double tt1 = MPI_Wtime();
  
  START_MASTER(threading)
  printf0("(total time = %f)\n\n",  tt1-tt0);
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)  
  return trace;
}



void hutchinson_diver_double_init( level_struct *l, struct Thread *threading ) {
  hutchinson_double_struct* h = &(l->h_double);

//PLAIN  
  h->buffer1 = NULL; //Solution
  h->buffer2 = NULL; //Source
//MLMC
  h->mlmc_b1 = NULL;
  h->mlmc_b1_double=NULL;
  
//BLOCK  
  h->block_buffer = NULL;
  h->X = NULL;            
  h->sample  = NULL;
  
  h->X_double = NULL;  
   
  SYNC_MASTER_TO_ALL(threading)
}

void hutchinson_diver_double_alloc( level_struct *l, struct Thread *threading ) {
  hutchinson_double_struct* h = &(l->h_double) ;

//For PLAIN hutchinson
  START_MASTER(threading)
  MALLOC(h->buffer1, complex_double, l->inner_vector_size ); //solution
  MALLOC(h->buffer2, complex_double, l->inner_vector_size ); // source
  ((vector_double *)threading->workspace)[0] = h->buffer1;
  ((vector_double *)threading->workspace)[1] = h->buffer2;

//For MLMC
  MALLOC( h->mlmc_b1, complex_double, l->inner_vector_size );   
  //MALLOC( h->mlmc_b1_double, complex_double, l->inner_vector_size );  
  MALLOC( h->mlmc_b1_double, complex_double, l->inner_vector_size );
  ((vector_double *)threading->workspace)[2] = h->mlmc_b1; 
  ((vector_double *)threading->workspace)[3] = h->mlmc_b1_double;

//for BLOCK
  MALLOC( h->block_buffer, complex_double, l->inner_vector_size/h->block_size );
  MALLOC( h->X, vector_double, h->block_size );
  MALLOC( h->X[0], complex_double, h->block_size* l->inner_vector_size);
  MALLOC( h->sample, complex_double, h->block_size*h->block_size*h->max_iters);
  ((vector_double *)threading->workspace)[4] = h->block_buffer;
  ((vector_double **)threading->workspace)[5] = h->X;  
  ((vector_double *)threading->workspace)[6] = h->X[0];
  ((vector_double *)threading->workspace)[7] = h->sample;

//for BLOCK mlmc  
  MALLOC( h->X_double, vector_double, h->block_size );//for mlmc
  MALLOC( h->X_double[0], complex_double, h->block_size* l->inner_vector_size);//for mlmc
  ((vector_double **)threading->workspace)[8] = h->X_double; //for mlmc
  ((vector_double *)threading->workspace)[9] = h->X_double[0]; //for mlmc 
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  h->buffer1 = ((vector_double *)threading->workspace)[0];
  h->buffer2 = ((vector_double *)threading->workspace)[1];
  h->mlmc_b1  = ((vector_double *)threading->workspace)[2];
  h->mlmc_b1_double = ((vector_double *)threading->workspace)[3];
  h->block_buffer =  ((vector_double *)threading->workspace)[4] ;
  h->X = ((vector_double **)threading->workspace)[5] ;
  h->X[0] = ((vector_double *)threading->workspace)[6] ;
  h->sample = ((vector_double *)threading->workspace)[7] ;
  h->X_double = ((vector_double **)threading->workspace)[8] ; //for mlmc
  h->X_double[0] = ((vector_double *)threading->workspace)[9] ; //for mlmc
    

  //----------- Setting the pointers to each column of X
  START_MASTER(threading)
  for(int i=0; i<l->h_double.block_size; i++)  l->h_double.X[i] = l->h_double.X[0]+i*l->inner_vector_size;
  for(int i=0; i<l->h_double.block_size; i++)  l->h_double.X_double[i] = l->h_double.X_double[0]+i*l->inner_vector_size;
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
}

void hutchinson_diver_double_free( level_struct *l, struct Thread *threading ) {
  hutchinson_double_struct* h = &(l->h_double) ;

  START_MASTER(threading)
  FREE( h->mlmc_b1, complex_double, l->inner_vector_size );   
  FREE( h->mlmc_b1_double, complex_double, l->inner_vector_size );   
  FREE(h->buffer1, complex_double, l->inner_vector_size ); //solution
  FREE(h->buffer2, complex_double, l->inner_vector_size ); // source
  
  FREE( h->block_buffer, complex_double, l->inner_vector_size/h->block_size );
  FREE( h->X[0], complex_double, h->block_size* l->inner_vector_size);
  FREE( h->X, vector_double, h->block_size );

  FREE( h->sample, complex_double, h->block_size*h->block_size*h->max_iters);

  FREE( h->X_double[0], complex_double, h->block_size* l->inner_vector_size);//for mlmc  
  FREE( h->X_double, vector_double, h->block_size );//for mlmc


  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
}




complex_double mlmc_hutchinson_diver_double( level_struct *l, struct Thread *threading ) {

  // FIXME : make this int a sort of input parameter ?
  int low_tol_restart_length = 5;

  START_MASTER(threading)
  g.avg_b1 = 0.0;
  g.avg_b2 = 0.0;
  g.avg_crst = 0.0;
  END_MASTER(threading)

  hutchinson_double_struct* h = &(l->h_double) ;
  complex_double rough_trace = h->rt;
  complex_double trace=0.0;

//--------------Setting the tolerance per level----------------
  int i;
  double est_tol = l->h_double.trace_tol;//  1E-6;
  double delta[g.num_levels];
  double tol[g.num_levels];
  double d0 = 0.8;
  delta[0] = d0; tol[0] = sqrt(delta[0]);
  //for(i=1 ;i<g.num_levels; i++){
  //  delta[i] = (1.0-d0)/(g.num_levels-1);
  //  tol[i] = sqrt(delta[i]);
  //}
  delta[1] = 0.20;
  delta[2] = 0.04;
  delta[3] = 0.01;
  tol[1] = sqrt(delta[1]);
  tol[2] = sqrt(delta[2]);
  tol[3] = sqrt(delta[3]);



//-------------------------------------------------------------  

//-----------------Setting variables in levels-----------------  
  int start, end, j, li;
  int nr_ests[g.num_levels];
  int counter[g.num_levels];
  double RMSD;
  double variance=0.0;  
  nr_ests[0]=l->h_double.max_iters; 
  for(i=1; i< g.num_levels; i++) nr_ests[i] = nr_ests[i-1]*1; //TODO: Change in every level?
  
  gmres_double_struct *ps[g.num_levels]; //float p_structs for coarser levels
  gmres_double_struct *ps_double[1];   //double p_struct for finest level
  level_struct *ls[g.num_levels];  
  ps_double[0] = &(g.p);
  ls[0] = l;

  level_struct *lx = l->next_level;
  double buff_tol[g.num_levels];
  buff_tol[0] =  ps_double[i]->tol; 
 /* START_MASTER(threading)
  if(g.my_rank==0) printf("LEVEL----------- \t %d  \t vector size: %d \t TOLERANCE: %e\n", 0, ls[0]->inner_vector_size, buff_tol[0]);
  END_MASTER(threading)*/

  for( i=1;i<g.num_levels;i++ ){
    SYNC_MASTER_TO_ALL(threading)
    ps[i] = &(lx->p_double); // p_double in each level
    ls[i] = lx;                  // level_struct of each level
    lx = lx->next_level;
    buff_tol[i]= ps[i]->tol;

    START_MASTER(threading)
    //if(g.my_rank==0) printf("LEVEL----------- \t %d  \t vector size: %d \tTOLERANCE: %.15f\n ", i, ls[i]->inner_vector_size, buff_tol[i]);
    END_MASTER(threading)         
  }

  // enforcing this at the coarsest level for ~10^-1 solves
  int tmp_length = ps[g.num_levels-1]->restart_length;
  ps[g.num_levels-1]->restart_length = low_tol_restart_length;

  complex_double es[g.num_levels];  //An estimate per level
  memset( es,0,g.num_levels*sizeof(complex_double) );
  complex_double tmpx = 0.0;
//---------------------------------------------------------------------  

//-----------------Hutchinson for all but coarsest level-----------------    
  for( li=0;li<(g.num_levels-1);li++ ){

    double tt0 = MPI_Wtime();

    variance=0.0;    counter[li]=0;
    complex_double sample[nr_ests[li]]; 
    compute_core_start_end( 0, ls[li]->inner_vector_size, &start, &end, l, threading );

    START_MASTER(threading)
    // if(g.my_rank==0) printf("LEVEL----------- \t %d  \t ests: %d  \t start: %d  \t end: %d \t vector size: %d \n ", li, nr_ests[li], start, ps_double[li]->v_end, ls[li]->inner_vector_size);
    END_MASTER(threading)
    
    START_MASTER(threading)
    printf0("starting level difference = %d\n\n", li);
    END_MASTER(threading)
    
    for(i=0;i<nr_ests[li]  ; i++){

      START_MASTER(threading) //Get random vector:
      if(li==0){
        vector_double_define_random_rademacher( ps_double[li]->b, 0, ls[li]->inner_vector_size, ls[li] );          
      }
      else
        vector_double_define_random_rademacher( ps[li]->b, 0, ls[li]->inner_vector_size, ls[li] );    
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)


      //-----------------Solve the system in current level and the next one----------------- 
      if(li==0){

        trans_double( l->sbuf_double[0], ps_double[li]->b, l->s_double.op.translation_table, l, threading );     
        restrict_double( ps[li+1]->b, l->sbuf_double[0], ls[li], threading ); // get rhs for next level.

        // this solve is at the second level
        double bfx1 = ps[li+1]->tol;
        ps[li+1]->tol = g.tol;
        double t0 = MPI_Wtime();
        int coarse_iters = fgmres_double( ps[li+1], ls[li+1], threading );               //solve in next level
        double t1 = MPI_Wtime();
        ps[li+1]->tol = bfx1;
        START_MASTER(threading)
        printf0("coarse iters = %d (time = %f)\n", coarse_iters, t1-t0);
        END_MASTER(threading)
          
        interpolate3_double( l->sbuf_double[1], ps[li+1]->x, ls[li], threading ); //      
        trans_back_double( h->mlmc_b1, l->sbuf_double[1], l->s_double.op.translation_table, l, threading );

        int fine_iters = fgmres_double( ps_double[li], ls[li], threading );   
        
        START_MASTER(threading)
        printf0("fine iters = %d\n", fine_iters);
        END_MASTER(threading)

      }
      else{
        restrict_double( ps[li+1]->b, ps[li]->b, ls[li], threading ); // get rhs for next level.
        ps[li+1]->tol = g.tol;
        g.coarse_tol = g.tol;
        double t0 = MPI_Wtime();

        //int coarse_iters = fgmres_double( ps[li+1], ls[li+1], threading );               //solve in next level

        gmres_double_struct *px = ps[li+1];
        level_struct *lx = ls[li+1];
        int start, end, coarsest_iters;
        compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);

        if( (li+2)==g.num_levels ){

          px->restart_length = tmp_length;

#ifdef POLYPREC
          px->preconditioner = px->polyprec_double.preconditioner;
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
#ifdef POLYPREC
          if ( lx->level==0 && lx->p_double.polyprec_double.update_lejas == 1 ) {
            if ( coarsest_iters >= 1.5*px->polyprec_double.d_poly ) {
              // re-construct Lejas
              re_construct_lejas_double( lx, threading );
            }
          }
#endif

          px->restart_length = low_tol_restart_length;

        } else {
          coarsest_iters = fgmres_double( px, lx, threading );
        }

        double t1 = MPI_Wtime();
        START_MASTER(threading)
        if( (li+1)==3 ){
          printf0("coarsest iters = %d (time = %f)\n", coarsest_iters, t1-t0);
        } else {
          printf0("coarse iters = %d (time = %f)\n", coarsest_iters, t1-t0);
        }
        END_MASTER(threading)

        ps[li+1]->tol = buff_tol[li+1];
        g.coarse_tol = buff_tol[li+1];
        interpolate3_double( h->mlmc_b1_double, ps[li+1]->x, ls[li], threading ); //   
        ps[li]->tol = g.tol;   
        t0 = MPI_Wtime();

        int fine_iters = fgmres_double( ps[li], ls[li], threading );

        //// checking relative residual of coarsest-level solves
        //double normx1 = global_norm_double( ps[li]->b,start,end,ls[li],threading );
        //apply_operator_double( ps[li]->w,ps[li]->x,ps[li],ls[li],threading ); // compute w = D*x
        //vector_double_minus( ps[li]->r,ps[li]->b,ps[li]->w,start,end,ls[li] ); // compute r = b - w
        //double normx2 = global_norm_double( ps[li]->r,start,end,ls[li],threading );
        //printf0("REL RESIDUAL SECOND LEVEL : %.16f\n", normx2/normx1);

        t1 = MPI_Wtime();
        START_MASTER(threading)
        printf0("fine iters = %d (time = %f)\n", fine_iters, t1-t0);
        END_MASTER(threading)
        ps[li+1]->tol = buff_tol[li+1];
      }     
      //------------------------------------------------------------------------------------
      
      //--------------Get the difference of the solutions and the corresponding sample--------------                        
      if(li==0){

        tmpx = global_inner_product_double( ps_double[li]->b, ps_double[li]->x, ps_double[li]->v_start, ps_double[li]->v_end, ls[li], threading );      
        printf0("first term = %.15f+i%.15f\n", CSPLIT(tmpx));
        tmpx = global_inner_product_double( ps_double[li]->b, h->mlmc_b1, ps_double[li]->v_start, ps_double[li]->v_end, ls[li], threading );      
        printf0("second term = %.15f+i%.15f\n", CSPLIT(tmpx));

        vector_double_minus( h->mlmc_b1, ps_double[li]->x, h->mlmc_b1, start, end, ls[li] );
        tmpx = global_inner_product_double( ps_double[li]->b, h->mlmc_b1, ps_double[li]->v_start, ps_double[li]->v_end, ls[li], threading );      
        sample[i]= tmpx;      
        es[li] += tmpx;
      }
      else{
        vector_double_minus( h->mlmc_b1_double, ps[li]->x, h->mlmc_b1_double, start, end, ls[li] );
        tmpx =  global_inner_product_double( ps[li]->b, h->mlmc_b1_double, ps[li]->v_start, ps[li]->v_end, ls[li], threading );      
        sample[i]= tmpx;      
        es[li] +=  tmpx;
      }
      //----------------------------------------------------------------------------------------------

      //--------------Get the Variance in the current level and use it in stop criteria---------------  
      variance = 0.0;
      for (j=0; j<i; j++) // compute the variance
        variance += (sample[j]- es[li]/(i+1)) *conj( sample[j]- es[li]/(i+1)); 
      variance /= (j+1);
      RMSD = sqrt(variance/(j+1)); //RMSD= sqrt(var+ bias²)
     
      START_MASTER(threading)
      //printf( "%d \t var %.15f \t Est %.15f  \t RMSD %.15f <  %.15f ?? \n ", i, variance, creal( es[li]/(i+1) ), RMSD,  cabs(rough_trace)*tol[li]*est_tol );
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)     
      counter[li]=i+1;
      if(i !=0 && RMSD< cabs(rough_trace)*tol[li]*est_tol && i>=l->h_double.min_iters-1){counter[li]=i+1; break;}
      //----------------------------------------------------------------------------------------------    
    }

    double tt1 = MPI_Wtime();

    START_MASTER(threading)
    printf0("ending difference = %d (total time = %f)\n\n", li, tt1-tt0);
    END_MASTER(threading)

  }
  
    START_MASTER(threading)
    printf0("starting coarsest level\n\n", li);
    END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)
//-----------------Hutchinson for just coarsest level-----------------       

  li = g.num_levels-1;
  
 /* START_MASTER(threading)
  if(g.my_rank==0)printf("LEVEL----------- \t %d  \t ests: %d  \t start: %d  \t end: %d \t vector size: %d \n ", li, nr_ests[li], start, ps[li]->v_end, ls[li]->inner_vector_size);
  END_MASTER(threading)*/
  
  complex_double sample[nr_ests[li]]; 
  variance = 0.0;
  counter[li]=0;

  double tt1c = MPI_Wtime();

  for(i=0;i<nr_ests[li]  ; i++){
    START_MASTER(threading)
    vector_double_define_random_rademacher( ps[li]->b, 0, ls[li]->inner_vector_size, ls[li] );    
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    
    ps[li]->tol = g.tol;
    g.coarse_tol = g.tol;
    double t0 = MPI_Wtime();

    //int coarsest_iters = fgmres_double( ps[li], ls[li], threading );

    gmres_double_struct *px = ps[li];
    level_struct *lx = ls[li];

    px->restart_length = tmp_length;

    int start, end, coarsest_iters;
    compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);
#ifdef POLYPREC
    px->preconditioner = px->polyprec_double.preconditioner;
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
#ifdef POLYPREC
    if ( lx->level==0 && lx->p_double.polyprec_double.update_lejas == 1 ) {
      if ( coarsest_iters >= 1.5*px->polyprec_double.d_poly ) {
        // re-construct Lejas
        re_construct_lejas_double( lx, threading );
      }
    }
#endif




    //// checking relative residual of coarsest-level solves
    //double normx1 = global_norm_double( px->b,start,end,lx,threading );
    //apply_operator_double( px->w,px->x,px,lx,threading ); // compute w = D*x
    //vector_double_minus( px->r,px->b,px->w,start,end,lx ); // compute r = b - w
    //double normx2 = global_norm_double( px->r,start,end,lx,threading );
    //printf0("REL RESIDUAL COARSEST LEVEL : %.16f\n", normx2/normx1);



    double t1 = MPI_Wtime();
    START_MASTER(threading)
    printf0("coarsest iters = %d (time = %f)\n", coarsest_iters, t1-t0);
    END_MASTER(threading)
    ps[li]->tol = buff_tol[li];
    g.coarse_tol = buff_tol[li];
    tmpx= global_inner_product_double( ps[li]->b, ps[li]->x, ps[li]->v_start, ps[li]->v_end, ls[li], threading );   
    sample[i]= tmpx;      
    es[li] += tmpx;

    variance = 0.0;
    for (j=0; j<i; j++) // compute the variance
      variance += (sample[j]- es[li]/(i+1)) *conj( sample[j]- es[li]/(i+1)); 
    variance /= (j+1);
    RMSD = sqrt(variance/(j+1)); //RMSD= sqrt(var+ bias²)
        
    START_MASTER(threading)
    //printf( "%d \t var %.15f \t Est %.15f  \t RMSD %.15f <  %.15f ?? \n ", i, variance, creal( es[li]/(i+1) ), RMSD,  cabs(rough_trace)*tol[li]*est_tol );
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    counter[li]=i+1;
    if(i !=0 && RMSD< cabs(rough_trace)*tol[li]*est_tol && i>=l->h_double.min_iters-1 ){counter[li]=i+1; break;}
  }

  double tt2c = MPI_Wtime();

    START_MASTER(threading)
    printf0("end of coarsest (total time = %f)\n\n", li, tt2c-tt1c);
    END_MASTER(threading)

  trace = 0.0;
  for( i=0;i<g.num_levels;i++ ){
    trace += es[i]/counter[i];    
  START_MASTER(threading)
  //if(g.my_rank==0)printf("Level:%d, ................................%f    %d\n", i, creal(es[i]),counter[i]);
  END_MASTER(threading)  

  }
  
  START_MASTER(threading)
  if(g.my_rank==0)
    //printf( "\n\n Trace: %.15f+i %.15f, \t ests in l_0: %d, \t ests in l_1: %d \n", CSPLIT( trace ), counter[0],counter[1],counter[2]);
      printf( "%.15f + i %.15f \t %d \t %d \n", CSPLIT( trace ), counter[0],counter[1],counter[2]);
  END_MASTER(threading)
  
  
  return trace;
  
}






complex_double mlmc_block_hutchinson_diver_double( level_struct *l, struct Thread *threading ) {

// FIXME : make this int a sort of input parameter ?
  int low_tol_restart_length = 5;

  START_MASTER(threading)
  g.avg_b1 = 0.0;
  g.avg_b2 = 0.0;
  g.avg_crst = 0.0;
  END_MASTER(threading)
  
  complex_double cov[g.num_levels][l->h_double.block_size*l->h_double.block_size]; 
    
  hutchinson_double_struct* h = &(l->h_double) ;
  complex_double rough_trace = h->rt;
  complex_double trace=0.0;
//--------------Setting the tolerance per level----------------
  int i;
  double est_tol = l->h_double.trace_tol;//  1E-6;
  double delta[g.num_levels];
  double tol[g.num_levels];
  double d0 = 0.75;
  delta[0] = d0; tol[0] = sqrt(delta[0]);
  //for(i=1 ;i<g.num_levels; i++){
  //  delta[i] = (1.0-d0)/(g.num_levels-1);
  //  tol[i] = sqrt(delta[i]);
  //}
  delta[1] = 0.20;
  delta[2] = 0.04;
  delta[3] = 0.01;
  tol[1] = sqrt(delta[1]);
  tol[2] = sqrt(delta[2]);
  tol[3] = sqrt(delta[3]);

//-------------------------------------------------------------  


//-----------------Setting variables in levels-----------------  
  int start, end, j, li;
  int nr_ests[g.num_levels];
  int counter[g.num_levels];
  double RMSD;
  nr_ests[0]=l->h_double.max_iters; 
  
  for(i=1; i< g.num_levels; i++) nr_ests[i] = nr_ests[i-1]*1; //TODO: Change in every level?
  
  gmres_double_struct *ps[g.num_levels]; //float p_structs for coarser levels
  //gmres_double_struct *ps_double[1];   //double p_struct for finest level
  level_struct *ls[g.num_levels];  
  ps[0] = &(g.p);
  ls[0] = l;
  
  level_struct *lx = l->next_level;
  double buff_tol[g.num_levels];
  buff_tol[0] =  ps[0]->tol; 
  START_MASTER(threading)
  if(g.my_rank==0) printf("LEVEL----------- \t %d  \t vector size: %ld \t TOLERANCE: %e\n", 0, ls[0]->inner_vector_size, buff_tol[0]);
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  for( i=1;i<g.num_levels;i++ ){
    SYNC_MASTER_TO_ALL(threading)
    ps[i] = &(lx->p_double); // p_PRECISION in each level
    ls[i] = lx;                  // level_struct of each level
    lx = lx->next_level;
    buff_tol[i]= ps[i]->tol;

    START_MASTER(threading)
    if(g.my_rank==0) printf("LEVEL----------- \t %d  \t vector size: %ld \tTOLERANCE: %f\n ", i, ls[i]->inner_vector_size, buff_tol[i]);
    END_MASTER(threading)
          
  }

  // enforcing this at the coarsest level for ~10^-1 solves
  int tmp_length = ps[g.num_levels-1]->restart_length;
  ps[g.num_levels-1]->restart_length = low_tol_restart_length;
  
  //-----------------------------BLOCK Variables----------------------------------------  
  int block_size=h->block_size, row=0, col=0, v=0;
  //vector_double* X = ls[0]->h_double.X;
  vector_double* X_double = ls[0]->h_double.X_double;
  vector_double sample = l->h_double.sample;  
  complex_double       es[g.num_levels][block_size*block_size];  //A BLOCK estimate per level
  complex_double variance[g.num_levels][block_size*block_size] ; 
  memset( es      ,0,g.num_levels*sizeof(complex_double)*block_size*block_size );
  memset( variance,0,g.num_levels*sizeof(complex_double)*block_size*block_size );
  complex_double tmpx =0.0;
  





//-----------------Hutchinson for all but coarsest level-----------------    
  for( li=0;li<(g.num_levels-1);li++ ){         
  
    double tt0 = MPI_Wtime();
    
    counter[li]=0;       l->h_double.total_variance=0.0;
    compute_core_start_end( 0, ls[li]->inner_vector_size, &start, &end, ls[li], threading );
    
    START_MASTER(threading)
    if(g.my_rank==0) printf("LEVEL----------- \t %d  \t ests: %d  \t start: %d  \t end: %d \t vector size: %ld \n \n", li, nr_ests[li], start, ps[li]->v_end, ls[li]->inner_vector_size);
    END_MASTER(threading)
       
    START_MASTER(threading)
    printf0("starting level difference = %d\n\n", li);
    END_MASTER(threading)
  //  exit(0);
    for(i=0;i<nr_ests[li] ; i++){
      START_MASTER(threading) //Get random vector:
      if(li==0){
        vector_double_define_random_rademacher( ps[li]->b, 0, ls[li]->inner_vector_size/block_size, ls[li] );

        //---------------------- Fill Big X----------------------------------------------      
        for(col=0; col<block_size; col++)
          for(row=col; row<ls[li]->inner_vector_size; row+=block_size)
            X_double[col][row] = ps[li]->b[row/12];                  
      }
      else{
        vector_double_define_random_rademacher( ps[li]->b, 0, ls[li]->inner_vector_size/block_size, ls[li] );
            
        //---------------------- Fill Big X----------------------------------------------      
        for(col=0; col<block_size; col++)
          for(row=col; row<ls[li]->inner_vector_size; row+=block_size)
            X_double[col][row] = ps[li]->b[row/12];    
      }    
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)
      
      //-----------------Solve the system in current level and the next one-----------------
      for(col=0; col<block_size; col++){
        if(li==0){
          //trans_double( l->sbuf_double[0], X[col], l->s_double.op.translation_table, l, threading );              
          restrict_double( ps[li+1]->b, X_double[col], ls[li], threading ); // get rhs for next level.
         
         
   
         
          // this solve is at the second level
          double bfx1 = ps[li+1]->tol;
          ps[li+1]->tol = g.tol;
          double t0 = MPI_Wtime();
          int coarse_iters = fgmres_double( ps[li+1], ls[li+1], threading );               //solve in next level
          double t1 = MPI_Wtime();
          ps[li+1]->tol = bfx1;
          START_MASTER(threading)
          printf0("coarse iters = %d (time = %f)\n", coarse_iters, t1-t0);
          END_MASTER(threading)
          
           
         
          
          interpolate3_double( h->mlmc_b1, ps[li+1]->x, ls[li], threading ); //      
          //trans_back_double( h->mlmc_b1, l->sbuf_double[1], l->s_double.op.translation_table, ls[li], threading ); //solution of next level in double
          
          vector_double_copy( ps[li]->b, X_double[col], ps[li]->v_start, ps[li]->v_end, ls[li] ); //Necesary for fgmres      
          int fine_iters = fgmres_double( ps[li], ls[li], threading );    
          START_MASTER(threading)
          //printf0("fine iters = %d\n", fine_iters);
          END_MASTER(threading)
        }     
        else{
          
          restrict_double( ps[li+1]->b, X_double[col], ls[li], threading ); // get rhs for next level.
          ps[li+1]->tol = g.tol;
          g.coarse_tol = g.tol;
          double t0 = MPI_Wtime();

          gmres_double_struct *px = ps[li+1];
          level_struct *lx = ls[li+1];
          int start, end, coarsest_iters;
          compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);
 
       
        
          coarsest_iters = fgmres_double( px, lx, threading );
        

        double t1 = MPI_Wtime();
        START_MASTER(threading)
        if( (li+2)==g.num_levels ){
          printf0("coarsest iters = %d (time = %f)\n", coarsest_iters, t1-t0);
        } else {
          printf0("coarse iters = %d (time = %f)\n", coarsest_iters, t1-t0);
        }
        END_MASTER(threading)

        ps[li+1]->tol = buff_tol[li+1];
        g.coarse_tol = buff_tol[li+1];
        interpolate3_double( h->mlmc_b1_double, ps[li+1]->x, ls[li], threading ); //   
        ps[li]->tol = g.tol;   
        t0 = MPI_Wtime();


        vector_double_copy( ps[li]->b, X_double[col], ps[li]->v_start, ps[li]->v_end, ls[li] ); //Necesary for fgmres      
        int fine_iters = fgmres_double( ps[li], ls[li], threading );

        t1 = MPI_Wtime();
        START_MASTER(threading)
        printf0("fine iters = %d (time = %f)\n", fine_iters, t1-t0);
        END_MASTER(threading)
        ps[li+1]->tol = buff_tol[li+1];  
         
      
        }
       
      //------------------------------------------------------------------------------------
      
      //--------------Get the difference of the solutions and the corresponding sample--------------                        
      if(li==0){     
        vector_double_minus( h->mlmc_b1, ps[li]->x, h->mlmc_b1, start, end, ls[li] );    
      }
      else{
        vector_double_minus( h->mlmc_b1_double, ps[li]->x, h->mlmc_b1_double, start, end, ls[li] );
      }
     
      //--------------Get the samples--------------  
      for(row=0; row<block_size; row++){       
        if(row==col){                  
          if(li==0){     
            tmpx = global_inner_product_double( X_double[row], h->mlmc_b1, ps[li]->v_start, ps[li]->v_end, ls[li], threading );      
          }
          else{
            tmpx =  global_inner_product_double(X_double[row], h->mlmc_b1_double, ps[li]->v_start, ps[li]->v_end, ls[li], threading );      
          }
        
          START_MASTER(threading)
          sample[col+row*block_size+ i*block_size*block_size ] = tmpx;
          es[li][col+row*block_size] += tmpx;
         
    
          END_MASTER(threading)      
          SYNC_MASTER_TO_ALL(threading)
        }//if statement     
      }//loop for row 
 
    }//loop for col    
    
      //---------------------Computing Covariance Matix for diag(BT)-------------------------
    
    START_MASTER(threading)
    l->h_double.total_variance=0.0;
    memset(cov, 0.0, sizeof(cov));
    int counter=0;
    for ( row=0; row< block_size*block_size; row+=block_size+1){ //over the diagonal, block_size elements  of BT 
      for ( col=0; col< block_size*block_size; col+=block_size+1){
        for (v=0; v<i; v++)
         cov[li][counter] += (sample[row+ v*block_size*block_size ] - es[li][row]/(i+1) )*
                               conj( sample[col+ v*block_size*block_size ] - es[li][col]/(i+1) ); 
                                        //conj(sample[j*block_size+ v*block_size*block_size ] - estimate[j*block_size]/(k+1) ); 
    
         cov[li][counter] /= (v+1);

         //l->h_PRECISION.total_variance += cov[i/block_size];
        //if(g.my_rank==0)printf("%d \t %d \t %f +i %f\n", i, counter,  CSPLIT(cov[li][counter]) ); 
        //if(g.my_rank==0)printf("Cov[ %d, %d] \t %f +i %f\n", i,j,  CSPLIT(cov[li][counter]) ); 
        counter++;        
      }
    }
    
    //-----------------Covariance for diag(gamma_5*BT)
    for(row=0; row< block_size; row++){
      for (col=0; col<block_size; col++){
       // if(row<block_size/2 && col< block_size/2 || row>block_size/2 && col> block_size/2 ){
l->h_double.total_variance += cov[li][row+col*block_size];
          //if(g.my_rank==0)printf("%d \t %d \t \t PLUS \t %f \n",  i, j, creal(cov[li][row/block_size] ) );
       // }
        //else{
        //  l->h_double.total_variance += cov[li][row/block_size];
          //if(g.my_rank==0)printf("%d \t %d \t \t minus \n",  i, j );
       // }
      }
    }
    
    
    END_MASTER(threading)      
    SYNC_MASTER_TO_ALL(threading)
    
    
    
    
    
    double RMSD = sqrt(l->h_double.total_variance)/(i+1);
    counter[li]=i+1;
    START_MASTER(threading) 
    if(g.my_rank==0)
    printf( "%d \t First var %f \t First Est %f  \t RMSD %f <  %f ?? \n ", i, creal(variance[li][0]), creal( es[li][0]/(i+1) ), RMSD,  cabs(rough_trace)*tol[li]*est_tol );
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    if(i !=0 && RMSD< cabs(rough_trace)*tol[li]*est_tol && i>=l->h_double.min_iters-1){counter[li]=i+1; break;}
    
   //----------------------------------------------------------------------------------------------    
   }// for i
   
   
    double tt1 = MPI_Wtime();

    START_MASTER(threading)
    printf0("ending difference = %d (total time = %f)\n\n", li, tt1-tt0);
    END_MASTER(threading)
 }// for li    


  START_MASTER(threading)
  printf0("starting coarsest level\n\n", li);
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)
//-----------------Hutchinson for just coarsest level-----------------       
  li = g.num_levels-1; 
  counter[li]=0; l->h_double.total_variance=0.0;
  START_MASTER(threading)
  if(g.my_rank==0)printf("LEVEL----------- \t %d  \t ests: %d  \t start: %d  \t end: %d \t vector size: %d \n\n ", li, nr_ests[li], start, ps[li]->v_end, ls[li]->inner_vector_size);
  END_MASTER(threading)
  
  double tt1c = MPI_Wtime();
 
  for(i=0;i<nr_ests[li]  ; i++){
    START_MASTER(threading)
    vector_double_define_random_rademacher( ps[li]->b, 0, ls[li]->inner_vector_size/block_size, ls[li] );
    //---------------------- Fill Big X----------------------------------------------      
    for(col=0; col<block_size; col++)
      for(row=col; row<ls[li]->inner_vector_size; row+=block_size)
        X_double[col][row] = ps[li]->b[row/12];                        
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    
    ps[li]->tol = g.tol;
    g.coarse_tol = g.tol;
    double t0 = MPI_Wtime();
    
    for(col=0; col<block_size; col++){   
      
      vector_double_copy( ps[li]->b, X_double[col], ps[li]->v_start, ps[li]->v_end, ls[li] );
      
      
      gmres_double_struct *px = ps[li];
      level_struct *lx = ls[li];

      px->restart_length = tmp_length;

      int start, end, coarsest_iters;
      compute_core_start_end(px->v_start, px->v_end, &start, &end, lx, threading);
#ifdef POLYPREC
      px->preconditioner = px->polyprec_double.preconditioner;
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
#ifdef POLYPREC
      if ( lx->level==0 && lx->p_double.polyprec_double.update_lejas == 1 ) {
        if ( coarsest_iters >= 1.5*px->polyprec_double.d_poly ) {
        // re-construct Lejas
          re_construct_lejas_double( lx, threading );
        }
      }
#endif
      
      
      
      
    double t1 = MPI_Wtime();
    START_MASTER(threading)
    printf0("coarsest iters = %d (time = %f)\n", coarsest_iters, t1-t0);
    END_MASTER(threading)
    ps[li]->tol = buff_tol[li];
    g.coarse_tol = buff_tol[li];
    for(row=0; row<block_size; row++){ 
      if(row==col){                        
        tmpx= global_inner_product_double( X_double[row], ps[li]->x, ps[li]->v_start, ps[li]->v_end, ls[li], threading );   
        
        START_MASTER(threading)
        sample[col+row*block_size+ i*block_size*block_size ] = tmpx;
        es[li][col+row*block_size] += tmpx;
      
        END_MASTER(threading)      
        SYNC_MASTER_TO_ALL(threading)
      }  
    }
   }
       
      //---------------------Computing Covariance Matix for diag(BT)-------------------------
    
    START_MASTER(threading)
    l->h_double.total_variance=0.0;
    memset(cov, 0.0, sizeof(cov));
    int counter=0;
    for ( row=0; row< block_size*block_size; row+=block_size+1){ //over the diagonal, block_size elements  of BT 
      for ( col=0; col< block_size*block_size; col+=block_size+1){
        for (v=0; v<i; v++)
         cov[li][counter] += (sample[row+ v*block_size*block_size ] - es[li][row]/(i+1) )*
                               conj( sample[col+ v*block_size*block_size ] - es[li][col]/(i+1) ); 
                                        //conj(sample[j*block_size+ v*block_size*block_size ] - estimate[j*block_size]/(k+1) ); 
    
         cov[li][counter] /= (v+1);

         //l->h_PRECISION.total_variance += cov[i/block_size];
        //if(g.my_rank==0)printf("%d \t %d \t %f +i %f\n", i, counter,  CSPLIT(cov[li][counter]) ); 
        //if(g.my_rank==0)printf("Cov[ %d, %d] \t %f +i %f\n", i,j,  CSPLIT(cov[li][counter]) ); 
        counter++;        
      }
    }
    
    //-----------------Covariance for diag(gamma_5*BT)
    for(row=0; row< block_size; row++){
      for (col=0; col<block_size; col++){
       // if(row<block_size/2 && col< block_size/2 || row>block_size/2 && col> block_size/2 ){
          l->h_double.total_variance += cov[li][row+col*block_size];
          //if(g.my_rank==0)printf("%d \t %d \t \t PLUS \t %f \n",  i, j, creal(cov[li][row/block_size] ) );
       // }
        //else{
        //  l->h_PRECISION.total_variance += cov[li][row/block_size];
          //if(g.my_rank==0)printf("%d \t %d \t \t minus \n",  i, j );
       // }
      }
    }
    
    
    END_MASTER(threading)      
    SYNC_MASTER_TO_ALL(threading)
   
    
    
    
    
    double RMSD = sqrt(l->h_double.total_variance)/(i+1);  
    counter[li]=i+1;
   
    START_MASTER(threading)
    if(g.my_rank==0)
    printf( "%d \t var %f \t First Est %f  \t RMSD %f <  %f ?? \n ", i, creal(l->h_double.total_variance), creal( es[li][0]/(i+1) ), RMSD,  cabs(rough_trace)*tol[li]*est_tol );
    END_MASTER(threading)
   
    if(i !=0 && RMSD< cabs(rough_trace)*tol[li]*est_tol && i>=l->h_double.min_iters-1){counter[li]=i+1; break;}

  }//for i//


 
    
  double tt2c = MPI_Wtime();
  START_MASTER(threading)
  printf0("end of coarsest (total time = %f)\n\n", li, tt2c-tt1c);
  END_MASTER(threading)


  START_MASTER(threading)
  complex_double     block_trace[block_size*block_size];
  complex_double  block_variance[block_size*block_size];
  memset( block_trace   ,0, sizeof(complex_double)*block_size*block_size );
  memset( block_variance,0, sizeof(complex_double)*block_size*block_size );
  //gathering the BLOCK estimates per level:
  for( li=0;li<(g.num_levels);li++ ){
    for(j=0; j<block_size; j++){  
      for (i=0; i< block_size; i++){
        block_trace   [j+i*block_size] +=       es[li][j+i*block_size]/ counter[li];
        block_variance[j+i*block_size] += variance[li][j+i*block_size];     
      }
    }
  }
  //temporal printing
  if(g.my_rank==0) printf( "\n\n------------------------- BLOCK TRACE-----------------\n");   
  for(j=0; j<block_size; j++){  
    for (i=0; i< block_size; i++){
      if(g.my_rank==0) printf( "%f\t", cabs( block_trace[j+i*block_size] )); 
      if(i==j) trace += block_trace[j+i*block_size];
    } 
    if(g.my_rank==0) printf( "\n");   
  }
  
  //temporal printing
 /* if(g.my_rank==0) printf( "\n\n------------------------- BLOCK Variance-----------------\n");   
  for(j=0; j<block_size; j++){  
    for (i=0; i< block_size; i++){
      if(g.my_rank==0) printf( "%f\t", cabs( block_variance[j+i*block_size] ));  
    } 
    if(g.my_rank==0) printf( "\n");   
  }*/
  
  END_MASTER(threading)
    
  START_MASTER(threading)
  if(g.my_rank==0)
    printf( "\n\n Trace: %f+i %f, \t ests in l_0: %d, \t ests in l_1: %d \n", CSPLIT( trace ), counter[0],counter[1]);
  END_MASTER(threading)

  
  return trace;
}















































void block_hutchinson_double( level_struct *l, struct Thread *threading , complex_double* estimate, 
                                 double* variance,
                                 complex_double* cov ){
  vector_double solution = l->h_double.buffer2; vector_double* X = l->h_double.X;
  vector_double buffer = l->h_double.block_buffer; vector_double sample = l->h_double.sample;
  complex_double rough_trace = l->h_double.rt;
  complex_double trace = l->h_double.trace;
  gmres_double_struct* p = &(g.p);
  int k, i, j, v, counter=0;
  int block_size=l->h_double.block_size, nr_ests=l->h_double.max_iters; 
  double trace_tol = l->h_double.trace_tol;
  

/*for(k=1; k<11; k++){
  vector_double_define_random_rademacher( l->h_double.buffer1, 0, l->inner_vector_size, l );    
  apply_operator_double(p->b, l->h_double.buffer1, p, l, threading ); // compute Ax=b
  solve_double( solution, p->b, l, threading); //get solution= A⁻¹x

   vector_double_minus( l->h_double.buffer2, solution, l->h_double.buffer1, p->v_start, p->v_start, l ); //  error = x_aprox - x
   double normx2 = global_norm_double( l->h_double.buffer2, p->v_start,p->v_end,l,threading );	    //norm
    printf("Error_NORM : %.16f\n", normx2);

//// checking relative residual of coarsest-level solves
        double normx1 = global_norm_double( p->b,p->v_start,p->v_end,l,threading ); //norm(b)
        apply_operator_double(p->w, solution, p, l, threading ); // compute w = D*x
        vector_double_minus( p->r,p->b,p->w,p->v_start,p->v_end,l ); // compute r = b - w
        normx2 = global_norm_double( p->r,p->v_start,p->v_end,l,threading );// norm(r)
        printf("REL RESIDUAL : %.16f\n", normx2/normx1); //relative residual
}

 
exit(0);
*/


  for( k=0; k<nr_ests; k++ ) { //Hutchinson Method loop
    START_MASTER(threading)
    // Initialize buffer with Rademacher
    vector_double_define_random_rademacher( buffer, 0, l->inner_vector_size/block_size, l );    
   //---------------------- Fill Big X----------------------------------------------
   //Zeros TODO: do a malloc or something
     for(j=0; j<block_size; j++)
      for(i=j; i<l->inner_vector_size; i++)
        X[j][i] = 0.0;  

    for(j=0; j<block_size; j++)
      for(i=j; i<l->inner_vector_size; i+=block_size)
        X[j][i] = buffer[i/12];
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    //-------------------------------------------------------------------------------
    

/*if(g.my_rank==0){  
    char a[100] ;    
    sprintf(a, "%s%d%s", "BIG",g.my_rank, ".txt");
    char fileSpec[strlen(a)+1];
    snprintf( fileSpec, sizeof( fileSpec ), "%s", a );
       
   FILE * fp;  
   fp = fopen (fileSpec, "w+");
      
   for(j=0; j<block_size; j++){
    for(i=0; i<l->inner_vector_size; i++){
       fprintf(fp, "%f\n",  creal(X[j][i]));
   //fprintf(fp, "%d \t %d\t %f ", j,i, creal(X[j][i]));
   }
      //fprintf(fp, "\n");
   }
   fclose(fp);  
}

exit(0);*/




    //-----------------------COMPUTING THE BLOCK TRACE SAMPLE------------------------- 
    for(j=0; j<block_size; j++){  
      solve_double( solution, X[j], l, threading); //get A⁻¹x 
              /*vector_double_copy( p->b, X[j], p->v_start, p->v_end, l );  
		          fgmres_double( p, l,  threading ); 
              vector_double_copy( solution, p->x, p->v_start, p->v_end, l );  */
    for (i=0; i< block_size; i++){
        if(i==j){
          complex_double tmpx = global_inner_product_double( X[i], solution, p->v_start, p->v_end, l, threading ); //compute x'A⁻¹x     
          START_MASTER(threading)
          sample[j+i*block_size+ k*block_size*block_size ] = tmpx;
          estimate[j+i*block_size] += sample[j+i*block_size+ k*block_size*block_size ];   
          
          END_MASTER(threading)      
          SYNC_MASTER_TO_ALL(threading)
        }
      }          
    }   
    //----------------------------------------------------------------------------    
    
    
    //---------------------Computing Covariance Matix for diag(BT)-------------------------
    
    START_MASTER(threading)
    l->h_double.total_variance=0.0;
    memset(cov, 0.0, sizeof(cov));
    int counter=0;
    for ( i=0; i< block_size*block_size; i+=block_size+1){ //over the diagonal, block_size elements  of BT 
      for ( j=0; j< block_size*block_size; j+=block_size+1){
        for (v=0; v<k; v++)
         cov[counter] += (sample[i+ v*block_size*block_size ] - estimate[i]/(k+1) )*
                               conj( sample[j+ v*block_size*block_size ] - estimate[j]/(k+1) ); 
                                        //conj(sample[j*block_size+ v*block_size*block_size ] - estimate[j*block_size]/(k+1) ); 
    
         cov[counter] /= (v+1);

         //l->h_double.total_variance += cov[i/block_size];
        //if(g.my_rank==0)printf("%d \t %d \t %f +i %f\n", k, counter,  CSPLIT(cov[counter]) ); 
        //if(g.my_rank==0)printf("Cov[ %d, %d] \t %f +i %f\n", i,j,  CSPLIT(cov[counter]) ); 
        counter++;        
      }
    }
    
    //-----------------Covariance for diag(gamma_5*BT)
    for(i=0; i< block_size; i++){
      for (j=0; j<block_size; j++){
        if(i<block_size/2 && j< block_size/2 || i>block_size/2 && j> block_size/2 ){
          l->h_double.total_variance += cov[i/block_size];
          //if(g.my_rank==0)printf("%d \t %d \t \t PLUS \n",  i, j );
        }
        else{
          l->h_double.total_variance +=cov[i/block_size];//-= cov[i/block_size];
          //if(g.my_rank==0)printf("%d \t %d \t \t minus \n",  i, j );
        }
      }
    }
    
    
    END_MASTER(threading)      
    SYNC_MASTER_TO_ALL(threading)
    
    
    
    
    
    
    
    //---------------------------------------------------------------------------- 
    double RMSD = sqrt(cabs(l->h_double.total_variance))/(k+1);
    counter=k;
    if(k !=0 && RMSD < cabs(rough_trace)*trace_tol && k> l->h_double.min_iters-1 ){counter=k; break; }
    
        complex_double tr=0.0;
        START_MASTER(threading) 

    for (i=0; i< block_size*block_size; i+=block_size+1){  //Compute gamma_5trace
      if(i<block_size*block_size/2) tr += estimate[i]/(k+1);//-= estimate[i]/(k+1);
      else tr += estimate[i]/(k+1);
            //if(g.my_rank==0)printf( "Trace: %f + i%f \t  i=%d \n ", CSPLIT(tr), i );
    }  
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    
    

    START_MASTER(threading) 
    /*if(g.my_rank==0)
    printf( "%d \t var %f \t RMSD  %f < %f ??\t current_trace %f \t rough_trace %f\n", k, cabs(l->h_double.total_variance), RMSD, cabs(rough_trace)*trace_tol, creal(tr), creal(rough_trace)  );*/
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    
    l->h_double.total_variance=0.0;  //Reset Total Variation
  } //LOOP Hutchinson
  
  
 
  
  //-----------------------COMPUTING TRACE -------------------------  
  START_MASTER(threading)
  for (i=0; i< block_size*block_size; i++)  {estimate[i] /=(counter+1); } //Average each estimate    
  //for (i=0; i< block_size*block_size; i+=block_size+1){  trace += estimate[i];} //Compute trace 
  
  for (i=0; i< block_size*block_size; i+=block_size+1){  //Compute gamma_5trace
    if(i<block_size*block_size/2) trace += estimate[i];
    else trace += estimate[i];//-= estimate[i];
      if(g.my_rank==0)printf( "Trace: %f + i%f \t  i=%d \n ", CSPLIT(trace), i );
  }  
  l->h_double.trace=trace;
  
  //---------------------------------------------------------------------------- 

  if(g.my_rank==0)printf( "Trace: %f + i%f \t using %d estimates \n ", CSPLIT(trace), counter+1  );
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

}



void block_hutchinson_driver_double( level_struct *l, struct Thread *threading ) {
  START_MASTER(threading)
  if(g.my_rank==0)printf("Inside block_hutchinson_driver(...) ...\n");
  END_MASTER(threading)
 //----------------------------INITIALIZE VARIABLES---------------------------------------------------- 
  complex_double estimate[l->h_double.block_size*l->h_double.block_size] ;
  double variance[l->h_double.block_size*l->h_double.block_size] ;
  complex_double cov[l->h_double.block_size*l->h_double.block_size]; 
  memset(estimate, 0.0, sizeof(estimate));
  memset(variance, 0.0, sizeof(variance));
  memset(cov, 0.0, sizeof(cov));
  
  int block_size= l->h_double.block_size;

  l->h_double.trace=0.0; l->h_double.rough_trace=0.0;
  l->h_double.total_variance=0.0; 
  
  int i;  
  

  /*//----------------------------Compute ROUGH block TRACE----------------------------------------------------  
  START_MASTER(threading)
  if(g.my_rank==0)printf("\t------computing ROUGH block trace ----------\n");
  END_MASTER(threading)

  l->h_double.max_iters = 10;  
  block_hutchinson_double(l, threading, estimate, variance, cov);

  START_MASTER(threading)
  if(g.my_rank==0)printf("\t------------- DONE--------------------\n\n");
  END_MASTER(threading)

  START_MASTER(threading)
  for (i=0; i< l->h_double.block_size*l->h_double.block_size; i+=l->h_double.block_size+1){
    if(i<l->h_double.block_size*l->h_double.block_size/2) l->h_double.rough_trace += estimate[i];
    else l->h_double.rough_trace -= estimate[i];  
   }
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)*/

//---------------------------Compute BLOCK TRACE------------------------------------------------ 
  START_MASTER(threading)
   if(g.my_rank==0) printf("\t-------computing the block trace -----------\n");
  END_MASTER(threading)
 
  l->h_double.total_variance=0.0; l->h_double.trace=0.0; 
  memset(estimate, 0.0, sizeof(estimate));
  memset(variance, 0.0, sizeof(variance)); memset(cov,0.0, sizeof(cov));
  
  block_hutchinson_double(l, threading, estimate, variance, cov);
  
  START_MASTER(threading)
  if(g.my_rank==0)printf("\t------------------DONE----------------\n\n");
  END_MASTER(threading)

   

  
  

//------------------------SAVE BLOCK TRACE INTO A FILE-----------------------------------
  START_MASTER(threading)
  int j;    
 if(g.my_rank==0){     
    char a[100] ;   
    sprintf(a, "%s%d%s", "BLOCK_TRACE",g.my_rank, ".txt");
    char fileSpec[strlen(a)+1];
    snprintf( fileSpec, sizeof( fileSpec ), "%s", a );  
    FILE * fp;  
    fp = fopen (fileSpec, "w+");
    
    sprintf(a, "%s%d%s", "Covar",g.my_rank, ".txt");
    char fileSpec_var[strlen(a)+1];
    snprintf( fileSpec_var, sizeof( fileSpec_var ), "%s", a );  
    FILE * fvar;  
    fvar = fopen (fileSpec_var, "w+");
    
    for (i=0; i< l->h_double.block_size; i++){
      for(j=0; j<l->h_double.block_size; j++){
        fprintf(fp, "%f\t", creal( estimate[i*block_size+j] ) );
        fprintf(fvar, "%f\t",  creal(   cov[i*block_size+j] ) );
        //fprintf(fp, "%f\t", creal(vec_BT[i*block_size+j]));
        //fprintf(fvar, "%f\t", creal(variance[i*block_size+j]));
       }
      fprintf(fp, "\n");
      fprintf(fvar, "\n");
    }
         
    fclose(fp);    fclose(fvar);  
  }     
  
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading) 
  
}










/*
 START_MASTER(threading)
  complex_PRECISION     block_trace[block_size*block_size];
  complex_PRECISION  block_variance[block_size*block_size];
  memset( block_trace   ,0, sizeof(complex_PRECISION)*block_size*block_size );
  memset( block_variance,0, sizeof(complex_PRECISION)*block_size*block_size );
  //gathering the BLOCK estimates per level:
  for( li=0;li<(g.num_levels);li++ ){
    for(j=0; j<block_size; j++){  
      for (i=0; i< block_size; i++){
        block_trace   [j+i*block_size] +=       es[li][j+i*block_size]/ counter[li];
        block_variance[j+i*block_size] += variance[li][j+i*block_size];     
      }
    }
  }
  //temporal printing
  if(g.my_rank==0) printf( "\n\n------------------------- BLOCK TRACE-----------------\n");   
  for(j=0; j<block_size; j++){  
    for (i=0; i< block_size; i++){
      if(g.my_rank==0) printf( "%f\t", cabs( block_trace[j+i*block_size] )); 
      if(i==j) trace += block_trace[j+i*block_size];
    } 
    if(g.my_rank==0) printf( "\n");   
  }
  
  //temporal printing
  if(g.my_rank==0) printf( "\n\n------------------------- BLOCK Variance-----------------\n");   
  for(j=0; j<block_size; j++){  
    for (i=0; i< block_size; i++){
      if(g.my_rank==0) printf( "%f\t", cabs( block_variance[j+i*block_size] ));  
    } 
    if(g.my_rank==0) printf( "\n");   
  }
  
  END_MASTER(threading)
*/

/*

// 1 (Radamacher vectors of size of level 0)
A1 - P1 * A2^{-1} * R1

// 2 (Radamacher vectors of size of level 1)
P1 * A2^{-1} * R1 - P1 * P2 * A3^{-1} * R2 * R1
A2^{-1} - P2 * A3^{-1} * R2

// 3
P1 * P2 * A3^{-1} * R2 * R1 - P1 * P2 * P3 * A4^{-1} * R3 * R2 * R1
A3^{-1} - P3 * A4^{-1} * R3

// 4
P1 * P2 * P3 * A4^{-1} * R3 * R2 * R1
A4^{-1}

*/
