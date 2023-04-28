
    #include "main.h"
    
    
  struct estimate {
    int counter; //required number of estimates.
    complex_double estimate;
  };


  struct sample {
    // required number of estimates
    int sample_size;
    // accumulated trace
    complex_double acc_trace;
    complex_double direct_trace;
  };


    
    
  void hutchinson_diver_double_init( level_struct *l, struct Thread *threading ) {
    hutchinson_double_struct* h = &(l->h_double);
    
    //MLMC
    h->mlmc_b1 = NULL;
    
    h->mlmc_b2 =NULL;
    
    h->rademacher_vector =NULL;
    
    SYNC_MASTER_TO_ALL(threading)
  }
    
  void hutchinson_diver_double_alloc( level_struct *l, struct Thread *threading ) {
    hutchinson_double_struct* h = &(l->h_double) ;
    
    //For MLMC
    PUBLIC_MALLOC( h->mlmc_b1, complex_double, l->inner_vector_size );
    PUBLIC_MALLOC( h->mlmc_b2, complex_double, l->inner_vector_size );
    
    PUBLIC_MALLOC( h->rademacher_vector, complex_double, l->inner_vector_size );
    
  }
    
  void hutchinson_diver_double_free( level_struct *l, struct Thread *threading ) {
    hutchinson_double_struct* h = &(l->h_double) ;
    
    PUBLIC_FREE( h->mlmc_b1, complex_double, l->inner_vector_size );   
    PUBLIC_FREE( h->mlmc_b2, complex_double, l->inner_vector_size );   
    PUBLIC_FREE( h->rademacher_vector, complex_double, l->inner_vector_size );   
  }
    

    
    



  // -------------------------------------------------------------------


  void rademacher_create( level_struct *l, hutchinson_double_struct* h, int type, struct Thread *threading );
  complex_double hutchinson_mlmc_difference( level_struct *l, hutchinson_double_struct* h, struct Thread *threading );
  complex_double hutchinson_split_intermediate( level_struct *l, hutchinson_double_struct* h, struct Thread *threading );
  complex_double hutchinson_deflated_split_orthogonal(level_struct *l, hutchinson_double_struct* h, struct Thread *threading);
  //complex_double hutchinson_split_orthogonal( level_struct *l, hutchinson_double_struct* h, struct Thread *threading );
  complex_double hutchinson_plain( level_struct *l, hutchinson_double_struct* h, struct Thread *threading );
  void apply_P_double( vector_double out, vector_double in, level_struct* l, struct Thread *threading );
  void apply_R_double( vector_double out, vector_double in, level_struct* l, struct Thread *threading );
  int apply_solver_double( level_struct* l, struct Thread *threading );
  gmres_double_struct* get_p_struct_double( level_struct* l );
  complex_double hutchinson_deflated_direct_term(level_struct *l, struct Thread *threading);

  // type : in case of 0 create Rademacher vectors at level l, in case of 1 create Rademacher vectors at level l->next_level
  struct sample hutchinson_blind_double( level_struct *l, hutchinson_double_struct* h, int type, struct Thread *threading ){

    int i, j;
    complex_double one_sample, variance, trace;
    double RMSD;
    struct sample estimate;

    // TODO : move this allocation to some init function
    complex_double* samples = (complex_double*) malloc( h->max_iters*sizeof(complex_double) );
    memset( samples, 0.0, h->max_iters*sizeof(complex_double) );

    estimate.acc_trace = 0.0;

    // stochastic part

    double t0 = MPI_Wtime();
    for( i=0; i<h->max_iters;i++ ){

      // 1. create Rademacher vector, stored in h->rademacher_vector
      rademacher_create( l, h, type, threading );

      // 2. apply the operator to the Rademacher vector
      // 3. dot product
      double t0 = MPI_Wtime();
      one_sample = h->hutch_compute_one_sample( l, h, threading );
      double t1 = MPI_Wtime();
      samples[i] = one_sample;
      
      // 4. compute estimated trace and variance, print something?
      estimate.acc_trace += one_sample;

      if( i!=0 ){
        variance = 0.0;
	estimate.sample_size = i+1;
        trace = estimate.acc_trace/estimate.sample_size;

        //trace = estimate.acc_trace/(i+1);
        for( j=0; j<i; j++ ){
          variance += conj(samples[j] - trace) * (samples[j] - trace);
        }
        variance = variance / j;
	    START_MASTER(threading);
        if(g.my_rank==0) printf( "%d\tVariance:\t%f\tTime:\t%f\n", i, creal(variance),t1-t0);
        END_MASTER(threading);
        RMSD = sqrt(creal(variance)/j);
        if( i > h->min_iters && RMSD < cabs(trace) * h->trace_tol * h->tol_per_level[l->depth]) break; 
      }
    }
    double t1 = MPI_Wtime();

    START_MASTER(threading);
    if(g.my_rank==0) printf( "%d\t \tvariance = %f+i%f \t t = %f, \n", i, CSPLIT(variance), t1-t0);// h->tol_per_level[l->depth]);
    END_MASTER(threading);
    estimate.sample_size = i;


    free(samples);

    return estimate;
  }


  complex_double hutchinson_driver_double( level_struct *l, struct Thread *threading ){
    
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "Trace computation via Hutchinson's method ...\n" );
    END_MASTER(thrading);
    
    int i;
    complex_double trace = 0.0;
    struct sample estimate;
    hutchinson_double_struct* h = &(l->h_double);
    level_struct* lx;
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\tfinest (and only) level ...\n" );
    END_MASTER(thrading);

    lx = l;
    // set the pointer to the finest-level Hutchinson estimator
    h->hutch_compute_one_sample = hutchinson_plain;
    
    double t_plain = MPI_Wtime();
    estimate = hutchinson_blind_double( lx, h, 0, threading );
    trace += estimate.acc_trace/estimate.sample_size;
    double t_plain1 = MPI_Wtime();
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\t... done\n" );
    if(g.my_rank==0)  printf( "Average solve time %f \n", (t_plain1-t_plain)/h->max_iters );
    END_MASTER(thrading);

    START_MASTER(threading);
    if(g.my_rank==0)  printf( "... done\n" );
    END_MASTER(thrading);
    
    //If deflation vectors are available
    if(g.trace_deflation_type[l->level] =! 3){
    trace += hutchinson_deflated_direct_term(l, threading);
    }
    
    return trace;
  }

  
  complex_double mlmc_hutchinson_driver_double( level_struct *l, struct Thread *threading ){
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "Trace computation via 'traditional' difference levels ...\n" );
    END_MASTER(thrading);

    int i;
    complex_double trace = 0.0;
    struct sample estimate;
    hutchinson_double_struct* h = &(l->h_double);
    level_struct* lx;

    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\tdifference levels ...\n" );
    END_MASTER(thrading);
    
    // for all but coarsest level
    lx = l;
    for( i=0; i<g.num_levels-1;i++ ){
      // set the pointer to the mlmc difference operator
      h->hutch_compute_one_sample = hutchinson_mlmc_difference;
      estimate = hutchinson_blind_double( lx, h, 0, threading );
      trace += estimate.acc_trace/estimate.sample_size;
      lx = lx->next_level;
    }
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\t... done\n" );
    END_MASTER(thrading);
    
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\tcoarsest level ...\n" );
    END_MASTER(thrading);
    
    // coarsest level
    // set the pointer to the coarsest-level Hutchinson estimator
    h->hutch_compute_one_sample = hutchinson_plain;
    estimate = hutchinson_blind_double( lx, h, 0, threading );
    trace += estimate.acc_trace/estimate.sample_size;
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\t... done\n" );
    END_MASTER(thrading);
    
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "... done\n" );
    END_MASTER(thrading);

    return trace;
  }


  complex_double split_mlmc_hutchinson_driver_double( level_struct *l, struct Thread *threading ){
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "Trace computation via split levels ...\n" );
    END_MASTER(thrading);

    int i;
    complex_double trace = 0.0;
    struct sample estimate;
    hutchinson_double_struct* h = &(l->h_double);
    level_struct* lx;
    
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\tfull-rank difference levels ...\n" );
    END_MASTER(thrading);
    
    // for all but coarsest level
    lx = l;
    for( i=0; i<g.num_levels-1;i++ ){  
      // set the pointer to the split intermediate operator
      h->hutch_compute_one_sample = hutchinson_split_intermediate;
      estimate = hutchinson_blind_double( lx, h, 1, threading );
      trace += estimate.acc_trace/estimate.sample_size;
      lx = lx->next_level;    
    }
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\t... done\n" );
    END_MASTER(thrading);

    START_MASTER(threading);
    if(g.my_rank==0) printf( "\torthogonalized difference levels ...\n" );
    END_MASTER(thrading);
    
    // for all but coarsest level
    lx = l;
    for( i=0; i<g.num_levels-1;i++ ){      
      // set the pointer to the split orthogonal operator
      h->hutch_compute_one_sample = hutchinson_split_orthogonal;//hutchinson_split_orthogonal;
      estimate = hutchinson_blind_double( lx, h, 0, threading );
      trace += estimate.acc_trace/estimate.sample_size;
      lx = lx->next_level; 
    }
    
    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\t... done\n" );
    END_MASTER(thrading);

    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\tcoarsest level ...\n" );
    END_MASTER(thrading);

    // coarsest level
    // set the pointer to the coarsest-level Hutchinson estimator
    h->hutch_compute_one_sample = hutchinson_plain;
    estimate = hutchinson_blind_double( lx, h, 0, threading );
    trace += estimate.acc_trace/estimate.sample_size;

    START_MASTER(threading);
    if(g.my_rank==0)  printf( "\t... done\n" );
    END_MASTER(thrading);

    return trace;
  }


  complex_double hutchinson_plain( level_struct *l, hutchinson_double_struct* h, struct Thread *threading ){

    {
      int start, end;
      gmres_double_struct* p = get_p_struct_double( l );
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

      vector_double_copy( p->b, h->rademacher_vector, start, end, l );
    }
    
    {
      apply_solver_double( l, threading );
    }
    // subtract the results and perform dot product
    {
      int start, end;
      gmres_double_struct* p = get_p_struct_double( l );
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      
      if(g.trace_deflation_type[l->level] =! 3){
        hutchinson_deflate_vector_double(p->x, l, threading); 
      }
       return global_inner_product_double( h->rademacher_vector, p->x, p->v_start, p->v_end, l, threading );   
    }
  }


  // the term tr( (I - P_{l} P_{l}^{H}) A_{l}^{-1} )
  complex_double hutchinson_split_orthogonal( level_struct *l, hutchinson_double_struct* h, struct Thread *threading ){

    // 1. project
    // 2. invert

    // FIRST TERM

    {
      int start, end;
      gmres_double_struct* p = get_p_struct_double( l );

      apply_R_double( h->mlmc_b2, h->rademacher_vector, l, threading );
      apply_P_double( h->mlmc_b1, h->mlmc_b2, l, threading );
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      vector_double_minus( h->mlmc_b1, h->rademacher_vector, h->mlmc_b1, start, end, l );

      vector_double_copy( p->b, h->mlmc_b1, start, end, l );
    }
    
    // SECOND TERM

    {
      apply_solver_double( l, threading );
    }

    // subtract the results and perform dot product
    {
      int start, end;
      gmres_double_struct* p = get_p_struct_double( l );
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      return global_inner_product_double( h->mlmc_b1, p->x, p->v_start, p->v_end, l, threading );   
    }
  }

  
  
  // the term tr( A_{l}^{-1} - P A_{l+1}^{-1} R )
  complex_double hutchinson_mlmc_difference( level_struct *l, hutchinson_double_struct* h, struct Thread *threading ){

    // FIRST TERM : result stored in p->x
    // apply A_{l}^{-1}

    { 
      int start, end;
      gmres_double_struct* p = get_p_struct_double( l );
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      vector_double_copy( p->b, h->rademacher_vector, start, end, l );
      // solution of this solve is in l->p_double.x
      apply_solver_double( l, threading );
    }
    // SECOND TERM : result stored in h->mlmc_b2
    // 1. Restrict
    // 2. invert
    // 3. Prolongate
    {
      gmres_double_struct* p = get_p_struct_double( l );
      
      apply_R_double( l->next_level->p_double.b, h->rademacher_vector, l, threading );
      // the input of this solve is l->next_level->p_double.x, the output l->next_level->p_double.b
      apply_solver_double( l->next_level, threading );
      apply_P_double( h->mlmc_b2, l->next_level->p_double.x, l, threading );
    }
  
    // subtract the results and perform dot product
    {
      int start, end;
      gmres_double_struct* p = get_p_struct_double( l);
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      vector_double_minus( h->mlmc_b1, p->x, h->mlmc_b2, start, end, l); 

      //if(l->depth ==0 && 1==0){  //Deflate from the left. (I-VV^*)Op
        hutchinson_deflate_vector_double(h->rademacher_vector, l, threading);
      //}
      return global_inner_product_double( h->rademacher_vector, h->mlmc_b1, p->v_start, p->v_end, l, threading );   
    }
  }



  // the term tr( R A_{l}^{-1} P - A_{l+1}^{-1} )
  complex_double hutchinson_split_intermediate( level_struct *l, hutchinson_double_struct* h, struct Thread *threading ){

    // FIRST TERM : result stored in h->mlmc_b1

    // 1. prolongate
    // 2. invert
    // 3. restrict
    {
      gmres_double_struct* p = get_p_struct_double( l );

      apply_P_double( p->b, h->rademacher_vector, l, threading );
      // the input of this solve is p->x, the output p->b
      apply_solver_double( l, threading );
      apply_R_double( h->mlmc_b1, p->x, l, threading );
    }

    // SECOND TERM : result stored in h->mlmc_b2

    // apply A_{l+1}^{-1}
    {
      int start, end;
      gmres_double_struct* p = get_p_struct_double( l->next_level );
      compute_core_start_end( 0, l->next_level->inner_vector_size, &start, &end, l->next_level, threading );
      vector_double_copy( p->b, h->rademacher_vector, start, end, l->next_level );
      // solution of this solve is in l->next_level->p_double.x
      apply_solver_double( l->next_level, threading );
      vector_double_copy( h->mlmc_b2, l->next_level->p_double.x, start, end, l->next_level );
    }

    // subtract the results and perform dot product
    {
      int start, end;
      gmres_double_struct* p = get_p_struct_double( l->next_level );
      compute_core_start_end( 0, l->next_level->inner_vector_size, &start, &end, l->next_level, threading );
      vector_double_minus( h->mlmc_b1, h->mlmc_b1, h->mlmc_b2, start, end, l->next_level ); // compute r = b - w
      return global_inner_product_double( h->rademacher_vector, h->mlmc_b1, p->v_start, p->v_end, l->next_level, threading );   
    }
  }
  

  void rademacher_create( level_struct *l, hutchinson_double_struct* h, int type, struct Thread *threading ){

    if( type==0 ){
      START_MASTER(threading)
      vector_double_define_random_rademacher( h->rademacher_vector, 0, l->inner_vector_size, l );
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)
    }
    else if( type==1 ){
      START_MASTER(threading)
      vector_double_define_random_rademacher( h->rademacher_vector, 0, l->next_level->inner_vector_size, l->next_level );
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)
    }
    else{ error0("Unknown value for type of Rademacher vector in relation to level of creation\n"); }
  }


  // apply the interpolation
  void apply_P_double( vector_double out, vector_double in, level_struct* l, struct Thread *threading ){

    if( l->depth==0 ){
      interpolate3_double( l->sbuf_double[0], in, l, threading );
      trans_back_double( out, l->sbuf_double[0], l->s_double.op.translation_table, l, threading );
    }
    else{
      interpolate3_double( out, in, l, threading );
    }
  }


  // apply the restriction
  void apply_R_double( vector_double out, vector_double in, level_struct* l, struct Thread *threading ){

    if( l->depth==0 ){
      trans_double( l->sbuf_double[0], in, l->s_double.op.translation_table, l, threading );     
      restrict_double( out, l->sbuf_double[0], l, threading );
    }
    else{
      restrict_double( out, in, l, threading );
    }
  }


  int apply_solver_double( level_struct* l, struct Thread *threading ){

    int nr_iters;
    double buff1, buff2;

    gmres_double_struct* p = get_p_struct_double( l );
    
    START_MASTER(threading)
    buff1 = p->tol;
    p->tol = g.tol;
    if( l->level==0 ){
      buff2 = g.coarse_tol;
      g.coarse_tol = g.tol;
    }
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)

    nr_iters = fgmres_double( p, l, threading );
    
    START_MASTER(threading);
    p->tol = buff1;
    if( l->level==0 ){
      g.coarse_tol = buff2;
    }
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)

    return nr_iters;
  }

  
  int apply_solver_powerit_double( level_struct* l, struct Thread *threading ){

    int nr_iters;
    double buff_coarsest_tol, buff_coarse_tol;
    
    gmres_double_struct* p = get_p_struct_double( l );
    
    START_MASTER(threading)
    buff_coarse_tol = p->tol;
    p->tol = l->powerit.tol_buff;
    if( l->level==0 ){
      buff_coarsest_tol = g.coarse_tol;
      g.coarse_tol = l->powerit.tol_buff;
    }
    END_MASTER(threading)
    SYNC_CORES(threading)

    nr_iters = fgmres_double( p, l, threading );

    START_MASTER(threading)
    p->tol = buff_coarse_tol;
    if( l->level==0 ){
      g.coarse_tol = buff_coarsest_tol;
    }
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)


    return nr_iters;
  }

  
  gmres_double_struct* get_p_struct_double( level_struct* l ){

    if( l->depth==0 ){
      return &(g.p);
    }
    else{
      return &(l->p_double);
    }
  }








    // direct term
  complex_double hutchinson_deflated_direct_term(level_struct *l, struct Thread *threading){
    double td0 = MPI_Wtime();
    int i, start, end;
    complex_double direct_trace = 0.0;
    gmres_double_struct* p = get_p_struct_double( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    // 0. create small matrix to store all dot products, let's call it small_T
    complex_double small_T[l->powerit.nr_vecs];
    
    vector_double* vecs_buff1;
    vector_double* vecs_buff2;

    vecs_buff1 = NULL;
    vecs_buff2 = NULL;
  
    PUBLIC_MALLOC( vecs_buff2, complex_double*, l->powerit.nr_vecs );
    PUBLIC_MALLOC( vecs_buff1, complex_double*, l->powerit.nr_vecs );
  
    vecs_buff1[0] = NULL;
    vecs_buff2[0] = NULL;
  
    PUBLIC_MALLOC( vecs_buff1[0], complex_double, l->powerit.nr_vecs*l->vector_size );
    PUBLIC_MALLOC( vecs_buff2[0], complex_double, l->powerit.nr_vecs*l->vector_size );
  
    START_MASTER(threading)
    for( i=1;i<l->powerit.nr_vecs;i++ ){
      vecs_buff1[i] = vecs_buff1[0] + i*l->vector_size;
      vecs_buff2[i] = vecs_buff2[0] + i*l->vector_size;
    }
    END_MASTER(threading)
    SYNC_CORES(threading)



    for( i=0; i<l->powerit.nr_vecs;i++ ){

      // 1. apply the operator on the ith deflation vector
      // TODO ...
      vector_double_copy(p->b, l->powerit.vecs[i], start, end, l);  
      apply_solver_double( l, threading );
      vector_double_copy(l->powerit.vecs_buff1, p->x, start, end, l);  
      
      // 2. dot product (do only the diagonal ones)
      // TODO ...
      small_T[i] = global_inner_product_double( l->powerit.vecs[i], l->powerit.vecs_buff1, p->v_start, p->v_end, l, threading );

      // 3. take trace of small_T, store in estimate.direct_trace
      // TODO ...
      direct_trace += small_T[i];
    }
    double td1 = MPI_Wtime();
    
    return direct_trace;
  }

  
  void hutchinson_deflate_vector_double(vector_double input, level_struct *l, struct Thread *threading ){
    int start, end;
    gmres_double_struct* p = get_p_struct_double( l);
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    complex_double aux[l->powerit.nr_vecs];
    
    for( int i = 0; i<l->powerit.nr_vecs; i++ ){
      aux[i] = global_inner_product_double(l->powerit.vecs[i], input, p->v_start, p->v_end, l, threading);	
    }
      
    vector_double_scale( l->powerit.vecs_buff1 , l->powerit.vecs[0], aux[0], start, end, l);
    for( int i = 1;  i< l->powerit.nr_vecs; i++ ){
      
      vector_double_copy(l->powerit.vecs_buff3, l->powerit.vecs_buff1, start, end, l);
      vector_double_scale( l->powerit.vecs_buff2, l->powerit.vecs[i], aux[i], start, end, l);
      vector_double_plus( l->powerit.vecs_buff1 , l->powerit.vecs_buff3 , l->powerit.vecs_buff2, start, end, l);

    }
    
    vector_double_minus(  input, input, l->powerit.vecs_buff1, start, end, l );
    
  }
