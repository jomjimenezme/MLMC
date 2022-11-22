
    #include "main.h"
    
    
    struct estimate {
        int counter; //required number of estimates.
        complex_double estimate;
    };    
    
    struct estimate compute_coarsest_trace(level_struct *ls, int nr_ests, double tol, gmres_double_struct *ps, hutchinson_double_struct* h, struct Thread *threading){
        
    
            
            //FIXME: set is as input parameter?
            int low_tol_restart_length = 5;
            int tmp_length = ps->restart_length;
            ps->restart_length = low_tol_restart_length;
            
            
            double buff_tol = ps->tol , RMSD = 0.0;
            complex_double es = 0.0;
            int j;
            complex_double rough_trace = h->rt;
            double est_tol = h->trace_tol;
            
        int li = g.num_levels-1;
        int start=0;
        
        START_MASTER(threading)
        printf0("starting coarsest level\n\n", li);
        END_MASTER(threading)
        SYNC_MASTER_TO_ALL(threading)
            struct estimate coarsest;
        
        START_MASTER(threading)
        if(g.my_rank==0)printf("LEVEL----------- \t %d  \t ests: %d  \t start: %d  \t end: %d \t vector size: %d \n ", li, nr_ests, start, ps->v_end, ls->inner_vector_size);
        END_MASTER(threading)
        
        complex_double sample[nr_ests]; 
        double variance = 0.0;
        int counter = 0; 
        
        double tt1c = MPI_Wtime();
        
        for(int i=0;i<nr_ests ; i++){
            START_MASTER(threading)
            vector_double_define_random_rademacher( ps->b, 0, ls->inner_vector_size, ls );    
            END_MASTER(threading)
            SYNC_MASTER_TO_ALL(threading)
            
            ps->tol = g.tol;
            g.coarse_tol = g.tol;
            double t0 = MPI_Wtime();
            
            
            //gmres_double_struct *ps = ps;
            //level_struct *ls = ls;
            
            ps->restart_length = tmp_length;
            
            int start, end, coarsest_iters;
            compute_core_start_end(ps->v_start, ps->v_end, &start, &end, ls, threading);
            #ifdef POLYPREC
            ps->preconditioner = ps->polyprec_double.preconditioner;
            #endif
            #ifdef GCRODR
            coarsest_iters = flgcrodr_double( ps, ls, threading );
            #else
            coarsest_iters = fgmres_double( ps, ls, threading );
            #endif
            START_MASTER(threading)
            g.avg_b1 += coarsest_iters;
            g.avg_b2 += 1;
            g.avg_crst = g.avg_b1/g.avg_b2;
            END_MASTER(threading)
            SYNC_MASTER_TO_ALL(threading)
            SYNC_CORES(threading)
            #ifdef POLYPREC
            if ( ls->level==0 && ls->p_double.polyprec_double.update_lejas == 1 ) {
                if ( coarsest_iters >= 1.5*ps->polyprec_double.d_poly ) {
                    // re-construct Lejas
                    re_construct_lejas_double( ls, threading );
                }
            }
            #endif
                        
            
            double t1 = MPI_Wtime();
            START_MASTER(threading)
            printf0("coarsest iters = %d (time = %f)\n", coarsest_iters, t1-t0);
            END_MASTER(threading)
            ps->tol = buff_tol;
            g.coarse_tol = buff_tol;
            complex_double tmpx= global_inner_product_double( ps->b, ps->x, ps->v_start, ps->v_end, ls, threading );   
            sample[i]= tmpx;      
            es += tmpx;
            
            variance = 0.0;
            for (j=0; j<i; j++) // compute the variance
                variance += (sample[j]- es/(i+1)) *conj( sample[j]- es/(i+1)); 
            variance /= (j+1);
            RMSD = sqrt(variance/(j+1)); //RMSD= sqrt(var+ bias²)
            
            START_MASTER(threading)
            if(g.my_rank==0)printf( "%d \t var %.15f \t Est %.15f  \t RMSD %.15f <  %.15f ?? \n ", i, creal(variance), creal( es/(i+1) ), RMSD,  cabs(rough_trace)*tol*est_tol );
            END_MASTER(threading)
            SYNC_MASTER_TO_ALL(threading)
            counter=i+1;
            if(i !=0 && RMSD< cabs(rough_trace)*tol*est_tol && i>=h->min_iters-1 ){counter=i+1; break;}
        }
        
        double tt2c = MPI_Wtime();
        
        START_MASTER(threading)
        printf0("end of coarsest (total time = %f)\n\n", li, tt2c-tt1c);
        END_MASTER(threading)
        
        coarsest.counter=counter;
        coarsest.estimate=es;
        return coarsest;
        
    }
    
    
    
    void hutchinson_diver_double_init( level_struct *l, struct Thread *threading ) {
        hutchinson_double_struct* h = &(l->h_double);
        
        //MLMC
        h->mlmc_b1 = NULL;
        
        
        SYNC_MASTER_TO_ALL(threading)
    }
    
    void hutchinson_diver_double_alloc( level_struct *l, struct Thread *threading ) {
        hutchinson_double_struct* h = &(l->h_double) ;
        
        //For MLMC
        PUBLIC_MALLOC( h->mlmc_b1, complex_double, l->inner_vector_size );   
        
    }
    
    void hutchinson_diver_double_free( level_struct *l, struct Thread *threading ) {
        hutchinson_double_struct* h = &(l->h_double) ;
        
        PUBLIC_FREE( h->mlmc_b1, complex_double, l->inner_vector_size );   
    }
    
    
    
    
    complex_double mlmc_hutchinson_diver_double( level_struct *l, struct Thread *threading ) {
        
        hutchinson_double_struct* h = &(l->h_double) ;
        complex_double rough_trace = h->rt;
        
        // FIXME : make this int a sort of input parameter ?
        int low_tol_restart_length = 5;
        
        START_MASTER(threading)
        g.avg_b1 = 0.0;     g.avg_b2 = 0.0;     g.avg_crst = 0.0;
        END_MASTER(threading)
        
        complex_double cov[12*12]; 
        double RMSD;
        
        complex_double trace=0.0; 
        //--------------Setting the tolerance per level----------------
        int i;
        double est_tol = l->h_double.trace_tol;//  1E-6;
        double delta[g.num_levels];
        double tol[g.num_levels];
        double d0 = 0.75;
        delta[0] = d0; tol[0] = sqrt(delta[0]);
        delta[1] = 0.20;    delta[2] = 0.04;    delta[3] = 0.01;
        tol[1] = sqrt(delta[1]);    tol[2] = sqrt(delta[2]);    tol[3] = sqrt(delta[3]);
        
        //-------------------------------------------------------------  
        
        
        //-----------------Setting variables in levels-----------------  
        int start, end, j, li;
        int nr_ests[g.num_levels];
        int counter[g.num_levels];
        
        nr_ests[0]=l->h_double.max_iters; 
        
        for(i=1; i< g.num_levels; i++) nr_ests[i] = nr_ests[i-1]*1; //TODO: Change in every level?
        
        gmres_double_struct *ps[g.num_levels]; //float p_structs for coarser levels
        level_struct *ls[g.num_levels];  
        ps[0] = &(g.p);
        ls[0] = l;
        
        level_struct *lx = l->next_level;
        double buff_tol[g.num_levels];
        buff_tol[0] =  ps[0]->tol; 
        START_MASTER(threading)
        if(g.my_rank==0) printf("LEVEL----------- \t %d  \t vector size: %d \t TOLERANCE: %e\n", 0, ls[0]->inner_vector_size, buff_tol[0]);
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
        
        complex_double variance = 0.0;
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
            // if(g.my_rank==0) printf("LEVEL----------- \t %d  \t ests: %d  \t start: %d  \t end: %d \t vector size: %d \n ", li, nr_ests[li], start, ps[li]->v_end, ls[li]->inner_vector_size);
            END_MASTER(threading)
            
            START_MASTER(threading)
            printf0("starting level difference = %d\n\n", li);
            END_MASTER(threading)
            
            for(i=0;i<nr_ests[li]  ; i++){
                //---------------------Get random vector-----------------------------------------------------------
                START_MASTER(threading) 
                vector_double_define_random_rademacher( ps[li]->b, 0, ls[li]->inner_vector_size, ls[li] );          
                END_MASTER(threading)
                SYNC_MASTER_TO_ALL(threading)
                
                
                //-----------------Solve the system in current level and the next one----------------- 
                if(li==0){                    
                    trans_double( l->sbuf_double[0], ps[li]->b, l->s_double.op.translation_table, l, threading );     
                    restrict_double( ps[li+1]->b, l->sbuf_double[0], ls[li], threading ); // get rhs for next level.
                    ps[li+1]->tol = g.tol;
                }
                else{
                    restrict_double( ps[li+1]->b, ps[li]->b, ls[li], threading ); // get rhs for next level.
                    ps[li+1]->tol = g.tol;
                    g.coarse_tol = g.tol;
                }
                double t0 = MPI_Wtime();
                level_struct *lx = ls[li+1];
                int coarsest_iters;
                
                if( (li+2)==g.num_levels ){
                    
                    ps[li+1]->restart_length = tmp_length;
                    
                    #ifdef POLYPREC
                    ps[li+1]->preconditioner = ps[li+1]->polyprec_double.preconditioner;
                    #endif
                    #ifdef GCRODR
                    coarsest_iters = flgcrodr_double( ps[li+1], lx, threading );
                    #else
                    coarsest_iters = fgmres_double( ps[li+1], lx, threading );
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
                        if ( coarsest_iters >= 1.5*ps[li+1]->polyprec_double.d_poly ) {
                            // re-construct Lejas
                            re_construct_lejas_double( lx, threading );
                        }
                    }
                    #endif
                    
                    ps[li+1]->restart_length = low_tol_restart_length;
                    
                } else {
                    coarsest_iters = fgmres_double( ps[li+1], lx, threading );
                }
                
                double t1 = MPI_Wtime();
                START_MASTER(threading)
                if( (li+1)==3 ){
                    printf0("coarsest iters = %d (time = %f)\n", coarsest_iters, t1-t0);
                } else {
                    printf0("coarse or fine iters = %d (time = %f)\n", coarsest_iters, t1-t0);
                }
                END_MASTER(threading)
                ps[li+1]->tol =buff_tol[li+1] ;
                
                
                //--------------------------------------Interpolate Solution----------------------------
                if(li==0){
                    interpolate3_double( l->sbuf_double[1], ps[li+1]->x, ls[li], threading ); //      
                    trans_back_double( h->mlmc_b1, l->sbuf_double[1], l->s_double.op.translation_table, l, threading );   
                }
                else{
                    interpolate3_double( h->mlmc_b1, ps[li+1]->x, ls[li], threading ); //   
                    g.coarse_tol = buff_tol[li+1];
                    ps[li]->tol = g.tol;  
                }
                //------------------------------------------------------------------------------------
                
                
                //--------------------------------------Solve Finer level----------------------------
                t0 = MPI_Wtime();
                
                int fine_iters = fgmres_double( ps[li], ls[li], threading );
                
                t1 = MPI_Wtime();
                START_MASTER(threading)
                printf0("fine iters = %d (time = %f)\n", fine_iters, t1-t0);
                END_MASTER(threading)
                //------------------------------------------------------------------------------------
                
                //--------------Get the difference of the solutions and the corresponding sample--------------                        
                vector_double_minus( h->mlmc_b1, ps[li]->x, h->mlmc_b1, start, end, ls[li] );
                tmpx =  global_inner_product_double( ps[li]->b, h->mlmc_b1, ps[li]->v_start, ps[li]->v_end, ls[li], threading );      
                sample[i]= tmpx;      
                es[li] +=  tmpx;
                //----------------------------------------------------------------------------------------------
                
                //--------------Get the Variance in the current level and use it in stop criteria---------------  
                variance = 0.0;
                for (j=0; j<i; j++) // compute the variance
                    variance += (sample[j]- es[li]/(i+1)) *conj( sample[j]- es[li]/(i+1)); 
                variance /= (j+1);
                RMSD = sqrt(variance/(j+1)); //RMSD= sqrt(var+ bias²)
                
                START_MASTER(threading)
                if(g.my_rank==0)printf( "%d \t var %.15f \t Est %.15f  \t RMSD %.15f <  %.15f ?? \n ", i, creal(variance), creal( es[li]/(i+1) ), RMSD,  cabs(rough_trace)*tol[li]*est_tol );
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
        
        
        //----------------------Hutchinson for just coarsest level-----------------               
        struct estimate coarsest;
        coarsest = compute_coarsest_trace(ls[li], nr_ests[li], tol[li], ps[li], h, threading);
        es[li]= coarsest.estimate;
        counter[li] = coarsest.counter;
        //--------------------------------------------------------------------------       

        
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
    
    
    
    
    
    
    
    
    
