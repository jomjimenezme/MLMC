
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
            
            /*double normx1 = global_norm_double( ps->b,start,end,ls,threading );
            apply_operator_double( ps->w,ps->x,ps,ls,threading ); // compute w = D*x
            vector_double_minus( ps->r,ps->b,ps->w,start,end,ls ); // compute r = b - w
            double normx2 = global_norm_double( ps->r,start,end,ls,threading );
            if(g.my_rank==0)printf("\n\t\t REL RESIDUAL  Coarsest: %.16f \t LEVEL= %d \n", normx2/normx1, li);*/
            
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
            counter+=1;
            if(i !=0 && RMSD< cabs(rough_trace)*tol*est_tol && i>=h->min_iters-1 ){counter+=1; break;}
        }
        
        double tt2c = MPI_Wtime();
        
        START_MASTER(threading)
        printf0("end of coarsest (total time = %f)\n\n", li, tt2c-tt1c);
        END_MASTER(threading)
        
        coarsest.counter=counter;
        coarsest.estimate=es;
        return coarsest;
        
    }
    
    
    struct estimate compute_difference_level_trace(level_struct *ls, int nr_ests, double tol, gmres_double_struct *ps, hutchinson_double_struct* h, struct Thread *threading){
        
        struct estimate difference;
        complex_double rough_trace = h->rt;
        double est_tol = h->trace_tol;//  1E-6;
        int li = ls->depth;
        complex_double sample[nr_ests]; 
        complex_double  es=0.0;
        double var =0.0, RMSD=0.0;
        int start, end, i, j, coarsest_iters;
        
        complex_double variance=0.0;    int counter=0;
        compute_core_start_end( 0, ls->inner_vector_size, &start, &end, ls, threading );
        //   for( int li=0;li<(g.num_levels-1);li++ ){
        
        level_struct *lx = ls->next_level;
        gmres_double_struct *px = &(lx->p_double);
        
        int tmp_length =h->tmp_length;
        int low_tol_restart_length = h-> low_tol_restart_length;
        
        complex_double buff_tol_coarser= ps->tol;
        
        double tt0 = MPI_Wtime();
        
        
        
        START_MASTER(threading)
        // if(g.my_rank==0) printf("LEVEL----------- \t %d  \t ests: %d  \t start: %d  \t end: %d \t vector size: %d \n ", li, nr_ests, start, ps->v_end, ls->inner_vector_size);
        END_MASTER(threading)
        
        START_MASTER(threading)
        printf0("starting level difference = %d\n\n", li);
        END_MASTER(threading)
        
        for(i=0;i<nr_ests  ; i++){
            //---------------------Get random vector-----------------------------------------------------------
            START_MASTER(threading) 
            vector_double_define_random_rademacher( ps->b, 0, ls->inner_vector_size, ls );          
            END_MASTER(threading)
            SYNC_MASTER_TO_ALL(threading)
            
            
            
            //-----------------Solve the system in current level and the next one----------------- 
            if(li==0){                    
                trans_double( ls->sbuf_double[0], ps->b, ls->s_double.op.translation_table, ls, threading );     
                restrict_double( px->b, ls->sbuf_double[0], ls, threading ); // get rhs for next level.
                px->tol = g.tol;
            }
            else{
                restrict_double( px->b, ps->b, ls, threading ); // get rhs for next level.
                px ->tol = g.tol;
                g.coarse_tol = g.tol;
            }
            
            double t0 = MPI_Wtime();
            
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
                printf0("coarse or fine iters = %d (time = %f)\n", coarsest_iters, t1-t0);
            }
            END_MASTER(threading)
            px->tol =buff_tol_coarser ;
            
            
            //--------------------------------------Interpolate Solution----------------------------
            if(li==0){
                interpolate3_double( ls->sbuf_double[1], px->x, ls, threading ); //      
                trans_back_double( h->mlmc_b1, ls->sbuf_double[1], ls->s_double.op.translation_table, ls, threading );   
            }
            else{
                interpolate3_double( h->mlmc_b1, px->x, ls, threading ); //   
                g.coarse_tol = buff_tol_coarser;
                ps->tol = g.tol;  
            }
            //------------------------------------------------------------------------------------
            
            
            //--------------------------------------Solve Finer level----------------------------
            t0 = MPI_Wtime();
            
            int fine_iters = fgmres_double( ps, ls, threading );
            
            t1 = MPI_Wtime();
            START_MASTER(threading)
            printf0("fine iters = %d (time = %f)\n", fine_iters, t1-t0);
            END_MASTER(threading)
            //------------------------------------------------------------------------------------
            
            //--------------Get the difference of the solutions and the corresponding sample--------------                        
            vector_double_minus( h->mlmc_b1, ps->x, h->mlmc_b1, start, end, ls );
            complex_double tmpx =  global_inner_product_double( ps->b, h->mlmc_b1, ps->v_start, ps->v_end, ls, threading );      
            sample[i]= tmpx;      
            es +=  tmpx;
            //----------------------------------------------------------------------------------------------
            
            //--------------Get the Variance in the current level and use it in stop criteria---------------  
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
            if(i !=0 && RMSD< cabs(rough_trace)*tol*est_tol && i>=h->min_iters-1){counter=i+1; break;}
            //----------------------------------------------------------------------------------------------    
        }
        
        double tt1 = MPI_Wtime();
        
        START_MASTER(threading)
        printf0("ending difference = %d (total time = %f)\n\n", li, tt1-tt0);
        END_MASTER(threading)
        
        difference.counter=counter;
        difference.estimate=es;
        return difference;
    }
    
    struct estimate compute_orthogonal_trace(level_struct *ls, int nr_ests, double tol, gmres_double_struct *ps, hutchinson_double_struct* h, struct Thread *threading){
        
        struct estimate difference;
        complex_double rough_trace = h->rt;
        double est_tol = h->trace_tol;//  
        int li = ls->depth;
        complex_double sample[nr_ests]; 
        complex_double  es=0.0;
        double var =0.0, RMSD=0.0;
        int start, end, i, j, iters;
        
        complex_double variance=0.0;    int counter=0;
        compute_core_start_end( 0, ls->inner_vector_size, &start, &end, ls, threading );
        
        level_struct *lx = ls->next_level;
        gmres_double_struct *px = &(lx->p_double);
        
        int tmp_length =h->tmp_length;
        int low_tol_restart_length = h-> low_tol_restart_length;
        
        
        double tt0 = MPI_Wtime();
        
        
        
        START_MASTER(threading)
        // if(g.my_rank==0) printf("LEVEL----------- \t %d  \t ests: %d  \t start: %d  \t end: %d \t vector size: %d \n ", li, nr_ests, start, ps->v_end, ls->inner_vector_size);
        END_MASTER(threading)
        
        START_MASTER(threading)
        printf0("starting orthogonal term at level  = %d\n\n", li);
        END_MASTER(threading)
        
        for(i=0;i<nr_ests  ; i++){
            //---------------------Get random vector-----------------------------------------------------------
            START_MASTER(threading) 
            vector_double_define_random_rademacher( ps->b, 0, ls->inner_vector_size, ls );          
            END_MASTER(threading)
            SYNC_MASTER_TO_ALL(threading)
            
            vector_double_copy(h->mlmc_b2, ps->b, start, end, ls);
            //-----------------Get projected Rademacher vector----------------- 
            if(li==0){                    
                trans_double( ls->sbuf_double[0], ps->b, ls->s_double.op.translation_table, ls, threading );     
                restrict_double( px->b, ls->sbuf_double[0], ls, threading ); // Restrict.
                
                interpolate3_double( ls->sbuf_double[1], px->b, ls, threading ); // Interpolate  
                trans_back_double( h->mlmc_b1, ls->sbuf_double[1], ls->s_double.op.translation_table, ls, threading );   
                
                
                vector_double_minus( ps->b, ps->b, h->mlmc_b1, start, end, ls );
            }
            else{
                restrict_double( px->b, ps->b, ls, threading ); // get rhs for next level.     
                interpolate3_double( h->mlmc_b1, px->b, ls, threading ); //   
                vector_double_minus( ps->b, ps->b, h->mlmc_b1, start, end, ls );
            }
            
            double t0 = MPI_Wtime();
            
            iters = fgmres_double( ps, ls, threading );
            
            /*double normx1 = global_norm_double( ps->b,start,end,ls,threading );
            apply_operator_double( ps->w,ps->x,ps,ls,threading ); // compute w = D*x
            vector_double_minus( ps->r,ps->b,ps->w,start,end,ls ); // compute r = b - w
            double normx2 = global_norm_double( ps->r,start,end,ls,threading );
            if(g.my_rank==0)printf("\n\t\t REL RESIDUAL  finer: %.16f \t LEVEL= %d \n", normx2/normx1, li);*/
            
            double t1 = MPI_Wtime();
            START_MASTER(threading)
                printf0("iters = %d (time = %f)\n", iters, t1-t0);
            END_MASTER(threading)
            SYNC_MASTER_TO_ALL(threading)


            //--------------Get the sample--------------                        
            complex_double tmpx =  global_inner_product_double( h->mlmc_b2, ps->x, ps->v_start, ps->v_end, ls, threading );      
            sample[i]= tmpx;      
            es +=  tmpx;
            //----------------------------------------------------------------------------------------------
            
            //--------------Get the Variance in the current level and use it in stop criteria---------------  
            variance = 0.0;
            for (j=0; j<i; j++) // compute the variance
                variance += (sample[j]- es/(i+1)) *conj( sample[j]- es/(i+1)); 
            variance /= (j+1);
            RMSD = sqrt(variance/(j+1)); //RMSD= sqrt(var+ bias²)
            
            START_MASTER(threading)
            if(g.my_rank==0)printf( "%d \t var %.15f \t Est %.15f  \t RMSD %.15f <  %.15f ?? \n ", i, creal(variance), creal( es/(i+1) ), RMSD,  cabs(rough_trace)*tol*est_tol );
            END_MASTER(threading)
            SYNC_MASTER_TO_ALL(threading)     
            counter+=1;
            if(i !=0 && RMSD< cabs(rough_trace)*tol*est_tol && i>=h->min_iters-1){counter+=1; break;}
            //----------------------------------------------------------------------------------------------    
        }
        
        double tt1 = MPI_Wtime();
        
        START_MASTER(threading)
        printf0("ending difference = %d (total time = %f)\n\n", li, tt1-tt0);
        END_MASTER(threading)
        
        difference.counter=counter;
        difference.estimate=es;
        return difference;
    }
    
    struct estimate compute_intermediate_trace(level_struct *ls, int nr_ests, double tol, gmres_double_struct *ps, hutchinson_double_struct* h, struct Thread *threading){
        
        struct estimate difference;
        complex_double rough_trace = h->rt;
        double est_tol = h->trace_tol*0.01;// 
        int li = ls->depth;
        complex_double sample[nr_ests]; 
        complex_double  es=0.0;
        double var =0.0, RMSD=0.0;
        int start, end, start_next, end_next, i, j, finer_iters, coarser_iters;
        
        complex_double variance=0.0;    int counter=0;
        compute_core_start_end( 0, ls->inner_vector_size, &start, &end, ls, threading );
        //   for( int li=0;li<(g.num_levels-1);li++ ){
        
        level_struct *lx = ls->next_level;
        gmres_double_struct *px = &(lx->p_double);
        compute_core_start_end( 0, lx->inner_vector_size, &start_next, &end_next, lx, threading );
        
       
        
       
        complex_double buff_tol_coarser= px->tol;
        
        
        double tt0 = MPI_Wtime();
        
        
        
        START_MASTER(threading)
        // if(g.my_rank==0) printf("LEVEL----------- \t %d  \t ests: %d  \t start: %d  \t end: %d \t vector size: %d \n ", li, nr_ests, start, ps->v_end, ls->inner_vector_size);
        END_MASTER(threading)
        
        START_MASTER(threading)
        printf0("starting level difference = %d\n\n", li);
        END_MASTER(threading)
        
        for(i=0;i<nr_ests  ; i++){
            ps->tol = g.tol;
            //---------------------Get random vector in the next level-----------------------------------------------------------
            START_MASTER(threading) 
            vector_double_define_random_rademacher( px->b, 0, lx->inner_vector_size, lx );          
            END_MASTER(threading)
            SYNC_MASTER_TO_ALL(threading)
            
            
            
            //-----------------Prolongate the rhs----------------- 
            if(li==0){              
                interpolate3_double( ls->sbuf_double[0], px->b, ls, threading ); //      
                trans_back_double( ps->b, ls->sbuf_double[0], ls->s_double.op.translation_table, ls, threading );   
                //px->tol = g.tol;
            }
            else{
                interpolate3_double( ps->b, px->b, ls, threading );
                /*px ->tol = g.tol;
                g.coarse_tol = g.tol;*/
            }
            
            double t0 = MPI_Wtime();
            
            //if( (li+2)==g.num_levels )
            //else
           
            //if(g.my_rank==0)printf("\n\t\tps->tol=  : %.16f \t ps->restart_lenght = %d", ps->tol, ps->restart_length);

            finer_iters = fgmres_double( ps, ls, threading );
            
            //// checking relative residual of solves
            /*double normx1 = global_norm_double( ps->b,start,end,ls,threading );
            apply_operator_double( ps->w,ps->x,ps,ls,threading ); // compute w = D*x
            vector_double_minus( ps->r,ps->b,ps->w,start,end,ls ); // compute r = b - w
            double normx2 = global_norm_double( ps->r,start,end,ls,threading );
            if(g.my_rank==0)printf("\n\t\t REL RESIDUAL  finer: %.16f \t LEVEL= %d \n", normx2/normx1, li);
            */
            double t1 = MPI_Wtime();
            START_MASTER(threading)
            
                printf0("finer iters = %d (time = %f)\n", coarser_iters, t1-t0);
        
            END_MASTER(threading)
            px->tol =buff_tol_coarser ;
            
            
            //--------------------------------------Restrict Solution----------------------------
            if(li==0){
                trans_double( ls->sbuf_double[1], ps->x, ls->s_double.op.translation_table, ls, threading );     
                restrict_double( h->mlmc_b1, ls->sbuf_double[1], ls, threading ); // get rhs for next level.                
            }
            else{
                restrict_double( h->mlmc_b1, ps->x, ls, threading ); // get rhs for next level.

            }
            //------------------------------------------------------------------------------------
            
            
            //--------------------------------------Solve Coarser level----------------------------
            t0 = MPI_Wtime();
            
            
            px ->tol = g.tol;
            g.coarse_tol = g.tol;
           // if(g.my_rank==0)printf("\n\t\tpx->tol=  : %.16f \t px->restart_lenght = %d ", px->tol, px->restart_length);
            
            int coarser_iters = fgmres_double( px, lx, threading );
            px->tol =buff_tol_coarser;
            g.coarse_tol = buff_tol_coarser;
            
           /* double normx3 = global_norm_double( px->b, start_next ,end_next ,lx ,threading );
            apply_operator_double( px->w, px->x, px, lx, threading ); // compute w = D*x
            vector_double_minus( px->r, px->b, px->w, start_next, end_next,lx ); // compute r = b - w
            double normx4 = global_norm_double( px->r, start_next, end_next, lx, threading );
            if(g.my_rank==0)printf("\n\t\tREL RESIDUAL  coarser: %.16f\n\n", normx4/normx3);*/
            
            t1 = MPI_Wtime();
            START_MASTER(threading)
            printf0("coarser iters = %d (time = %f)\n", coarser_iters, t1-t0);
            END_MASTER(threading)
            //------------------------------------------------------------------------------------
            
            //--------------Get the difference of the solutions and the corresponding sample--------------                        
            vector_double_minus( h->mlmc_b1, h->mlmc_b1, px->x, start_next, end_next, lx );
            complex_double tmpx =  global_inner_product_double( px->b, h->mlmc_b1, px->v_start, px->v_end, lx, threading );      
            sample[i]= tmpx;      
            es +=  tmpx;
            //----------------------------------------------------------------------------------------------
            
            //--------------Get the Variance in the current level and use it in stop criteria---------------  
            variance = 0.0;
            for (j=0; j<i; j++) // compute the variance
                variance += (sample[j]- es/(i+1)) *conj( sample[j]- es/(i+1)); 
            variance /= (j+1);
            RMSD = sqrt(variance/(j+1)); //RMSD= sqrt(var+ bias²)
            
            START_MASTER(threading)
            if(g.my_rank==0)printf( "%d \t var %.15f \t Est %.15f  \t RMSD %.15f <  %.15f ?? \n ", i, creal(variance), creal( es/(i+1) ), RMSD,  cabs(rough_trace)*tol*est_tol );
            END_MASTER(threading)
            SYNC_MASTER_TO_ALL(threading)     
            counter+=1;
            if(i !=0 && RMSD< cabs(rough_trace)*tol*est_tol && i>=h->min_iters-1){counter+=i; break;}
            //----------------------------------------------------------------------------------------------    
        }
        
        double tt1 = MPI_Wtime();
        
        START_MASTER(threading)
        printf0("ending difference = %d (total time = %f)\n\n", li, tt1-tt0);
        END_MASTER(threading)
        
        difference.counter=counter;
        difference.estimate=es;
        return difference;
    }
    
    
    void hutchinson_diver_double_init( level_struct *l, struct Thread *threading ) {
        hutchinson_double_struct* h = &(l->h_double);
        
        //MLMC
        h->mlmc_b1 = NULL;
        
        h->mlmc_b2 =NULL;
        
        SYNC_MASTER_TO_ALL(threading)
    }
    
    void hutchinson_diver_double_alloc( level_struct *l, struct Thread *threading ) {
        hutchinson_double_struct* h = &(l->h_double) ;
        
        //For MLMC
        PUBLIC_MALLOC( h->mlmc_b1, complex_double, l->inner_vector_size );   
        PUBLIC_MALLOC( h->mlmc_b2, complex_double, l->inner_vector_size );   
        
    }
    
    void hutchinson_diver_double_free( level_struct *l, struct Thread *threading ) {
        hutchinson_double_struct* h = &(l->h_double) ;
        
        PUBLIC_FREE( h->mlmc_b1, complex_double, l->inner_vector_size );   
    }
    
    complex_double mlmc_hutchinson_diver_double( level_struct *l, struct Thread *threading ) {
        
        hutchinson_double_struct* h = &(l->h_double) ;
        
        // FIXME : make this int a sort of input parameter ?
        h->low_tol_restart_length =5;
        int low_tol_restart_length = h-> low_tol_restart_length;
        
        
        
        START_MASTER(threading)
        g.avg_b1 = 0.0;     g.avg_b2 = 0.0;     g.avg_crst = 0.0;
        END_MASTER(threading)
        
        
        complex_double trace=0.0; 
        //--------------Setting the tolerance per level----------------
        int i;
        double delta[g.num_levels];
        double tol[g.num_levels];
        double d0 = 0.75;
        delta[0] = d0; tol[0] = sqrt(delta[0]);
        delta[1] = 0.20;    delta[2] = 0.04;    delta[3] = 0.01;
        tol[1] = sqrt(delta[1]);    tol[2] = sqrt(delta[2]);    tol[3] = sqrt(delta[3]);
        
        //-------------------------------------------------------------  
        
        
        //-----------------Setting variables in levels-----------------  
        int li;
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
        h->tmp_length = ps[g.num_levels-1]->restart_length;
        ps[g.num_levels-1]->restart_length = low_tol_restart_length;
        
        complex_double es[g.num_levels];  //An estimate per level
        memset( es,0,g.num_levels*sizeof(complex_double) );
        //---------------------------------------------------------------------  
        //-----------------Hutchinson for all but coarsest level-----------------    
        for( li=0;li<(g.num_levels-1);li++ ){
            struct estimate difference = compute_difference_level_trace(ls[li], nr_ests[li], tol[li], ps[li], h, threading);
            es[li]= difference.estimate;
            counter[li] = difference.counter;    
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
            printf( "%.15f + i %.15f \t %d \t %d \t %d \n", CSPLIT( trace ), counter[0],counter[1],counter[2]);
        END_MASTER(threading)
        
        
        return trace;
        
    }
    
    
    
    complex_double split_mlmc_hutchinson_diver_double( level_struct *l, struct Thread *threading ) {
        
        hutchinson_double_struct* h = &(l->h_double) ;
        
        // FIXME : make this int a sort of input parameter ?
        h->low_tol_restart_length =5;
        int low_tol_restart_length = h-> low_tol_restart_length;
        
        
        
        START_MASTER(threading)
        g.avg_b1 = 0.0;     g.avg_b2 = 0.0;     g.avg_crst = 0.0;
        END_MASTER(threading)
        
        
        complex_double trace=0.0; 
        //--------------Setting the tolerance per level----------------
        int i;
        double delta[g.num_levels];
        double tol[g.num_levels];
        double d0 = 0.75;
        delta[0] = d0; tol[0] = sqrt(delta[0]);
        delta[1] = 0.20;    delta[2] = 0.04;    delta[3] = 0.01;
        tol[1] = sqrt(delta[1]);    tol[2] = sqrt(delta[2]);    tol[3] = sqrt(delta[3]);
        
        //-------------------------------------------------------------  
        
        
        //-----------------Setting variables in levels-----------------  
        int li;
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
        h->tmp_length = ps[g.num_levels-1]->restart_length;
        ps[g.num_levels-1]->restart_length = low_tol_restart_length;
        
        complex_double es[g.num_levels];  //An estimate per level
        memset( es,0,g.num_levels*sizeof(complex_double) );
        //---------------------------------------------------------------------  
        trace = 0.0;
        
        //-----------------Hutchinson for the orthogonally deflated term-----------------    
        for( li=0;li<(g.num_levels-1);li++ ){
            struct estimate orthogonal = compute_orthogonal_trace(ls[li], nr_ests[li], tol[li], ps[li], h, threading);
            /*es[li]+= orthogonal.estimate;
            counter[li] += orthogonal.counter;    */
            trace += orthogonal.estimate/orthogonal.counter;
        }
        //---------------------------------------------------------------------  
        
        //-----------------Hutchinson for the intermediate term-----------------    
        
        for( li=0;li<(g.num_levels-1);li++ ){
               struct estimate intermediate = compute_intermediate_trace(ls[li], nr_ests[li], tol[li], ps[li], h, threading);
                    /*es[li]+= intermediate.estimate;
                    counter[li] = intermediate.counter;    */
                    trace += intermediate.estimate/intermediate.counter;

        }
        //---------------------------------------------------------------------  
        
        
        //----------------------Hutchinson for coarsest level-----------------               
        li = g.num_levels-1;
        struct estimate coarsest;
        coarsest = compute_coarsest_trace(ls[li], nr_ests[li], tol[li], ps[li], h, threading);
        /*es[li]= coarsest.estimate;
        counter[li] = coarsest.counter;*/
        trace += coarsest.estimate/coarsest.counter;

        //--------------------------------------------------------------------------       
        
        
      /*  for( i=0;i<g.num_levels;i++ ){
            trace += es[i]/counter[i];    
            START_MASTER(threading)
            //if(g.my_rank==0)printf("Level:%d, ................................%f    %d\n", i, creal(es[i]),counter[i]);
            END_MASTER(threading)  
            
        }*/
        
        START_MASTER(threading)
        if(g.my_rank==0)
            //printf( "\n\n Trace: %.15f+i %.15f, \t ests in l_0: %d, \t ests in l_1: %d \n", CSPLIT( trace ), counter[0],counter[1],counter[2]);
            printf( "%.15f + i %.15f \t %d \t %d \t %d \n", CSPLIT( trace ), counter[0],counter[1],counter[2]);
        END_MASTER(threading)
        
        
        return trace;
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
