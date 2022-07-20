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

#ifdef POLYPREC


/*-----------------------------------------------*/

void print_matrix_double(complex_double* A, int mv, int mh )
{
  int i,j;

  // printf("\n\n");
  // for (i=0; i < mv; i++)
  // {
  //     for(j=0; j < mh; j++)
  //     {
  //             fprintf(stdout, "%6.6f +i%6.6f\t", creal(A[i*mh + j]), cimag(A[i*mh+j]));
  //     }
  //     fprintf(stdout, "\n");
  // }
  // printf("--\n");

  printf("\n\n");
  for (i=0; i < mh; i++)
  {
    for(j=0; j < mv; j++)
    {
      fprintf(stdout, "%6.6f +i%6.6f\t", creal(A[j*mh + i]), cimag(A[j*mh+i]));
    }
    fprintf(stdout, "\n");
  }
  printf("--\n");
  printf("\n\n");
}


void print_vector_double( char* desc, vector_double w, int n)
{
  int j;
  printf0( "\n %s\n", desc );
  for( j = 0; j < n; j++ ) printf0( " (%6.6f,%6.6f)", creal(w[j]), cimag(w[j]) );
  printf0( "\n" );
}


void harmonic_ritz_double( gmres_double_struct *p )
{
  int i, j, d;
  complex_double h_dd;

  d = p->polyprec_double.d_poly;
  h_dd = p->polyprec_double.Hc[d-1][d];
  memset(p->polyprec_double.dirctslvr.b, 0.0, sizeof(complex_double)*(d-1));
  p->polyprec_double.dirctslvr.b[d-1] = 1.;

  for (i=0; i<d; i++)
    for (j=0; j<d; j++)
      p->polyprec_double.Hcc[i*d + j ] = conj(p->polyprec_double.Hc[j][i]);

  p->polyprec_double.dirctslvr.dirctslvr_double(&p->polyprec_double.dirctslvr);

  for (i=0; i<d; i++)
    p->polyprec_double.Hc[d-1][i] += h_dd*h_dd*p->polyprec_double.dirctslvr.x[i];
    
  p->polyprec_double.eigslvr.eigslvr_double(&p->polyprec_double.eigslvr);
}


/*-----------------------------------------------*/

void leja_ordering_double( gmres_double_struct *p )
{

  int i, j, ii, d_poly;
  int max_j, exchange_cols;
  complex_double tmp, leja;

  complex_double** L;
  complex_double* col_prods;

  if (g.low_level_meas == 1) {
    print_vector_double("Eigenvalues", p->polyprec_double.h_ritz, p->polyprec_double.eigslvr.N);
  }

  d_poly = p->polyprec_double.d_poly;
  L = p->polyprec_double.L;
  col_prods = p->polyprec_double.col_prods;

  // Create a matrix made of n+1 rows, each row is x (all rows equal).
  for (i=0; i<d_poly+1; i++ )
    memcpy( L[i], p->polyprec_double.h_ritz, sizeof(complex_double)*(d_poly) );

  leja = 0; 

  for (i=0; i < d_poly-1; i++)
  {
    for (j=i; j<d_poly; j++ ) 
      L[i][j] = cabs( L[i][j] - leja );

    for (j = i; j < d_poly; j++)
    {
      col_prods[j] = 1.;
      for (ii = 0; ii <= i; ii++)
        col_prods[j] *= L[ii][j];
    }
        
    exchange_cols = 0;
    max_j = i;
    for (j=i+1; j<d_poly; j++ )
    {
      if ( creal(col_prods[j]) > creal(col_prods[max_j]) )
      {
        max_j = j; 
        exchange_cols = 1;
      }
    }
        
    if (exchange_cols)
    {
      for (ii=0; ii<d_poly+1; ii++ )
      {
        tmp = L[ii][i];
        L[ii][i] = L[ii][max_j];
        L[ii][max_j] = tmp;
      } 
    }

    leja = L[d_poly][i];

  }

  memcpy( p->polyprec_double.lejas, p->polyprec_double.L[d_poly], sizeof(complex_double)*(d_poly) );
  if (g.low_level_meas == 1) {
    print_vector_double("lejas", p->polyprec_double.lejas, d_poly);
  }
}



void update_lejas_double( gmres_double_struct *p, level_struct *l, struct Thread *threading )
{
  int start, end;
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  vector_double random_rhs, buff0;
  random_rhs = p->polyprec_double.random_rhs;
  double buff3;
  vector_double buff4;

  // vector_double_define(random_rhs, , start, end, l);

  int buff1, buff2;

  buff0 = p->b;
  buff2 = p->num_restart;
  buff1 = p->restart_length;
  buff3 = p->tol;
  buff4 = p->x;

  START_MASTER(threading)
  p->b = random_rhs;
  p->num_restart = 1;
  p->restart_length = p->polyprec_double.d_poly;
  p->preconditioner = NULL;
  p->tol = 1E-20;
  p->x = p->polyprec_double.xtmp;
  END_MASTER(threading)

  int fgmres_itersx;

  // TODO: add while loop to 
  START_MASTER(threading)
  l->dup_H = 1;
  END_MASTER(threading)

  START_LOCKED_MASTER(threading)
  vector_double_define_random( random_rhs, p->v_start, p->v_end, l );
  END_LOCKED_MASTER(threading)

  fgmres_itersx = fgmres_double(p, l, threading);

  START_MASTER(threading)
  l->dup_H = 0;
  END_MASTER(threading)

  START_MASTER(threading)
  p->b = buff0;
  p->num_restart = buff2;
  p->restart_length = buff1;
  p->tol = buff3;
  p->x = buff4;
  END_MASTER(threading)

  if ( fgmres_itersx == p->polyprec_double.d_poly ) {
    START_MASTER(threading)
    p->polyprec_double.preconditioner = p->polyprec_double.preconditioner_bare;
    l->p_double.polyprec_double.update_lejas = 0;
    END_MASTER(threading)

    printf0( "Updated polynomial preconditioner!\n" );

    SYNC_MASTER_TO_ALL(threading);
    SYNC_CORES(threading);

  } else { return; }

  START_MASTER(threading)
  harmonic_ritz_double(p);
  leja_ordering_double(p);
  END_MASTER(threading)

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)
}


void re_construct_lejas_double( level_struct *l, struct Thread *threading ) {

  update_lejas_double(&(l->p_double), l, threading);
}


void apply_polyprec_double( vector_double phi, vector_double Dphi, vector_double eta,
                               int res, level_struct *l, struct Thread *threading )
{

  START_MASTER(threading)
  //printf0("entering polyprec ...\n");
  END_MASTER(threading)

  int i, start, end;

  compute_core_start_end(l->p_double.v_start, l->p_double.v_end, &start, &end, l, threading);

  int d_poly = l->p_double.polyprec_double.d_poly;
  vector_double accum_prod = l->p_double.polyprec_double.accum_prod;
  vector_double product = l->p_double.polyprec_double.product;
  vector_double temp = l->p_double.polyprec_double.temp;
  vector_double lejas = l->p_double.polyprec_double.lejas;

  vector_double_copy( product, eta, start, end, l );
  vector_double_define(accum_prod, 0.0, start, end, l);

  vector_double_saxpy(accum_prod, accum_prod, product, 1./lejas[0], start, end, l);
  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)
  for (i = 1; i < d_poly; i++)
  {
#ifdef PERS_COMMS
    g.pers_comms_id2 = l->p_double.restart_length + g.pers_comms_nrZxs;
    g.use_pers_comms1 = 1;
#endif
    apply_operator_double(temp, product, &l->p_double, l, threading);
#ifdef PERS_COMMS
    g.pers_comms_id2 = -1;
    g.use_pers_comms1 = 0;
#endif

    vector_double_saxpy(product, product, temp, -1./lejas[i-1], start, end, l);
    vector_double_saxpy(accum_prod, accum_prod, product, 1./lejas[i], start, end, l);
  }

  vector_double_copy( phi, accum_prod, start, end, l );

  SYNC_MASTER_TO_ALL(threading)
  SYNC_CORES(threading)

  START_MASTER(threading)
  //printf0("exiting polyprec ...\n");
  END_MASTER(threading)

}

#endif
