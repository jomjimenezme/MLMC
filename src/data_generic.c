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

// vector storage for double precision
void vector_double_define( vector_double phi, complex_double value, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_double_START( _SET );
  if ( phi != NULL ) {
    int i;
    for ( i=start; i<end; i++ )
      phi[i] = value;
  } else {
    error0("Error in \"vector_double_define\": pointer is null\n");
  }
  if(thread == 0 && start != end)
    PROF_double_STOP( _SET, 1 );
}


void vector_double_define_random( vector_double phi, int start, int end, level_struct *l ) {

  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_double_START( _SET );
  if ( phi != NULL ) {
    int i;
    for ( i=start; i<end; i++ )
      phi[i] = (double)(((double)rand()/(double)RAND_MAX))-0.5 + ( (double)((double)rand()/(double)RAND_MAX)-0.5)*_Complex_I;
  } else {
    error0("Error in \"vector_double_define_random\": pointer is null\n");
  }
  if(thread == 0 && start != end)
    PROF_double_STOP( _SET, 1 );
}

void vector_double_define_random_rademacher( vector_double phi, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_double_START( _SET );
  if ( phi != NULL ) {
    int i;
    for ( i=start; i<end; i++ )
      if(   (double)((double)rand()<(double)RAND_MAX/2.0)   ) phi[i]=  (double) (-1);
      else phi[i]= (double)(1);
  } else {
    error0("Error in \"vector_double_define_random\": pointer is null\n");
  }
  if(thread == 0 && start != end)
  PROF_double_STOP( _SET, 1 );
  
}
