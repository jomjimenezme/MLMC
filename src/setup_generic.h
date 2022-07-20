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

#ifndef SETUP_double_HEADER
  #define SETUP_double_HEADER

  struct Thread;
  
  void coarse_grid_correction_double_setup( level_struct *l, struct Thread *threading );
  void coarse_grid_correction_double_free( level_struct *l );
  void interpolation_double_define( vector_double* V, level_struct *l, struct Thread *threading );
  void iterative_double_setup( int setup_iter, level_struct *l, struct Thread *threading );
  void re_setup_double( level_struct *l, struct Thread *threading );
  void inv_iter_inv_fcycle_double( int setup_iter, level_struct *l, struct Thread *threading );
  
#endif

