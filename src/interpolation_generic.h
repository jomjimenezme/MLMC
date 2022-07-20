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

#ifndef INTERPOLATION_double_HEADER
  #define INTERPOLATION_double_HEADER

  struct Thread;
  
  void interpolation_double_alloc( level_struct *l );
  void interpolation_double_free( level_struct *l );
  void interpolation_double_dummy_alloc( level_struct *l );
  void interpolation_double_dummy_free( level_struct *l );
  
  void interpolate_double( vector_double phi, vector_double phi_c, level_struct *l, struct Thread *threading );
  void interpolate3_double( vector_double phi, vector_double phi_c, level_struct *l, struct Thread *threading );
  void restrict_double( vector_double phi_c, vector_double phi, level_struct *l, struct Thread *threading );
  
  void define_interpolation_double_operator( complex_double **interpolation, level_struct *l, struct Thread *threading );
#endif

