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

#ifndef COARSENING_double_HEADER
  #define COARSENING_double_HEADER
  
  void interpolation_double_struct_init( interpolation_double_struct *is );
  void coarsening_index_table_double_alloc( interpolation_double_struct *is, level_struct *l );
  void coarsening_index_table_double_free( interpolation_double_struct *is, level_struct *l );
  void coarsening_index_table_double_define( interpolation_double_struct *is, schwarz_double_struct *s,
                                                level_struct *l );
  
#endif
