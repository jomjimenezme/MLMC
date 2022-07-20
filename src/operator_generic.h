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

#ifndef OPERATOR_double_HEADER
  #define OPERATOR_double_HEADER

  struct Thread;
  
  void operator_double_init( operator_double_struct *op );
  void operator_double_alloc( operator_double_struct *op, const int type, level_struct *l );
  void operator_double_define( operator_double_struct *op, level_struct *l );
  void operator_double_free( operator_double_struct *op, const int type, level_struct *l );

  void operator_double_set_couplings( operator_double_struct *op, level_struct *l );
  void operator_double_set_self_couplings( operator_double_struct *op, level_struct *l );
  void operator_double_set_neighbor_couplings( operator_double_struct *op, level_struct *l );
  
  void operator_double_test_routine( operator_double_struct *op, level_struct *l, struct Thread *threading );
  
#endif
