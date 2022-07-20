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

#ifndef BLOCK_JACOBI_double_HEADER
  #define BLOCK_JACOBI_double_HEADER

  void block_jacobi_double_struct_init( gmres_double_struct *p );

  void block_jacobi_double_struct_alloc( int m, int n, long int vl, double tol, const int type, const int prec_kind,
                                            void (*precond)(), void (*eval_op)(), gmres_double_struct *p, level_struct *l );

  void block_jacobi_double_struct_free( gmres_double_struct *p, level_struct *l );

  void block_jacobi_apply_double( vector_double out, vector_double in, gmres_double_struct *p, level_struct *l, struct Thread *threading );

  //void block_jacobi_restore_from_buffer_double( vector_double out, gmres_double_struct *p, level_struct *l, struct Thread *threading );

  void block_jacobi_update_double( vector_double in, vector_double out, gmres_double_struct *p, level_struct *l, struct Thread *threading );

#endif
