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
 
#ifndef ODDEVEN_double_HEADER
  #define ODDEVEN_double_HEADER

struct Thread;

  void hopping_term_double( vector_double eta, vector_double phi, operator_double_struct *op,
                               const int amount, level_struct *l, struct Thread *threading );
  
  void oddeven_setup_double( operator_double_struct *in, level_struct *l );
  void oddeven_free_double( level_struct *l );
  
  void oddeven_to_serial_double( vector_double out, vector_double in, level_struct *l, struct Thread *threading );
  void serial_to_oddeven_double( vector_double out, vector_double in, level_struct *l, struct Thread *threading );
  
  void oddeven_to_block_double( vector_double out, vector_double in, level_struct *l, struct Thread *threading );
  void block_to_oddeven_double( vector_double out, vector_double in, level_struct *l, struct Thread *threading );
  
  void block_hopping_term_double( vector_double eta, vector_double phi,
                                     int start, int amount, schwarz_double_struct *s, level_struct *l, struct Thread *threading );
  void block_n_hopping_term_double( vector_double eta, vector_double phi, 
                                       int start, int amount, schwarz_double_struct *s, level_struct *l, struct Thread *threading );
  void block_diag_oo_inv_double( vector_double eta, vector_double phi,
                                    int start, schwarz_double_struct *s, level_struct *l, struct Thread *threading );
  void block_diag_oo_double( vector_double eta, vector_double phi,
                                int start, schwarz_double_struct *s, level_struct *l, struct Thread *threading );
  void block_diag_ee_double( vector_double eta, vector_double phi,
                                int start, schwarz_double_struct *s, level_struct *l, struct Thread *threading );
  
  void apply_schur_complement_double( vector_double out, vector_double in, operator_double_struct *op, level_struct *l, struct Thread *threading );
  void solve_oddeven_double( gmres_double_struct *p, operator_double_struct *op, level_struct *l, struct Thread *threading );
  void g5D_apply_schur_complement_double( vector_double out, vector_double in, operator_double_struct *op, level_struct *l, struct Thread *threading );
  void g5D_solve_oddeven_double( gmres_double_struct *p, operator_double_struct *op, level_struct *l, struct Thread *threading );
  
  void schwarz_double_oddeven_setup( schwarz_double_struct *s, level_struct *l );
  
  void apply_block_schur_complement_double( vector_double out, vector_double in, int start,
                                               schwarz_double_struct *s, level_struct *l, struct Thread *threading );
  void block_solve_oddeven_double( vector_double phi, vector_double r, vector_double latest_iter,
                                      int start, schwarz_double_struct *s, level_struct *l, struct Thread *threading );
  
  void oddeven_double_test( level_struct *l );
  void block_oddeven_double_test( level_struct *l, struct Thread *threading );
  
#endif
