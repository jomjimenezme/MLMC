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

#ifndef LOCAL_POLYPREC_double_HEADER
  #define LOCAL_POLYPREC_double_HEADER


  void local_set_ghost_double( vector_double phi, const int mu, const int dir,
                                  comm_double_struct *c, const int amount, level_struct *l );

  void coarse_local_n_hopping_term_double( vector_double out, vector_double in, operator_double_struct *op,
                                              const int amount, level_struct *l, struct Thread *threading );

  void coarse_local_hopping_term_double( vector_double out, vector_double in, operator_double_struct *op,
                                            const int amount, level_struct *l, struct Thread *threading );

  void coarse_local_apply_schur_complement_double( vector_double out, vector_double in,
                                                      operator_double_struct *op, level_struct *l,
                                                      struct Thread *threading );

  void local_apply_polyprec_double( vector_double phi, vector_double Dphi, vector_double eta,
                                       int res, level_struct *l, struct Thread *threading );

  void local_re_construct_lejas_double( level_struct *l, struct Thread *threading );

  void local_fgmres_double_struct_init( local_gmres_double_struct *p );

  void local_fgmres_double_struct_alloc( int m, int n, long int vl, double tol, const int type, const int prec_kind,
                                            void (*precond)(), void (*eval_op)(), local_gmres_double_struct *p, level_struct *l );

  void local_fgmres_double_struct_free( local_gmres_double_struct *p, level_struct *l );

  int local_fgmres_double( local_gmres_double_struct *p, level_struct *l, struct Thread *threading );

  int local_arnoldi_step_double( vector_double *V, vector_double *Z, vector_double w,
                                    complex_double **H, complex_double* buffer, int j, void (*prec)(),
                                    local_gmres_double_struct *p, level_struct *l, struct Thread *threading );

  void local_process_multi_inner_product_double( int count, complex_double *results, vector_double *phi, vector_double psi,
                                                    int start, int end, level_struct *l, struct Thread *threading );

  void local_qr_update_double( complex_double **H, complex_double *s,
                                  complex_double *c, complex_double *gamma, int j,
                                  level_struct *l, struct Thread *threading );

  void coarse_local_diag_oo_inv_double( vector_double y, vector_double x, operator_double_struct *op, 
                                           level_struct *l, struct Thread *threading );

  void coarse_local_diag_ee_double( vector_double y, vector_double x, operator_double_struct *op, level_struct *l, struct Thread *threading );

  void local_harmonic_ritz_double( local_gmres_double_struct *p );
  void local_leja_ordering_double( local_gmres_double_struct *p );
  void local_update_lejas_double( local_gmres_double_struct *p, level_struct *l, struct Thread *threading );
  void local_re_construct_lejas_double( level_struct *l, struct Thread *threading );
  void local_apply_polyprec_double( vector_double phi, vector_double Dphi, vector_double eta,
                                       int res, level_struct *l, struct Thread *threading );

#endif
