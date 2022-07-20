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

#ifndef LINSOLVE_double_HEADER
  #define LINSOLVE_double_HEADER

  struct Thread;
  
  void fgmres_double_struct_init( gmres_double_struct *p );
  void fgmres_double_struct_alloc( int m, int n, long int vl, double tol, const int type, const int prec_kind,
                                      void (*precond)(), void (*eval_op)(), gmres_double_struct *p, level_struct* l );
  void fgmres_double_struct_free( gmres_double_struct *p, level_struct *l );
  
  int fgmres_double( gmres_double_struct *p, level_struct *l, struct Thread *threading );
  void fgcr_double( gmres_double_struct *p, level_struct *l );
  void cgn_double( gmres_double_struct *p, level_struct *l, struct Thread *threading );
  void bicgstab_double( gmres_double_struct *ps, level_struct *l, struct Thread *threading );
  void local_minres_double( vector_double phi, vector_double eta, vector_double latest_iter,
                               int start, schwarz_double_struct *s, level_struct *l, struct Thread *threading );
  int arnoldi_step_double( vector_double *V, vector_double *Z, vector_double w,
                              complex_double **H, complex_double* buffer, int j, void (*prec)(),
                              gmres_double_struct *p, level_struct *l, struct Thread *threading );
  void qr_update_double( complex_double **H, complex_double *s,
                            complex_double *c, complex_double *gamma, int j,
                            level_struct *l, struct Thread *threading );
  void compute_solution_double( vector_double x, vector_double *V, complex_double *y, complex_double *gamma,
                                   complex_double **H, int j, int ol, gmres_double_struct *p, level_struct *l, struct Thread *threading );
  
#endif
