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

#ifndef SCHWARZ_double_HEADER
  #define SCHWARZ_double_HEADER

struct Thread;

  void block_double_boundary_op( vector_double eta, vector_double phi, int k,
                                    schwarz_double_struct *s, level_struct *l );
  void n_block_double_boundary_op( vector_double eta, vector_double phi, int k,
                                      schwarz_double_struct *s, level_struct *l );
  void coarse_block_double_boundary_op( vector_double eta, vector_double phi,
                                           int k, schwarz_double_struct *s, level_struct *l );
  void n_coarse_block_double_boundary_op( vector_double eta, vector_double phi,
                                             int k, schwarz_double_struct *s, level_struct *l );
  
  void smoother_double_def( level_struct *l );
  void smoother_double_free( level_struct *l );
  
  void schwarz_double_init( schwarz_double_struct *s, level_struct *l );
  void schwarz_double_alloc( schwarz_double_struct *s, level_struct *l );
  void schwarz_layout_double_define( schwarz_double_struct *s, level_struct *l );
  void schwarz_double_boundary_update( schwarz_double_struct *s, level_struct *l );
  
  void schwarz_double_setup( schwarz_double_struct *s, operator_double_struct *op_in, level_struct *l );
  
  void schwarz_double_def( schwarz_double_struct *s, operator_double_struct *op, level_struct *l );
  void schwarz_double_free( schwarz_double_struct *s, level_struct *l );
  
  void additive_schwarz_double( vector_double phi, vector_double D_phi, vector_double eta, const int cycles, int res, 
                                   schwarz_double_struct *s, level_struct *l, struct Thread *threading );
  void red_black_schwarz_double( vector_double phi, vector_double D_phi, vector_double eta, const int cycles, int res,
                                    schwarz_double_struct *s, level_struct *l, struct Thread *threading );
  void schwarz_double( vector_double phi, vector_double D_phi, vector_double eta, const int cycles, int res, 
                          schwarz_double_struct *s, level_struct *l, struct Thread *threading );
  void sixteen_color_schwarz_double( vector_double phi, vector_double D_phi, vector_double eta, const int cycles, int res, 
                                        schwarz_double_struct *s, level_struct *l, struct Thread *threading );
  
  void trans_double( vector_double out, vector_double in, int *tt, level_struct *l, struct Thread *threading );
  void trans_back_double( vector_double out, vector_double in, int *tt, level_struct *l, struct Thread *threading );
  
  void schwarz_double_mvm_testfun( schwarz_double_struct *s, level_struct *l, struct Thread *threading );
  
  static inline int connect_link_double( int t, int z, int y, int x, int mu, int dir, int *dt, int *it, 
                                            schwarz_double_struct *s, level_struct *l ) {
    
    int coord[4];
    coord[T]=t; coord[Z]=z; coord[Y]=y; coord[X]=x;
    coord[mu]+=dir;
    if ( l->global_splitting[mu] > 1 ) {
      return site_mod_index( coord[T], coord[Z], coord[Y], coord[X], dt, it );
    } else {
      coord[mu] = (coord[mu]+l->local_lattice[mu])%l->local_lattice[mu];
      return site_index( coord[T], coord[Z], coord[Y], coord[X], dt, it );
    }
  }

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_double
static inline void set_double_D_vectorized( double *out1, double *out2, complex_double *in ) {
  // out1: column major, out2: row major
  for ( int i=0; i<3; i++ ) { // column
    for ( int j=0; j<3; j++ ) { // row
      out1[8*i  +j] = creal_double(in[3*j+i]);
      out1[8*i+4+j] = cimag_double(in[3*j+i]);
      out2[8*i  +j] = creal_double(in[j+3*i]);
      out2[8*i+4+j] = cimag_double(in[j+3*i]);
    }
    out1[8*i+3] = 0.0;
    out1[8*i+7] = 0.0;
    out2[8*i+3] = 0.0;
    out2[8*i+7] = 0.0;
  }
}
#endif

#endif 
