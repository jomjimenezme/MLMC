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

#ifndef COARSE_OPERATOR_double_HEADER
  #define COARSE_OPERATOR_double_HEADER

  #include "blas_vectorized.h"

  struct Thread;
  
  void coarse_operator_double_alloc( level_struct *l );
  void coarse_operator_double_free( level_struct *l );
  void coarse_operator_double_free_vectorized( operator_double_struct *op, level_struct *l );
  void coarse_operator_double_setup( vector_double *V, level_struct *l );
  void coarse_operator_double_setup_finalize( level_struct *l, struct Thread *threading );
  void coarse_operator_double_set_couplings( operator_double_struct *op, level_struct *l, struct Thread *threading );
  void coarse_operator_double_set_self_couplings( operator_double_struct *op, level_struct *l, struct Thread *threading );
  void coarse_operator_double_set_neighbor_couplings( operator_double_struct *op, level_struct *l, struct Thread *threading );
  
  void set_coarse_self_coupling_double( vector_double buffer1, vector_double buffer2,
                                           vector_double *V, const int n, level_struct *l );
  void set_coarse_neighbor_coupling_double( vector_double buffer1, vector_double buffer2,
                                               vector_double *V, const int mu, const int n, level_struct *l );

  void coarse_self_couplings_double( vector_double eta, vector_double phi,
                                        operator_double_struct *op, int start, int end, level_struct *l );
  void coarse_spinwise_self_couplings_double( vector_double eta1, vector_double eta2, vector_double phi, 
                                                 config_double clover, int length, level_struct *l );
  
  void coarse_gamma5_double( vector_double eta, vector_double phi, int start, int end, level_struct *l );
  void coarse_tau1_gamma5_double( vector_double eta, vector_double phi, int start, int end, level_struct *l );
  void apply_coarse_operator_double( vector_double eta, vector_double phi,
                                        operator_double_struct *op, level_struct *l, struct Thread *threading );
  void g5D_apply_coarse_operator_double( vector_double eta, vector_double phi, operator_double_struct *op,
                                            level_struct *l, struct Thread *threading );
  void apply_coarse_operator_dagger_double( vector_double eta, vector_double phi,
                                               operator_double_struct *op, level_struct *l, struct Thread *threading );
  void coarse_block_operator_double( vector_double eta, vector_double phi, int start,
                                        schwarz_double_struct *s, level_struct *l, struct Thread *threading );
  void coarse_aggregate_self_couplings_double( vector_double eta1, vector_double eta2, vector_double phi, 
                                                  schwarz_double_struct *s, level_struct *l );

  void coarse_aggregate_neighbor_couplings_double( vector_double eta1, vector_double eta2, vector_double phi, const int mu, schwarz_double_struct *s, level_struct *l );

  void set_block_diagonal_double( vector_double spin_0_1, vector_double spin_2_3, vector_double *V, const int n, config_double block, level_struct *l );

  void coarse_aggregate_block_diagonal_double( vector_double eta1, vector_double eta2, vector_double phi, config_double block, level_struct *l );
 
  void coarse_operator_double_test_routine( level_struct *l, struct Thread *threading );
  
  // eta += D*phi, D stored columnwise
  static inline void mv_double( const vector_double eta, const complex_double *D,
                                   const vector_double phi, const register int n ) {
    register int i, j, k=0;

    for ( i=0; i<n; i++ )
      for ( j=0; j<n; j++, k++ )
        eta[j] += D[k]*phi[i];
  }

  // eta -= D*phi, D stored columnwise
  static inline void nmv_double( const vector_double eta, const complex_double *D,
                                    const vector_double phi, const register int n ) {
    register int i, j, k=0;
    
    for ( i=0; i<n; i++ )
      for ( j=0; j<n; j++, k++ )
        eta[j] -= D[k]*phi[i];
  }

  // eta += D^Dagger*phi, D stored columnwise
  static inline void mvh_double( const vector_double eta, const complex_double *D,
                                    const vector_double phi, const register int n ) {
    register int i, j, k=0;    

    for ( i=0; i<n; i++ )
      for ( j=0; j<n; j++, k++ )
        eta[i] += conj_double(D[k])*phi[j];
  }

  // eta -= D^Dagger*phi, D stored columnwise
  static inline void nmvh_double( const vector_double eta, const complex_double *D,
                                     const vector_double phi, const register int n ) {
    register int i, j, k=0; 

    for ( i=0; i<n; i++ )
      for ( j=0; j<n; j++, k++ )
        eta[i] -= conj_double(D[k])*phi[j];
  }

  // eta = D*phi, D hermitian and stored columnwise packed
  static inline void mvp_double( const vector_double eta, const complex_double *D,
                                    const vector_double phi, const register int n ) {
    register int i, j, k;

    eta[0] = D[0]*phi[0];
    for ( i=1, k=1; i<n; i++ ) {
      eta[i] = conj_double(D[k])*phi[0];
      eta[0] += D[k]*phi[i]; k++;
      for ( j=1; j<i; j++, k++ ) {
        eta[j] += D[k]*phi[i];
        eta[i] += conj_double(D[k])*phi[j];
      }
      eta[i] += D[k]*phi[i]; k++;
    }
  }

  // eta += D*phi, D hermitian and stored columnwise packed
  static inline void pmvp_double( const vector_double eta, const complex_double *D,
                                    const vector_double phi, const register int n ) {
    register int i, j, k;

    eta[0] += D[0]*phi[0];
    for ( i=1, k=1; i<n; i++ ) {
      eta[i] += conj_double(D[k])*phi[0];
      eta[0] += D[k]*phi[i]; k++;
      for ( j=1; j<i; j++, k++ ) {
        eta[j] += D[k]*phi[i];
        eta[i] += conj_double(D[k])*phi[j];
      }
      eta[i] += D[k]*phi[i]; k++;
    }
  }

  // eta += D*phi, D hermitian and stored columnwise packed
  static inline void mmvp_double( const vector_double eta, const complex_double *D,
                                     const vector_double phi, const register int n ) {
    register int i, j, k;

    eta[0] -= D[0]*phi[0];
    for ( i=1, k=1; i<n; i++ ) {
      eta[i] -= conj_double(D[k])*phi[0];
      eta[0] -= D[k]*phi[i]; k++;
      for ( j=1; j<i; j++, k++ ) {
        eta[j] -= D[k]*phi[i];
        eta[i] -= conj_double(D[k])*phi[j];
      }
      eta[i] -= D[k]*phi[i]; k++;
    }
  }

  // eta += D*phi, D anti-hermitian and stored columnwise packed
  static inline void pamvp_double( const vector_double eta, const complex_double *D,
                                    const vector_double phi, const register int n ) {
    register int i, j, k;

    eta[0] += D[0]*phi[0];
    for ( i=1, k=1; i<n; i++ ) {
      eta[i] -= conj_double(D[k])*phi[0];
      eta[0] += D[k]*phi[i]; k++;
      for ( j=1; j<i; j++, k++ ) {
        eta[j] += D[k]*phi[i];
        eta[i] -= conj_double(D[k])*phi[j];
      }
      eta[i] += D[k]*phi[i]; k++;
    }
  }
  
  // eta -= D*phi, D anti-hermitian and stored columnwise packed
  static inline void mamvp_double( const vector_double eta, const complex_double *D,
                                    const vector_double phi, const register int n ) {
    register int i, j, k;

    eta[0] -= D[0]*phi[0];
    for ( i=1, k=1; i<n; i++ ) {
      eta[i] += conj_double(D[k])*phi[0];
      eta[0] -= D[k]*phi[i]; k++;
      for ( j=1; j<i; j++, k++ ) {
        eta[j] -= D[k]*phi[i];
        eta[i] += conj_double(D[k])*phi[j];
      }
      eta[i] -= D[k]*phi[i]; k++;
    }
  }

  static inline void coarse_self_couplings_clover_double( vector_double eta, vector_double phi,
                                                             config_double clover, int length, level_struct *l ) {
    
    int site_var = l->num_lattice_site_var,
      num_eig_vect = l->num_parent_eig_vect,
      clover_step_size1 = (num_eig_vect * (num_eig_vect+1))/2,
      clover_step_size2 = SQUARE(num_eig_vect);
    config_double clover_pt = clover;
    vector_double phi_pt=phi, eta_pt=eta, phi_end_pt=phi+length;
    // U(x) = [ A B      , A=A*, D=D*, C = -B*
    //          C D ]
    // storage order: upper triangle of A, upper triangle of D, B, columnwise
    // diagonal coupling
#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      while ( phi_pt < phi_end_pt ) {
        // A
        mvp_double( eta_pt, clover_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect;//1
        phi_pt += num_eig_vect;//1
        mvp_double( eta_pt, clover_pt, phi_pt, num_eig_vect );
        // D
        eta_pt += num_eig_vect;//2
        phi_pt += num_eig_vect;//2
        clover_pt += clover_step_size1; 
        mvp_double( eta_pt, clover_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect;//3
        phi_pt += num_eig_vect;//3
        mvp_double( eta_pt, clover_pt, phi_pt, num_eig_vect );
        // C = -B*
        eta_pt -= num_eig_vect;//2
        phi_pt -= 3*num_eig_vect;//0
        clover_pt += clover_step_size1;
        nmvh_double( eta_pt, clover_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect;//3
        phi_pt += num_eig_vect;//1
        nmvh_double( eta_pt, clover_pt, phi_pt, num_eig_vect );
        // B
        eta_pt -= 3*num_eig_vect;//0
        phi_pt += num_eig_vect;//2
        mv_double( eta_pt, clover_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect;//1
        phi_pt += num_eig_vect;//3
        mv_double( eta_pt, clover_pt, phi_pt, num_eig_vect );
        eta_pt += 3*num_eig_vect;//4
        phi_pt += num_eig_vect;//4
        clover_pt += clover_step_size2;
      }
    } else
#endif
      while ( phi_pt < phi_end_pt ) {
        // A
        mvp_double( eta_pt, clover_pt, phi_pt, num_eig_vect );
        clover_pt += clover_step_size1; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        // D
        mvp_double( eta_pt, clover_pt, phi_pt, num_eig_vect );
        clover_pt += clover_step_size1; phi_pt -= num_eig_vect;
        // C = -B*
        nmvh_double( eta_pt, clover_pt, phi_pt, num_eig_vect );
        phi_pt += num_eig_vect; eta_pt -= num_eig_vect;
        // B
        mv_double( eta_pt, clover_pt, phi_pt, num_eig_vect );
        clover_pt += clover_step_size2; phi_pt += num_eig_vect; eta_pt += site_var;
      }
  }

  static inline void coarse_add_block_diagonal_double( vector_double eta, vector_double phi,
                                                        config_double block, int length, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
      block_step_size = (num_eig_vect * (num_eig_vect+1))/2;
    config_double block_pt = block;
    vector_double phi_pt=phi, eta_pt=eta, phi_end_pt=phi+length;
    // U(x) = [ A 0      , A=A*, D=D* diag. excluded
    //          0 D ]
    // storage order: upper triangle of A, upper triangle of D, columnwise
    // diagonal coupling
#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      while ( phi_pt < phi_end_pt ) {
        // A
        pmvp_double( eta_pt, block_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        mmvp_double( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        // D
        pmvp_double( eta_pt, block_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        mmvp_double( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
      }
    } else
#endif
      while ( phi_pt < phi_end_pt ) {
        // A
        pmvp_double( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        // D
        pmvp_double( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
      }
  }

  static inline void coarse_add_anti_block_diagonal_double( vector_double eta, vector_double phi,
                                                               config_double block, int length, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
      block_step_size = (num_eig_vect * (num_eig_vect+1))/2;
    config_double block_pt = block;
    vector_double phi_pt=phi, eta_pt=eta, phi_end_pt=phi+length;
    // U(x) = [ A 0      , A=-A*, D=-D* diag. excluded
    //          0 D ]
    // storage order: upper triangle of A, upper triangle of D, columnwise
    // diagonal coupling
#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      while ( phi_pt < phi_end_pt ) {
        // A
        pamvp_double( eta_pt, block_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        mamvp_double( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        // D
        pamvp_double( eta_pt, block_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        mamvp_double( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
      }
    } else
#endif
      while ( phi_pt < phi_end_pt ) {
        // A
        pamvp_double( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        // D
        pamvp_double( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
      }
  }

  static inline void coarse_add_doublet_coupling_double( vector_double eta, vector_double phi,
                                                          config_double block, int length, level_struct *l ) {
    
#ifdef HAVE_TM1p1
    int num_eig_vect = l->num_parent_eig_vect,
      block_step_size = (num_eig_vect * (num_eig_vect+1))/2;
    config_double block_pt = block;
    vector_double phi_pt=phi, eta_pt=eta, phi_end_pt=phi+length;
    // U(x) = [ 0 A      , A=-A*, D=-D* diag. excluded
    //          D 0 ]
    // storage order: upper triangle of A, upper triangle of D, columnwise
    // diagonal coupling
    
    while ( phi_pt < phi_end_pt ) {
      // A
      pamvp_double( eta_pt, block_pt, phi_pt+num_eig_vect, num_eig_vect );
      pamvp_double( eta_pt+num_eig_vect, block_pt, phi_pt, num_eig_vect );
      block_pt += block_step_size; eta_pt += 2*num_eig_vect; phi_pt += 2*num_eig_vect;
      // D
      pamvp_double( eta_pt, block_pt, phi_pt+num_eig_vect, num_eig_vect );
      pamvp_double( eta_pt+num_eig_vect, block_pt, phi_pt, num_eig_vect );
      block_pt += block_step_size; eta_pt += 2*num_eig_vect; phi_pt += 2*num_eig_vect;
    }
#else
    warning0("coarse_add_doublet_coupling_double called without HAVE_TM1p1 defined.\n");
    return;
#endif
}
  
  static inline void coarse_hopp_double( vector_double eta, vector_double phi,
                                            config_double D, level_struct *l ) {
  
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      // A  
      nmv_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//1
      nmv_double( eta, D, phi, num_eig_vect );
      // C
      eta += num_eig_vect;//2
      phi -= num_eig_vect;//0
      D += num_eig_vect2;
      nmv_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//1
      nmv_double( eta, D, phi, num_eig_vect );
      // B
      eta -= 3*num_eig_vect;//0
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      nmv_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//3
      nmv_double( eta, D, phi, num_eig_vect );
      // D
      eta += num_eig_vect;//2
      phi -= num_eig_vect;//2
      D += num_eig_vect2;
      nmv_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//3
      nmv_double( eta, D, phi, num_eig_vect );
    } else {
#endif
      // A  
      nmv_double( eta, D, phi, num_eig_vect );
      // C
      eta += num_eig_vect;
      D += num_eig_vect2;
      nmv_double( eta, D, phi, num_eig_vect );
      // B
      phi += num_eig_vect;
      eta -= num_eig_vect;
      D += num_eig_vect2;
      nmv_double( eta, D, phi, num_eig_vect );
      // D
      eta += num_eig_vect;
      D += num_eig_vect2;
      nmv_double( eta, D, phi, num_eig_vect );
#ifdef HAVE_TM1p1
    }
#endif
  }


  static inline void coarse_daggered_hopp_double( vector_double eta, vector_double phi,
                                                     config_double D, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      // A* 
      nmvh_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//1
      nmvh_double( eta, D, phi, num_eig_vect );
      // -C*
      eta -= num_eig_vect;//0
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      mvh_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//3
      mvh_double( eta, D, phi, num_eig_vect );
      // -B*
      eta += num_eig_vect;//2
      phi -= 3*num_eig_vect;//0
      D += num_eig_vect2;
      mvh_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//1
      mvh_double( eta, D, phi, num_eig_vect );
      // D*
      eta -= num_eig_vect;//2
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      nmvh_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//3
      nmvh_double( eta, D, phi, num_eig_vect );
    } else {
#endif
      // A* 
      nmvh_double( eta, D, phi, num_eig_vect );
      // -C*
      phi += num_eig_vect;
      D += num_eig_vect2;
      mvh_double( eta, D, phi, num_eig_vect );
      // -B*
      eta += num_eig_vect;
      phi -= num_eig_vect;
      D += num_eig_vect2;
      mvh_double( eta, D, phi, num_eig_vect );
      // D*
      phi += num_eig_vect;
      D += num_eig_vect2;
      nmvh_double( eta, D, phi, num_eig_vect );
#ifdef HAVE_TM1p1
    }
#endif
  }
  
  static inline void coarse_n_hopp_double( vector_double eta, vector_double phi,
                                              config_double D, level_struct *l ) {
  
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      // A  
      mv_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//1
      mv_double( eta, D, phi, num_eig_vect );
      // C
      eta += num_eig_vect;//2
      phi -= num_eig_vect;//0
      D += num_eig_vect2;
      mv_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//1
      mv_double( eta, D, phi, num_eig_vect );
      // B
      eta -= 3*num_eig_vect;//0
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      mv_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//3
      mv_double( eta, D, phi, num_eig_vect );
      // D
      eta += num_eig_vect;//2
      phi -= num_eig_vect;//2
      D += num_eig_vect2;
      mv_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//3
      mv_double( eta, D, phi, num_eig_vect );
    } else {
#endif
      // A  
      mv_double( eta, D, phi, num_eig_vect );
      // C
      eta += num_eig_vect;
      D += num_eig_vect2;
      mv_double( eta, D, phi, num_eig_vect );
      // B
      phi += num_eig_vect;
      eta -= num_eig_vect;
      D += num_eig_vect2;
      mv_double( eta, D, phi, num_eig_vect );
      // D
      eta += num_eig_vect;
      D += num_eig_vect2;
      mv_double( eta, D, phi, num_eig_vect );
#ifdef HAVE_TM1p1
    }
#endif
  }

  static inline void coarse_n_daggered_hopp_double( vector_double eta, vector_double phi,
                                                       config_double D, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      // A* 
      mvh_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//1
      mvh_double( eta, D, phi, num_eig_vect );
      // -C*
      eta -= num_eig_vect;//0
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      nmvh_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//3
      nmvh_double( eta, D, phi, num_eig_vect );
      // -B*
      eta += num_eig_vect;//2
      phi -= 3*num_eig_vect;//0
      D += num_eig_vect2;
      nmvh_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//1
      nmvh_double( eta, D, phi, num_eig_vect );
      // D*
      eta -= num_eig_vect;//2
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      mvh_double( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//3
      mvh_double( eta, D, phi, num_eig_vect );
    } else {
#endif
      // A* 
      mvh_double( eta, D, phi, num_eig_vect );
      // -C*
      phi += num_eig_vect;
      D += num_eig_vect2;
      nmvh_double( eta, D, phi, num_eig_vect );
      // -B*
      eta += num_eig_vect;
      phi -= num_eig_vect;
      D += num_eig_vect2;
      nmvh_double( eta, D, phi, num_eig_vect );
      // D*
      phi += num_eig_vect;
      D += num_eig_vect2;
      mvh_double( eta, D, phi, num_eig_vect );
#ifdef HAVE_TM1p1
    }
#endif
  }

  static inline void coarse_spinwise_hopp_double( vector_double eta1, vector_double eta2, 
                                                     vector_double phi, config_double D, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

    // A  
    mv_double( eta1, D, phi, num_eig_vect );
    // C
    eta1 += num_eig_vect;
    D += num_eig_vect2;
    mv_double( eta1, D, phi, num_eig_vect );
    // B
    phi += num_eig_vect;
    D += num_eig_vect2;
    mv_double( eta2, D, phi, num_eig_vect );
    // D
    eta2 += num_eig_vect;
    D += num_eig_vect2;
    mv_double( eta2, D, phi, num_eig_vect );
  }


  static inline void coarse_spinwise_daggered_hopp_double( vector_double eta1, vector_double eta2,
                                                              vector_double phi, config_double D, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);  
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

    // A* 
    mvh_double( eta1, D, phi, num_eig_vect );
    // -C*
    phi += num_eig_vect;
    D += num_eig_vect2;
    nmvh_double( eta2, D, phi, num_eig_vect );
    // -B*
    eta1 += num_eig_vect;
    phi -= num_eig_vect;
    D += num_eig_vect2;
    nmvh_double( eta1, D, phi, num_eig_vect );
    // D*
    eta2 += num_eig_vect;
    phi += num_eig_vect;
    D += num_eig_vect2;
    mvh_double( eta2, D, phi, num_eig_vect );
  }

  static inline void coarse_spinwise_n_hopp_double( vector_double eta1, vector_double eta2,
                                                       vector_double phi, config_double D, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

    // A  
    nmv_double( eta1, D, phi, num_eig_vect );
    // C
    eta1 += num_eig_vect;
    D += num_eig_vect2;
    nmv_double( eta1, D, phi, num_eig_vect );
    // B
    phi += num_eig_vect;
    D += num_eig_vect2;
    nmv_double( eta2, D, phi, num_eig_vect );
    // D
    eta2 += num_eig_vect;
    D += num_eig_vect2;
    nmv_double( eta2, D, phi, num_eig_vect );
  }


  static inline void coarse_spinwise_n_daggered_hopp_double( vector_double eta1, vector_double eta2,
                                                                vector_double phi, config_double D, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);  
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

    // A* 
    nmvh_double( eta1, D, phi, num_eig_vect );
    // -C*
    phi += num_eig_vect;
    D += num_eig_vect2;
    mvh_double( eta2, D, phi, num_eig_vect );
    // -B*
    eta1 += num_eig_vect;
    phi -= num_eig_vect;
    D += num_eig_vect2;
    mvh_double( eta1, D, phi, num_eig_vect );
    // D*
    eta2 += num_eig_vect;
    phi += num_eig_vect;
    D += num_eig_vect2;
    nmvh_double( eta2, D, phi, num_eig_vect );
  }

#endif
