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

#ifndef GHOST_double_HEADER
  #define GHOST_double_HEADER
    
  void negative_sendrecv_double( vector_double phi, const int mu, comm_double_struct *c, level_struct *l );
  
  // as negative_sendrecv_double, but for count vectors stored in phi in vector-fused data layout
  // buffer must be big enough to hold the surface data for count vectors (in one direction)
  void negative_sendrecv_double_vectorized( complex_double *phi, const int mu, comm_double_struct *c, level_struct *l, int count, complex_double *buffer );
  void negative_wait_double( const int mu, comm_double_struct *c, level_struct *l );
  
  void ghost_alloc_double( int buffer_size, comm_double_struct *c, level_struct *l );
  void ghost_free_double( comm_double_struct *c, level_struct *l );
  void ghost_sendrecv_init_double( const int type, comm_double_struct *c, level_struct *l );
  void ghost_sendrecv_double( vector_double phi, const int mu, const int dir,
                                 comm_double_struct *c, const int amount, level_struct *l );
  void ghost_wait_double( vector_double phi, const int mu, const int dir,
                             comm_double_struct *c, const int amount, level_struct *l );
  
  void ghost_update_double( vector_double phi, const int mu, const int dir, comm_double_struct *c, level_struct *l );
  void ghost_update_wait_double( vector_double phi, const int mu, const int dir, comm_double_struct *c, level_struct *l );

#ifdef PERS_COMMS
  void pers_comms_init_double( vector_double phi, const int mu, const int dir, const int pers_comms_id1, const int pers_comms_id2,
                                  comm_double_struct *c, const int amount, level_struct *l );
  void pers_comms_open_double( level_struct *l );
  void pers_comms_free_double( vector_double phi, const int mu, const int dir, const int pers_comms_id1, const int pers_comms_id2,
                                  comm_double_struct *c, const int amount, level_struct *l );
  void pers_comms_close_double( level_struct *l );
#endif

#endif
