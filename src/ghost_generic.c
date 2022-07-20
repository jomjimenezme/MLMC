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

#include "main.h"

void negative_sendrecv_double( vector_double phi, const int mu, comm_double_struct *c, level_struct *l ) {
  // send dir = -1
  if( l->global_splitting[mu] > 1 ) {    
    
    int i, j, num_boundary_sites = c->num_boundary_sites[2*mu+1], boundary_start,
      *boundary_table = c->boundary_table[2*mu+1], n = l->num_lattice_site_var;
    
    vector_double buffer, tmp_pt, buffer_pt;
    
    boundary_start = l->num_inner_lattice_sites;
    for ( i=0; i<mu; i++ )
      boundary_start += c->num_boundary_sites[2*i];

    buffer = l->vbuf_double[8]+n*(boundary_start-l->num_inner_lattice_sites);
    buffer_pt = buffer;
    
    for ( i=0; i<num_boundary_sites; i++ ) {
      tmp_pt = phi + n*boundary_table[i];
      for ( j=0; j<n; j++, buffer_pt++, tmp_pt++ )
        *buffer_pt = *tmp_pt;
    }
    
    MPI_Irecv( phi+n*boundary_start, n*num_boundary_sites, MPI_COMPLEX_double,
               l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(c->rreqs[2*mu+1]) );
    MPI_Isend( buffer, n*num_boundary_sites, MPI_COMPLEX_double,
               l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(c->sreqs[2*mu+1]) );
  }
}


void negative_sendrecv_double_vectorized( complex_double *phi, const int mu, comm_double_struct *c,
                                             level_struct *l, int count, complex_double *buffer ) {
  // send dir = -1
  if( l->global_splitting[mu] > 1 ) {

    int i, j, num_boundary_sites = c->num_boundary_sites[2*mu+1], boundary_start,
        *boundary_table = c->boundary_table[2*mu+1], n = l->num_lattice_site_var;

    complex_double *tmp_pt;
    complex_double *buffer_pt;

    boundary_start = l->num_inner_lattice_sites;
    for ( i=0; i<mu; i++ )
      boundary_start += c->num_boundary_sites[2*i];

    buffer_pt = buffer;

    for ( i=0; i<num_boundary_sites; i++ ) {
      tmp_pt = phi + count*n*boundary_table[i];
      for ( j=0; j<count*n; j++, buffer_pt++, tmp_pt++ )
        *buffer_pt = *tmp_pt;
    }

    MPI_Irecv( phi+count*n*boundary_start, count*n*num_boundary_sites, MPI_COMPLEX_double,
               l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(c->rreqs[2*mu+1]) );
    MPI_Isend( buffer, count*n*num_boundary_sites, MPI_COMPLEX_double,
               l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(c->sreqs[2*mu+1]) );
  }
}


void negative_wait_double( const int mu, comm_double_struct *c, level_struct *l ) {
 
  if( l->global_splitting[mu] > 1 ) {
    MPI_Wait( &(c->sreqs[2*mu+1]), MPI_STATUS_IGNORE );
    MPI_Wait( &(c->rreqs[2*mu+1]), MPI_STATUS_IGNORE );
  }
}


void ghost_alloc_double( int buffer_size, comm_double_struct *c, level_struct *l ) {
  
  int mu, nu, factor=1;
  
  if ( l->depth > 0 ) {
    c->offset = l->num_lattice_site_var;
  } else {
    c->offset = l->num_lattice_site_var/2;
    if ( g.method < 5 )
      factor = 2;
  }

#ifdef HAVE_TM1p1
  factor *= 2;
#endif
  
  if ( buffer_size <= 0 ) {
    c->comm_start[0] = c->offset*l->num_inner_lattice_sites;
    c->comm_start[1] = c->offset*l->num_inner_lattice_sites;
    for ( mu=0; mu<4; mu++ ) {
      if ( mu > 0 ) {
        c->comm_start[2*mu] = c->comm_start[2*(mu-1)] + buffer_size;
        c->comm_start[2*mu+1] = c->comm_start[2*(mu-1)+1] + buffer_size;
      }
      buffer_size = c->offset;
      for ( nu=0; nu<4; nu++ ) {
        if ( nu != mu ) {
          buffer_size *= l->local_lattice[nu];
        }
      }
      c->length[2*mu] = buffer_size;
      c->length[2*mu+1] = buffer_size;
      c->max_length[mu] = factor*buffer_size;
      MALLOC( c->buffer[2*mu], complex_double, factor*buffer_size );
      MALLOC( c->buffer[2*mu+1], complex_double, factor*buffer_size );
      c->in_use[2*mu] = 0;
      c->in_use[2*mu+1] = 0;
    }
  } else {
    for ( mu=0; mu<4; mu++ ) {
      c->max_length[mu] = buffer_size;
#ifdef HAVE_TM1p1
      MALLOC( c->buffer[2*mu], complex_double, 2*buffer_size );
      MALLOC( c->buffer[2*mu+1], complex_double, 2*buffer_size );
#else
      MALLOC( c->buffer[2*mu], complex_double, buffer_size );
      MALLOC( c->buffer[2*mu+1], complex_double, buffer_size );
#endif
    }
  }
  
  if ( l->vbuf_double[8] == NULL ) {
#ifdef HAVE_TM1p1
    MALLOC( l->vbuf_double[8], complex_double, 2*l->vector_size );
#else
    MALLOC( l->vbuf_double[8], complex_double, l->vector_size );
#endif
  }
}


void ghost_free_double( comm_double_struct *c, level_struct *l ) {
  
  int mu;
  
  for ( mu=0; mu<4; mu++ ) {    
    FREE( c->buffer[2*mu], complex_double, c->max_length[mu] );
    FREE( c->buffer[2*mu+1], complex_double, c->max_length[mu] );
  }
  
  if ( l->vbuf_double[8] != NULL ) {
#ifdef HAVE_TM1p1
    FREE( l->vbuf_double[8], complex_double, 2*l->vector_size );
#else
    FREE( l->vbuf_double[8], complex_double, l->vector_size );
#endif
  }
}


void ghost_sendrecv_init_double( const int type, comm_double_struct *c, level_struct *l ) {
  
  int mu; 
  
  if ( type == _COARSE_GLOBAL ) {
    c->comm = 1;
    for ( mu=0; mu<4; mu++ ) {
      ASSERT( c->in_use[2*mu] == 0 );
      ASSERT( c->in_use[2*mu+1] == 0 );
    }
  }
}


void ghost_sendrecv_double( vector_double phi, const int mu, const int dir,
                               comm_double_struct *c, const int amount, level_struct *l ) {

  //gmres_double_struct *p = &(l->p_double);

  // does not allow sending in both directions at the same time
  if( l->global_splitting[mu] > 1 ) {
    
    int i, j, *table=NULL, mu_dir = 2*mu-MIN(dir,0), offset = c->offset,
        length[2] = {0,0}, comm_start = 0, table_start = 0;
    vector_double buffer, phi_pt;
    
    if ( amount == _FULL_SYSTEM ) {
      length[0] = (c->num_boundary_sites[2*mu])*offset;
      length[1] = (c->num_boundary_sites[2*mu+1])*offset;
      comm_start = c->comm_start[mu_dir];
      table_start = 0;
    } else if ( amount == _EVEN_SITES ) {
      length[0] = c->num_even_boundary_sites[2*mu]*offset;
      length[1] = c->num_even_boundary_sites[2*mu+1]*offset;
      comm_start = c->comm_start[mu_dir];
      table_start = 0;
    } else if ( amount == _ODD_SITES ) {
      length[0] = c->num_odd_boundary_sites[2*mu]*offset;
      length[1] = c->num_odd_boundary_sites[2*mu+1]*offset;
      comm_start = c->comm_start[mu_dir]+c->num_even_boundary_sites[mu_dir]*offset;
      table_start = c->num_even_boundary_sites[mu_dir];
    }
    
#ifdef HAVE_TM1p1
    if ( g.n_flavours == 2 ) {
      length[0] *= 2;
      length[1] *= 2;
      comm_start *= 2;
      offset *= 2;
    }
#endif

    ASSERT( c->in_use[mu_dir] == 0 );
    c->in_use[mu_dir] = 1;
    
    if ( MAX(length[0],length[1]) > c->max_length[mu] ) {
      printf("CAUTION: my_rank: %d, not enough comm buffer\n", g.my_rank ); fflush(0);
      ghost_free_double( c, l );
      ghost_alloc_double( MAX(length[0],length[1]), c, l );
    }
    
    buffer = (vector_double)c->buffer[mu_dir];

    // dir = senddir
    if ( dir == 1 ) {
      // data to be communicated is stored serially in the vector phi
      // recv target is a buffer
      // afterwards (in ghost_wait) the data has to be distributed onto the correct sites
      // touching the respective boundary in -mu direction
      
      phi_pt = phi + comm_start;
      if ( length[1] > 0 ) {
        PROF_double_START( _OP_COMM );
#ifdef PERS_COMMS
        if( l->level==0 && g.use_pers_comms1==1 && g.use_pers_comms2==1 ){
          if( g.pers_comms_id1==0 ){
            MPI_Start( &(g.pers_comms_recvrs_plus[2*mu][g.pers_comms_id2]) );
          }else{
            MPI_Start( &(g.pers_comms_recvrs_minus[2*mu][g.pers_comms_id2]) );
          }
        } else {
          MPI_Irecv( buffer, length[1], MPI_COMPLEX_double,
                     l->neighbor_rank[2*mu+1], 2*mu, g.comm_cart, &(c->rreqs[2*mu]) );
        }
#else
        MPI_Irecv( buffer, length[1], MPI_COMPLEX_double,
                   l->neighbor_rank[2*mu+1], 2*mu, g.comm_cart, &(c->rreqs[2*mu]) );
#endif
        PROF_double_STOP( _OP_COMM, 1 );
      }
      if ( length[0] > 0 ) {
        PROF_double_START( _OP_COMM );
#ifdef PERS_COMMS
        if( l->level==0 && g.use_pers_comms1==1 && g.use_pers_comms2==1 ){
          if( g.pers_comms_id1==0 ){
            MPI_Start( &(g.pers_comms_sendrs_plus[2*mu][g.pers_comms_id2]) );
          }else{
            MPI_Start( &(g.pers_comms_sendrs_minus[2*mu][g.pers_comms_id2]) );
          }
        } else {
          MPI_Isend( phi_pt, length[0], MPI_COMPLEX_double,
                     l->neighbor_rank[2*mu], 2*mu, g.comm_cart, &(c->sreqs[2*mu]) );
        }
#else
        MPI_Isend( phi_pt, length[0], MPI_COMPLEX_double,
                   l->neighbor_rank[2*mu], 2*mu, g.comm_cart, &(c->sreqs[2*mu]) );
#endif
        PROF_double_STOP( _OP_COMM, 0 );
      }
      
      
    } else if ( dir == -1 ) {
      // data to be communicated is stored on the sites touching the boundary in -mu direction
      // this data is gathered in a buffer in the correct ordering
      // which is required on the boundary of the vector phi
      int num_boundary_sites = length[1]/offset;
      
      table = c->boundary_table[2*mu+1]+table_start;
      for ( j=0; j<num_boundary_sites; j++ ) {
        phi_pt = phi + table[j]*offset;
        for ( i=0; i<offset; i++ ) {
          buffer[i] = phi_pt[i];
        }
        buffer += offset;
      }
      
      buffer = (vector_double)c->buffer[mu_dir];      
      phi_pt = phi + comm_start;
      
      if ( length[0] > 0 ) {
        PROF_double_START( _OP_COMM );
#ifdef PERS_COMMS
        if( l->level==0 && g.use_pers_comms1==1 && g.use_pers_comms2==1 ){
          if( g.pers_comms_id1==0 ){
            MPI_Start( &(g.pers_comms_recvrs_plus[2*mu+1][g.pers_comms_id2]) );
          }else{
            MPI_Start( &(g.pers_comms_recvrs_minus[2*mu+1][g.pers_comms_id2]) );
          }
        } else {
          MPI_Irecv( phi_pt, length[0], MPI_COMPLEX_double,
                     l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(c->rreqs[2*mu+1]) );
        }
#else
        MPI_Irecv( phi_pt, length[0], MPI_COMPLEX_double,
                   l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(c->rreqs[2*mu+1]) );
#endif
        PROF_double_STOP( _OP_COMM, 1 );
      }
      if ( length[1] > 0 ) {
        PROF_double_START( _OP_COMM );
#ifdef PERS_COMMS
        if( l->level==0 && g.use_pers_comms1==1 && g.use_pers_comms2==1 ){
          if( g.pers_comms_id1==0 ){
            MPI_Start( &(g.pers_comms_sendrs_plus[2*mu+1][g.pers_comms_id2]) );
          }else{
            MPI_Start( &(g.pers_comms_sendrs_minus[2*mu+1][g.pers_comms_id2]) );
          }
        } else {
          MPI_Isend( buffer, length[1], MPI_COMPLEX_double,
                     l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(c->sreqs[2*mu+1]) );
        }
#else
        MPI_Isend( buffer, length[1], MPI_COMPLEX_double,
                   l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(c->sreqs[2*mu+1]) );
#endif
        PROF_double_STOP( _OP_COMM, 0 );
      }
      
    } else ASSERT( dir == 1 || dir == -1 );

    //if( l->level==0 && g.use_pers_comms1==1 && g.use_pers_comms2==1 ){
    //  printf0("(SR) (%d,%d) mu=%d,dir=%d,sgn=%d,id=%d,oe=%d,comm_start=%d,g.in_pt=%p,g.out_pt=%p\n", g.use_pers_comms1,g.use_pers_comms2,mu,dir,g.pers_comms_id1,g.pers_comms_id2,amount,comm_start,p->pers_comms_ins[g.pers_comms_id2]+comm_start,p->pers_comms_outs[g.pers_comms_id2]+comm_start);
    //}

  }
}


void ghost_wait_double( vector_double phi, const int mu, const int dir,
                           comm_double_struct *c, const int amount, level_struct *l ) {
  
  if( l->global_splitting[mu] > 1 ) {    
    int mu_dir = 2*mu-MIN(dir,0);
    int i, j, *table, offset = c->offset, length[2]={0,0}, table_start = 0;
    vector_double buffer, phi_pt;

#ifdef HAVE_TM1p1
    if ( g.n_flavours == 2 )
      offset *= 2;
#endif
      
    if ( amount == _FULL_SYSTEM ) {
      length[0] = (c->num_boundary_sites[2*mu])*offset;
      length[1] = (c->num_boundary_sites[2*mu+1])*offset;
      table_start = 0;
    } else if ( amount == _EVEN_SITES ) {
      length[0] = c->num_even_boundary_sites[2*mu]*offset;
      length[1] = c->num_even_boundary_sites[2*mu+1]*offset;
      table_start = 0;
    } else if ( amount == _ODD_SITES ) {
      length[0] = c->num_odd_boundary_sites[2*mu]*offset;
      length[1] = c->num_odd_boundary_sites[2*mu+1]*offset;
      table_start = c->num_even_boundary_sites[mu_dir];
    }

    ASSERT( c->in_use[mu_dir] == 1 );
    
    if ( dir == 1 ) {
      
      int num_boundary_sites = length[0]/offset;
      
      buffer = (vector_double)c->buffer[mu_dir];      
      table = c->boundary_table[2*mu+1] + table_start;

      if ( length[0] > 0 ) {
        PROF_double_START( _OP_IDLE );
#ifdef PERS_COMMS
        if( l->level==0 && g.use_pers_comms1==1 && g.use_pers_comms2==1 ){
          if( g.pers_comms_id1==0 ){
            MPI_Wait( &(g.pers_comms_sendrs_plus[2*mu][g.pers_comms_id2]),MPI_STATUS_IGNORE );
          }else{
            MPI_Wait( &(g.pers_comms_sendrs_minus[2*mu][g.pers_comms_id2]),MPI_STATUS_IGNORE );
          }
        } else {
          MPI_Wait( &(c->sreqs[2*mu]), MPI_STATUS_IGNORE );
        }
#else
        MPI_Wait( &(c->sreqs[2*mu]), MPI_STATUS_IGNORE );
#endif
        PROF_double_STOP( _OP_IDLE, 0 );
      }
      if ( length[1] > 0 ) {
        PROF_double_START( _OP_IDLE );
#ifdef PERS_COMMS
        if( l->level==0 && g.use_pers_comms1==1 && g.use_pers_comms2==1 ){
          if( g.pers_comms_id1==0 ){
            MPI_Wait( &(g.pers_comms_recvrs_plus[2*mu][g.pers_comms_id2]),MPI_STATUS_IGNORE );
          }else{
            MPI_Wait( &(g.pers_comms_recvrs_minus[2*mu][g.pers_comms_id2]),MPI_STATUS_IGNORE );
          }
        } else {
          MPI_Wait( &(c->rreqs[2*mu]), MPI_STATUS_IGNORE );
        }
#else
        MPI_Wait( &(c->rreqs[2*mu]), MPI_STATUS_IGNORE );
#endif
        PROF_double_STOP( _OP_IDLE, 1 );
      }
      
      if ( l->depth == 0 ) {
        for ( j=0; j<num_boundary_sites; j++ ) {
          phi_pt = phi + table[j]*offset;
          
          for ( i=0; i<offset; i++ )
            phi_pt[i] = buffer[i];
          
          buffer += offset;
        }
      } else {
        for ( j=0; j<num_boundary_sites; j++ ) {
          phi_pt = phi + table[j]*offset;
          
          for ( i=0; i<offset; i++ )
            phi_pt[i] += buffer[i];
          
          buffer += offset;
        }
      }
    } else if ( dir == -1 ) {
      
      if ( length[1] > 0 ) {
        PROF_double_START( _OP_IDLE );
#ifdef PERS_COMMS
        if( l->level==0 && g.use_pers_comms1==1 && g.use_pers_comms2==1 ){
          if( g.pers_comms_id1==0 ){
            MPI_Wait( &(g.pers_comms_sendrs_plus[2*mu+1][g.pers_comms_id2]),MPI_STATUS_IGNORE );
          }else{
            MPI_Wait( &(g.pers_comms_sendrs_minus[2*mu+1][g.pers_comms_id2]),MPI_STATUS_IGNORE );
          }
        } else {
          MPI_Wait( &(c->sreqs[2*mu+1]), MPI_STATUS_IGNORE );
        }
#else
        MPI_Wait( &(c->sreqs[2*mu+1]), MPI_STATUS_IGNORE );
#endif
        PROF_double_STOP( _OP_IDLE, 0 );
      }
      if ( length[0] > 0 ) {
        PROF_double_START( _OP_IDLE );
#ifdef PERS_COMMS
        if( l->level==0 && g.use_pers_comms1==1 && g.use_pers_comms2==1 ){
          if( g.pers_comms_id1==0 ){
            MPI_Wait( &(g.pers_comms_recvrs_plus[2*mu+1][g.pers_comms_id2]),MPI_STATUS_IGNORE );
          }else{
            MPI_Wait( &(g.pers_comms_recvrs_minus[2*mu+1][g.pers_comms_id2]),MPI_STATUS_IGNORE );
          }
        } else {
          MPI_Wait( &(c->rreqs[2*mu+1]), MPI_STATUS_IGNORE );
        }
#else
        MPI_Wait( &(c->rreqs[2*mu+1]), MPI_STATUS_IGNORE );
#endif
        PROF_double_STOP( _OP_IDLE, 1 );
      }
      
    } else ASSERT( dir == 1 || dir == -1 );

    /*
    TODO : 
    		1. the problem at hand right now is that at the beginning we don't have the polynomial
           preconditioner, but we later do ... so, to fix this, open communications channels for both
           p->V and p->Za (when pipelining and preconditioner) or for p->V and p->Z (when no pipelining and yes preconditioner)

        2. NEXT ISSUE : the send/recv buffers in persistent comm.s should be set plus comm_start
    */

    /*
    if( l->level==0 && g.use_pers_comms1==1 && g.use_pers_comms2==1 ){
      //gmres_double_struct *p = &(l->p_double);
      printf0("\n");
      //printf0("%p\n", phi+comm_start);
      //printf0("%p\n", p->pers_comms_ins[g.pers_comms_id2]+comm_start);
      printf0("dir=%d\n", dir);
      printf0("mu=%d\n", mu);
      printf0("sgn=%d\n", g.pers_comms_id1);
      printf0("id=%d\n", g.pers_comms_id2);
      //printf0("comm_start=%d\n", comm_start);
      printf0("oe=%d\n", amount);
      printf0("\n");
      //printf0("\n\nexiting ... (%d,%p,%p) (%p) (%p) \n\n\n", g.pers_comms_id2, p->pers_comms_ins[g.pers_comms_id2],p->op->buffer[0], phi, p->Za[0]);
      MPI_Barrier( MPI_COMM_WORLD );
      exit(0);
    }
    */

    //if( l->level==0 && g.use_pers_comms1==1 && g.use_pers_comms2==1 ){
    //  printf0("(W) (%d,%d) mu=%d,dir=%d,sgn=%d,id=%d,oe=%d\n", g.use_pers_comms1,g.use_pers_comms2,mu,dir,g.pers_comms_id1,g.pers_comms_id2,amount);
    //}

    c->in_use[mu_dir] = 0;
  }
}


void ghost_update_double( vector_double phi, const int mu, const int dir, comm_double_struct *c, level_struct *l ) {

  if( l->global_splitting[mu] > 1 ) {
    int i, j, mu_dir = 2*mu-MIN(dir,0), nu, inv_mu_dir = 2*mu+1+MIN(dir,0), length, *table=NULL,
        comm_start, num_boundary_sites, site_var;
    vector_double buffer, recv_pt, phi_pt;
    
    site_var = l->num_lattice_site_var;
    length = c->num_boundary_sites[mu_dir]*l->num_lattice_site_var;
    num_boundary_sites = c->num_boundary_sites[mu_dir];
    buffer = (vector_double)c->buffer[mu_dir];
    
    if ( dir == -1 )
      comm_start = l->vector_size;
    else
      comm_start = l->inner_vector_size;
    for ( nu=0; nu<mu; nu++ ) {
      comm_start += c->num_boundary_sites[2*nu]*l->num_lattice_site_var;
    }
    
    ASSERT( c->in_use[mu_dir] == 0 );
    c->in_use[mu_dir] = 1;
    
    recv_pt = phi + comm_start;
    if ( length > 0 ) {
      PROF_double_START( _OP_COMM );
      MPI_Irecv( recv_pt, length, MPI_COMPLEX_double,
                 l->neighbor_rank[mu_dir], mu_dir, g.comm_cart, &(c->rreqs[mu_dir]) );
      PROF_double_STOP( _OP_COMM, 1 );
    }
    
    table = c->boundary_table[inv_mu_dir];
    for ( j=0; j<num_boundary_sites; j++ ) {
      phi_pt = phi + table[j]*site_var;
      
      for ( i=0; i<site_var; i++ ) {
        buffer[i] = phi_pt[i];
      }
      buffer += site_var;
    }
    buffer = (vector_double)c->buffer[mu_dir];
    
    if ( length > 0 ) {
      PROF_double_START( _OP_COMM );
      MPI_Isend( buffer, length, MPI_COMPLEX_double,
                 l->neighbor_rank[inv_mu_dir], mu_dir, g.comm_cart, &(c->sreqs[mu_dir]) );
      PROF_double_STOP( _OP_COMM, 0 );
    }
  }
}


void ghost_update_wait_double( vector_double phi, const int mu, const int dir, comm_double_struct *c, level_struct *l ) {
  
  if( l->global_splitting[mu] > 1 ) {
    int mu_dir = 2*mu-MIN(dir,0), length = c->num_boundary_sites[mu_dir]*l->num_lattice_site_var;
    
    ASSERT( c->in_use[mu_dir] == 1 );
      
    if ( length > 0 ) {
      PROF_double_START( _OP_IDLE );
      MPI_Wait( &(c->sreqs[mu_dir]), MPI_STATUS_IGNORE );
      MPI_Wait( &(c->rreqs[mu_dir]), MPI_STATUS_IGNORE );
      PROF_double_STOP( _OP_IDLE, 1 );
    }
    c->in_use[mu_dir] = 0;
  }
}





#ifdef PERS_COMMS

void pers_comms_open_double( level_struct *l ){

  g.use_pers_comms1=0;
  g.use_pers_comms2=0;

  gmres_double_struct *p = &(l->p_double);

  // FIXME : remove the factor of 4 in the next MALLOCs?

  int mx = p->restart_length;
  // pp_nr is the number of vector associated to POLYPREC
  int nrZxs,pp_nr;
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  if( p->preconditioner==NULL ){
    nrZxs = 0;
    pp_nr = 0;
    MALLOC( p->pers_comms_ins,vector_double,mx );
    MALLOC( p->pers_comms_outs,vector_double,mx );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_recvrs_plus[ix],MPI_Request,mx*4 );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_sendrs_plus[ix],MPI_Request,mx*4 );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_recvrs_minus[ix],MPI_Request,mx*4 );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_sendrs_minus[ix],MPI_Request,mx*4 );
  } else {
    // both p->V and p->Za
    nrZxs = g.pers_comms_nrZas;
    pp_nr = 1;
    MALLOC( p->pers_comms_ins,vector_double,mx+nrZxs+pp_nr );
    MALLOC( p->pers_comms_outs,vector_double,mx+nrZxs+pp_nr );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_recvrs_plus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_sendrs_plus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_recvrs_minus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_sendrs_minus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
  }
#else
  if( p->preconditioner==NULL ){
    nrZxs = 0;
    pp_nr = 0;
    MALLOC( p->pers_comms_ins,vector_double,mx );
    MALLOC( p->pers_comms_outs,vector_double,mx );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_recvrs_plus[ix],MPI_Request,mx*4 );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_sendrs_plus[ix],MPI_Request,mx*4 );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_recvrs_minus[ix],MPI_Request,mx*4 );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_sendrs_minus[ix],MPI_Request,mx*4 );
  } else {
    // both p->V and p->Z
    nrZxs = g.pers_comms_nrZs;
    pp_nr = 1;
    MALLOC( p->pers_comms_ins,vector_double,mx+nrZxs+pp_nr );
    MALLOC( p->pers_comms_outs,vector_double,mx+nrZxs+pp_nr );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_recvrs_plus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_sendrs_plus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_recvrs_minus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
    for(int ix=0;ix<8;ix++) MALLOC( g.pers_comms_sendrs_minus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
  }
#endif

  g.pers_comms_nrZxs = nrZxs;

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  if( p->preconditioner==NULL ){
    for( int ix=0;ix<mx;ix++ ) p->pers_comms_ins[ix] = p->V[ix];
    for( int ix=0;ix<mx;ix++ ) p->pers_comms_outs[ix] = p->block_jacobi_double.xtmp;
  }
  else{
    //int nrZas = g.pers_comms_nrZas;
    for( int ix=0;ix<mx;ix++ ) p->pers_comms_ins[ix] = p->V[ix];
    for( int ix=0;ix<mx;ix++ ) p->pers_comms_outs[ix] = p->block_jacobi_double.xtmp;
    for( int ix=0;ix<nrZxs;ix++ ) p->pers_comms_ins[mx+ix] = p->Za[ix];
    for( int ix=0;ix<nrZxs;ix++ ) p->pers_comms_outs[mx+ix] = p->block_jacobi_double.xtmp;
    for( int ix=0;ix<pp_nr;ix++ ) p->pers_comms_ins[mx+nrZxs+ix] = l->p_double.polyprec_double.product;
    for( int ix=0;ix<pp_nr;ix++ ) p->pers_comms_outs[mx+nrZxs+ix] = p->block_jacobi_double.xtmp;
  }
#else
  if( p->preconditioner==NULL ){
    for( int ix=0;ix<mx;ix++ ) p->pers_comms_ins[ix] = p->V[ix];
    for( int ix=0;ix<mx;ix++ ) p->pers_comms_outs[ix] = p->block_jacobi_double.xtmp;
  }
  else{
    //int nrZs = g.pers_comms_nrZs;
    for( int ix=0;ix<mx;ix++ ) p->pers_comms_ins[ix] = p->V[ix];
    for( int ix=0;ix<mx;ix++ ) p->pers_comms_outs[ix] = p->block_jacobi_double.xtmp;
    for( int ix=0;ix<nrZxs;ix++ ) p->pers_comms_ins[mx+ix] = p->Z[ix];
    for( int ix=0;ix<nrZxs;ix++ ) p->pers_comms_outs[mx+ix] = p->block_jacobi_double.xtmp;
    for( int ix=0;ix<pp_nr;ix++ ) p->pers_comms_ins[mx+nrZxs+ix] = l->p_double.polyprec_double.product;
    for( int ix=0;ix<pp_nr;ix++ ) p->pers_comms_outs[mx+nrZxs+ix] = p->block_jacobi_double.xtmp;
  }
#endif

  int minus_dir_param,plus_dir_param;

  // 0 : plus (odd!) hopping, send and recv
  for( int mu=0;mu<4;mu++ ){
    int dir;
    //+1
    dir=1;
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
    for(int ix=0;ix<mx+nrZxs+pp_nr;ix++) {
      pers_comms_init_double( p->op->buffer[0], mu, dir, 0, ix, &(p->op->c), plus_dir_param, l );
    }
    //-1
    dir=-1;
    for(int ix=0;ix<mx+nrZxs+pp_nr;ix++) {
      pers_comms_init_double( p->pers_comms_ins[ix], mu, dir, 0, ix, &(p->op->c), minus_dir_param, l );
    }
  }

  // 1 : minus (even!) hopping, send and recv
  for( int mu=0;mu<4;mu++ ){
    int dir;
    //+1
    dir=1;
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
    for(int ix=0;ix<mx+nrZxs+pp_nr;ix++) {
      pers_comms_init_double( p->pers_comms_outs[ix], mu, dir, 1, ix, &(p->op->c), plus_dir_param, l );
    }
    //-1
    dir=-1;
    for(int ix=0;ix<mx+nrZxs+pp_nr;ix++) {
      pers_comms_init_double( p->op->buffer[1], mu, dir, 1, ix, &(p->op->c), minus_dir_param, l );
    }
  }

  //printf0("\n");

/*
#ifdef GCRODR
  //int mx = p->restart_length;
  //int kx = p->gcrodr_double.k;

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  if( p->preconditioner==NULL ){
    // nrZxs should be zero here anyways
    for( int ix=0;ix<kx;ix++ ) p->pers_comms_ins[mx+nrZxs+ix] = p->gcrodr_double.C[ix];
    for( int ix=0;ix<kx;ix++ ) p->pers_comms_outs[mx+nrZxs+ix] = p->block_jacobi_double.xtmp;
  }
  else{
    for( int ix=0;ix<kx;ix++ ) p->pers_comms_ins[mx+nrZxs+ix] = p->gcrodr_double.C[ix];
    for( int ix=0;ix<kx;ix++ ) p->pers_comms_outs[mx+nrZxs+ix] = p->block_jacobi_double.xtmp;
    for( int ix=0;ix<kx;ix++ ) p->pers_comms_ins[mx+nrZxs+kx+ix] = p->gcrodr_double.PC[ix];
    for( int ix=0;ix<kx;ix++ ) p->pers_comms_outs[mx+nrZxs+kx+ix] = p->block_jacobi_double.xtmp;
  }

  // 0 : plus (odd!) hopping, send and recv
  for( int mu=0;mu<4;mu++ ){
    int dir;
    //+1
    dir=1;
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
    for(int ix=mx+nrZxs;ix<mx+nrZxs+kx*2;ix++) {
      pers_comms_init_double( p->op->buffer[0], mu, dir, 0, ix, &(p->op->c), plus_dir_param, l );
    }
    //-1
    dir=-1;
    for(int ix=mx+nrZxs;ix<mx+nrZxs+kx*2;ix++) {
      pers_comms_init_double( p->pers_comms_ins[ix], mu, dir, 0, ix, &(p->op->c), minus_dir_param, l );
    }
  }

  // 1 : minus (even!) hopping, send and recv
  for( int mu=0;mu<4;mu++ ){
    int dir;
    //+1
    dir=1;
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
    for(int ix=mx+nrZxs;ix<mx+nrZxs+kx*2;ix++) {
      pers_comms_init_double( p->pers_comms_outs[ix], mu, dir, 1, ix, &(p->op->c), plus_dir_param, l );
    }
    //-1
    dir=-1;
    for(int ix=mx+nrZxs;ix<mx+nrZxs+kx*2;ix++) {
      pers_comms_init_double( p->op->buffer[1], mu, dir, 1, ix, &(p->op->c), minus_dir_param, l );
    }
  }
#else
    //if( p->preconditioner==NULL ){
    //  for( int ix=0;ix<mx;ix++ ) g.pers_comms_ins[ix] = p->V[ix];
    //}
    //else{
    //  for( int ix=0;ix<mx;ix++ ) g.pers_comms_ins[ix] = p->Z[ix];
    //}
#endif

#else

    //vector_double product = l->p_double.polyprec_double.product;
    //vector_double temp = l->p_double.polyprec_double.temp;

#endif
*/

}


void pers_comms_init_double( vector_double phi, const int mu, const int dir, const int pers_comms_id1, const int pers_comms_id2,
                                comm_double_struct *c, const int amount, level_struct *l ) {
  // does not allow sending in both directions at the same time
  if( l->global_splitting[mu] > 1 ) {

    int mu_dir = 2*mu-MIN(dir,0), offset = c->offset,
        length[2] = {0,0}, comm_start = 0;
    vector_double buffer, phi_pt;
    
    if ( amount == _FULL_SYSTEM ) {
      length[0] = (c->num_boundary_sites[2*mu])*offset;
      length[1] = (c->num_boundary_sites[2*mu+1])*offset;
      comm_start = c->comm_start[mu_dir];
    } else if ( amount == _EVEN_SITES ) {
      length[0] = c->num_even_boundary_sites[2*mu]*offset;
      length[1] = c->num_even_boundary_sites[2*mu+1]*offset;
      comm_start = c->comm_start[mu_dir];
    } else if ( amount == _ODD_SITES ) {
      length[0] = c->num_odd_boundary_sites[2*mu]*offset;
      length[1] = c->num_odd_boundary_sites[2*mu+1]*offset;
      comm_start = c->comm_start[mu_dir]+c->num_even_boundary_sites[mu_dir]*offset;
    }
    
#ifdef HAVE_TM1p1
    if ( g.n_flavours == 2 ) {
      length[0] *= 2;
      length[1] *= 2;
      comm_start *= 2;
      offset *= 2;
    }
#endif

    buffer = (vector_double)c->buffer[mu_dir];
    
    // dir = senddir
    if ( dir == 1 ) {
      phi_pt = phi + comm_start;

      if ( length[1] > 0 ) {

        // 0 is plus hopping, 1 is minus hopping
        if( pers_comms_id1==0 ){
          MPI_Recv_init( buffer, length[1], MPI_COMPLEX_double, l->neighbor_rank[2*mu+1], 2*mu, g.comm_cart, &(g.pers_comms_recvrs_plus[2*mu][pers_comms_id2]) );
        } else {
          MPI_Recv_init( buffer, length[1], MPI_COMPLEX_double, l->neighbor_rank[2*mu+1], 2*mu, g.comm_cart, &(g.pers_comms_recvrs_minus[2*mu][pers_comms_id2]) );
        }

      }
      if ( length[0] > 0 ) {

        // 0 is plus hopping, 1 is minus hopping
        if( pers_comms_id1==0 ){
          MPI_Send_init( phi_pt, length[0], MPI_COMPLEX_double, l->neighbor_rank[2*mu], 2*mu, g.comm_cart, &(g.pers_comms_sendrs_plus[2*mu][pers_comms_id2]) );
        } else {
          MPI_Send_init( phi_pt, length[0], MPI_COMPLEX_double, l->neighbor_rank[2*mu], 2*mu, g.comm_cart, &(g.pers_comms_sendrs_minus[2*mu][pers_comms_id2]) );
        }

      }

    } else if ( dir == -1 ) {
      buffer = (vector_double)c->buffer[mu_dir];      
      phi_pt = phi + comm_start;
      
      if ( length[0] > 0 ) {

        // 0 is plus hopping, 1 is minus hopping
        if( pers_comms_id1==0 ){
          MPI_Recv_init( phi_pt, length[0], MPI_COMPLEX_double, l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(g.pers_comms_recvrs_plus[2*mu+1][pers_comms_id2]) );
        } else {
          MPI_Recv_init( phi_pt, length[0], MPI_COMPLEX_double, l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(g.pers_comms_recvrs_minus[2*mu+1][pers_comms_id2]) );
        }

      }
      if ( length[1] > 0 ) {

        // 0 is plus hopping, 1 is minus hopping
        if( pers_comms_id1==0 ){
          MPI_Send_init( buffer, length[1], MPI_COMPLEX_double, l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(g.pers_comms_sendrs_plus[2*mu+1][pers_comms_id2]) );
        } else {
          MPI_Send_init( buffer, length[1], MPI_COMPLEX_double, l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(g.pers_comms_sendrs_minus[2*mu+1][pers_comms_id2]) );
        }

      }

    } else ASSERT( dir == 1 || dir == -1 );

    //printf0("mu=%d \tdir=%d \tid=%d \tphi_pt=%p sgn=%d oe=%d (even=%d) comm_start=%d, reqobj_p1=%p, reqobj_m1=%p\n", mu, dir, pers_comms_id2, phi+comm_start, pers_comms_id1, amount, _EVEN_SITES, comm_start, &(g.pers_comms_sendrs_minus[2*mu][pers_comms_id2]), &(g.pers_comms_sendrs_minus[2*mu+1][pers_comms_id2]));

  }
}


void pers_comms_close_double( level_struct *l ){

  g.use_pers_comms1=0;
  g.use_pers_comms2=0;

  gmres_double_struct *p = &(l->p_double);

  int mx = p->restart_length;
  // pp_nr is the number of vector associated to POLYPREC
  int nrZxs,pp_nr;
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  if( p->preconditioner==NULL ){
    nrZxs = 0;
    pp_nr = 0;
  } else {
    // both p->V and p->Za
    nrZxs = g.pers_comms_nrZas;
    pp_nr = 1;
  }
#else
  if( p->preconditioner==NULL ){
    nrZxs = 0;
    pp_nr = 0;
  } else {
    // both p->V and p->Z
    nrZxs = g.pers_comms_nrZs;
    pp_nr = 1;
  }
#endif

  g.pers_comms_nrZxs = nrZxs;

  int minus_dir_param,plus_dir_param;

  // 0 : plus (odd!) hopping, send and recv
  for( int mu=0;mu<4;mu++ ){
    int dir;
    //+1
    dir=1;
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
    for(int ix=0;ix<mx+nrZxs+pp_nr;ix++) {
      pers_comms_free_double( p->op->buffer[0], mu, dir, 0, ix, &(p->op->c), plus_dir_param, l );
    }
    //-1
    dir=-1;
    for(int ix=0;ix<mx+nrZxs+pp_nr;ix++) {
      pers_comms_free_double( p->pers_comms_ins[ix], mu, dir, 0, ix, &(p->op->c), minus_dir_param, l );
    }
  }

  // 1 : minus (even!) hopping, send and recv
  for( int mu=0;mu<4;mu++ ){
    int dir;
    //+1
    dir=1;
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
    for(int ix=0;ix<mx+nrZxs+pp_nr;ix++) {
      pers_comms_free_double( p->pers_comms_outs[ix], mu, dir, 1, ix, &(p->op->c), plus_dir_param, l );
    }
    //-1
    dir=-1;
    for(int ix=0;ix<mx+nrZxs+pp_nr;ix++) {
      pers_comms_free_double( p->op->buffer[1], mu, dir, 1, ix, &(p->op->c), minus_dir_param, l );
    }
  }

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  if( p->preconditioner==NULL ){
    nrZxs = 0;
    pp_nr = 0;
    FREE( p->pers_comms_ins,vector_double,mx );
    FREE( p->pers_comms_outs,vector_double,mx );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_recvrs_plus[ix],MPI_Request,mx*4 );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_sendrs_plus[ix],MPI_Request,mx*4 );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_recvrs_minus[ix],MPI_Request,mx*4 );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_sendrs_minus[ix],MPI_Request,mx*4 );
  } else {
    // both p->V and p->Za
    nrZxs = g.pers_comms_nrZas;
    pp_nr = 1;
    FREE( p->pers_comms_ins,vector_double,mx+nrZxs+pp_nr );
    FREE( p->pers_comms_outs,vector_double,mx+nrZxs+pp_nr );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_recvrs_plus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_sendrs_plus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_recvrs_minus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_sendrs_minus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
  }
#else
  if( p->preconditioner==NULL ){
    nrZxs = 0;
    pp_nr = 0;
    FREE( p->pers_comms_ins,vector_double,mx );
    FREE( p->pers_comms_outs,vector_double,mx );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_recvrs_plus[ix],MPI_Request,mx*4 );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_sendrs_plus[ix],MPI_Request,mx*4 );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_recvrs_minus[ix],MPI_Request,mx*4 );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_sendrs_minus[ix],MPI_Request,mx*4 );
  } else {
    // both p->V and p->Z
    nrZxs = g.pers_comms_nrZs;
    pp_nr = 1;
    FREE( p->pers_comms_ins,vector_double,mx+nrZxs+pp_nr );
    FREE( p->pers_comms_outs,vector_double,mx+nrZxs+pp_nr );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_recvrs_plus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_sendrs_plus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_recvrs_minus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
    for(int ix=0;ix<8;ix++) FREE( g.pers_comms_sendrs_minus[ix],MPI_Request,(mx+nrZxs+pp_nr)*4 );
  }
#endif
}


void pers_comms_free_double( vector_double phi, const int mu, const int dir, const int pers_comms_id1, const int pers_comms_id2,
                                comm_double_struct *c, const int amount, level_struct *l ) {
  // does not allow sending in both directions at the same time
  if( l->global_splitting[mu] > 1 ) {

    int mu_dir = 2*mu-MIN(dir,0), offset = c->offset,
        length[2] = {0,0}, comm_start = 0;
    
    if ( amount == _FULL_SYSTEM ) {
      length[0] = (c->num_boundary_sites[2*mu])*offset;
      length[1] = (c->num_boundary_sites[2*mu+1])*offset;
      comm_start = c->comm_start[mu_dir];
    } else if ( amount == _EVEN_SITES ) {
      length[0] = c->num_even_boundary_sites[2*mu]*offset;
      length[1] = c->num_even_boundary_sites[2*mu+1]*offset;
      comm_start = c->comm_start[mu_dir];
    } else if ( amount == _ODD_SITES ) {
      length[0] = c->num_odd_boundary_sites[2*mu]*offset;
      length[1] = c->num_odd_boundary_sites[2*mu+1]*offset;
      comm_start = c->comm_start[mu_dir]+c->num_even_boundary_sites[mu_dir]*offset;
    }
    
#ifdef HAVE_TM1p1
    if ( g.n_flavours == 2 ) {
      length[0] *= 2;
      length[1] *= 2;
      comm_start *= 2;
      offset *= 2;
    }
#endif
    
    // dir = senddir
    if ( dir == 1 ) {
      if ( length[1] > 0 ) {
        // 0 is plus hopping, 1 is minus hopping
        if( pers_comms_id1==0 ){
          MPI_Request_free( &(g.pers_comms_recvrs_plus[2*mu][pers_comms_id2]) );
        } else {
          MPI_Request_free( &(g.pers_comms_recvrs_minus[2*mu][pers_comms_id2]) );
        }
      }
      if ( length[0] > 0 ) {
        // 0 is plus hopping, 1 is minus hopping
        if( pers_comms_id1==0 ){
          MPI_Request_free( &(g.pers_comms_sendrs_plus[2*mu][pers_comms_id2]) );
        } else {
          MPI_Request_free( &(g.pers_comms_sendrs_minus[2*mu][pers_comms_id2]) );
        }
      }
    } else if ( dir == -1 ) {
      if ( length[0] > 0 ) {
        // 0 is plus hopping, 1 is minus hopping
        if( pers_comms_id1==0 ){
          MPI_Request_free( &(g.pers_comms_recvrs_plus[2*mu+1][pers_comms_id2]) );
        } else {
          MPI_Request_free( &(g.pers_comms_recvrs_minus[2*mu+1][pers_comms_id2]) );
        }
      }
      if ( length[1] > 0 ) {
        // 0 is plus hopping, 1 is minus hopping
        if( pers_comms_id1==0 ){
          MPI_Request_free( &(g.pers_comms_sendrs_plus[2*mu+1][pers_comms_id2]) );
        } else {
          MPI_Request_free( &(g.pers_comms_sendrs_minus[2*mu+1][pers_comms_id2]) );
        }
      }

    } else ASSERT( dir == 1 || dir == -1 );

  }
}

#endif
