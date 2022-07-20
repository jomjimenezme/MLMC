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

#ifdef MUMPS_ADDS

#include "mumps_double.h"


void mumps_setup_double(level_struct *l, struct Thread *threading){

  START_MASTER(threading)

  double t0,t1;
  t0 = MPI_Wtime();

  gmres_double_struct* px = &(l->p_double);
  operator_double_struct* op = px->op; //&(l->op_double);//
  config_double clover_pt = l->p_double.op->clover;

  int num_eig_vect = l->num_parent_eig_vect;
  int site_var = l->num_lattice_site_var,
      clover_step_size1 = (num_eig_vect * (num_eig_vect+1))/2,
      clover_step_size2 = SQUARE(num_eig_vect);

  int nr_nodes = l->num_inner_lattice_sites;
  int i, j, k; // k = index in matrix
  int c, r;	// col no., row no.
  int skip = 8 * SQUARE(site_var);	//skip number of elements in Blockrow in large matrix (for self coupl. only skip = 0, else 8 * SQUARE(site_var))
  int j_start = nr_nodes * g.my_rank * site_var, i_start = nr_nodes * g.my_rank * site_var;

  for (j = 0, k = 0; j < nr_nodes; j++){
    for (i = 0; i < SQUARE(site_var); i++, k++){
      *(px->mumps_Is +k) = i_start + j * site_var + (int)(i/site_var);	// col indices
      *(px->mumps_Js +k) = j_start + j * site_var + (i % site_var); 	// row indices
    }
    k += skip;
  }

//      nr_nodes * g.my_rank 			= nodes per process * p_id
//	nr_nodes * (g.my_rank +1)		= first node of next process
  for (j = 0; j < nr_nodes; j++){
    // A store column-wise
    for (k = 0, r = 0; r < num_eig_vect; r++, k++){
      for (c = 0; c < r; c++, k++){
        px->mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] = *(clover_pt + k);
        px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] = conj_double(*(clover_pt + k));
      }
      px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] = *(clover_pt + k);
    }

    //remove this line as well, as soon k++ is removed
    clover_pt += clover_step_size1;

	// D store column-wise
    for (k = 0, r = num_eig_vect; r < 2*num_eig_vect; r++, k++){
      for (c = num_eig_vect; c < r; c++, k++){
        px->mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] = *(clover_pt + k);
        px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] = conj_double(*(clover_pt + k));
      }
      px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] = *(clover_pt + k);
    }

    //remove this line as well, as soon k++ is removed
    clover_pt += clover_step_size1;

    // C store column-wise
    for (r = num_eig_vect, k = 0; r < 2*num_eig_vect; r++){
      for (c = 0; c < num_eig_vect; c++, k++){
        px->mumps_vals[j*9*SQUARE(site_var) + (r * site_var) + c] = -1.0*(conj_double(*(clover_pt + k)));
        //px->mumps_vals[j*9*SQUARE(site_var) + (r +  num_eig_vect) * site_var + c] =  -1.0*(conj_double(*(clover_pt + k)));
      }
    }

    //no clover_pt correction / change this once k++ is removed
    // B store column-wise / transposed from former storage
    for (r = 0, k = 0; r < num_eig_vect; r++){
      for (c = 0; c < num_eig_vect; c++, k++){
        px->mumps_vals[j*9*SQUARE(site_var) + (c * site_var) + r + num_eig_vect] = *(clover_pt + k);
      }
    }
    clover_pt += clover_step_size2;
    // skipping num for hopping terms in mumps_vals
  }  // end for loop over the blocks

#ifdef HAVE_TM
  // #########################################################################
  // twisted mass-term:
  int block_step_size = (num_eig_vect * (num_eig_vect+1))/2;
  config_double tm_block_pt = op->tm_term;	//TODO add  + start * block_step_size *2

  for (j = 0; j < nr_nodes; j++){
    // A store column-wise
    for (k = 0, r = 0; r < num_eig_vect; r++, k++){
      for (c = 0; c < r; c++, k++){
        px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] += *(tm_block_pt + k);
        px->mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] -= conj_double(*(tm_block_pt + k));
      }
      px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] += *(tm_block_pt + k);
    }

    //remove this line as well, as soon k++ is removed
    tm_block_pt += block_step_size;

    // D store column-wise
    for (k = 0, r = num_eig_vect; r < 2*num_eig_vect; r++, k++){
      for (c = num_eig_vect; c < r; c++, k++){
        px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] += *(tm_block_pt + k);
        px->mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] -= conj_double(*(tm_block_pt + k));
      }
      px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] += *(tm_block_pt + k);
    }

    //remove this line as well, as soon k++ is removed
    tm_block_pt += block_step_size;
  }  // end for loop over the blocks
#endif	//tm term

  // #########################################################################
  // #2 hopping-term:
  // #########################################################################
 
  /*	vals = [[self_coupling of site 1][T-_coupling site 1][T+_coupling site 1][Z-_coupling site 1][Z+_coupling site 1] .... [X-_coupling site 1][X+_coupling site 1]
  [self_coup site 2][T-_coup site 2]....[X+_coup site N]]
  each of the inner [] contain num_link_var elements -> to store 1 block row in matrix (entire coupling of one site) we need 9 * num_link_var elements

  vals =   [[self, T-, T+, Z-, Z+, Y-, Y+, X-, X+][self, T-, T+, Z-, Z+, Y-, Y+, X-, X+]....]
  */

  //START_NO_HYPERTHREADS(threading)
  {

  int index,
      num_4link_var=4*4*l->num_parent_eig_vect*l->num_parent_eig_vect,
      num_link_var=4*l->num_parent_eig_vect*l->num_parent_eig_vect,
      start=0;

  //int core_start;
  //int core_end;
  //compute_core_start_end_custom(start, start+nr_nodes, &core_start, &core_end, l, threading, 1);
  int core_start = start;
  int core_end = start+nr_nodes;

  int comm_nr[4] = {0, 0, 0, 0}; 	//number elements to communicate in all directions
  int dir;
  int node;

  // FIXME (?) : USE c->num_boundary_sites[2*mu+1] instead
  // calculate the number of elements to communicate:
  for (dir = T; dir <= X; dir++){
    for (node = core_start; node < core_end; node++){
      if (op->neighbor_table[5*node+1+dir] >= nr_nodes){
        comm_nr[dir]++;
      }
    }
  }

  //allocate memory for buffers:
  int *buff_i_send[4], *buff_i_recv[4];				// will contain node nr. in reciever processors domain
  complex_double *buff_d_send[4], *buff_d_recv[4];		// will contain mu- coupling

  int buffer_i_pt, buffer_d_pt;
  MPI_Request r;
  MPI_Status s;

  int i_start = g.my_rank * l->num_inner_lattice_sites * site_var, 
      j_start = g.my_rank * l->num_inner_lattice_sites * site_var; 	//contains own global row and col. start indices
  int neighbors_j_start; //contains column start index of neighboring process
  
  //printf("r: %d, \tT: %d, \tZ: %d, \tY: %d, \tX: %d\n", g.my_rank, comm_nr[T], comm_nr[Z], comm_nr[Y], comm_nr[X]);

  int *boundary_table;// = op->c.boundary_table[...];
  int bt_index;

  int num_site_var=site_var;
  
  for (dir = T; dir <= X; dir++){
    boundary_table = op->c.boundary_table[2*dir];
    buffer_i_pt = 0;

    buff_i_send[dir] = NULL;
    buff_i_recv[dir] = NULL;
    MALLOC(buff_i_send[dir], int, 2 * comm_nr[dir]);
    MALLOC(buff_i_recv[dir], int, 2 * comm_nr[dir]);

    buff_d_send[dir] = NULL;
    buff_d_recv[dir] = NULL;
    MALLOC(buff_d_send[dir], complex_double, num_link_var * comm_nr[dir]);
    MALLOC(buff_d_recv[dir], complex_double, num_link_var * comm_nr[dir]);

    memset(buff_i_send[dir], 0, 2 * comm_nr[dir] * sizeof(int));
    memset(buff_i_recv[dir], 0, 2 * comm_nr[dir] * sizeof(int));
    memset(buff_d_send[dir], 0, num_link_var * comm_nr[dir] * sizeof(complex_double));
    memset(buff_d_recv[dir], 0, num_link_var * comm_nr[dir] * sizeof(complex_double));

    bt_index = 0;
    neighbors_j_start = l->neighbor_rank[2*dir] * l->num_inner_lattice_sites * site_var;
    for (node = core_start; node < core_end; node ++){
      index = 5 * node;

      // make mu+ couplings as usual (Values + Row indices aka. Is)
      // A
      /*
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
	printf0("k :%d, neighbor_table: %d, index: %d, n4link: %d, dir: %d, nlink: %d\n", k, op->neighbor_table[index], index, num_4link_var, dir, num_link_var);
        //printf0("k :%4d, D: %+-f%+-fi\n", k, CSPLIT(-1.0 * *(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + k)));
      }
      MPI_Finalize();
      exit(0);
      */
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        //find correct block row		skip self coupl.		find pos of mu+ coupl. ("2*mu +1" due to structure of vals[self, T-, T+, Z-, Z+...]
        //l->p_double.mumps_vals[9 * num_link_var)*op->neighbor_table[index] + num_link_var + 		(2*dir + 1)*num_link_var + k]
        *(l->p_double.mumps_vals + 	(9 * num_link_var)*op->neighbor_table[index] + num_link_var + 		(2*dir + 1)*num_link_var + k) = 
			-1.0 * *(op->D + 	num_4link_var*op->neighbor_table[index] + 	dir*num_link_var + k);
	//					find correct block row				start of mu- coupling
        *(l->p_double.mumps_Is +	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		(2*dir +1)*num_link_var + k) = 
			i_start +	num_site_var * op->neighbor_table[index] + 		k%((int)(num_site_var*0.5));
        //	proc start		block row start						fast changing index
      }

      // C
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_double.mumps_vals +	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 	(2*dir+1)*num_link_var + 	1 * (int)SQUARE(num_site_var/2) + k) = 
			-1.0 * *(op->D + 	num_4link_var*op->neighbor_table[index] + 	dir*num_link_var + 	1 * (int)SQUARE(num_site_var/2) + k);
        *(l->p_double.mumps_Is + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 	(2*dir +1)*num_link_var + 	1 * (int)SQUARE(num_site_var/2) + k) = 
			i_start +	num_site_var * op->neighbor_table[index] + 		k%(int)(num_site_var*0.5) + 	(int)(num_site_var*0.5);
      }

      // B
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_double.mumps_vals + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		(2*dir+1)*num_link_var + 	2 * (int)SQUARE(num_site_var/2) + k) = 
			-1.0 * *(op->D +	num_4link_var*op->neighbor_table[index] + 	dir*num_link_var + 	2 * (int)SQUARE(num_site_var/2) + k); 	//columwise in D
        *(l->p_double.mumps_Is + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 	(2*dir+1)*num_link_var + 	2 * (int)SQUARE(num_site_var/2) + k) = 
			i_start +	num_site_var * op->neighbor_table[index] + 		k%(int)(num_site_var*0.5) + 	0;
      }

      // D
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_double.mumps_vals + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 	(2*dir+1)*num_link_var + 	3 * (int)SQUARE(num_site_var/2) + k) = 
			-1.0 * *(op->D + 	num_4link_var*op->neighbor_table[index] + 	dir*num_link_var + 	3 * (int)SQUARE(num_site_var/2) + k); 	//columwise in D
        *(l->p_double.mumps_Is + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 	(2*dir +1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) = 
			i_start +	num_site_var * op->neighbor_table[index] + 		k%(int)(num_site_var*0.5) + 	(int)(num_site_var*0.5);
      }

      // DO THE COL INDICES aka. Js
      //FIXME (?) : maybe remove comm_nr > 0 ?
      if (comm_nr[dir] > 0 && op->neighbor_table[index+1+dir] >= l->num_inner_lattice_sites){
        for (k = 0; k < SQUARE(num_site_var/2); k ++){
          //A
          *(l->p_double.mumps_Js +	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		(2*dir +1)*num_link_var + k) = 
			neighbors_j_start + 	num_site_var * boundary_table[bt_index] +	k/((int)(num_site_var*0.5));
          //C
	  *(l->p_double.mumps_Js + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		(2*dir +1)*num_link_var + 	1 * (int)SQUARE(num_site_var/2) + k) = 
			neighbors_j_start + 	num_site_var * boundary_table[bt_index] +	k/((int)(num_site_var*0.5));
          //B
          *(l->p_double.mumps_Js + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		(2*dir+1)*num_link_var + 	2 * (int)SQUARE(num_site_var/2) + k) = 
			neighbors_j_start + 	num_site_var * boundary_table[bt_index] +	k/((int)(num_site_var*0.5)) +	(int)(num_site_var*0.5);
          //D
          *(l->p_double.mumps_Js + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		(2*dir + 1)*num_link_var + 	3 * (int)SQUARE(num_site_var/2) + k) = 
			neighbors_j_start + 	num_site_var * boundary_table[bt_index] +	k/((int)(num_site_var*0.5)) + 	(int)(num_site_var*0.5);

        }  //  boundary_table[op->neighbor_table[index +1 + dir] % l->num_inner_lattice_sites]

	bt_index++;
      } else {	//my neighbor is on same processor
        for (k = 0; k < SQUARE(num_site_var/2); k ++){
          //A
          *(l->p_double.mumps_Js +	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		(2*dir +1)*num_link_var + k) = 
          		j_start +	num_site_var * op->neighbor_table[index +1 + dir] +	k/((int)(num_site_var*0.5));
          //C 
          *(l->p_double.mumps_Js + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		(2*dir +1)*num_link_var + 	1 * (int)SQUARE(num_site_var/2) + k) = 
			j_start + 	num_site_var * op->neighbor_table[index +1 + dir] +	k/(int)(num_site_var*0.5) + 	0;
          //B
          *(l->p_double.mumps_Js + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		(2*dir+1)*num_link_var + 	2 * (int)SQUARE(num_site_var/2) + k) = 
			j_start +	num_site_var * op->neighbor_table[index +1 + dir] +	k/(int)(num_site_var*0.5) + 	(int)(num_site_var*0.5);
          //D
          *(l->p_double.mumps_Js + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		(2*dir+1)*num_link_var + 	3 * (int)SQUARE(num_site_var/2) + k) = 
			j_start +	num_site_var * op->neighbor_table[index +1 + dir] +	k/(int)(num_site_var*0.5) + 	(int)(num_site_var*0.5);
          //		p_start		find start of T- coupling 		slow changing index
        }
      }

      // check whether dir neighbor is in halo

      if (op->neighbor_table[index+1+dir] >= l->num_inner_lattice_sites) {
        //printf("(proc=%d) i am node %d, coupling to node %d\n", g.my_rank, op->neighbor_table[index], op->neighbor_table[index+1+dir]);

        // write mu- coupling to buffer
	// also write global target and source site no. to buffer
	// send both buffers

	buffer_d_pt = buffer_i_pt * num_link_var;

        //write mu- couplings:
	//  A* #################################
        for (k = 0; k < SQUARE(site_var/2); k ++){
          *(buff_d_send[dir] + buffer_d_pt + k) = 
				-1.0 * conj_double(*(op->D	+ num_4link_var*op->neighbor_table[index] + 	dir*num_link_var + 	k));
        }
	// -C* #################################
        for (k = 0; k < SQUARE(site_var/2); k ++){
          *(buff_d_send[dir] + buffer_d_pt + 1 * (int)SQUARE(site_var/2) +	k) = 
				 1.0 * conj_double(*(op->D 	+ num_4link_var*op->neighbor_table[index] + 	dir*num_link_var + 	1*(int)SQUARE(site_var/2) + k));
        }
	// -B* #################################
        for (k = 0; k < SQUARE(site_var/2); k ++){
          *(buff_d_send[dir] + buffer_d_pt + 2 * (int)SQUARE(site_var/2) + 	k) = 
				 1.0 * conj_double(*(op->D 	+ num_4link_var*op->neighbor_table[index] + 	dir*num_link_var + 	2 * (int)SQUARE(site_var/2) + k));
        }
	// D* #################################
        for (k = 0; k < SQUARE(site_var/2); k ++){
          *(buff_d_send[dir] + buffer_d_pt + 3 * (int)SQUARE(site_var/2) +	k) = 
				-1.0 * conj_double(*(op->D	+ num_4link_var*op->neighbor_table[index] + 	dir*num_link_var + 	3 * (int)SQUARE(site_var/2) + k));
        }

	//write i to buffer:
        *(buff_i_send[dir] + 2 * buffer_i_pt) = boundary_table[op->neighbor_table[index + 1 + dir] % comm_nr[dir]];
        *(buff_i_send[dir] + 2 * buffer_i_pt + 1) = op->neighbor_table[index];
 	buffer_i_pt++;
      }
    }	//loop over nodes

    //send both buffers:
    /*
    MPI_Isend(cont void *buf, int count, MPI_Datatype dt, int dest, int tag, MPI_Comm comm, MPI_Request *req);
    */
    if (comm_nr[dir] > 0){
      MPI_Isend(buff_d_send[dir], comm_nr[dir] * num_link_var, MPI_COMPLEX_double, l->neighbor_rank[2*dir], dir, g.comm_cart, &r);
      MPI_Isend(buff_i_send[dir], 2 * comm_nr[dir], MPI_INT, l->neighbor_rank[2*dir], dir, g.comm_cart, &r);
      //printf("process %d send message in direction %d to process %d\n", g.my_rank, dir, l->neighbor_rank[2*dir]);
      // HOW TO FIND NEIGHBOR?
      //int l.neighbor_rank[8] contains ranks of neighbors
      //in the order [T+ T- Z+ Z- ...]
      MPI_Barrier(MPI_COMM_WORLD);
      // ensures buffer-sending order: T, Z, Y, X
      // FIXME remove barrier. already replaced by tag = dir
    }
  }  // loop over directions

  // mu- couplings
  for (dir = T; dir <= X; dir++){
    if (comm_nr[dir] > 0){
      //there is stuff to communicate in direction dir

      //MPI_Recv(void *buf, int count, MPI_Datatype dt, int source, int tag, MPI_Comm comm, MPI_Status *status);
      //		len					[T+, T-, Z+, Z-...]
      MPI_Recv(buff_d_recv[dir], num_link_var * comm_nr[dir], MPI_COMPLEX_double, l->neighbor_rank[2*dir+1], dir, g.comm_cart, &s);
      MPI_Recv((buff_i_recv[dir]), 2 * comm_nr[dir], MPI_INT, l->neighbor_rank[2*dir+1], dir, g.comm_cart, &s);
      // tag = dir
      //printf("process %d recieved message from process %d\n", g.my_rank, l->neighbor_rank[2*dir+1]);

      neighbors_j_start = l->neighbor_rank[2*dir+1] * l->num_inner_lattice_sites * site_var;	
      // copy buffer content to mumps_vals
      for (buffer_i_pt = 0; buffer_i_pt < comm_nr[dir]; buffer_i_pt++){
        buffer_d_pt = num_link_var * buffer_i_pt;

	// A* #################################
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(l->p_double.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1]	+ num_link_var + 	2*dir*num_link_var + k) = 
			*(buff_d_recv[dir] + buffer_d_pt + k);
          // straight copy (orders are set before send)
	  *(l->p_double.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] 	+ num_link_var +	2*dir*num_link_var + k) = 
			i_start +	num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + 		k/(int)(num_site_var*0.5);
	  *(l->p_double.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] 	+ num_link_var +	2*dir*num_link_var + k) = 
			neighbors_j_start +	num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] +	k%(int)(num_site_var*0.5);
        }
	// -C* #################################
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(l->p_double.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1]+ num_link_var + 	2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 1 * SQUARE((int)(num_site_var*0.5)) + k);
          // straight copy (orders are set before send)
	  *(l->p_double.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] 	+ num_link_var +	2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start +	num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + 		k/(int)(num_site_var*0.5);
	  *(l->p_double.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] 	+ num_link_var +	2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			neighbors_j_start +	num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] +	k%(int)(num_site_var*0.5) + num_site_var*0.5;
        }
	// -B* #################################
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(l->p_double.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1]+ num_link_var + 	2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 2 * SQUARE((int)(num_site_var*0.5)) + k);
          // straight copy (orders are set before send)
	  *(l->p_double.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] 	+ num_link_var +	2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start +	num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + 		k/(int)(num_site_var*0.5) + num_site_var*0.5;
	  *(l->p_double.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] 	+ num_link_var +	2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			neighbors_j_start +	num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] +	k%(int)(num_site_var*0.5);
        } 
	// D* #################################
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(l->p_double.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1]+ num_link_var + 	2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 3 * SQUARE((int)(num_site_var*0.5)) + k);
          // straight copy (orders are set before send)
	  *(l->p_double.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] 	+ num_link_var +	2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start +	num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + 		k/(int)(num_site_var*0.5) + num_site_var*0.5;
	  *(l->p_double.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] 	+ num_link_var +	2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			neighbors_j_start +	num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] +	k%(int)(num_site_var*0.5) + num_site_var*0.5;
        }

        // message comes from process l->neighbor_rank[2*dir+1]
        // col nr. in sparse matrix:  	l->neighbor_rank[2*dir+1] * 
        //				l->num_inner_lattice_sites * 
        //				site_var
      }
    }	//end if (comm_nr[dir] > 0)

    // regular mu- coupling for all nodes except communicated ones
    buffer_i_pt = 0;
    for (i = 0; i < core_end; i++){	//loop over lattice sites
      if (comm_nr[dir] > 0){
        while (i == *(buff_i_recv[dir] + 2 * buffer_i_pt + 1)){
          // skip this node because it was already communicated
          i++;
          buffer_i_pt++;
          if (i >= core_end) break;
        }
	if (i >= core_end) break;
      }
      index = 5 * i;
      // regular mu- coupling
      // A* ##################################
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        // skip to proc start	find correct block row				skip self coupl.	find pos of T- coupl.
        *(l->p_double.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 	 	2*dir*num_link_var + k) = 
			-1.0 * conj_double(*(op->D + 	num_4link_var*op->neighbor_table[index] + 	dir*num_link_var + k));
	//				    proc start		find correct block row				start of T- coupling

        *(l->p_double.mumps_Is + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		2*dir*num_link_var + k) = 
			i_start +	num_site_var * op->neighbor_table[index + 1 + dir] + 		k/(int)(num_site_var*0.5);
	//	proc start		block row start						fast changing index

        *(l->p_double.mumps_Js + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		2*dir*num_link_var + k) = 
			j_start +	num_site_var * op->neighbor_table[index] +	k%(int)(num_site_var*0.5);
	//		no p_start	find start of T- coupling 				slow changing index	
      }
      // -C* ##################################
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_double.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 	 	2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			1.0 * conj_double(*(op->D + 	num_4link_var*op->neighbor_table[index] + 	dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k));
        *(l->p_double.mumps_Is + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start +	num_site_var * op->neighbor_table[index + 1 + dir] + 		k/(int)(num_site_var*0.5);
        *(l->p_double.mumps_Js + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			j_start +	num_site_var * op->neighbor_table[index] +	k%(int)(num_site_var*0.5) + num_site_var*0.5;
      }
      // -B* #################################
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_double.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 	 	2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			1.0 * conj_double(*(op->D + 	num_4link_var*op->neighbor_table[index] + 	dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k));
        *(l->p_double.mumps_Is + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start +	num_site_var * op->neighbor_table[index + 1 + dir] + 		k/(int)(num_site_var*0.5) + num_site_var*0.5;
        *(l->p_double.mumps_Js + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			j_start +	num_site_var * op->neighbor_table[index] +	k%(int)(num_site_var*0.5);
      }
      // D* #################################
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_double.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 	 	2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			-1.0 * conj_double(*(op->D + 	num_4link_var*op->neighbor_table[index] + 	dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k));
        *(l->p_double.mumps_Is + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start +	num_site_var * op->neighbor_table[index + 1 + dir] + 		k/(int)(num_site_var*0.5) + num_site_var*0.5;
        *(l->p_double.mumps_Js + 	(9 * num_link_var)*op->neighbor_table[index] + 	num_link_var + 		2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			j_start +	num_site_var * op->neighbor_table[index] +	k%(int)(num_site_var*0.5) + num_site_var*0.5;
      }

    }	//loop over nodes
  }	//loop over directions

  }
  //END_NO_HYPERTHREADS(threading)

  // increase global indices by 1 to match fortran indexing.
  // spmv doesn't work anymore!!!

  int nnz_loc = SQUARE(site_var) * nr_nodes *9;
  for (i = 0; i < nnz_loc; i++){	//increase indices by one to match fortran indexing
    *(l->p_double.mumps_Js + i ) = *(l->p_double.mumps_Js + i ) +1;
    *(l->p_double.mumps_Is + i ) = *(l->p_double.mumps_Is + i ) +1;
  }

  t1 = MPI_Wtime();

  printf0("MUMPS pre-setup time (seconds) : %f\n",t1-t0);

  END_MASTER(threading)
  SYNC_CORES(threading)
}



void mumps_solve_double( vector_double phi, vector_double Dphi, vector_double eta,
                            int res, level_struct *lx, struct Thread *threading )
{

  //printf0("Applying MUMPS solver ... \n");

  START_MASTER(threading)

  gmres_double_struct* px = &(lx->p_double);

  int i;
  //int start, end, i;
  //compute_core_start_end_custom(0, lx->num_inner_lattice_sites, &start, &end, lx, threading, 1);
  //int start = 0;
  //int end = lx->num_inner_lattice_sites;

  // ######### SET UP RHS #############
  int rhs_len = lx->p_double.v_end-lx->p_double.v_start;
  for (i = 0; i < rhs_len; i++){	//set the rhs-indices to global values
    *(px->mumps_irhs_loc + i) = g.my_rank * rhs_len + i+1;		//+1 due to fortran indexing
  }

  vector_double_copy(px->mumps_rhs_loc, eta, px->v_start, px->v_end, lx );

  // centralized solution, definitely mention it!
  if (g.my_rank == 0){
    // FIXME : do some sort of casting here, to avoid warnings at compile-time
    g.mumps_id.rhs = px->mumps_SOL;
  }

  g.mumps_id.job = 3;		// solve
  zmumps_c(&(g.mumps_id));

  int send_count = (lx->p_double.v_end-lx->p_double.v_start);
  MPI_Scatter(px->mumps_SOL, send_count, MPI_COMPLEX_double, phi, send_count, MPI_COMPLEX_double, 0, MPI_COMM_WORLD);	//scatter again to have px->x filled with mumps' solution

  END_MASTER(threading)
  SYNC_CORES(threading);
}

#endif
