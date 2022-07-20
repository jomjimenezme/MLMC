#ifndef MUMPS_double_HEADER
  #define MUMPS_double_HEADER

  #define ICNTL(I) icntl[(I) -1]	//macro according to docu //bridges from fortran indices to c

  void mumps_setup_double(level_struct *l, struct Thread *threading);

  void mumps_solve_double( vector_double phi, vector_double Dphi, vector_double eta,
                              int res, level_struct *l, struct Thread *threading );


#endif
