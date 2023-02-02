#ifndef HUTCHINSON_double_HEADER
  #define HUTCHINSON_double_HEADER


  struct Thread;

  complex_double mlmc_hutchinson_diver_double( level_struct *l, struct Thread *threading );

  void hutchinson_diver_double_init( level_struct *l, struct Thread *threading );
  void hutchinson_diver_double_alloc( level_struct *l, struct Thread *threading );
  void hutchinson_diver_double_free( level_struct *l, struct Thread *threading );
  
  complex_double hutchinson_driver_double( level_struct *l, struct Thread *threading );
  complex_double mlmc_hutchinson_driver_double( level_struct *l, struct Thread *threading );
  complex_double split_mlmc_hutchinson_driver_double( level_struct *l, struct Thread *threading );

  // different operators for Hutchinson and block power iteration
  complex_double hutchinson_split_orthogonal( level_struct *l, hutchinson_double_struct* h, struct Thread *threading );
  
  
  void apply_P_double( vector_double out, vector_double in, level_struct* l, struct Thread *threading );
  void apply_R_double( vector_double out, vector_double in, level_struct* l, struct Thread *threading );

  int apply_solver_powerit_double( level_struct* l, struct Thread *threading );
#endif
