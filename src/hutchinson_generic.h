#ifndef HUTCHINSON_double_HEADER
  #define HUTCHINSON_double_HEADER


  struct Thread;
  void block_hutchinson_driver_double( level_struct *l, struct Thread *threading );
  complex_double hutchinson_driver_double( level_struct *l, struct Thread *threading );
  void mlmc_hutchinson_rough_trace_double( level_struct *l, struct Thread *threading);
  complex_double mlmc_hutchinson_diver_double( level_struct *l, struct Thread *threading );
	
	complex_double mlmc_block_hutchinson_diver_double( level_struct *l, struct Thread *threading );
	
	
	void hutchinson_diver_double_init( level_struct *l, struct Thread *threading );
	void hutchinson_diver_double_alloc( level_struct *l, struct Thread *threading );
	void hutchinson_diver_double_free( level_struct *l, struct Thread *threading );

#endif
