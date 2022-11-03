#ifndef HUTCHINSON_double_HEADER
  #define HUTCHINSON_double_HEADER


  struct Thread;
  
    complex_double mlmc_hutchinson_diver_double( level_struct *l, struct Thread *threading );
		
	
	void hutchinson_diver_double_init( level_struct *l, struct Thread *threading );
	void hutchinson_diver_double_alloc( level_struct *l, struct Thread *threading );
	void hutchinson_diver_double_free( level_struct *l, struct Thread *threading );

#endif
