#ifndef POWERIT_double_HEADER
  #define POWERIT_double_HEADER

  void block_powerit_double_init_and_alloc( int spec_type, int op_id, int depth_bp_op, int nr_vecs, int nr_bpi_cycles, double bp_tol, level_struct* l, struct Thread* threading );
  void block_powerit_double( int op_id, int depth_bp_op, level_struct *l, struct Thread* threading );
  void block_powerit_double_free( level_struct* l, struct Thread* threading );

  void block_powerit_driver_double( level_struct* l, struct Thread* threading );
  
  void blind_bp_op_double_apply( op_id, lx, threading );
  void powerit_non_diff_op( level_struct *l, int i, struct Thread *threading );
  void powerit_diff_op( level_struct *l, int i, struct Thread *threading );
  void powerit_split_op( level_struct *l, int i, struct Thread *threading );
  
  // this functional aplies anoperator to a bp vector 
  void (*apply_to_one_vector)();

#endif
