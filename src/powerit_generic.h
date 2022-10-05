#ifndef POWERIT_double_HEADER
  #define POWERIT_double_HEADER

  void block_powerit_double_init_and_alloc( char* spec_type, char* op_name, int depth_bp_op, int nr_vecs, int nr_bpi_cycles, double bp_tol, level_struct* l, struct Thread* threading );
  void block_powerit_double( char* op_name, int depth_bp_op, level_struct *l, struct Thread *threading );
  void block_powerit_double_free( char* op_name, int depth_bp_op, level_struct* l, struct Thread* threading );

#endif
