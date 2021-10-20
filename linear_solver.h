
#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

void forward_elimination( int R, int C, float (*out)[C]);

void back_substitution(int R, int C, float * out, const float (*in)[C] );

void naive_gauss(int R, int C, float (*fw_elim_mat)[C], float * sol_vector );

void tdma_forward_elimination(int R, int C, float(*mat)[C] );

#endif
