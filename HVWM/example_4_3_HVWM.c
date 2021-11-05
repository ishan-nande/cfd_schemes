/*
*   Example 4.3 from
*   "An introduction to computationalf fluid dynamics,
*    the finite volume method", H.Versteeg,W.Malaskara, 2nd edition
*    L = length of the rod in metres 
*    N = number of control volumes 
*    delta_x = distance between each control volume centre
*    T_B = fixed temperature at end B of the rod in Celsius
*    T_amb = temperature surrounding the rod in Celsius
*    a_W = coefficient for west of node P
*    a_E = coefficient fo eastt of node P
*    S_U/S_P = linearized source term portions
*    n_square = constant term from hP/ka in 1/m2
*/

#include<stdio.h>

#include"lis.h"

#define N 5
#define L 1.0
#define T_B 100
#define n_square 25
#define T_amb 20

LIS_INT main(int argc, char *argv[])
{ 
    LIS_Comm comm = LIS_COMM_WORLD;
    LIS_VECTOR b,x;
    LIS_MATRIX A;
    LIS_SOLVER solver;
    
    lis_initialize(&argc, &argv);

    lis_matrix_create(comm, &A);
    lis_matrix_set_size(A,0, (LIS_INT)N);

    lis_vector_create(comm, &b);
    lis_vector_set_size(b,0,(LIS_INT)N);
    
   float delta_x = L/N; 

   float a_W,a_E;
   a_W = 1/delta_x;
   a_E = a_W; //east and west coefficients have the same formulation for all nodes 
   
   //common heat source linearization term for nodes 2 to N-1  
   float S_u = n_square * T_amb * delta_x;
   float S_P = - n_square * delta_x;

   //S_u and S_P terms obtained at node 1 from values at B
   float S_u_B = (2*T_B)/delta_x + ( n_square * T_amb * delta_x); 
   float S_P_B = -2/delta_x - (n_square*delta_x);

   //S_u and S_P terms obtained at node 5 from insulation boundary condition
   float S_u_zero_q = n_square * T_amb * delta_x;
   float S_P_zero_q = -n_square * delta_x;

    //build the system
    for(LIS_INT i=0; i<N; i++)
    {
        if(i==0)//boundary condition on B
        {
            lis_matrix_set_value(LIS_INS_VALUE,i,i,a_E-S_P_B,A);
            lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-a_E,A);
            lis_vector_set_value(LIS_INS_VALUE,i, S_u_B,b);  
        }
        else if(i==4)//zero heat flux boundary condition on other end of rod 
        {
            lis_matrix_set_value(LIS_INS_VALUE,i,i,a_W-S_P_zero_q,A);
            lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-a_W,A);
            lis_vector_set_value(LIS_INS_VALUE,i, S_u_zero_q,b);            
        }
        else//nodes from i+1 to N-1
        {
            lis_matrix_set_value(LIS_INS_VALUE,i,i,a_W+a_E-S_P,A);//a_P formulation
            lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-a_W,A);
            lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-a_E,A);
            lis_vector_set_value(LIS_INS_VALUE,i, S_u,b);      
        }
    }

    lis_matrix_set_type(A, LIS_MATRIX_CSR);
    lis_matrix_assemble(A);

    //solver setup 
    lis_solver_create(&solver);
    lis_solver_set_option("-i gs -p none -print all -maxiter[1000]", solver);
    lis_vector_duplicate(A,&x);
    lis_solve(A,b,x,solver);
    //lis_solver_get_iter(solver,&iter);
    //lis_printf(comm,"number of iterations = %D\n",iter);
          
    //data output 
    lis_output_matrix(A, LIS_FMT_MM, "./data_output/matrx_example_4_3_HVWM");
    lis_output_vector(b, LIS_FMT_MM, "./data_output/b_vector_example_4_3_HVWM");
    lis_output_vector(x, LIS_FMT_MM, "./data_output/x_vector_example_4_3_HVWM");

    lis_matrix_destroy(A);
    lis_vector_destroy(b);
    lis_vector_destroy(x);
    lis_solver_destroy(solver);

    lis_finalize();
     
    return 0;
}