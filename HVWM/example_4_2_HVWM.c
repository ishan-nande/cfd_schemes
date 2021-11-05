/*
*   Example 4.2 from
*   "An introduction to computationalf fluid dynamics,
*    the finite volume method", H.Versteeg,W.Malaskara, 2nd edition
*    k = thermal conductvity in W/(mK)
*    L = length of the rod in metres 
*    N = number of control volumes 
*    delta_x = distance between each control volume centre
*    T = temperature in Kelvin
*    T_A = fixed temperature at end A of the rod in Celsius
*    T_B = fixed temperature at end B of the rod
*    a_W = coefficient for west of node P
*    a_E = coefficient fo eastt of node P
*    S_U/S_P = linearized source term portions
*     q = heat generation source in W/m^3
*/

#include<stdio.h>

#include"lis.h"

#define N 5
#define L 0.02
#define T_A 100
#define T_B 200
#define k 0.5 //in W/mK
#define Area 1 //in m2
#define q 1000*1000 //in W/m3

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
   a_W = (k*Area)/delta_x;
   a_E = a_W; //east and west coefficients have the same formulation for all nodes 
   
   //common heat source linearization term for nodes 2 to N-1 
   float S_u = q * Area * delta_x;

   //S_u and S_P terms obtained at node 1 from values at A
   float S_u_A = (2*k*Area*T_A)/delta_x + (q*Area*delta_x); 
   float S_P_A = -(2*k*Area)/delta_x;

   //S_u and S_P terms obtained at node 5 from values at B
   float S_u_B = (2*k*Area*T_B)/delta_x + (q*Area*delta_x); 
   float S_P_B = -(2*k*Area)/delta_x ;

    //build the system
    for(LIS_INT i=0; i<N; i++)
    {
        if(i==0)//boundary condition on A
        {
            lis_matrix_set_value(LIS_INS_VALUE,i,i,a_E-S_P_A,A);
            lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-a_E,A);
            lis_vector_set_value(LIS_INS_VALUE,i,S_u_A,b);
        }
        else if(i==4)//boundary condition on B
        {
            lis_matrix_set_value(LIS_INS_VALUE,i,i,a_W-S_P_B,A);
            lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-a_W,A);
            lis_vector_set_value(LIS_INS_VALUE,i,S_u_B,b);            
        }
        else//internal nodes from 2 to N-1
        {
            lis_matrix_set_value(LIS_INS_VALUE,i,i,a_W+a_E,A);
            lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-a_W,A);
            lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-a_E,A);
            lis_vector_set_value(LIS_INS_VALUE,i,S_u,b);            
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
    lis_output_matrix(A, LIS_FMT_MM, "./data_output/matrx_example_4_2_HVWM");
    lis_output_vector(b, LIS_FMT_MM, "./data_output/b_vector_example_4_2_HVWM");
    lis_output_vector(x, LIS_FMT_MM, "./data_output/x_vector_example_4_2_HVWM");


    lis_matrix_destroy(A);
    lis_vector_destroy(b);
    lis_vector_destroy(x);
    lis_solver_destroy(solver);

    lis_finalize();    
   
    return 0;
}