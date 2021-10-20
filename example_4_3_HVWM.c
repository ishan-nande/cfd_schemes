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
#include<stdlib.h>
#include<string.h>

#include"linear_solver.h"
#include"common_utils.h"

#define N 5
#define L 1.0
#define T_B 100
#define n_square 25
#define T_amb 20

int main(int argc, char *argv[])
{ 

   float delta_x = L/N; 

   float a_W,a_E;
   a_W = 1/delta_x;
   a_E = a_W;   //east and west coefficients have the same formulation for all nodes 
   
   /*
   *  common heat source linearization term for nodes 2 to N-1 
   */  
   float S_u = n_square * T_amb * delta_x;
   float S_P = - n_square * delta_x;

   /*
   * S_u and S_P terms obtained at node 1 from values at B
   */
   float S_u_B = (2*T_B)/delta_x + ( n_square * T_amb * delta_x); 
   float S_P_B = -2/delta_x - (n_square*delta_x);

   /*
   * S_u and S_P terms obtained at node 5 from insulation boundary condition
   */
   float S_u_zero_q = n_square * T_amb * delta_x;
   float S_P_zero_q = -n_square * delta_x;

   //allocate memory for system coefficients and constants 
   float (*eq_system)[N+1];
   eq_system = (float (*)[N+1])calloc( N*(N+1), sizeof(float));
   if( eq_system==NULL)
    {
        puts("Mem alloc for equation system failed");
        exit(EXIT_FAILURE);
    }

    /*
    * build the discretized system matrix 
    */
   for(int i=0; i<N; i++)
   {
       if(i==0) //take into account boundary condition on B
       {
           eq_system[i][i] = a_E - S_P_B;
           eq_system[i][i+1] = -a_E;
           eq_system[i][N] = S_u_B;
       }
       else if(i==4) //take into account zero heat flux boundary condition on other end of rod 
       {
            eq_system[i][i] = a_W - S_P_zero_q;
            eq_system[i][i-1] = -a_W;
            eq_system[i][N] = S_u_zero_q;          
       }
       else //nodes from i+1 to N-1
       {
           eq_system[i][i] = a_W + a_E - S_P; //a_P forrmulation
           eq_system[i][i-1] = -a_W;
           eq_system[i][i+1] = -a_E;
           eq_system[i][N] = S_u;

       }
   }

    
    //allocate memory for forward eliminated matrix
    float (*ptr)[N+1];
    ptr=( float (*)[N+1])calloc(N*(N+1),sizeof(float));
    if( ptr==NULL)
    {
        puts("Mem alloc for ptr failed");
        exit(EXIT_FAILURE);
    }

    //copy elements of the discretized system to the forward elimination matrix
    float (*fw_elim_mat)[N+1] = memcpy( ptr, eq_system, sizeof(float)*N*(N+1) );
    
    //allocate memory for solution vector x
    float * x;
    x = (float *)calloc( N, sizeof(float) );
    if( x==NULL)
    {
        puts("Mem alloc for x vector failed");
        exit(EXIT_FAILURE);
    }

    tdma_forward_elimination(N, N+1, fw_elim_mat);
    back_substitution(N, N+1, x, fw_elim_mat);
    //print_array(N, N+1, fw_elim_mat);
    print_vector(N, x);

    free(fw_elim_mat);
    free(x);
    free(eq_system);
    return 0;
}