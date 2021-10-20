/*
*   Example 4.1 from
*   "An introduction to computationalf fluid dynamics,
*    the finite volume method", H.Versteeg,W.Malaskara, 2nd edition
*    k = thermal conductvity in W/(mK)
*    L = length of the rod in metres 
*    N = number of control volumes 
*    delta_x = distance between each control volume centre
*    T = temperature in Kelvin
*    T_A = fixed temperature at end A of the rod in Celsius
*    T_B = fixed temperature at end B of the rod in Celsius
*    a_W = coefficient for west of node P
*    a_E = coefficient fo eastt of node P
*    S_U/S_P = linearized source term portions
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"linear_solver.h"
#include"common_utils.h"

#define N 5 
#define L 0.5
#define T_A 100
#define T_B 500
#define k 1000 //in W/mK
#define A 1e-2 //in m2

int main(int argc, char *argv[])
{ 

   float delta_x = L/N; 

   float a_W,a_E;
   a_W = (k*A)/delta_x;
   a_E = a_W;   //east and west coefficients have the same formulation
   
   float S_P = (-2*k*A)/delta_x;
   float S_u_A = (2*k*A*T_A)/delta_x; //S_u term generated from node A bc
   float S_u_B = (2*k*A*T_B)/delta_x; //S_u term generated from node B bc

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
       if(i==0) //boundary condition on A
       {
           eq_system[i][i] = a_E - S_P;
           eq_system[i][i+1] = -a_E;
           eq_system[i][N] = S_u_A;
       }
       else if(i==4) //boundary condition on B 
       {
            eq_system[i][i] = a_W - S_P;
            eq_system[i][i-1] = -a_W;
            eq_system[i][N] = S_u_B;          
       }
       else //nodes from i+1 to N-1
       {
           eq_system[i][i] = a_W + a_E; //a_P forrmulation
           eq_system[i][i-1] = -a_W;
           eq_system[i][i+1] = -a_E;
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