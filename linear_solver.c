#include<stdio.h>
#include<stdlib.h>

#include"linear_solver.h"
#include"common_utils.h"


//forward elimination of a system
void forward_elimination( int R, int C, float (*out)[C])
{
     for(int k=0; k<R; k++)               //k tracks the row containg the pivot element  
     {
         for(int i=k+1; i<R; i++)         //i tracks the row on which normalization will be performed
         {                                     
             float factor = out[i][k]/out[k][k]; //normalization factor, out[k][k] is the current pivot 
            
             for(int j=k; j<C; j++)      //j tracks the column elements for pivot element row and the normalization row  
            {
                 if( j==(C-1))  //operate on constant value of the system
                 {
                    out[i][j] = out[i][j]-factor*out[k][j]; //column operations on b
                 }
                 else
                 {
                    out[i][j] = out[i][j]-factor*out[k][j]; //column operations on a
                 }
             }
         }
     }      
}

//backward substitution for x vector solution
void back_substitution(int R, int C, float * out, const float (*in)[C] )
{
    //last variable  calculation
    out[R-1] = in[R-1][C-1] / in[R-1][C-2];

    for(int i=R-2; i>=0; i--)   //row tracker 
    {
        float sum = in[i][C-1];

        for(int j=R-1; j>=i+1; j--)
        {
            sum = sum - in[i][j]*out[j];    //column operations 
        }

        out[i] = sum / in[i][i];    //i_th variable calculation
    }
}

//naive Gauss elimination 
void naive_gauss(int R, int C, float (*fw_elim_mat)[C], float * sol_vector )
{
    forward_elimination(R, C, fw_elim_mat);
    back_substitution(R, C, sol_vector, fw_elim_mat);
}
 
//TDMA forward elimination
void tdma_forward_elimination(int R, int C, float(*mat)[C] )
{
    /*
    *   operations start from row 1
    */
    for(int i=1; i<R; i++) 
    {
        mat[i][i] = mat[i][i] - (mat[i][i-1]*mat[i-1][i])/mat[i-1][i-1]; //set diagonal elements 
        mat[i][C-1]  = mat[i][C-1] - (mat[i-1][C-1]*mat[i][i-1])/mat[i-1][i-1]; //set constants in the last column
        mat[i][i-1] = 0; //set lower diagonal term to zero
    }
}