#include<stdio.h>
#include"common_utils.h"

//print contents of the array
void print_array(int r, int c, const float (*ar)[c])
{
    for(int i=0; i<r; i++)
    {
        for(int j=0; j<c; j++)
        {
            printf("%.4f ", ar[i][j]);
        }
        printf("\n");
    }
}

//print vector contents 
void print_vector( int r, const float * ar)
{
    for(int i=0; i<r; i++)
        printf("%.6f ", ar[i]);
}

//copy contents of an array memory block into another array memory block
void copy_array(int r, int c, float const source[r][c], float (*dest)[c])
{
    for(int i=0;i<r;i++)
    {
        for(int j=0;j<c; j++)
        {
            dest[i][j] = source[i][j];
        }
    }
}
//write system matrix and solution vector to files
void write_to_file( int r, int c, const float * sol_vector, const float (*sys_aug_mat)[c])
{
    
}

