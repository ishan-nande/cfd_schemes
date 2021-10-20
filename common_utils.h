#ifndef COMMON_UTILS_H
#define COMMON_UTILS_H

void print_array(int r, int c, const float (*ar)[c]);

void print_vector(int , const float * vec );

void copy_array(int r, int c, float const source[r][c], float (*dest)[c]);

#endif
