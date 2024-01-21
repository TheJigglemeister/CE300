#ifndef GENERALCODES_H_INCLUDED
#define GENERALCODES_H_INCLUDED
double LineInterp(double t,         /**interpolate at this value**/
                   int n,           /**number of data points**/
                   double *x,       /**using these points**/
                   double *y);
double** NewMatrix(int row, int column); /**create matrix a[row][column]**/
double* NewVector(int elems); /**create vector a[elems]**/
double* GJLinSolv(double**A, int size, double*B); /**gauss-jordan elimination, return x from Ax = B**/

int CountData(FILE* fp); /**count number of lines in file; POTENTIAL CAUSE FOR LONG RUNTIMES**/
int AllocateLoad(FILE* fp, double*u, double *px, double *py, int nData); /**construct load vector; for verification of original RBS-DEM; POTENTIAL CAUSE FOR LONG RUNTIMES**/





#endif // GENERALCODES_H_INCLUDED
