#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>


double LineInterp(double t,         /**interpolate at this value**/
                   int n,           /**number of data points**/
                   double *x,       /**using these points**/
                   double *y)
{
    int i;

    if((t < x[0] || (t > x[n-1])))
    {
        /**if t is outside of x, cannot interpolate**/
        printf("Cannot interpolate; out of bounds\n");
        return 0;
    }

    for(i=0;i<n;i++)
    {
        if(fabs(x[i]-t)<0.00001)
        {
            /**if t is sufficiently close to a given datapoint x,
            return corresponding ydata instead of interpolating**/
            return y[i];
        }
        if(x[i]>t)
        {
            double a = (t - x[i])  /   (x[i-1]-x[i]);
            double b = (t - x[i-1])/   (x[i]-x[i-1]);

            /**manipulated equation for interpolation**/
            return y[i-1]*a+y[i]*b;
        }
    }
    return 0;

}

double** NewMatrix(int row, int column) /**create matrix a[row][column]**/
{
    int i;
    double** a = calloc(row,sizeof(double*));
    for(i=0;i<row;i++) a[i]=calloc(column,sizeof(double));

    return a;
}

double* NewVector(int elems) /**create vector a[elems]**/
{
    double* a = calloc(elems,sizeof(double));

    return a;
}

double* GJLinSolv(double**A, int size, double*B) /**gauss-jordan elimination, return x from Ax = B**/
{
    int i, j, k;
    double *x= NewVector(size), **C= NewMatrix(size, size+1), f;

    /**take the augmented matrix form C = [A|B]**/
    for(i=0;i<size;i++)
    {
        for(j=0;j<size;j++)
            {
                C[i][j] = A[i][j];
            }
            C[i][size]=B[i];
    }

    /**apply gauss-jordan elimination**/
    for(i=0;i<size;i++)
    {
        if( fabsf(C[i][i]) < 0.0001)
        {
            printf("\nZero in diagonal found. Cannot perform Gauss-Jordan Elimintation.\n");
			exit(5);
        }

        for(j=0;j<size;j++)
        {
            if(i!=j)
            {
                f = C[j][i]/C[i][i];

                for(k=0;k<size+1;k++)
                {
                    C[j][k] = C[j][k] - (f*C[i][k]);
                }
            }
        }
    }

//    printf("\n\n");
    /**obtaining the solution by extracting B from [A|B]**/
    for(i=0;i<size;i++)
    {
        x[i]= C[i][size]/C[i][i];

    }

    return x;
}

int CountData(FILE* fp)/**count number of lines in file; POTENTIAL CAUSE FOR LONG RUNTIMES**/
{
    int count=1;   /**start at 1 to account for n[0]**/
    char c;

    /**count number of lines**/
    for (c = getc(fp); !feof(fp); c = getc(fp))
    {
        if (c == '\n')
        {
            count++;
        }
    }
    rewind(fp);

    return count;
}

int AllocateLoad(FILE* fp, double*u, double *px, double *py, int nData)
/**construct load vector; for verification of original RBS-DEM; POTENTIAL CAUSE FOR LONG RUNTIMES**/
{
    int i=0, count=0;

    double a;
//    fscanf(fp,"%lf", &a);
    for(i=0;i<nData;i++)
    {
        fscanf(fp,"%lf %lf %lf",&u[i],&px[i], &py[i]);    //per time basis
        count++;
    }
    rewind(fp);

    return count;
}

