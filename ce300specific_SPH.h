#ifndef CE300SPECIFIC_H_INCLUDED
#define CE300SPECIFIC_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "generalcodes_SPH.h"


#define Length  0.876  /**specimen length(height) in m**/
#define Mass  97.467   /**specimen mass in kg**/

double length, mass;

double l(double L, int i, int j);      /**l_ij**/

double F(double t, int n, double *tLoad, double *pLoad);         /**load interpolation**/

double* SPH();

double* df_SmallAngleApprox(double t,
                            int nComp,  /**number of components**/
                            double *y,  /**y is of size "order" (6 x no. of elems)**/
                            int nData, double *px, double *py, double *tx, double *ty/**load vector; number of load data, load values, assigned time interval values**/
                            , double*FxTip, double*FyTip, double*MzTip); /**SPH particle loads**/
double** RK4(
double * (*df)(double t,int nComp, double *y,
                            int nData, double *px,
                            double *py, double *tx, double *ty, double*FxTip, double*FyTip, double*MzTip),   /**function derivative**/
int nComp,                                                         /**number of components**/
double a, double b,					                               /**left and right end points**/
double * ya,							                           /**initial condition**/
int nStep,									                       /**Number of steps**/
double SPHx,                                                       /**x-pos of fluid**/
double SPHy,                                                       /**y-pos of fluid**/
int nData, double *px, double *py, double *tx, double *ty, double*FxTip, double*FyTip, double*MzTip)         /**load vector; number of load data, load values, assigned time interval values**/
;
#endif // CE300SPECIFIC_H_INCLUDED
