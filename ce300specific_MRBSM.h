#ifndef CE300SPECIFIC_H_INCLUDED
#define CE300SPECIFIC_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "generalcodes_MRBSM.h"


#define Length  0.876  /**specimen length(height) in m**/
#define Mass  97.467   /**specimen mass in kg**/

double length, mass;

typedef struct rigidbody /**rigid body information**/
{
    double dt;   /** time step **/
    char* fname; /** output file name **/

    double Kx;
    double Ky;
    double Kz;


    int n;       /** number of rigid body components**/
    double l;    /** rigid body length **/
    double m;    /** mass total **/

}
rigidbody;

typedef struct appliedforce /**information on current system state and integration algorithm**/
{
    int n;       /** number of rigid body components**/

    double l;    /** rigid body length **/
    double m;    /** mass total **/

    /** forces **/
    double *Fx;
    double *Fy;
    double *Mz;
}
appliedforce;

void default_mrbsm(rigidbody* R) /** default stuff for rigid body **/
;

void default_load(appliedforce* F, int nData) /** default stuff for load **/
;

double l(double L, int i, int j);      /**l_ij**/

double F(double t, int n, double *tLoad, double *pLoad);         /**load interpolation**/

double* SPH();

double* df_SmallAngleApprox(double t,
                            int nComp,  /**number of components**/
                            double *y,  /**y is of size "order" (6 x no. of elems)**/
                            int nData, double *px, double *py, double *tx, double *ty, double*FxTip, double*FyTip, double*MzTip)/**load vector; number of load data, load values, assigned time interval values**/
;
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

double** LeapFrog(
double * (*df)(double t,int nComp, double *y,
                            int nData, double *px,
                            double *py, double *tx, double *ty, double*FxTip, double*FyTip, double*MzTip),   /**function derivative**/
int nComp,                                                         /**number of components**/
double a, double b,					                               /**left and right end points (boundaries of time; basically duration)**/
double * ya,							                           /**initial condition**/
int nStep,									                       /**Number of steps**/
double SPHx,                                                       /**x-pos of fluid**/
double SPHy,                                                       /**y-pos of fluid**/
int nData, double *px, double *py, double *tx, double *ty, double*FxTip, double*FyTip, double*MzTip) /**load vector; number of load data, load values, assigned time interval values**/
;

double* df_general(double t,
                            int nComp,  /**number of components**/
                            double *y,  /**y is of size "order" (6 x no. of elems)**/
                            int nData, double *px, double *py, double *tx, double *ty, double*FxTip, double*FyTip, double*MzTip)/**load vector; number of load data, load values, assigned time interval values**/
;

double** LeapFrogStart(
double * (*df)(),   /**function derivative; currently solves for acceleration**/
int nComp,                                                         /**number of components**/
double a,
double * ya,							                           /**initial condition**/
//int nData, double *px, double *py, double *tx, double *ty, double*FxTip, double*FyTip, double*MzTip
appliedforce* forceinfo) /**load vector; number of load data, load values, assigned time interval values**/
;

double** LeapFrogStep(
double * (*df)(),   /**function derivative**/
int nComp,                                                         /**number of components**/
double a,
double * ya,							                           /**initial condition**/
int nData, double *px, double *py, double *tx, double *ty, double*FxTip, double*FyTip, double*MzTip) /**load vector; number of load data, load values, assigned time interval values**/
;


double* df_SmallAngleApproxNew(double t,
                            int nComp,  /**number of components**/
                            double *y,  /**y is of size "order" (6 x no. of elems)**/
                            int nData, double *px, double *py, double *tx, double *ty, rigidbody *RigidBody
                            )/**load vector; number of load data, load values, assigned time interval values**/
;
#endif // CE300SPECIFIC_H_INCLUDED
