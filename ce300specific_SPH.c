#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "generalcodes_SPH.h"
#include "ce300specific_SPH.h"

#define SQR(x) (x*x)
#define pi() acos(-1)

double F(double t, int n, double *tLoad, double *pLoad)         /**load interpolation**/
{
    return LineInterp(t,n,tLoad,pLoad);
//    return .0100;
}

double* df_SmallAngleApprox(double t,
                            int nComp,  /**number of components**/
                            double *y,  /**y is of size "order" (6 x no. of elems)**/
                            int nData, double *px, double *py, double *tx, double *ty, double*FxTip, double*FyTip, double*MzTip)/**load vector; number of load data, load values, assigned time interval values**/
{
    int order, nEqn;
    order = 3*nComp;
    nEqn = 2*order;
    double *theta;

    double *dy=NewVector(nEqn);               /**answer vector**/

    int h=0, i=0, j=0, k=0;             /**counters**/

    double **M = NewMatrix(order,order) /**mass matrix**/,
            *R = NewVector(order)       /**load vector**/,
            *u                          /**horizontal displacement**/,
            *v                          /**vertical displacement**/,
            *udot                       /**horizontal velocity**/,
            *vdot                       /**vertical velocity**/,
            *theta_b                    /**theta_b is the relative rotation of the previous element**/,
            *theta_dot                  /**angular velocity**/;

    double *DY;
    /**spring constants for rubber; assuming homogeneous & isotropic, Kx=Ky=Kz=K**/

    /**for validation purposes, sample is Wight & Sozen 1973, No. 25.033 (east)**/
    /**Kx and Ky are at yielding, Kz is at cracking*/
    double Kx=2.575E+07, Ky=1.193E+09, Kz=2.866E+07;

    /**kinematic parameters; split initial guess y into components**/
    u           = &y[0];
    v           = &y[nComp];
    theta       = &y[2*nComp];
    udot        = theta+nComp;
    vdot        = udot+nComp;
    theta_dot   = vdot+nComp;

    theta_b=NewVector(nComp);
    /**theta_b[i] is the rotation of the (i-1)th component; basically theta_ij; necessary for some mass submatrices**/
    theta_b[0]=0;
    for (i=1;i<nComp;i++) theta_b[i]=theta[i-1];

    double **Px_U=NewMatrix(nComp,nComp), **Px_V=NewMatrix(nComp,nComp), **Px_Theta=NewMatrix(nComp,nComp), /**mass submatrices**/
           **Py_U=NewMatrix(nComp,nComp), **Py_V=NewMatrix(nComp,nComp), **Py_Theta=NewMatrix(nComp,nComp),
           **Mz_U=NewMatrix(nComp,nComp), **Mz_V=NewMatrix(nComp,nComp), **Mz_Theta=NewMatrix(nComp,nComp);

    mass = Mass/nComp;
    length = Length/nComp;

        /**Submatrices**/
        /**mass_Px_U      A**/
        for (i=0;i<nComp;i++){
            for (j=0;j<nComp;j++)
            {
                if (i>=j)   {for (k=i;k<nComp;k++)  Px_U[i][j]+=mass;}
                else        {for (k=j;k<nComp;k++)  Px_U[i][j]+=mass;}
            }
        }
        /**mass_Px_Theta  B**/
        for (i=0;i<nComp;i++){
            for (j=0;j<nComp;j++){
                if (i>=j)   {for (k=i;k<nComp;k++)   Px_Theta[i][j]+=mass*l(length,k,j);}
                else        {for (k=j;k<nComp;k++)   Px_Theta[i][j]+=mass*l(length,k,j);}
            }
        }
        /**mass_Mz_U      C**/
        for (i=0;i<nComp;i++){
            for (j=0;j<nComp;j++){
                if(i>=j){for (h=i;h<nComp;h++)  {for (k=h;k<nComp;k++)  Mz_U[i][j]+=mass*l(length,k,h);}}

                else    {for (h=i;h<j;h++)      {for (k=j;k<nComp;k++)  Mz_U[i][j]+=mass*l(length,k,h);}
                         for (h=j;h<nComp;h++)  {for (k=h;k<nComp;k++)  Mz_U[i][j]+=mass*l(length,k,h);}}
                             }
                         }
        /**mass_Mz_Theta  D**/
        for (i=0;i<nComp;i++){
            for (j=0;j<nComp;j++)
            {
                if      (i>j)
                {
                    for (h=i;h<nComp;h++)  {
                        for (k=h;k<nComp;k++)  {
                            Mz_Theta[i][j] += mass*l(length,k,h)*l(length,k,j);
                                           }
                                       }
                }
                else if (i==j)
                {
                    for (h=i;h<nComp;h++)  {
                        for (k=h;k<nComp;k++)  {
                            Mz_Theta[i][j] += mass*l(length,k,h)*l(length,k,j);
                                           }
                                       }
                    Mz_Theta[i][j] += mass*(length*length)/12.0;   //I[j]
                }
                else
                {
                    for (h=j;h<nComp;h++)  {
                        for (k=h;k<nComp;k++)  {
                            Mz_Theta[i][j] += mass*l(length,k,h)*l(length,k,j);
                                           }
                                       }

                    for (h=i;h<j;h++)  {
                        for (k=j;k<nComp;k++)  {
                            Mz_Theta[i][j] += mass*l(length,k,h)*l(length,k,j);
                                           }
                                       }
                    Mz_Theta[i][j] += mass*(length*length)/12.0;   //I[j]
                }
            }
        }

        /**mass_Px_V      Q**/
        for (i=0;i<nComp;i++){
            for (j=0;j<nComp;j++)
            {
                if (i>=j)   {for (k=i;k<nComp;k++)   Px_V[i][j]+=mass*theta_b[i];}
                else        {for (k=j;k<nComp;k++)   Px_V[i][j]+=mass*theta_b[i];}
            }
            }
        /**mass_Py_U      W**/
        for (i=0;i<nComp;i++){
            for (j=0;j<nComp;j++)
            {
                if (i>=j)   {for (k=i;k<nComp;k++)   Py_U[i][j]+=mass*theta_b[k];}
                else        {for (k=j;k<nComp;k++)   Py_U[i][j]+=mass*theta_b[k];}
            }
            }
        /**mass_Py_V      E**/
        for (i=0;i<nComp;i++){
            for (j=0;j<nComp;j++)
            {
                if (i>=j)   {for (k=i;k<nComp;k++){Py_V[i][j]+=mass;}}
                else        {for (k=j;k<nComp;k++){Py_V[i][j]+=mass;}}
            }
        }
        /**mass_Py_Theta  R**/
        for (i=0;i<nComp;i++){
            for (j=0;j<nComp;j++)
            {
                if (i>=j)
                {
                    for (k=i;k<nComp;k++)
                    {
                        Py_Theta[i][j]+=l(length,k,j)*mass*theta_b[k];
                    }
                }
                else
                {
                    for (k=j;k<nComp;k++)
                    {
                        Py_Theta[i][j]+=l(length,k,i)*mass*theta_b[k];
                    }

                }
            }
        }
        /**mass_Mz_V      T**/
        for (i=0;i<nComp;i++){
            for (j=0;j<nComp;j++){

                if (i>=j){
                        for (h=i;h<nComp;h++)
                        {
                            for (k=i;k<nComp;k++)
                            {
                                if (k==h && nComp!=1)             {Mz_V[i][j] += theta[k]*l(length,k,h)*mass;}
                                else if (k!=h || nComp==1)        {Mz_V[i][j] += theta_b[k]*l(length,k,h)*mass;}
                            }
                        }
                    }
                else{
                        for (k=j;k<nComp;k++)
                        {
                            if (k==i)   {Mz_V[i][j] += theta[k]*l(length,k,i)*mass;}
                            else        {Mz_V[i][j] += theta_b[k]*l(length,k,i)*mass;}
                        }


                        for (k=i+1;k<nComp;k++)
                        {
                            for (h=j;h<nComp;h++)
                            {
                                if (k==h)   {Mz_V[i][j] += theta[k]*l(length,k,h)*mass;}
                                else        {Mz_V[i][j] += theta_b[k]*l(length,k,h)*mass;}
                            }
                        }
                    }
            }
        }

        /**Modular mass matrix construction**/
        for (i=0;i<nComp;i++){
            for (j=0;j<nComp;j++){
                M[i][j]                 =    Px_U[i][j];        /**A**/
                M[i][j+nComp]           =    -Px_V[i][j];            /**Q; set to negative to apply MRBSM-SPH orientation**/
                M[i][j+2*nComp]         =    Px_Theta[i][j];    /**B**/

                M[i+nComp][j]           =    -Py_U[i][j];            /**W; set to negative to apply MRBSM-SPH orientation**/
                M[i+nComp][j+nComp]     =    Py_V[i][j];             /**E; does not change**/
                M[i+nComp][j+2*nComp]   =    -Py_Theta[i][j];        /**R; set to negative to apply MRBSM-SPH orientation**/

                M[i+2*nComp][j]         =    Mz_U[i][j];        /**C**/
                M[i+2*nComp][j+nComp]   =    -Mz_V[i][j];            /**T; set to negative to apply MRBSM-SPH orientation**/
                M[i+2*nComp][j+2*nComp] =    Mz_Theta[i][j];    /**D**/

            }
        }

        /**net force vector R, MRBSM only (original)**/
        for (i=0;i<nComp;i++){
            R[i]        = F(t, nData, tx, px)
                        - F(t, nData, ty, py)*theta[i]  /**doesn't change from RBS-DEM to MRBSM-SPH**/
                        - Kx*u[i];     /**Px - Py*theta - Kx**/

            R[i+nComp]  = F(t, nData, ty, py)
                        - F(t, nData, tx, px)*theta[i]  /**doesn't change from RBS-DEM to MRBSM-SPH**/
                        - Ky*v[i];     /**Py - Px*theta - Ky**/

            for (k=i;k<nComp;k++){
            R[i+2*nComp] += F(t, nData, tx, px)*length
                          - F(t, nData, ty, py)*length*theta[k] /**doesn't change from RBS-DEM to MRBSM-SPH**/
                          - Kz*theta[k];     /**PxL - PyL*theta - Kz**/
            }
        }

        /**net force vector R, MRSBM-SPH (ce300)**/
//        for(i=0;i<nComp;i++)
//        {
//            R[i] = FxTip[i];
//            R[i+nComp] = FyTip[i];
//            R[i+2*nComp] = MzTip[i];
//        }

        /**debug: manually check for when a 0 in the diagonal occurs**/
//        FILE*CheckMass;
//        CheckMass = fopen("MassMatrix.csv","a");
//
//        TestCount++;
//        fprintf(CheckMass,"\n%d\n", TestCount);
//        for(i=0;i<order;i++)
//        {
//            for(j=0;j<order;j++)
//            {
//                fprintf(CheckMass,"%lf,",M[i][j]);
//            }fprintf(CheckMass,"\n");
//        }
//        for(i=0;i<order;i++) {if(M[i][i]==0){printf("%lf\n", M[i][i]);}}
        /**can confirm, no zeros in diagonal, but submatrices Q, W, R and T go undefined after a while**/

        /**direct solution via gauss-jordan elimination**/
        DY = GJLinSolv(M,order,R);

        /**assembly technique that enables use of RK4 1st order ODE**/
        for(i=0;i<nComp;i++)
        {
            dy[i]         = udot[i];            /**du**/
            dy[i+nComp]   = vdot[i];            /**dv**/
            dy[i+2*nComp] = theta_dot[i];       /**dtheta**/
            dy[i+3*nComp] = DY[i];              /**ddu**/
            dy[i+4*nComp] = DY[i+nComp];        /**ddv**/
            dy[i+5*nComp] = DY[i+2*nComp];      /**ddtheta**/
        }


        /**free the mass submatrices**/
        free(Px_U); free(Px_V); free(Px_Theta);
        free(Py_U); free(Py_V); free(Py_Theta);
        free(Mz_U); free(Mz_V); free(Mz_Theta);

        /**since dy is the only tensor needed, free everything else**/
        free(M);
        free(theta_b);
        free(R);
        free(DY);

        return dy;
}

/**matrix RK4; algorithm to be used in the study

 * Classical Runge-Kutta method (RK4) for systems of ODEs
 *      y<k+1> = y<k> + (h/6)*(k1 + 2k2 + 2k3 + k4)
 * that returns: Matrix with the following columns
 *   | t  |  y1  | ... | yn | y1dot | ...  | yndot |
 **/
double** RK4(
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
{
    int order, nEqn;
    order = 3*nComp;
    nEqn = 2*order;

    int i, j;                                                   /**counters**/
	double h = (b-a)/((double)(nStep));                         /**interval width**/
	double **m = NewMatrix(nStep+1, nEqn+1);                    /**answer matrix; t, y1, y2, ..., yn, y1dot, y2dot, ..., yndot**/
                                                                /**size+1 accounts for column for t, as well as row t=0**/
	double ti, *yk = calloc(nEqn, sizeof(double));              /**ti - first column, yk - RK evaluations of y**/
	double *k1, *k2, *k3, *k4;                                  /**RK4 parameters**/
    double *ynew, *yi;	                                        /**pointers to answer matrix**/

	m[0][0] = a;                                                /**answer matrix, first vector takes points a to b**/
	memcpy(&m[0][1], ya, nEqn*sizeof(double));                  /**copy initial conditions to answer matrix**/

//    m[0] = ya;

	for(i=0; i<nStep; i++)  //this loop might have to be removed altogether, then have the SPH perform repeated function calls to RK4
	{
		ti = m[i][0];
		m[i+1][0] = a + (i+1)*h;
//		printf("%lf\n", m[i][7]);

		yi = &m[i][1];
		k1 = df(ti,nComp, yi, nData,px,py,tx,ty,FxTip, FyTip, MzTip);
//        printf("RK step %d - k1 done\n", i);

        for(j=0; j<nEqn; j++) yk[j] = yi[j] + (0.5*h)*k1[j];
		k2 = df(ti+(0.5*h),nComp, yk, nData,px,py,tx,ty,FxTip, FyTip, MzTip);
//        printf("RK step %d - k2 done\n", i);

        for(j=0; j<nEqn; j++) yk[j] = yi[j] + (0.5*h)*k2[j];
		k3 = df(ti+(0.5*h),nComp, yk, nData,px,py,tx,ty,FxTip, FyTip, MzTip);
//        printf("RK step %d - k3 done\n", i);

        for(j=0; j<nEqn; j++) yk[j] = yi[j] + h*k3[j];
		k4 = df(ti+h,nComp, yk, nData,px,py,tx,ty,FxTip, FyTip, MzTip);
//        printf("RK step %d - k4 done\n", i);

		ynew = &m[i+1][1];
		for(j=0; j<nEqn; j++)
        {
            ynew[j] = yi[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
//            printf("%lf\t", ynew[i]);
        }
//        printf("RK step %d - y_new done\n", i);
        printf("t = %lf\n", ti);
        free(k1);
        free(k2);
        free(k3);
        free(k4);

	}

	free(yk);

	return m;
}
