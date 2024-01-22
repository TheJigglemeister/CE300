#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <stdio.h>

#include ".\libraries - mrbsm\ce300specific_MRBSM.h"

double** main_mrbsm()
/**change to a function that asks for fluid params, sys params, and then returns centroid positions**/
//double**MRBSM_Centroid(double*Xcentroid, double*Ycentroid)  //centroid matrix for MRBSM; initial position as input
{
    int h=0, i=0, j=0, k=0,  /**counters**/
    nComp,                   /**number of components**/
    order,                   /**with/without axial deformations; to be used for sizing**/
    nEqn,                    /**number of ODEs = N x 3 elems x 2 (main eqn and state eqn for RK4) **/
    nStep;                   /**number of RK steps**/

    /**necessary additional return variable**/
    printf("Modified Rigid Body-Spring Method, Small-Angle Approximation Formulation\n\n");

    nComp = 6;  //0 in diagonal found when nComp > 5; to be (a: solved / b: ignored and restricted)
    if (nComp <= 0 )
    {
        printf("\n\nInvalid component count. Program terminating . . .\n\n");   exit (999);
    }

    printf("Using a %d-component model . . . \n", nComp);

    order = 3*nComp;            /**for DOFs: u, v, theta**/
    nEqn = 2*order;             /**for each DOF, 2 equations; base and derivative**/
    //number of equations = 6 x number of components

    double *t, *y0, **resp, SPHx=0, SPHy=0;

    FILE*fp=fopen_s("loads.txt","r");

    if(fp==NULL)
    {
        printf("\nInput file not found. Program terminated.\n");
        exit (1);
    }

    int nData = CountData(fp);    /**first, count number of data, close file, reopen file*/
    printf("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n%d\nAAAAAAAAAAAAA\n\n",nData);

    double *uData = NewVector(nData) /**store hor. disp. of moving end here**/
         , *px = NewVector(nData)    /**after knowing data count, allocate memory for loads*/
         , *py = NewVector(nData)
         , *tx, *ty;                 /**placeholder; anticipates X(tx) and Y(ty) having different distributions in time**/

    int count = AllocateLoad(fp, uData, px, py, nData);  /**re-read input data, this time storing loads, while counting number of lines*/
    fclose(fp);

    double duration=20,  /**total loading duration**/
           f;            /**necessary parameter for allocation of t[i],
                         to be assigned with corresponding px[i] and py[i]**/
//           ,T;
    f=duration/nData;
//    T=nData/duration;
    t=NewVector(nData);

    for(i=0;i<nData;i++)
    {
        //time vector assigned
        t[i]=i*f;
        //loads converted from kN to N
        px[i]*=1000;
        py[i]*=1000;

        printf("\nuData[%d]\t=\t%lf\t\tpx[%d]\t=\t%lf\t\tpy[%d]\t=\t%lf",i,uData[i],i,px[i],i,py[i]);

    }
    nStep=t[nData-1]/0.0001;    /**arbitrarily decent amount of RK steps**/

    y0=NewVector(nEqn);
    FILE*outputfile;
    outputfile = fopen_s("results.csv","w");

    //for some reason, cannot create csv file after this line
    resp = LeapFrog(
               df_SmallAngleApprox,     /**function**/
               nComp,                   /**number of components**/
               0,t[nData-1],            /**endpoints**/
               y0,                      /**initial conditions**/
               nStep,                   /**number of RK steps**/
               SPHx,SPHy,               /**position of fluid particles**/
               nData,                   /**number of loading data**/
               px,py,t,t,px,px,px);              /**induced loads, application times
                                        currently expects px and py to occur simultaneously,
                                        but code is written in anticipation of
                                        non-simultaneous px and py**/


    double **u = NewMatrix(nData,nComp),
           **v = NewMatrix(nData,nComp),
           **theta = NewMatrix(nData,nComp),            //element xdisp, ydisp, rotation
           **x = NewMatrix(nData,nComp),
           **y = NewMatrix(nData,nComp),
           **phi = NewMatrix(nData,nComp);              //global xdisp, ydisp, rotation of each elem

    int nA, nB,  //indices for extracting solutions from resp matrix
        ppdt = nStep/nData;
    for(i=0;i<nData;i++)
    {
        for(j=0;j<nComp;j++)
        {
            nA = i*ppdt;
            nB = 1+j;
            u[i][j] = resp[nA][nB];

            nB += nComp;
            v[i][j] = resp[nA][nB];

            nB += nComp;
            theta[i][j] = resp[nA][nB];
        }
    }

    //phi is the total rotation of an element with respect to the initial axis;
    //top is rotation of topmost element
    double  *top = NewVector(nData), dTop;    i=0;    length = Length/nComp;  mass = Mass/nComp;
    for(k=0;k<nData;k++)
    {
        phi[k][i] = 0;
        top[k] = 0;
        for(i=0;i<nComp;i++)
        {
            top[k] += theta[k][i];
            for(j=0;j<=i;j++)
            phi[k][i] += theta[k][j];
        }
    }

    //calculation of global displacements of last element
    for(k=0;k<nData;k++)
    {
        for(i=0;i<nComp;i++)
        {
            x[k][i] = 0;    y[k][i] = 0;
            for(j=0;j<=i;j++)
            {
                //x and y coords of centroid of ith element
                /**CRUCIAL; CENTROIDS TO BE EXPORTED TO SPH300 ALGO**/
                x[k][i] += l(length,i,j)*sin(phi[k][j]) + u[k][j]*cos(phi[k][j-1]) + v[k][j]*sin(phi[k][j-1]) ;
                y[k][i] += l(length,i,j)*cos(phi[k][j]) - u[k][j]*sin(phi[k][j-1]) + v[k][j]*cos(phi[k][j-1]) ;
            }

        }
    }

    //exporting of output to a csv file
    fprintf(outputfile,"t,PEER_dTop,PEER_P,dTop,");
    for(j=0;j<nComp;j++){fprintf(outputfile,"u%d,", j+1);}
    for(j=0;j<nComp;j++){fprintf(outputfile,"v%d,", j+1);}
    for(j=0;j<nComp;j++){fprintf(outputfile,"theta%d,", j+1);}

    for(i=0;i<nData;i++)
    {
        fprintf(outputfile,"\n%lf,%lf,%lf,", t[i],uData[i],px[i]/1000);
        dTop = 1000 * (x[i][nComp-1] + (length/2)*top[i]);
        fprintf(outputfile,"%lf,", dTop);

        for(j=0;j<nComp;j++){fprintf(outputfile,"%lf,", u[i][j]);}
        for(j=0;j<nComp;j++){fprintf(outputfile,"%lf,", v[i][j]);}
        for(j=0;j<nComp;j++){fprintf(outputfile,"%lf,", theta[i][j]);}

    }

    printf("\n\nAnalysis done. Output can be found in \"results.csv\".\n\n");

    return 0;
    //return centroids;
}
