#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



#include <stdio.h>
#include "main-sph.h"
#include "main-mrbsm.h"


void main_combo()
{

}

int main()
{
    int i, j, k, l;

goto yeet;
/**SPH ALGORITHM**/
{
    /**SPH INPUT PARAMETERS**/
    int argc;   char** argv;
    double t_sph=0;
    sim_param_t params;

    /**initialize**/
    /**check if SPH parameters are available; if error, check sources**/
    if (get_params(argc, argv, &params) != 0)
    exit(-1);
    sim_state_t* state = init_particles(&params);
    FILE* fp = fopen(params.fname, "w");

    /**data headers + initial step t0**/
    write_header(fp, state->n);
    write_frame_data(fp, state->n, state->x, NULL);

    /**leapfrog starting line, different from whole routine (startup mechanism)**/
    compute_accel(state, &params);
    leapfrog_start(state, params.dt, &params);
    check_state(state);

    /**SPH ROUTINE**/
    for(t_sph=params.dt;t_sph<4.0;t_sph+=params.dt)
    {
        main_sph(argc, argv, params, state, t_sph);     //contains leapfrog, collision, damping, contact detection
        printf("SPH -- ti = %f\n", t_sph);              //print t_i in cmd screen to show that code is running
        write_frame_data(fp,state->n,state->x,NULL);    //print output in a .csv file
    }

}

/**MRBSM ALGORITHM**/
{
    rigidbody*R;
    R = calloc(sizeof(rigidbody),1);
    default_mrbsm(R);

    int nComp=5, nData;
    double**resp;
    resp = calloc(sizeof(double*),3);
    for(i=0;i<3;i++)
    {
        resp[i]=calloc(sizeof(double),nComp);
    }
        /** store load info; future: remove/replace with "reactive" data storing, from SPH output**/
    FILE*loadinput=fopen("loads.txt","r");
    fscanf(loadinput,"%d",&nData);  //record the first line of the input file, store as data count
    appliedforce*F;
    F=calloc(sizeof(appliedforce),1);
    default_load(F, nData);

    double *t, *nothing=calloc(sizeof(double),nData);
//    printf("%d\n\n", nData);

    //for verification if load data is properly read
//    AllocateLoad(loadinput,nothing,F->Fx,F->Fy,nData);


    double   *x = calloc(sizeof(double),nComp*6);     //
    double  *vh = calloc(sizeof(double),nComp*6);     //
    double *ddx = calloc(sizeof(double),nComp*6);     //ddx = a = (F-Kx)/M


    double duration = 20.0, freq;
    freq = duration/nData;
    t = calloc(sizeof(double), nData);

    for(i=0;i<nData;i++)
    {
//        t[i]= i*freq;
//        F->Fx[i]*=1000;
//        F->Fy[i]*=1000;
//        printf("\n%d\t%lf\t%lf\t%lf\t%lf",i, t[i], nothing[i], F->Fx[i], F->Fy[i]);
    }

    int nStep;
    nStep=t[nData-1]/0.0001;    /**arbitrarily decent amount of RK steps**/


//    printf("\n\n%lf", t[nData-1]);
    double currtime=0.0, dt=0.01;
    int nEqn = 6*nComp;

    double *k1, *k2, *yk,*yi, *ynew;
    main_mrbsm();


    double test;
    i=0;

//    for(currtime=0;currtime<t[nData-2];currtime+=dt)
//    {
//        if(currtime < t[nData-1])
//        {
//            k1 = df_SmallAngleApproxNew(currtime+0.5*dt,nComp, x, nData,F->Fx,F->Fy,t,t,R);
//            for(j=nEqn/2.0;j<nEqn;j++){
//            yk[j] = yi[j] + k1[j]*0.5*dt;
//        }   //option to store velocity values here (midsteps)
//
//            //positions, updated at fullstep
//            k2 = df_SmallAngleApproxNew(currtime+dt,nComp, x, nData,F->Fx,F->Fy,t,t,R);
//            for(j=0;j<nEqn/2.0;j++){
//            ynew[j] = yi[j] + k2[j]*dt;
//        }
//
//            //velocities, re-updated at another halfstep
//            for(j=nEqn/2.0;j<nEqn;j++){
//            ynew[j] = yi[j] + k1[j]*0.5*dt;
//        }   //option to store velocity values here (fullsteps)
//        free(k1);
//        free(k2);
//        }
//    }

//    while(currtime < t[nData-1])    //"-2" instead of "-1" guarantees simulation is always within bounds, albeit last datapoint never gets used
//    {
//        currtime+=dt;
//            k1 = df_SmallAngleApproxNew(currtime+0.5*dt,nComp, x, nData,F->Fx,F->Fy,t,t,R);
//            for(j=nEqn/2.0;j<nEqn;j++){
//            yk[j] = yi[j] + k1[j]*0.5*dt;
//        }   //option to store velocity values here (midsteps)
//
//            //positions, updated at fullstep
//            k2 = df_SmallAngleApproxNew(currtime+dt,nComp, x, nData,F->Fx,F->Fy,t,t,R);
//            for(j=0;j<nEqn/2.0;j++){
//            ynew[j] = yi[j] + k2[j]*dt;
//        }
//
//            //velocities, re-updated at another halfstep
//            for(j=nEqn/2.0;j<nEqn;j++){
//            ynew[j] = yi[j] + k1[j]*0.5*dt;
//        }   //option to store velocity values here (fullsteps)
//        free(k1);
//        free(k2);
//
//    }
//    printf("\n\n%lf\t\t%lf",currtime, t[nData-1]);
return 1235;
//return 5;
//
//double A[5] = {2,4,8,16,32};
//double B[5] = {1,2,3,4,5};
//double C;
//
//
//printf("\n\n\nC = %lf\n\n",C);


}


yeet:

    return 0;
}

void trash()
{
    int i;

    /**MRBSM INPUT PARAMETERS**/
    double mass, length, Kx, Ky, Kz;
    //change this such that input is configuration; loop with time
    /**MRBSM ROUTINE**/

    double**resp;

    rigidbody*R;
    R=calloc(sizeof(rigidbody),1);
    default_mrbsm(R);

    appliedforce*F;
    F=calloc(sizeof(appliedforce),1);

    int nComp = 5, nData;

    resp = calloc(sizeof(double*),3);   //X,V,A of x,y,z

    double SPHx, SPHy,*px,*py,*t;


    for(i=0;i<3;i++){
    resp[i] = calloc(sizeof(double),nComp);
    }

    /** store load info; future: remove/replace with "reactive" data storing, from SPH output**/
    FILE*loadinput=fopen("loads.txt","r");
    fscanf(loadinput,"%d",&nData);  //record the first line of the input file, store as data count
    default_load(F, nData);

    double *nothing=calloc(sizeof(double),nData);
    printf("%d\n\n", nData);

    //for verification if load data is properly read
    AllocateLoad(loadinput,nothing,F->Fx,F->Fy,nData);


    double   *x = calloc(sizeof(double),nComp*6);     //
    double  *vh = calloc(sizeof(double),nComp*6);     //
    double *ddx = calloc(sizeof(double),nComp*6);     //ddx = a = (F-Kx)/M


    double tmax = 20.0, freq;
    freq = tmax/nData;
    t = calloc(sizeof(double), nData);


    for(i=0;i<nData;i++)
    {
        t[i]= i*freq;
    }
    int nStep;
    nStep=t[nData-1]/0.0001;    /**arbitrarily decent amount of RK steps**/


    for(i=0;i<nData;i++)
    {
//        printf("\n%d\t%lf\t%lf\t%lf\t%lf",i, t[i], nothing[i], F->Fx[i], F->Fy[i]);
    }

    printf("\n\n");
    double count, dt=0.01;

    for(count=0;count<tmax;count+=dt)
    {
//        printf("\n%d\t%lf\t%lf\t%lf\t%lf",i, dt, nothing[i], F->Fx[i], F->Fy[i]);
//        ddx = df_SmallAngleApprox(t[i],nComp,x,nData,F->Fx,F->Fy,t,t);
//        printf("%lf\n",ddx[0]);
    }


    /**leapfrog starting line, different from whole routine (startup mechanism)**/


return 0;

}
