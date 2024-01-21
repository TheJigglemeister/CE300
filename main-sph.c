/**CE 300 - SPH code, as implemented by David Bindel (Fall 2011) - Application of Parallel Computers (CS 5200)**/
/**http://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf**/
/**for 2D, as of Sep 2021**/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>


//#include (up 1 level, include \MRBSM\MRBSM.h)
#include ".\libraries - mrbsm\ce300specific_MRBSM.h"

#define VERSION_TAG "SPHView00"

typedef struct sim_param_t /**parameters describing the simulation**/
{
    char* fname; /* File name */
    int nframes; /* Number of frames */
    int npframe; /* Steps per frame */
    double h; /* Particle size */
    double dt; /* Time step */
    double rho0; /* Reference density */
    double k; /* Bulk modulus */
    double mu; /* Viscosity */
    double g; /* Gravity strength */
}
sim_param_t;

typedef struct sim_state_t /**information on current system state and integration algorithm**/
{
    int n; /* Number of particles */
    double mass; /* Particle mass */
    double* restrict rho; /* Densities */

    /**horizontal - even, vertical - odd**/
    double* restrict x; /* Positions */
    double* restrict vh; /* Velocities (half step) */
    double* restrict v; /* Velocities (full step) */
    double* restrict a; /* Acceleration */
}
sim_state_t;

sim_state_t* alloc_state(int n)
{
    sim_state_t* s = (sim_state_t*) calloc(1, sizeof(sim_state_t));
    s->n   =  n;    /**number of particles**/
    s->rho = (double*) calloc(  n, sizeof(double));

    /**combines horizontal (even) and vertical (odd) components of motion into one vector**/
    s->x   = (double*) calloc(2*n, sizeof(double));
    s->vh  = (double*) calloc(2*n, sizeof(double));
    s->v   = (double*) calloc(2*n, sizeof(double));
    s->a   = (double*) calloc(2*n, sizeof(double));
    return s;
}

void free_state(sim_state_t* s)
{
    free(s->a);
    free(s->v);
    free(s->vh);
    free(s->x);
    free(s->rho);
    free(s);
}

typedef int (*domain_fun_t)(double, double); /**1 = occupied by fluid, 0 elsewhere**/

static void default_params(sim_param_t* params)
/**SPH particle parameters**/
/**validation params are sufficient for actual simulation**/
{
    params->fname = "output.csv";
////    params->nframes = 400;
//    params->nframes = 5e4;  //validation
////    params->npframe = 100;
//    params->npframe = 1;    //validation
////    params->dt = 1e-4;
//    params->dt = 5e-5;   //validation ; [s]
////    params->h = 5e-2;
//    params->h = 5e-2;   //validation
//    params->rho0 = 1000;    //water ; [kg/m^3]
////    params->k = 1000;
//    params->k = 1000;   //validation; bulk modulus; normal 2.1e9 N/m^2 --> 2.1e3 MN/m^2 or 2.1e3 N/mm^2?
//    /**MAY 30, 2022: this is NOT the bulk modulus; this is from pi = k(rhoi-rho0) such that cs >> vi; units seem to be m^2/s^2**/
//    /**Ideal Gas Law: k = rho*(dP/drho) --> pi = k(rhoi-rho0)/rhoi **/
//    params->mu = .1;    //validation ; [N-s/m^2]
//    params->g = 9.8;    // [m/s]

    params->nframes = 400;
    params->npframe = 1;
    params->dt = 1e-3;
    params->h = 5e-2;
    params->rho0 = 1000;
    params->k = 1e3;
    params->mu = 0.1;
    params->g = 9.8;
}


void compute_density(sim_state_t* s, sim_param_t* params)
/**direct application of SPH formula, Wpoly6**/
{
    int n = s->n;
    double* restrict rho = s->rho;
    const double* restrict x = s->x;
    double h = params->h;
    double h2 = h*h;
    double h8 = ( h2*h2 )*( h2*h2 );
    double C = 4 * s->mass / M_PI / h8;
    memset(rho, 0, n*sizeof(double));

    for (int i = 0; i < n; i++)
    {
        rho[i] += 4 * s->mass / M_PI / h2;
        for (int j = i+1; j < n; j++)
        {
            double dx = x[2*i+0]-x[2*j+0];
            double dy = x[2*i+1]-x[2*j+1];
            double r2 = dx*dx + dy*dy;
            double z = h2-r2;

            if (z > 0) /**check for neighbors, basically**/
            {
            double rho_ij = C*z*z*z;
            rho[i] += rho_ij;
            rho[j] += rho_ij;
            }
        }
    }
}

void compute_accel(sim_state_t* state, sim_param_t* params)
/**computed via (sum of interacting forces/density) + gravity**/
{
    // Unpack basic parameters
    const double h = params->h;
    const double rho0 = params->rho0;
    const double k = params->k;
    const double mu = params->mu;
    const double g = params->g;
    const double mass = state->mass;
    const double h2 = h*h;

    // Unpack system state
    const double* restrict rho = state->rho;
    const double* restrict x = state->x;
    const double* restrict v = state->v;
    double* restrict a = state->a;
    int n = state->n;

    // Compute density and color
    compute_density(state, params);

    // Start with gravity and surface forces
    for (int i = 0; i < n; i++)
    {
        a[2*i+0] = 0;
        a[2*i+1] = -g;
    }


    // Constants for interaction term
    double C0 = mass / M_PI / ( (h2)*(h2) );
    double Cp = 15*k;
    double Cv = -40*mu;

    // Now compute interaction forces
    for (int i = 0; i < n; i++)
    {
        const double rhoi = rho[i];
        for (int j = i+1; j < n; j++)
        {
            double dx = x[2*i+0]-x[2*j+0];
            double dy = x[2*i+1]-x[2*j+1];
            double r2 = dx*dx + dy*dy;

            if (r2 < h2)
            {
                const double rhoj = rho[j];
                double q = sqrt(r2)/h;
                double u = 1-q;
                double w0 = C0 * u/rhoi/rhoj;    //incorporating rho_i to convert force to acceleration
                double wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
                double wv = w0 * Cv;
                double dvx = v[2*i+0]-v[2*j+0];
                double dvy = v[2*i+1]-v[2*j+1];
                a[2*i+0] += (wp*dx + wv*dvx);
                a[2*i+1] += (wp*dy + wv*dvy);
                a[2*j+0] -= (wp*dx + wv*dvx);
                a[2*j+1] -= (wp*dy + wv*dvy);
            }
        }
    }
}

static int damp_reflect(int which, double barrier, double* x, double* v, double* vh)
/**e < 1; velocity damped after inelastic collision**/
/**vertical barrier: which == 0, horizontal barrier: which == 1**/
/**just does stuff when particles hit a boundary**/
{
    /**vertical: which = 1**/
    /**horizontal: which = 0**/
    // Coefficient of resitiution
    const double DAMP = 0.75;

    // Ignore degenerate cases
    if (fabs(v[which]) == 0)
    return;   //arbitrary, specific return value for non-execution

    // Scale back the distance traveled based on time from collision
    double tbounce = (x[which]-barrier)/v[which];
    x[0] -= v[0]*(1-DAMP)*tbounce;
    x[1] -= v[1]*(1-DAMP)*tbounce;

    // Reflect the position and velocity
    x[which] = 2*barrier-x[which];
    v[which] = -v[which];
    vh[which] = -vh[which];

    // Damp the velocities
    v[0] *= DAMP; vh[0] *= DAMP;
    v[1] *= DAMP; vh[1] *= DAMP;
}

static void reflect_bc(sim_state_t* s
//                       ,mrbsm*R
                       ,sim_param_t*params)

/**boundary condition; since boundaries aren't modeled, function is made to enforce BCs**/
/**ACTUALLY dictates the shape of the boundary; EDIT THIS TO INCORPORATE MRBSM**/
{
    // Boundaries of the computational domain
    /**unit box confinement: 0<x<1 && 0<y<1**/
    const double XMIN = 0.0;
    const double XMAX = 1.0;
    const double YMIN = 0.0;
    const double YMAX = 1.0;
    /**validation, Koshizuka and Oka (1996)**/
//    const double XMIN = 0.0;
//    const double XMAX = 4.0;
//    const double YMIN = 0.0;
//    const double YMAX = 3.0;
    /**open box: 0<x<0.1m && 0<y<inf, but x changes on a per-rigidbody basis**/
//    double XMIN = 0;
//    double YMIN = 0;
//    double XMAX = 0.1;        //is actually moving boundary; should restate this
//    double YMAX = 1000000.0;    //absurdly high to denote free surface
    /**anticipated: open box, moving boundaries**/
//    double**XMAX = (double**)calloc(R->nComp, sizeof(double*)); //centroid of rigid blocks
//    double**YMAX = (double**)calloc(R->nComp, sizeof(double*)); //centroid of rigid blocks
//
//    for (i=0;i<nComp;i++)
//    {
//          XMAX[i] = (double*)calloc(nSubs, sizeof(double));
//          YMAX[i] = (double*)calloc(nSubs, sizeof(double));
//        //again, re-discretize rigid blocks into finite rigid body subdivision, performing geometric measurements again
//    }
//    subdivide(params->h,R,XMAX,YMAX);  //new upper bounds go into XMAX and YMAX matrices

    double* restrict vh = s->vh;
    double* restrict v = s->v;
    double* restrict x = s->x;
    int n = s->n;

    for (int i = 0; i < n; i++, x += 2, v += 2, vh += 2)
    {   //if x[0/1] exceeds min/max, it means collision with boundary occured; initiate damp_reflect
        if (x[0] < XMIN) damp_reflect(0, XMIN, x, v, vh);
        if (x[0] > XMAX) damp_reflect(0, XMAX, x, v, vh);   //expected: changing values because deformable boundary
        if (x[1] < YMIN) damp_reflect(1, YMIN, x, v, vh);
        if (x[1] > YMAX) damp_reflect(1, YMAX, x, v, vh);   //expected: never to happen because open top

//        probably something like:
//        if(x[0]<XMIN)       damp_reflect(0,XMIN,x,v,vh);  //same as original
//        int switchX =1, switchY = 1, j;
//        while(switchX == 0 && switchY == 0) //start at 0; if either x or y boundary is hit, end loop
//        {
//        //arrange MRBSM bounds in descending order, then check x>XMAX; first XMAX hit will force the loop to end,
//        //denoting ONE point of contact
//        for(j=(nSubs-1);j>=0;j--)
//          {
//              if(x[0]>XMAX[j])    switchX = damp_reflect(0,XMAX[j],x,v,vh); //required; actual deformation of the boundary is mostly in the x direction
//              if(x[1]>YMAX[j])    switchY = damp_reflect(1,YMAX[j],x,v,vh); //??; may not be too important? some orders of magnitude smaller than x deformation
//          }
//        }
//
//            if(x[1]<YMIN)       damp_reflect(1,YMIN,x,v,vh);  //same as original
    }

}

void leapfrog_step(sim_state_t* s, double dt
//                   ,mrbsm* R
                   , sim_param_t*param
                   )
/**preferred method due to simplicity, high rel. acc, conservative for particle methods**/
{
    const double* restrict a = s->a;
    double* restrict vh = s->vh;
    double* restrict v = s->v;
    double* restrict x = s->x;
    int n = s->n;

    for (int i=0;i<2*n;i++) vh[i] += a[i] * dt;
    for (int i=0;i<2*n;i++) v[i] = vh[i] + a[i] * dt / 2;
    for (int i=0;i<2*n;i++) x[i] += vh[i] * dt;

//    reflect_bc(s);
    reflect_bc(s,
//               R,
               param);
}

void leapfrog_start(sim_state_t* s, double dt
//                    ,mrbsm* R
                    , sim_param_t *param
                    )
/**start gets a different algorithm because both x and v are at t=0**/
{
    const double* restrict a = s->a;
    double* restrict vh = s->vh;
    double* restrict v = s->v;
    double* restrict x = s->x;
    int n = s->n;

    for (int i=0;i<2*n;i++) vh[i] = v[i] + a[i] * dt / 2;
    for (int i=0;i<2*n;i++) v[i] += a[i] * dt;
    for (int i=0;i<2*n;i++) x[i] += vh[i] * dt;

//    reflect_bc(s);
    reflect_bc(s,
//               R,
               param);
}

/**sample indicator functions; shape of fluid**/
int box_indicator(double x, double y)
{
    return (x < 0.5) && (y < 0.5);    //original
//    return (x < 1.0) && (y < 2.0);     //validation, Koshizuka and Oka (1996)
//    return (x < 1.0) && (y < 1.0);      //ce300 config
}

int circ_indicator(double x, double y)
{
    double dx = (x-0.5);
    double dy = (y-0.3);
    double r2 = dx*dx + dy*dy;

    return (r2 < 0.25*0.25);
}

sim_state_t* place_particles(sim_param_t* param, domain_fun_t indicatef)
/**fills region indicated by "indicatef" argument; assigns ONLY position of fluid;
currently: unit box**/
/**loop bounds dependent on container size**/
{
    double h = param->h;
    double hh = h/1.3; /**mesh smaller than particle radius; allow particles to overlap a bit**/

    double x, y;
    FILE*XYCheck;
    XYCheck = fopen("bindel_XY.csv", "w");

    // Count mesh points that fall in indicated region.
    int count = 0;
    for (x = 0; x < 1; x += hh){
        for (y = 0; y < 1; y += hh){
            fprintf(XYCheck,"%lf,%lf\n",x,y);
            count += indicatef(x,y);
            }}
//            printf("count = %d\n", count);  exit(3);
    // Populate the particle data structure
    sim_state_t* s = alloc_state(count);

    int p = 0;

    for (x=0;x<1;x+=hh)
    {
        for (y=0;y<1;y+=hh)
        {
            if (indicatef(x,y))
            {
                s->x[2*p+0] = x;
                s->x[2*p+1] = y;
                s->v[2*p+0] = 0;
                s->v[2*p+1] = 0;
                p++;
            }
        }
    }
return s;
}

void normalize_mass(sim_state_t* s, sim_param_t* param)
/**calculate for mass such that density is roughly consistent**/
{
    s->mass = 1;
    compute_density(s, param);
    double rho0 = param->rho0;
    double rho2s = 0;
    double rhos = 0;

    for (int i=0;i<s->n;i++)
    {
        rho2s += (s->rho[i])*(s->rho[i]);
        rhos += s->rho[i];
    }

    s->mass *= ( rho0*rhos / rho2s );
}

sim_state_t* init_particles(sim_param_t* param)/**initialize; place particles + normalize mass**/
{
    sim_state_t* s = place_particles(param, box_indicator);
    normalize_mass(s, param);
    return s;
}

void check_state(sim_state_t* s)
/**check if particles are within bounds; REPLACE TO REFLECT MRBSM DISPLACEMENTS**/
{
    for (int i=0;i<s->n;i++)
    {
        double xi = s->x[2*i];
        double yi = s->x[2*i+1];
        assert( xi >= 0 || xi <= 1 );
        assert( yi >= 0 || yi <= 1 );

        /**validation, Koshizuka & Oka 1996**/
//        assert( xi >= 0.0 || xi <= 4.0 );
//        assert( yi >= 0.0 || yi <= 3.0 );
    }
}

static void print_usage()
{
    sim_param_t param;
    default_params(&param);
    fprintf(stderr,
        "nbody\n"
        "\t-h: print this message\n"
        "\t-o: output file name (%s)\n" //a
        "\t-F: number of frames (%d)\n" //b
        "\t-f: steps per frame (%d)\n" //c
        "\t-t: time step (%e)\n" //d
        "\t-s: particle size (%e)\n" //e
        "\t-d: reference density (%g)\n" //f
        "\t-k: bulk modulus (%g)\n" //g
        "\t-v: dynamic viscosity (%g)\n" //h
        "\t-g: gravitational strength (%g)\n", //i
    param.fname, //a
    param.nframes, //b
    param.npframe, //c
    param.dt, //d
    param.h, //e
    param.rho0, //f
    param.k, //g
    param.mu, //h
    param.g); //i
}

int get_params(int argc, char** argv, sim_param_t* params)
{
    extern char* optarg;
    const char* optstring = "ho:F:f:t:s:d:k:v:g:";
    int c;

    #define get_int_arg(c, field) \
    case c: params->field = atoi(optarg); break
    #define get_flt_arg(c, field) \
    case c: params->field = (double) atof(optarg); break



    default_params(params);
    while ((c = getopt(argc, argv, optstring)) != -1)
    {
        switch (c)
        {
            case 'h':
            print_usage();
            return -1;

            case 'o':
            strcpy(params->fname = malloc(strlen(optarg)+1), optarg);

            break;
            get_int_arg('F', nframes);
            get_int_arg('f', npframe);
            get_flt_arg('t', dt);
            get_flt_arg('s', h);
            get_flt_arg('d', rho0);
            get_flt_arg('k', k);
            get_flt_arg('v', mu);
            get_flt_arg('g', g);

            default:
            fprintf(stderr, "Unknown option\n");
            return -1;
        }
    }
    return 0;
}

void write_header(FILE* fp, int n)
{
    fprintf(fp, "%s,%d,%d\n", VERSION_TAG, n, 1);
}

void write_frame_data(FILE* fp, int n, double* x, int* c)
{
    for (int i=0;i<n;i++)
    {
        double xi = *x++;
        double yi = *x++;
        int   ci = c ? *c++ : 0;


        fprintf(fp, "%e,%e,%d\n", xi, yi, ci);
    }
}

void split(double*R, //subject vector
           double*x, //xcomp
           double*y, //ycomp
           int n    //number of elements of split
           )
{
    int i;
    for(i=0;i<n;i++)
    {
        x[i] = R[2*i];
        y[i] = R[2*i+1];
    }

}

void sortascend (double*R,  //subject vector
                 double*r,  //result vector
                 int n     //number of elements
                )
{
    int i,j;
    double hold;

    /**copy elements into new vector**/
    for(i=0;i<n;i++)
    {
        r[i]=R[i];
    }

    /**sort new vector**/
    for(i=0;i<n;i++)
    {
        for(j=i+1;j<n;j++)
        {
            if(r[i]>r[j])
            {
                hold=r[i];
                r[i]=r[j];
                r[j]=hold;
            }
        }
    }
}

void flip(double*x,   //subject vector
          double*y, //storage vector
          int n   //number of elements
          )
/**prelim soln [1] to flip y-axis; store x0->xN to a different vector y; xN->y0 to x0->yN;
>keeps original vector as is, stores flipped version in a different vector;
>also separates horizontal and vertical components**/
{
    int i;

    double*u, *v;
    u = (double*)calloc(n, sizeof(double));
    v = (double*)calloc(n, sizeof(double));

    /**split x into u and v**/
    for(i=0;i<n;i++)
    {
        u[i]    =   x[i*2];
        v[i]    =   x[(i*2)+1];
    }

    /**allocate u and v to y, in reverse order**/
    for(i=0;i<n;i++)
    {
        y[(2*i)]      =   u[(n-i)-1];
        y[(2*i)+1]    =   v[(n-i)-1];
    }

    free(u);
    free(v);
}

double** main_sph(int argc, char** argv, sim_param_t params, sim_state_t* state, double t)
{

    int i, count=0;


    FILE* fp = fopen(params.fname, "a");

    int nframes = params.nframes;
    int npframe = params.npframe;
    double dt = params.dt;
    int n = state->n;


//    write_header(fp, state->n);
//    write_frame_data(fp, state->n, state->x, NULL);
//    compute_accel(state, &params);
//    leapfrog_start(state, dt, &params);
//    check_state(state);



//    for (int frame = 1; frame < nframes; frame++)
//    {
        for (int i = 0; i < npframe; i++)
        {
            compute_accel(state, &params);
            leapfrog_step(state, dt, &params);
            check_state(state);
        }
//        write_frame_data(fp, n, state->x, NULL);
//        count++;
//        printf("ti = %f\n", count*dt);
//    }
    double** ANSWER=calloc(sizeof(double*),2);

    ANSWER[0]=calloc(sizeof(double),state->n+1);
    ANSWER[1]=calloc(sizeof(double),state->n+1);


    split(state->x, ANSWER[0],ANSWER[1],state->n);

    return ANSWER;

}
