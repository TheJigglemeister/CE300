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
#include ".\libraries - sph\generalcodes_SPH.h"

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

sim_state_t* alloc_state(int n);

void free_state(sim_state_t* s);
typedef int (*domain_fun_t)(double, double); /**1 = occupied by fluid, 0 elsewhere**/

static void default_params(sim_param_t* params)
/**SPH particle parameters**/
/**validation params are sufficient for actual simulation**/
{
}


void compute_density(sim_state_t* s, sim_param_t* params);

void compute_accel(sim_state_t* state, sim_param_t* params);
static int damp_reflect(int which, double barrier, double* x, double* v, double* vh);
static void reflect_bc(sim_state_t* s
//                       ,mrbsm*R
                       ,sim_param_t*params)
;

void leapfrog_step(sim_state_t* s, double dt
//                   ,mrbsm* R
                   , sim_param_t*param
                   );
void leapfrog_start(sim_state_t* s, double dt
//                    ,mrbsm* R
                    , sim_param_t *param
                    );
/**sample indicator functions; shape of fluid**/
int box_indicator(double x, double y);
int circ_indicator(double x, double y);
sim_state_t* place_particles(sim_param_t* param, domain_fun_t indicatef);
void normalize_mass(sim_state_t* s, sim_param_t* param);
sim_state_t* init_particles(sim_param_t* param);
void check_state(sim_state_t* s);
static void print_usage();
int get_params(int argc, char** argv, sim_param_t* params);
void write_header(FILE* fp, int n);
void write_frame_data(FILE* fp, int n, double* x, int* c);
void split(double*R, //subject vector
           double*x, //xcomp
           double*y, //ycomp
           int n    //number of elements of split
           );
void sortascend (double*R,  //subject vector
                 double*r,  //result vector
                 int n     //number of elements
                );
void flip(double*x,   //subject vector
          double*y, //storage vector
          int n   //number of elements
          );

double** main_sph(int argc, char** argv
                  , sim_param_t params, sim_state_t* state, double t
                  );
