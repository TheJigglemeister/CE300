#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "generalcodes_SPH.h"
#include "ce300specific_SPH.h"

#define SQR(x) (x*x)
#define pi() acos(-1)

/**matrix RK4; algorithm to be used in the study

 * Classical Runge-Kutta method (RK4) for systems of ODEs
 *      y<k+1> = y<k> + (h/6)*(k1 + 2k2 + 2k3 + k4)
 * that returns: Matrix with the following columns
 *   | t  |  y1  | ... | yn | y1dot | ...  | yndot |
 **/