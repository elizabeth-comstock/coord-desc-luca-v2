# include <cstdlib>
# include <cstdio>
# include <cmath>
# include <cstring>

# include "vec.h"

/* Optimization based on H.H. Rosenbrock, "An Authomatic Method for Finding the
   Greatest or Least Value of a Function", The Computer Journal, Vol 3(3), 1960,
   pp 175-184

   Written by Luca di Mare 25 Jul 2020 17:02:44
 */

/* Update 6 August 2020: adding valley detection.
 */

// Rosenbrock function
double fun( double *x )
{
    double u,v;
    u= x[1]-x[0]*x[0];
    v= 1.- x[0];
    return 100*u*u+ v*v;
}

/** Line search
    x0                  Initial position
    l                   Search direction
    y                   Step to local minimum
    f                   Estimated local minimum
 **/
void lsrch( int dim, double *x0, double *l, double &y, double &f, int &ismin, int &results )
{
    double fP,fd,fB,fA,f0,ft;
    double gradA, gradB, gradt;
    double delta=1.e-4;
    double xt[N];
    double xd[N];
    double xA[N];
    double xB[N];
    double xP[N];

    int lsn = 0;

    /* A line search contains the following steps:
       1. Determine search direction, positive or negative. This is performed
          by taking a small step to the left or right of the current position
          and evaluation the function value.
       2. Perform the actual line search, calculating the step length.
     */

    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       1. For each vector, determine search direction
          fmav function linearly adds vectors
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    // values declared here because they reset with each iteration
    bool dirfound = 0;  // direction finding success flag
    int sdn = 1;        // hard iteration limit so no going out of control
    f0 = fun(x0);     // original function value

    // DEBUG
    // lsn += 1;
    printf("***ORTHOGONAL VECTOR***\n");
    printf("v  = %f, %f\t\n\n",l[0],l[1]);
    printf("***DETERMINE SEARCH DIRECTION***\n");
    printf("c0 = %f, %f \tF0 = %f \n", x0[0],x0[1],f0);

    // actually determine search direction
    while(dirfound == 0 && sdn <= 4)
    {
        // take a small step from x0 and write coordinates to xt
        fmav(dim, x0, delta, l, xt);
        ft = fun(xt); // new function value

        // DEBUG
        // lsn += 1;
        printf("ct = %f, %f \tFt = %f \t\t", xt[0],xt[1],ft);

        if (ft < f0){
            // calc gradient at original pt, taking descent direction as positive
            gradA = (ft-f0)/delta;
            // We are moving in a direction of descent! Terminate the loop.
            dirfound = 1;
            printf("grad0 = %f \n", gradA);
            printf("Dir trial %i success! \td = %f \n\n", sdn,delta);
        } else {
            // reverse direction and halve, and check other side
            printf("Dir trial %i fail! \td = %f \n", sdn,delta);
            delta *= -0.1;      // decreased from -0.5 to make more performant
        }

        sdn++;
    }

    /* Minima detection: if 4 iterations of direction search have
       failed at this point, we're probably at a minima along the axis.
     */
    if(dirfound == 0){
        printf("Minima found along test direction! \n\n");
        ismin = 1;
        results = 0;
        f = f0;
        y = (1.e-6) * delta;
        return;
    }

    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       2. Perform line search, calculating step length
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    printf("***PERFORM LINE SEARCH***\n");
    // First, attempt regula falsi method. If that fails, revert back to previous method.

    // STEP LENGTHS
    double eB = 5000 * delta;
    double esuccess;
    double e0 = 0;
    double eA = e0;
    double eP, eQ;

    // ROSENBROCK SLEDGEHAMMER
    double a = 3;
    double b = 0.5;     // values recommended by paper
    int succount = 0;
    int failcount = 0;

        /* %%%%%%%%%%%%%%%%%%%%%%%%%%%
            REGULA FALSI, FIVE TRIALS
           %%%%%%%%%%%%%%%%%%%%%%%%%%% */
    fmav( dim, x0, 0, l, xA );  // duplicate x0 to xA and keep x0 as original
    fmav( dim, x0, eB, l, xB ); // set a point x1 some distance away
    fA = f0;
    fB = fun(xB);
    lsn += 1;

    for (int i=1; i<=2; i++)
    {
        // DEBUG
        printf("A = %f, %f \tFA = %f\n",xA[0],xA[1],fA);
        printf("B = %f, %f \tFB = %f   \teB = %f   \n",xB[0],xB[1],fB,eB);
        // interpol( dim, x0, x1, -f0, f1, xt );   // find trial estimated minima point
        // treat the descent direction as an axis, calculate step size, then calculate coordinates and function value
        eP = ipstep( eA, eB, -fA, fB );
        fmav( dim, x0, eP, l, xP );
        fP = fun(xP);
        // DEBUG
        printf("P = %f, %f \tFP = %f    \teP = %f\t",xP[0],xP[1],fP,eP);

        // evaluate gradient at x1
        fmav( dim, xB, delta, l, xd );  // offset x1 by delta
        fd = fun(xd);
        gradB = (fd-fB)/delta;

        /* Equation 15 can't be implemented on vectors!
           Instead we take A and B to be the step length along the searchdir.
         */
        eQ = (fB-fA+eA*gradA-eB*gradB)/(gradA-gradB);

        lsn += 2;
        printf("eQ = %f\tlsn = %i\n\n",eQ,lsn);

        // now narrow the bracket
        if (eP < eQ) {
            eA = eP;
            eB = eQ;
        } else {
            eA = eQ;
            eB = eP;
        }

        /* and update coordinates
           Q (x1) MUST BE UPDATED BEFORE P (x0) BECAUSE x0 IS OVERWRITTEN
         */
        fmav( dim, x0, eA, l, xA );
        fmav( dim, x0, eB, l, xB );

        fA = fun(xA);
        fB = fun(xB);

        esuccess = eP;
        lsn += 2;

        /*
        // if regula falsi is successful, skip the loop!
        if( ft < f1 && ft < f0 ){
            esuccess = e * (f0 / (f0+f1));  // recalculate successful step size taken for output

            printf("Regula falsi successful! \n");
        } else {
            printf("Regula falsi failed! Reverting to previous method...\n\n");
            f1 = f0;        // use f1 as holder variable for most successful trial
        }
        */
    }
    /*
    // if regula falsi fails, fall back to Rosenbrock sledgehammer
    while((succount < 1 || failcount < 5) && lsn <= 500)
    {
        // do a trial
        fmav( dim, x0, e, l, xt );
        ft = fun(xt);
        // DEBUG
        printf("ct = %f, %f \tFt = %f   \te = %f   \t",xt[0],xt[1],ft,e);

        if (ft <= f1) {
            f1 = ft;            // update function value to lowest found
            esuccess = e;
            e = e*a;
            succount += 1;
            failcount = 0;      // we want one success and five fails in succession
        } else {
            e = e*b;
            failcount += 1;
        }
        printf("%i %i %i \n",succount,failcount,lsn);
        lsn++;
    }
    */
    f = fP;         // write new function value
    y = esuccess;
    results = lsn;
    printf("\n");   // OCD

    return;
}

int main()
{
    int j,h,k;          // iteration variables
    int dim=N;          // number of dimensions of the problem
    int outloops=30;    // set number of outer loops here

    double v[N][N]={ {1,0},{0,1} }; // orthogonal vectors, {0,1) and (1,0) for Rosenbrock
    double y[N]={0,0};
    double f[N]={1,1};
    double a[N][N]={0,0};

    double x[N]= {-1.2,1};          // starting point, (-1.2,1) for Rosenbrock

    /* DEBUG: results matrix to hold number of function evaluations
       First column holds the iteration number.
       Second  and third column hold the number of evaluations for line search, each direction.
       Fourth column holds the total number of evaluations for the iteration.
     */
    int results[outloops][4];

    // DEBUG
    printf( "Starting coordinates: c = %f, %f \n", x[0],x[1] );
    printf( "Starting value:       f = %f \n\n", fun(x) );

    for( k=0;k<outloops;k++ )     // outer loop, through each SET of orthogonal vectors
    {
        // DEBUG
        printf( "***STEP %i*** \n\n", k );

        int ismin[N]={0,0};   // adds minima detection functionality

        // Section 3.1 of Rosenbrock's paper with modified line search
        for(h=0; h<2; h++)
        {
            printf( "Direction %i \n\n" ,h );
            /* Feed into the line search function:
               dim  number of dimensions
               x    starting position
               v[h] the given orthogonal vector
               y[h] step to local minimum - this is a scalar variable
               f[h] estimated local minimum - this is a scalar variable
             */
            lsrch( dim, x, v[h], y[h], f[h], ismin[h], results[k][h+1] );
            printf("***LINE SEARCH COMPLETE! RESULTS:***\n");
            printf( "Step = %f, %f \tF = %f \n\n", y[0],y[1],f[h] );

            /* Here we are checking which direction gave us the smallest value.
               We will use this for two operations:
               1. update starting point
               2. generate new set of orthogonal directions
               hmin is the direction that yields the minimum, fmin is that minimum.
             */
            if( fun(x) > f[h] ){ fmav( dim, x,y[h],v[h],   x); }
        }

        // DEBUG
        results[k][0] = k + 1;
        results[k][3] = results[k][1] + results[k][2];

        // Generate new set of orthogonal directions, see Equation 8 in the paper.
        for( h=0;h<dim;h++ )
        {
            memset( a[h],0,N*sizeof(double) );
            for( j=h;j<dim;j++ ){ fmav( dim, a[h],y[j],v[j], a[h] ); };
        }

        // Equation 9
        for( h=0;h<dim;h++ )
        {
            memcpy( v[h],a[h],N*sizeof(double) );
            for( j=0;j<h;j++ ){ gsot( dim, v[h], v[j] ); }
            norml( dim, v[h] );
        }
    }

    // DEBUG: print results array
    printf("***NUMBER OF FUNCTION EVALUATIONS FOR LINE SEARCH IN EACH ITERATION AND DIRECTION***\n");
    printf("Step 1\tDir 0\tDir 1\tTotal\n");
    for( h=0;h<outloops;h++ )
    {
        for( j=0;j<4;j++ ){printf("%i\t", results[h][j]);}
        printf(" \n");
    }

    exit(0);
}
