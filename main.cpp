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
    double f1,f0,ft;
    double delta=1.e-6;
    double x1[N], xt[N];
    double zeros[N] = {0,0};
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
    f0 = fun( x0 );     // original function value
    lsn += 1;

    // DEBUG
    printf("***ORTHOGONAL VECTOR***\n");
    printf("v  = %f, %f\t\n\n",l[0],l[1]);
    printf("***DETERMINE SEARCH DIRECTION***\n");
    printf("c0 = %f, %f \tF0 = %f \n\n", x0[0],x0[1],f0);

    // actually determine search direction
    while(dirfound == 0 && sdn <= 2)
    {
        // take a small step from x0 and write coordinates to x1
        fmav( dim, x0, delta,l, x1 );
        f1 = fun( x1 ); // new function value
        lsn += 1;
        // DEBUG
        printf("c1 = %f, %f \tF1 = %f \t\t", x1[0],x1[1],f1);

        if (f1 < f0){
            // We are moving in a direction of descent! Terminate the loop.
            printf("Dir trial %i success! \td = %f \n\n",sdn,delta);
            dirfound = 1;
        } else {
            // reverse direction and halve, and check other side
            printf("Dir trial %i fail! \td = %f \n",sdn,delta);
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
    /* Try a step of arbitrary length e, pointing the right way.

       UPDATE: simple modification to prevent overshoot. Instead of:
           double e = 0.5 * stepdir; (Version A)
       We have:
           double e = 5000 * delta;  (Version B)
       So that step size scales with results from previous direction search.
       Furthermore eliminates need to use holder variable stepdir.
       Version iteration:  A    B    C
       Step 4  Dir 1      13   13   15
       Step 26 Dir 1      37   24   16
       Step 29 Dir 0      36   21   16
       Feel free to implement gradient-based method later on.
     */
    double e = 5000 * delta;
    double esuccess = e;        // holder variable for successful step length
    double a = 3;
    double b = 0.5;             // values recommended by paper
    int succount = 0;
    int failcount = 0;
    f1 = f0;                    // now use f1 as a holder variable for the most successful trial
    // calculate actual step length
    while((succount < 1 || failcount < 3) && lsn <= 500)
    {
        // do a trial
        fmav( dim, x0, e, l, xt );
        ft = fun(xt);
        lsn += 1;
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
    }

    f = f1;         // write new function value
    y = esuccess;
    results = lsn;
    printf("\n");   // OCD

    return;
}

int main()
{
    int j,h,k;  // iteration variables
    int dim=N;  // number of dimensions of the problem

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
    int results[30][4];
    int funceval = 0;

    // DEBUG
    printf( "Starting coordinates: c = %f, %f \n", x[0],x[1] );
    printf( "Starting value:       f = %f \n\n", fun(x) );

    FILE *fi;
    fi= fopen( "points.dat","w" );

    for( k=0;k<30;k++ )     // outer loop, through each SET of orthogonal vectors
    {
        // DEBUG
        printf( "***STEP %i*** \n\n", k );

        fprintf( fi,"% 12.5e %12.5e \n", x[0],x[1] );

        int    hmin=-1;
        double fmin=1e16;
        int    ismin[N]={0,0};   // adds minima detection functionality

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
        funceval += results[k][3];

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
    printf("***NUMBER OF FUNCTION EVALUATIONS FOR EACH ITERATION AND DIRECTION***\n");
    printf("Step 1\tDir 0\tDir 1\tTotal\n");
    for( h=0;h<30;h++ )
    {
        for( j=0;j<4;j++ ){printf("%i\t", results[h][j]);}
        printf(" \n");
    }

    printf(" \n");
    printf("%i function evaluations\n\n", funceval);

    fclose( fi );

    exit(0);
}
