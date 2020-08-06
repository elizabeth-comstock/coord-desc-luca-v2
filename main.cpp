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
void lsrch( int dim, double *x0, double *l, double &y, double &f )
{
    double f1,f0,ft;
    double delta=1.e-4;
    double x1[N];
    double xt[N];
    double zeros[N] = {0,0};

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
    // int stepdir = 1;    // either 1 or -1, multiply to e to set direction
    // actually determine search direction
    while(dirfound == 0 && sdn <= 20)
    {
        // take a small step from x0 and write coordinates to x1
        fmav( dim, x0, delta,l, x1 );

        // calculate function value for each of the points
        f0 = fun( x0 );
        f1 = fun( x1 );

        // DEBUG
        printf("***ORTHOGONAL VECTOR***\n");
        printf("v  = %f, %f\t\n\n",l[0],l[1]);
        printf("***DETERMINE SEARCH DIRECTION***\n");
        printf("c0 = %f, %f \tF0 = %f \n", x0[0],x0[1],f0);
        printf("c1 = %f, %f \tF1 = %f \n", x1[0],x1[1],f1);

        if (f1 < f0){
            // We are moving in a direction of descent! Terminate the loop.
            printf("Dir trial %i success! \t\td = %f \n\n",sdn,delta);
            dirfound = 1;
        } else {
            // reverse direction and halve, and check other side
            printf("Dir trial %i fail! \td = %f \n\n",sdn,delta);
            delta *= -0.1;      // decreased from -0.5 to make more performant
            // stepdir *= -1;
        }

        sdn++;
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
    // double e = 0.5 * stepdir;
    double e = 5000 * delta;
    double esuccess = e;        // holder variable for successful step length
    double a = 3;
    double b = 0.5;             // values recommended by paper
    int lsn = 1;                // hard iteration limit so no going out of control
    int succount = 0;
    int failcount = 0;
    f1 = f0;                    // now use f1 as a holder variable for the most successful trial
    // calculate actual step length
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

//fmav( dim, zeros, esuccess, l, &y );    // write step length
    f = f1;         // write new function value
    y = esuccess;
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

    // DEBUG
    printf( "Starting coordinates: c = %f, %f \n", x[0],x[1] );
    printf( "Starting value:       f = %f \n\n", fun(x) );

    for( k=0;k<30;k++ )     // outer loop, through each SET of orthogonal vectors
    {
        // DEBUG
        printf( "***STEP %i*** \n\n", k );

        int    hmin=-1;
        double fmin=1e16;

        // Section 3.1 of Rosenbrock's paper with modified line search
        // printf( "Search directions\n" );

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
            lsrch( dim, x,v[h], y[h], f[h] );
            printf("***LINE SEARCH COMPLETE! RESULTS:***\n");
            printf( "Step = %f, %f \tF = %f \n\n", y[0],y[1],f[h] );
            // printf( "F  = %f \tStep  = %f, %f \n\n", f[h],y[0],y[1] );

            // Here we a checking which direction gave us the smallest
            // value. We will use this for two operations
            // 1) update starting point (see line 194
            // 2) generate new set of orthogonal directions
            // hmin is the direction that yields the minimum, fmin
            // is that minimum
             if( fun(x) > f[h] ){ fmav( dim, x,y[h],v[h],   x); }

        }

// Generate new set of orthogonal directions, see Equation 8 in the paper.
        for( h=0;h<dim;h++ )
        {
            memset( a[h],0,N*sizeof(double) );
            for( j=h;j<dim;j++ ){ fmav( dim, a[h],y[j],v[j], a[h] ); };
        }

// Update starting point
        //fmav( dim, x,1.,a[0], x );

/* Equation 9 */
        for( h=0;h<dim;h++ )
        {
            memcpy( v[h],a[h],N*sizeof(double) );
            for( j=0;j<h;j++ ){ gsot( dim, v[h], v[j] ); }
            norml( dim, v[h] );
        }
    }

    exit(0);
}
