#ifndef VEC_H_INCLUDED
#define VEC_H_INCLUDED

#define N 2

/* Linear combination of vectors in dim dimensions
   For each element in vectors x and y, computes x + a*y
   and writes result to z.
 */
inline void fmav( int dim, double *x, double a, double *y, double *z )
{
    for( int i=0;i<dim;i++ )
    {
        z[i]= x[i]+ a*y[i];
    }
}

/* Linear interpolation function for regula falsi.
   There are two versions:
   interpol() provides the actual coordinates
   ipstep() provides the required step length
 */
inline void interpol( int dim, double *x1, double *x2, double b1, double b2, double *x )
{
    for( int i=0;i<dim;i++ )
    {
        x[i]= ( b1*x2[i] - b2*x1[i] )/( b1 - b2 );
    }
}

inline double ipstep( double e1, double e2, double f1, double f2 )
{
    return ( f1*e2 - f2*e1 )/( f1 - f2 );
}

/* Compares function value at four reference points for line search */
inline char least( double A, double B, double P )
{
    if ( A<B && A<P ){
        return 'A';
    } else if ( B<A && B<P ){
        return 'B';
    } else {
        return 'P';
    }
}

/* Normalize vector to unit norm */
inline void norml( int dim, double *x )
{
    double u=0;
    for( int i=0;i<dim;i++ )
    {
        u+= x[i]*x[i];
    }
    u= sqrt(u);
    u= 1./u;
    for( int i=0;i<dim;i++ )
    {
        x[i]*= u;
    }
    return;
}

/* Step in Graham Schmidt orthogonalization */
inline void gsot( int dim, double *x, double *z )
{
    double u=0;
    for( int i=0;i<dim;i++ )
    {
        u+= x[i]*z[i];
    }
    for( int i=0;i<dim;i++ )
    {
        x[i]-= u*z[i];
    }
    return;
}

#endif // VEC_H_INCLUDED
