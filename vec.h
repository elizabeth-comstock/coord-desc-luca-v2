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
