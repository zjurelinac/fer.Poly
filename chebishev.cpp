/***************************************************************************************************

    File "chebishev.cpp"
    --------------------------------------------------------------------------------------------
    Function approximation using Chebyshev polynomials of the first kind.

    Details:
        Approximation of function cos x on interval [-1, 1] with Chebyshev polynomials upto
        degree 8. Algorithm hereby implemented can be found at:
            http://mathworld.wolfram.com/ChebyshevApproximationFormula.html

    Author:
        zvonimir.jurelinac@fer.hr

***************************************************************************************************/
#define _USE_MATH_DEFINES   // neccessary to include the definition of PI from cmath

#include <cmath>
#include <cstdio>

// Uses the Poly class from the previous problem
#include "poly.h"

int main(){
    Poly cheb[ 10 ];
    double c[ 9 ];
    double p1[] = { 1 }, p2[] = { 0, 1 }, x, y;

    cheb[ 0 ] = Poly( 0, p1 );
    cheb[ 1 ] = Poly( 1, p2 );

    // Generating Chebyshev polynomials from their recursive relation
    for( int i = 2; i <= 8; ++i )
        cheb[ i ] = 2*cheb[ i-1 ].shift( 1 ) - cheb[ i-2 ];

    // Approximating function cos( x ) via Chebyshev polynomials
    for( int i = 0; i <= 8; ++i ){
        /*  c_i = 2/N * sum( from k=1 to N ) of f( x_k )*T_i( x_k ), where x_k is the k-th zero of
            N-th Chebyshev polynomial, ie. when N = 9( ie. number of Chebyshev polynomials used ),
            f = cos x, c_i = 2/9 * sum of cos( cos ( x_k ) ) * cos( i*x_k ) */

        c[ i ] = 0.0;
        for( int k = 1; k <= 9; ++k ){
            double xk = M_PI * ( k - 0.5 ) / 9.0;
            c[ i ] += cos( cos( xk ) ) * cos( i * xk );
        }
        c[ i ] /= 4.5;
    }

    for( int i = 0; i <= 8; ++i )
        printf( "%.5lf  ", c[ i ] );

    printf( "\nEnter the value of x from [-1, 1 ] to evaluate cos at: " );
    scanf( "%lf", &x );

    y = - c[ 0 ]/2;
    for( int i = 0; i <= 8; ++i )
        y += c[ i ] * cheb[ i ]( x );

    printf( "%.10lf\n", y );

    printf( "Difference from the standard function: %.12lf\n", y - cos( x ) );

    return 0;
}
