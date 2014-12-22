/***************************************************************************************************

    File "poly.cpp"
    --------------------------------------------------------------------------------------------
    Contains implementation details

    class Polynomial:

        Implements all the basic functionalities of a mathematical polynomial,
        including overloaded operators for a more natural usage on symbolic level.

    Author:
        zvonimir.jurelinac@fer.hr

***************************************************************************************************/

#include <cstdio>
#include <cmath>

#include <algorithm>
#include <vector>
#include <string>
#include <complex>

#include "poly.h"

using namespace std;



/* Poly constructors */

Poly::Poly(){
    N = -1;
    a = NULL;
}

Poly::Poly( int n ) : N( n ) {
    a = new double[ N+1 ];
    for( int i = 0; i <= N; ++i )
        a[ i ] = 0;
    name = "_";
}

Poly::Poly( int n, double *cs ) : N( n ) {
    a = new double[ N+1 ];
    for( int i = 0; i <= N; ++i )
        a[ i ] = cs[ i ];
    name = "_";
}

Poly::Poly( const Poly &P ){
    N = P.N;
    a = new double[ N+1 ];
    for( int i = 0; i <= N; ++i )
        a[ i ] = P[ i ];
    name = P.name;
}


/* Poly simple and unary operators*/

double& Poly::operator []( const int i ){
    return a[ i ];
}

double const& Poly::operator []( const int i ) const{
    return a[ i ];
}

Poly Poly::operator -(){
    Poly P = Poly( *this );
    for( int i = 0; i <= P.N; ++i )
        P[ i ] = - P[ i ];
    return P;
}

// Horner's method for evaluation of polynomials
double Poly::operator()( double x ){
    double r = 0;
    for( int i = N; i >= 0; --i )
        r = r*x + a[ i ];
    return r;
}

bool Poly::operator ==( double d ){
    for( int i = N; i > 0; --i )
        if( fabs( a[ i ] ) > EPS ) return false;
    return a[ 0 ] == d;
}

/* Polynomial arithmetics */

Poly operator +( const Poly P, const Poly Q ){
    int M = max( P.N, Q.N );
    Poly nP = Poly( M );

    for( int i = 0; i <= P.N; ++i )
        nP.a[ i ] += P[ i ];

    for( int i = 0; i <= Q.N; ++i )
        nP.a[ i ] += Q[ i ];


    return nP;
}

Poly operator -( const Poly P, const Poly Q ){
    int M = max( P.N, Q.N );
    Poly nP = Poly( M );

    for( int i = 0; i <= P.N; ++i )
        nP.a[ i ] += P[ i ];

    for( int i = 0; i <= Q.N; ++i )
        nP.a[ i ] -= Q[ i ];

    nP.reduce();

    return nP;
}

Poly operator *( const Poly P, const Poly Q ){
    Poly nP = Poly( P.N + Q.N );

    for( int i = 0; i <= P.N; ++i )
        for( int j = 0; j <= Q.N; ++j )
            nP[ i + j ] += P[ i ] * Q[ j ];

    nP.reduce();

    return nP;
}

Poly operator *( const double q, const Poly P ){
    Poly nP = Poly( P );
    for( int i = 0; i <= P.N; ++i )
        nP[ i ] *= q;

    return nP;
}

// Symbolic division of two polynomials, using long division algorithm
Poly operator /( const Poly P, const Poly Q ){
    if( P.N < Q.N ) return Poly( 0 );

    int M = P.N - Q.N;
    Poly Pp = Poly( P ), Rr = Poly( M );

    for( int i = 0; i <= Pp.N; ++i )
        Pp[ i ] /= Q[ Q.N ];
    Poly Qq = Q.normalized();

    while( Pp.N >= Qq.N ){
        Poly Tt = Poly( 0 );

        Tt[ 0 ] = Pp[ Pp.N ];
        Tt = Tt.shift( Pp.N - Qq.N );

        Rr = Rr + Tt;
        Pp = Pp - Tt*Qq;
    }

    Rr.reduce();

    return Rr;
}

// Finding the reminder of division
Poly operator %( const Poly P, const Poly Q ){
    return P - ( P / Q ) * Q;
}

/* Poly miscelaneous functions */

void Poly::print() const{
    printf( "%s( x ) := ", name.c_str() );
    for( int i = N; i >= 2; --i )
        if( a[ i ] != 0 )
            printf( "%+.2lf * x^%d ", a[ i ], i );

    if( fabs( a[ 1 ] ) > EPS )
        printf( "%+.2lf * x ", a[ 1 ] );

    if( fabs( a[ 0 ] ) > EPS )
        printf( "%+.2lf", a[ 0 ] );
    puts( "" );
}

void Poly::dump() const{
    printf( "poly <%s>; degree( %d ); coeffs = [ ", name.c_str(), N );
    for( int i = N; i >= 0; --i )
        printf( "%.16lf, ", a[ i ] );
    printf( "]\n" );
}

void Poly::reduce(){
    int i;
    for( i = N; i >= 0 && fabs( a[ i ] ) < EPS; --i );
    N = i;
}

Poly Poly::normalized() const{
    if( N < 0 ) return *this;
    Poly nP = Poly( N );
    for( int i = 0; i < N; ++i )
        nP[ i ] = a[ i ] / a[ N ];
    nP[ N ] = 1;

    return nP;
}

Poly Poly::shift( int k ){
    Poly nP = Poly( N + k );
    for( int i = 0; i <= N; ++i )
        nP[ i + k ] = a[ i ];

    return nP;
}

bool Poly::isNumber(){
    for( int i = N; i > 0; --i )
        if( fabs( a[ i ] ) > EPS ) return false;
    return true;
}

/* Polynomial mathematical functions */

Poly Poly::differentiate(){
    Poly nP = Poly( N-1 );
    for( int i = 1; i <= N; ++i )
        nP[ i-1 ] = i * a[ i ];

    return nP;
}


Poly gcd( Poly a, Poly b ){
    if( b.N < 0 || b == 0 )
        return a.normalized();

    return gcd( b, a % b );
}


/* Main part */
vector< double > Poly::roots(){
    Poly P = Poly( *this );
    Poly dP = P.differentiate();

    vector< double > troots, roots;

    // Eliminate multiple roots

    P = ( P / gcd( P, dP ) ).normalized();

    dP = P.differentiate();

    double bound = 0;
    for( int i = 0; i < N; ++i )
        bound = max( bound, fabs( P[ i ] ) );

    bound += 1;

    double delta = 2*bound/Intervals, x, tx;
    for( x = -bound; x < bound; x += delta ){
        // if there is a sign change in this interval
        if( P( x )*P( x + delta ) < 0 ){
            tx = x;
            // look for a root
            while( fabs( P( tx ) ) > EPS )
                tx -= P( tx ) / dP( tx );


            troots.push_back( tx );
            Poly tp = Poly( 1 );
            tp[ 0 ] = -tx;
            tp[ 1 ] = 1;
            P = P / tp;
        }
    }

    P = Poly( *this );
    P.name = "P";

    int s = troots.size();
    for( int i = 0; i < s; ++i ){
        Poly tp = Poly( 1 );
        tp[ 0 ] = -troots[ i ];
        tp[ 1 ] = 1;

        do{
            P = P / tp;
            roots.push_back( troots[ i ] );
        } while( fabs( P( troots[ i ] ) ) < EPS2 );
    }

    return roots;
}


vector< Poly > Poly::factorize(){
    Poly P = Poly( *this ).normalized(), Q, T, r1, r2;
    vector< Poly > V;

    double as[ 3 ], u, v, c, d, g, h, D, err, nu, nv;
    int iter;

    while( P.N > 2 ){
        u = 1;
        v = 1;
        as[ 2 ] = 1;

        iter = 0;

        do{
            //printf( "%lf %lf\n", u, v );
            as[ 0 ] = v;
            as[ 1 ] = u;
            as[ 2 ] = 1;

            T = Poly( 2, as );

            Q = P / T;
            r1 = P % T;

            r2 = Q % T;

            c = r1[ 1 ];
            d = r1[ 0 ];

            g = r2[ 1 ];
            h = r2[ 0 ];

            D = 1.0 / ( v*g*g + h*( h - u*g ) );

            nu = u - ( g*d - c*h ) * D;
            nv = v - ( d*( g*u - h ) - c*g*v ) * D;

            err = max( fabs( u - nu ), fabs( v - nv ) );

            u = nu;
            v = nv;

            ++iter;

       } while( err > EPS3 && iter < MAX_ITER );

       if( iter == MAX_ITER ) // Cannot converge to the solution from the initial condition
            break;

        as[ 0 ] = v;
        as[ 1 ] = u;
        as[ 2 ] = 1;

        T = Poly( 2, as );

        V.push_back( T );

        P = P / T;
    }
    V.push_back( P );
    return V;
}


vector< Complex > Poly::complexRoots(){
    vector< Complex > roots;
    vector< Poly > qs = this->factorize();
    Complex b, c, x1, x2, D;
    Poly Q;

    for( int i = 0; i < qs.size(); ++i ){
        Q = qs[ i ];
        if( Q.N > 2 ) continue;
        b = Q[ 1 ] / 2.0;
        c = Q[ 0 ];
        D = b*b - c;
        x1 = -b + sqrt( D );
        x2 = -b - sqrt( D );
        roots.push_back( x1 );
        roots.push_back( x2 );
    }

    return roots;
}
