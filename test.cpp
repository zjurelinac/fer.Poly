/***************************************************************************************************

    File "test.cpp"
    --------------------------------------------------------------------------------------------
    Contains tests of polynomial functionality,
    particulary of finding roots by Newton's method.

    Test case 1:
        (x-1)^2*(x+2)^3*(x^2+x+2)*(x-5)
            = x^8 - 18x^6 - 36x^5 - 7x^4 + 44x^3 + 80x^2 + 16x - 80

    Test case 2:
        x^10 - x^3 + 1

    Author:
        zvonimir.jurelinac@fer.hr

***************************************************************************************************/

#include <cstdio>

#include "poly.h"

int main(){
    double cs1[] = { -80, 16, 80, 44, -7, -36, -18, 0, 1 },
           cs2[] = { 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1 };

    Poly F( 8, cs1 );

    vector< double > rs1 = F.roots();

    puts( "Testing polynomial F(x) = (x-1)^2*(x+2)^3*(x-5):" );
    for( int i = 0; i < rs1.size(); ++i )
        printf( "Root %d = %.8lf\n", i+1, rs1[ i ] );

    if( rs1.size() == 0 )
        puts( "The polynomial F has no real roots." );



    Poly G( 10, cs2 );

    vector< double > rs2 = G.roots();

    puts( "Testing polynomial G(x) = x^10 - x^3 + 1:" );
    for( int i = 0; i < rs2.size(); ++i )
        printf( "Root %d = %.8lf\n", i+1, rs2[ i ] );

    if( rs2.size() == 0 )
        puts( "The polynomial G has no real roots." );

    puts( "Quadratic factors of polynomial F are:" );
    vector< Poly > factors = F.factorize();
    for( int i = 0; i < factors.size(); ++i )
        factors[ i ].print();

    // Finds all the quadratic factors of
    puts( "Quadratic factors of polynomial G are:" );
    factors = G.factorize();
    for( int i = 0; i < factors.size(); ++i )
        factors[ i ].print();

    vector< Complex > rs = F.complexRoots();
    for( int i = 0; i < rs.size(); ++i ){
        printf( "%.4lf  %+.4lf i\n", rs[ i ].real(), rs[ i ].imag() );
    }

    return 0;
}
