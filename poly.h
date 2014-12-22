/***************************************************************************************************

    File "poly.h"
    --------------------------------------------------------------------------------------------
    Class definition - fields and methods

    class Polynomial:

        Implements all the basic functionalities of a mathematical polynomial,
        including overloaded operators for a more natural usage on symbolic level.

    Author:
        zvonimir.jurelinac@fer.hr

***************************************************************************************************/

#include <string>
#include <vector>
#include <complex>

using std::string;
using std::vector;
using std::complex;

// Define a default complex number type - double precision
typedef complex< double > Complex;

// Arithmetic operations precision
const double EPS = 1e-12;

// Root-finding precision
const double EPS2 = 1e-10;

// Bairstow precision
const double EPS3 = 5e-3;

// Maximum number of iterations for Bairstow method
const int MAX_ITER = 100;

// Number of intervals of subdivision when using Newton's method
const int Intervals = 1000;


struct Poly{
    // Polynomial degree
    int N;

    // Array of coefficients
    double *a;

    // Polynomial name
    string name;

    // Empty constructor
    Poly();

    // Construct a zero polynomial with a given degree ( all coefficients are 0 )
    Poly( int degree );

    // Construct a polynomial from a list of coefficients
    Poly( int degree, double *coefficients );

    // Copy constructor - construct a polynomial from another one, with all the same coefficients
    Poly( const Poly &p );

    // Overloaded binary + operator, returns a sum of two polynomials
    friend Poly operator +( const Poly P, const Poly Q );

    // Overloaded binary - operator, returns a difference between two polynomials
    friend Poly operator -( const Poly P, const Poly Q );

    // Overloaded * operator, returns a result of polynomial multiplication of two polynomials
    friend Poly operator *( const Poly P, const Poly Q );

    // Overloaded * operator, returns a polynomial multiplied by a constant
    friend Poly operator *( const double q, const Poly P );

    // Overloaded / operator, returns a quotient of polynomial division
    friend Poly operator /( const Poly P, const Poly Q );

    // Overloaded % operator, returns a remainder of polynomial division
    friend Poly operator %( const Poly P, const Poly Q );

    // Overloaded unary - operator, returns a polynomial with all its coefficient negated
    Poly operator -();

    // Utility operator for testing whether a polynomial is equal to a constant d
    bool operator ==( double d );

    // Overloaded array indexing operator, returns i-th coefficient of a polynomial
    double &operator []( const int i );

    // Overloaded array indexing operator, returns i-th coefficient of a polynomial,
    // used with constant( non-modifiable ) polynomials
    double const &operator []( const int i ) const;

    // Polynomial evaluation at point x using Horner's method
    double operator()( double x );

    // Calculates a greatest common divisor of two polynomials using Euclid algorithm
    friend Poly gcd( Poly a, Poly b );

    // Returns a list of real roots of a polynomial, which are obtained using Newton's method
    vector< double > roots();

    // Utility function, tests whether polynomial is of degree 0
    bool isNumber();

    // Utility function, reduces the degree of a polynomial if it's leading coefficients are == 0,
    // used inside subtraction and division
    void reduce();

    // Returns a normalized polynomial so that it's leading coefficient is 1
    Poly normalized() const;

    // Return derivation P' of a polynomial P
    Poly differentiate();

    // Utility function, return a polynomial multiplied by x^k - all it's coefficients shifted k
    // places to the right
    Poly shift( int k );

    // Display a polynomial on screen
    void print() const;

    // Debug function, displays memory contents of the class etc.
    void dump() const;

    // Factorize polynomial using Bairstow method
    vector< Poly > factorize();

    // Obtain all roots (including complex) of polynomial
    vector< Complex > complexRoots();
};
