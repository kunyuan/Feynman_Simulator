//---------------------------------------------------------------------------
// ©Nikolai V. Shokhirev, 2004-2008  <nikolai@shokhirev.com>
// http://www.shokhirev.com/nikolai.html
// Variant 2 - based on compound assignment
//---------------------------------------------------------------------------
// Based on:
//   Danny Kalev. Overloading Operator + the Right Way,
//   http://gethelp.devx.com/techtips/cpp_pro/10min/2001/august/10min0801.asp
//   ClÈment Lavoillotte. Yet another complex number class in C++
//   http://clement.lavoillotte.info/cpp_complex
//   Donnie Pinkston. C++ Operator Overloading Guidelines.
//   Mayukh Bose. C++ Operator Overloading Tutorial.
//   http://www.mayukhbose.com/tutorials/overloading/
//
#ifndef Complex_h
#define Complex_h

#include <iosfwd>
#include <iomanip>
#include <math.h>
#include "utility.h"

// Complex number class
class Complex {
public:
    Complex(); // Default constructor
    Complex(real re, real im = 0.0);
    Complex(const Complex& c); // Copy constructor

    real Re; // real part
    real Im; // imaginary part

    Complex& operator=(const Complex&);
    Complex& operator=(const real&);

    // define the compound assignment operators first
    Complex& operator+=(const Complex&);
    Complex& operator+=(const real&);
    Complex& operator-=(const Complex&);
    Complex& operator-=(const real&);
    Complex& operator*=(const Complex&);
    Complex& operator*=(const real&);
    Complex& operator/=(const Complex&);
    Complex& operator/=(const real&);

    friend std::ostream& operator<<(std::ostream&, const Complex&);
    friend std::istream& operator>>(std::istream&, Complex&);
};

inline bool Equal(const Complex& c1, const Complex& c2, real eps = eps0)
{
    return (fabs(c1.Re - c2.Re) < eps && fabs(c1.Im - c2.Im) < eps);
}

inline bool Equal(const Complex& c1, real r, real i, real eps = eps0)
{
    return (fabs(c1.Re - r) < eps && fabs(c1.Im - i) < eps);
}

inline bool IsZero(const Complex& c)
{
    return (c.Re == 0.0) && (c.Im == 0.0);
}

// Nonmember operators (to allow implicit conversion of the left operand)
Complex operator+(const Complex&, const Complex&);
Complex operator+(const Complex&); // unary + operator
Complex operator-(const Complex&, const Complex&);
Complex operator-(const Complex&); // unary - operator
Complex operator*(const Complex&, const Complex&);
Complex operator/(const Complex&, const Complex&);

Complex operator+(const Complex&, real);
Complex operator+(real, const Complex&);
Complex operator-(const Complex&, real);
Complex operator-(real, const Complex&);
Complex operator*(real, const Complex&);
Complex operator*(const Complex&, real);
Complex operator/(const Complex&, real);
// Complex library

// mod2 = Re*Re + Im*Im
inline real mod2(const Complex& c)
{
    return (c.Re * c.Re + c.Im * c.Im);
}

// sqrt(Re*Re + Im*Im)
inline real mod(const Complex& c)
{
    return sqrt(mod2(c));
}

inline real arg(const Complex& c)
{
    return c.Im != 0.0 ? atan2(c.Im, c.Re) : 0.0;
}

inline Complex phase(const Complex& c)
{
    if (IsZero(c))
        return 0.0;
    return c / mod(c);
}

inline Complex exp(const Complex& c)
{
    return (exp(c.Re) * Complex(cos(c.Im), sin(c.Im)));
}

// Return the principal branch of the square root (non-negative real part).
inline Complex sqrt(const Complex& x)
{
    real mag = mod(x);
    if (mag == 0.0)
        return Complex(0.0, 0.0);
    else if (x.Re > 0) {
        real re = sqrt(0.5 * (mag + x.Re));
        return Complex(re, 0.5 * x.Im / re);
    }
    else {
        real im = sqrt(0.5 * (mag - x.Re));
        if (x.Im < 0)
            im = -im;
        return Complex(0.5 * x.Im / im, im);
    }
}

inline Complex polar(real r, real t)
{
    return Complex(r * cos(t), r * sin(t));
}

// Complex exponentiation
inline Complex pow(const Complex& z, const Complex& w)
{
    real u = w.Re;
    real v = w.Im;
    if (IsZero(z))
        return IsZero(w) ? 1.0 : 0.0;
    real logr = 0.5 * log(mod2(z));
    real th = (arg(z));
    real phi = logr * v + th * u;
    return exp(logr * u - th * v) * Complex(cos(phi), sin(phi));
}

inline Complex pow(const Complex& z, real u)
{
    if (IsZero(0.0))
        return u == 0.0 ? 1.0 : 0.0;
    real logr = 0.5 * log(mod2(z));
    real theta = u * arg(z);
    return exp(logr * u) * Complex(cos(theta), sin(theta));
}

#endif
