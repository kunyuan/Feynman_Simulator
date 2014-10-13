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

#include <string>
#include <iostream>
#include <iomanip>
#include <math.h>
//#include <limits>
#include "utility.h"

// Complex number class
class Complex {
  public:
    Complex(); // Default constructor
    Complex(real re, real im = 0.0);
    Complex(const Complex &c); // Copy constructor

    ~Complex(); // Destructor

    real Re; // real part
    real Im; // imaginary part

    Complex &operator=(const Complex &);
    Complex &operator=(const real &);

    // define the compound assignment operators first
    Complex &operator+=(const Complex &);
    Complex &operator+=(const real &);
    Complex &operator-=(const Complex &);
    Complex &operator-=(const real &);
    Complex &operator*=(const Complex &);
    Complex &operator*=(const real &);
    Complex &operator/=(const Complex &);
    Complex &operator/=(const real &);

    friend std::ostream &operator<<(std::ostream &, const Complex &);
    friend std::istream &operator>>(std::istream &, Complex &);

    friend bool IsZero(const Complex &c);
};

// Nonmember operators (to allow implicit conversion of the left operand)
Complex operator+(const Complex &, const Complex &);
Complex operator+(const Complex &); // unary + operator
Complex operator-(const Complex &, const Complex &);
Complex operator-(const Complex &); // unary - operator
Complex operator*(const Complex &, const Complex &);
Complex operator/(const Complex &, const Complex &);

// Complex library

// mod2 = Re*Re + Im*Im
real mod2(const Complex &c);

// sqrt(Re*Re + Im*Im)
real mod(const Complex &c);

bool Equal(Complex c1, Complex c2, real eps = eps0);

bool Equal(Complex c1, real r, real i, real eps = eps0);

bool IsZero(const Complex &c);

Complex exp(const Complex &c);

std::string str(Complex c);
// More complex functions ...

#endif
