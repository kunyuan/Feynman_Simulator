//
//  dyson.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/4/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__dyson__
#define __Feynman_Simulator__dyson__

#include "utility/convention.h"
#include "utility/complex.h"

namespace weight {
class Weight;
class G;
class W;
class Sigma;
class Polar;
}
namespace para {
class Parameter;
}

namespace calc {
class Calculator {

  public:
    real Beta;
    int Vol;
    int TauBin;
    weight::G *G;
    weight::W *W;
    weight::Sigma *Sigma;
    weight::Polar *Polar;
    bool BuildNew(para::Parameter &, weight::Weight &);
    void DeriveG();
    void DeriveW();

    //  protected:
};

void MatrixInverse(Complex *, int);
void MatrixMultiply(Complex *, Complex *, int);
void MatrixMultiply(Complex *, Complex *, int, int);
int TestCalculator();
}

#endif /* defined(__Feynman_Simulator__dyson__) */
