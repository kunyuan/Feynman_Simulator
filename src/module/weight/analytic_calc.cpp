//
//  analytic.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/10/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_inherit.h"

namespace weight {

void G::InitialWithBare()
{
    DeltaTWeight = 0.0;
    SmoothWeight = 0.0;
    BareWeight = 0.0;
    Complex mu = Complex(0.0, PI / 2.0 / _Beta);
    int spin_down = SpinIndex(DOWN, DOWN);
    int spin_up = SpinIndex(UP, UP);
    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (_Lat.IsLocal(sub))
            continue;
        int coor = 0;
        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = exp(mu * BinToTau(tau)) / Complex(1.0, 1.0);
            BareWeight[spin_down][sub][coor][tau] = weight;
            BareWeight[spin_up][sub][coor][tau] = weight;
        }

        BareWeight = 0.0;
        //TODO: add bare G initialization
    }
}

void W::InitialWithBare()
{
    DeltaTWeight = 0.0;
    SmoothWeight = 0.0;
    BareWeight = 0.0;
    //TODO: add bare W initialization
}
}
