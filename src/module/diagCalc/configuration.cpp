//
// Created by yuan huang on 10/11/17.
//
#include "diag_calculator.h"
#include "utility/dictionary.h"
#include "lattice/lattice.h"
#include <iostream>
#include <array>

using namespace diagCalc;

void Conf::Initialize(){
    Order = 2;
    MeasureGLine = true;

    for (int i = 0; i < Order*2; i++){
        R_list[i] = Site(0, {0,0});
        Tau_list[i] = 0.0;
    }
}

