//
//  dyson.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/4/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "dyson.h"
#include "module/parameter/parameter.h"
#include "module/observable/weight.h"

using namespace dyson;
using namespace weight;

bool Dyson::BuildNew(para::ParaDyson &para, Weight &weight)
{
    Beta = para.Beta;
    G = weight.G;
    W = weight.W;
    Sigma = weight.Sigma;
    Polar = weight.Polar;
    return true;
}
