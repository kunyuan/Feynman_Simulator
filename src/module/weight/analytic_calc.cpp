//
//  analytic.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/10/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_inherit.h"
#include <vector>
using namespace weight;

void G::_InitialBare()
{
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
    }
}

void G::StartWithBare()
{
    DeltaTWeight = 0.0;
    SmoothWeight = BareWeight;
}

void W::_InitialBare()
{
    BareWeight = 0.0;
    int Lx = _Lat.Size[0], Ly = _Lat.Size[1];
    assert(Lx > 1 && Ly > 1 && D == 2);
    int spinindex = SpinIndex(UP, UP);
    int sublatA2B = _Lat.Sublat2Index(0, 1);
    vector<initializer_list<int>> coord;
    coord.push_back({0, 0});
    coord.push_back({0, Ly - 1});
    coord.push_back({Lx - 1, 0});
    coord.push_back({Lx - 1, Ly - 1});
    for (auto e : coord)
        BareWeight[spinindex][sublatA2B][_Lat.Vec2Index(e)] = 1.0;
    int sublatB2A = _Lat.Sublat2Index(1, 0);
    coord.clear();
    coord.push_back({0, 0});
    coord.push_back({0, 1});
    coord.push_back({1, 0});
    coord.push_back({1, 1});
    for (auto e : coord)
        BareWeight[spinindex][sublatB2A][_Lat.Vec2Index(e)] = 1.0;
}

void W::StartWithBare()
{
    DeltaTWeight = BareWeight;
    SmoothWeight = 0.0;
}
