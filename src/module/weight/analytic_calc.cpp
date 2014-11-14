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

void G::_InitialBareWeight()
{
    BareWeight = 0.0;
    Complex mu = Complex(0.0, PI / 2.0 / _Beta);
    int spin_down = SpinIndex(DOWN, DOWN);
    int spin_up = SpinIndex(UP, UP);
    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (_Lat.IsOnSameSubLat(sub))
            continue;
        int coor = 0;
        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = exp(mu * BinToTau(tau)) / Complex(1.0, 1.0);
            BareWeight[spin_down][sub][coor][tau] = weight;
            BareWeight[spin_up][sub][coor][tau] = weight;
        }
    }
}

void G::InitialWithBare()
{
    _InitialBareWeight();
    DeltaTWeight = 0.0;
    SmoothWeight = BareWeight;
}

void G::SetTest()
{
    //Spin independent; dr==0; exp(i*tau)
    SmoothWeight = 0.0;
    int spin_down = SpinIndex(DOWN, DOWN);
    int spin_up = SpinIndex(UP, UP);
    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (!_Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _Lat.Vec2Index({0,0});
        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = exp(Complex(0.0, BinToTau(tau)));
            SmoothWeight[spin_down][sub][coor][tau] = weight;
            SmoothWeight[spin_up][sub][coor][tau] = weight;
        }
    }
    
    DeltaTWeight = 0.0;
    BareWeight = 0.0;
}

void G::InitialWithDiagCounter()
{
    SmoothWeight = 0.0;
    int spin_down = SpinIndex(DOWN, DOWN);
    int spin_up = SpinIndex(UP, UP);
    
    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (!_Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _Lat.Vec2Index({0,0});
        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = Complex(1.0, 0.0);
            SmoothWeight[spin_down][sub][coor][tau] = weight;
            SmoothWeight[spin_up][sub][coor][tau] = weight;
        }
    }
    
    DeltaTWeight=0.0;
    BareWeight=0.0;
}

void W::_InitialBareWeight()
{
    BareWeight = 0.0;
    int Lx = _Lat.Size[0], Ly = _Lat.Size[1];
    assert(Lx > 1 && Ly > 1 && D == 2);
    int spinindex = SpinIndex(UP, //InOfW/InOfVertex
                              UP, //InOfW/OutOfVertex
                              UP, //OutOfW/InOfVertex
                              UP);//OutOfW/OutOfVertex

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
    
    int index = 0;
    for(int sinin=0; sinin<2; sinin++)
        for(int sinout =0; sinout<2; sinout++)
            for(int soutin=0; soutin<2; soutin++)
                for(int soutout=0; soutout<2; soutout++)
                {
                    index = SpinIndex(spin(sinin), spin(sinout), spin(soutin), spin(soutout));
                    BareWeight[index] = BareWeight[spinindex];
                    if(sinin==sinout && soutin==soutout){
                        if(sinin!=soutin)
                            BareWeight[index] *=-1.0;
                    }else if(sinin+soutin==sinout+soutout){
                        BareWeight[index] *= 2.0;
                        
                    }
                }
}

void W::InitialWithBare()
{
    _InitialBareWeight();
    DeltaTWeight = BareWeight;
    SmoothWeight = 0.0;
}

void W::SetTest()
{
    //spin conserved; dr==0; exp(-i*tau)
    SmoothWeight = 0.0;
    int Lx = _Lat.Size[0], Ly = _Lat.Size[1];
    assert(Lx > 1 && Ly > 1 && D == 2);
    int spin_up = SpinIndex(UP, //InOfW/InOfVertex
                            UP, //InOfW/OutOfVertex
                            UP, //OutOfW/InOfVertex
                            UP);//OutOfW/OutOfVertex

    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (!_Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _Lat.Vec2Index({0,0});
        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = exp(Complex(0.0, -BinToTau(tau)));
            SmoothWeight[spin_up][sub][coor][tau] = weight;
        }
    }
    
    int spinindex ;
    for(int sinin=0; sinin<2; sinin++)
        for(int sinout =0; sinout<2; sinout++)
            for(int soutin=0; soutin<2; soutin++)
                for(int soutout=0; soutout<2; soutout++)
                {
                    spinindex = SpinIndex(spin(sinin), spin(sinout), spin(soutin), spin(soutout));
                    SmoothWeight[spinindex] = SmoothWeight[spin_up];
                    if(sinin==sinout && soutin==soutout){
                        if(sinin!=soutin)
                            SmoothWeight[spinindex] *=-1.0;
                    }else if(sinin+soutin==sinout+soutout){
                        SmoothWeight[spinindex] *= 2.0;
                        
                    }
                }

    DeltaTWeight = 0.0;
    BareWeight = 0.0;
}

void W::InitialWithDiagCounter()
{
    //spin==UP,UP,UP,UP; dr==0; independent of tau
    SmoothWeight = 0.0;
    int Lx = _Lat.Size[0], Ly = _Lat.Size[1];
    assert(Lx > 1 && Ly > 1 && D == 2);
    int spin_up = SpinIndex(UP, //InOfW/InOfVertex
                            UP, //InOfW/OutOfVertex
                            UP, //OutOfW/InOfVertex
                            UP);//OutOfW/OutOfVertex

    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (!_Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _Lat.Vec2Index({0,0});
        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = Complex(1.0, 0.0);
            SmoothWeight[spin_up][sub][coor][tau] = weight;
        }
    }
    
    DeltaTWeight=0.0;
    BareWeight=0.0;
}

