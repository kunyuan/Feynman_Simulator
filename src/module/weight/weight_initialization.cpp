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

void G::Initial(model Model_)
{
    _Model = Model_;
    switch (_Model) {
        case TEST:
            _InitialTest();
            break;
        case DIAGRAMCOUNTER:
            _InitialDiagCounter();
            break;
        case J1J2:
            _InitialBareSpin();
            break;
        case HUBBARD:
            _InitialBareHubbard();
            break;
    }
}

void G::_InitialTest()
{
    //Spin independent; dr==0; exp(i*tau)
    SmoothWeight = 0.0;
    MeasureWeight = 0.0;
    int spin_down = SpinIndex(DOWN, DOWN);
    int spin_up = SpinIndex(UP, UP);
    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (!_Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _Lat.Vec2Index({0, 0});

        MeasureWeight[spin_down][sub][coor] = Complex(1.0, 0.0);
        MeasureWeight[spin_up][sub][coor] = Complex(1.0, 0.0);

        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = exp(Complex(0.0, BinToTau(tau)));
            SmoothWeight[spin_down][sub][coor][tau] = weight;
            SmoothWeight[spin_up][sub][coor][tau] = weight;
        }
    }

    DeltaTWeight = 0.0;
    BareWeight = 0.0;
}

void G::_InitialDiagCounter()
{
    SmoothWeight = 0.0;
    MeasureWeight = 0.0;
    int spin_down = SpinIndex(DOWN, DOWN);
    int spin_up = SpinIndex(UP, UP);

    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (!_Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _Lat.Vec2Index({0, 0});
        MeasureWeight[spin_up][sub][coor] = Complex(1.0, 0.0);
        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = Complex(1.0, 0.0);
            SmoothWeight[spin_down][sub][coor][tau] = weight;
            SmoothWeight[spin_up][sub][coor][tau] = weight;
        }
    }

    DeltaTWeight = 0.0;
    BareWeight = 0.0;
}

void G::_InitialBareSpin()
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
    DeltaTWeight = 0.0;
    SmoothWeight = BareWeight;
    MeasureWeight = Complex(1.0, 0.0);
}

void G::_InitialBareHubbard()
{
    BareWeight = 0.0;
    //Initialize in momentum space first
    //    int spin_down = SpinIndex(DOWN, DOWN);
    //    int spin_up = SpinIndex(UP, UP);
    //    for (int sub = 0; sub < _Shape[SUB]; sub++) {
    //        if (_Lat.IsOnSameSubLat(sub))
    //            continue;
    //        int coor = 0;
    //        for (int tau = 0; tau < _Shape[TAU]; tau++) {
    //            //            BareWeight[spin_down][sub][coor][tau] = weight;
    //            //            BareWeight[spin_up][sub][coor][tau] = weight;
    //        }
    //    }

    DeltaTWeight = 0.0;
    SmoothWeight = BareWeight;
}

// interaction

void W::Initial(model Model_)
{
    _Model = Model_;
    switch (_Model) {
        case TEST:
            _InitialTest();
            break;
        case DIAGRAMCOUNTER:
            _InitialDiagCounter();
            break;
        case J1J2:
            _InitialBareJ1J2();
            break;
        case HUBBARD:
            _InitialBareHubbard();
            break;
    }
}

void W::_InitialTest()
{
    //spin conserved; dr==0; exp(-i*tau)
    DeltaTWeight = 0.0;
    BareWeight = 0.0;

    SmoothWeight = 0.0;
    MeasureWeight = 0.0;
    int Lx = _Lat.Size[0], Ly = _Lat.Size[1];
    ASSERT_ALLWAYS(Lx > 1 && Ly > 1, "System size should be bigger than 1!");
    int spin_up = SpinIndex(UP,  //InOfW/InOfVertex
                            UP,  //InOfW/OutOfVertex
                            UP,  //OutOfW/InOfVertex
                            UP); //OutOfW/OutOfVertex

    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (!_Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _Lat.Vec2Index({0, 0});

        MeasureWeight[spin_up][sub][coor] = Complex(1.0, 0.0);

        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = exp(Complex(0.0, -BinToTau(tau)));
            SmoothWeight[spin_up][sub][coor][tau] = weight;
        }
    }

    for (auto i : GetSpinIndexVector_FourSpinsFileter(UpUp2UpUp))
        SmoothWeight[i] = SmoothWeight[spin_up];

    for (auto i : GetSpinIndexVector_FourSpinsFileter(UpDown2UpDown)) {
        SmoothWeight[i] = SmoothWeight[spin_up];
        SmoothWeight[i] *= -1.0;
    }
    for (auto i : GetSpinIndexVector_FourSpinsFileter(UpDown2DownUp)) {
        SmoothWeight[i] = SmoothWeight[spin_up];
        SmoothWeight[i] *= 2.0;
    }
}

void W::_InitialDiagCounter()
{
    //spin==UP,UP,UP,UP; dr==0; independent of tau
    SmoothWeight = 0.0;
    MeasureWeight = 0.0;
    int Lx = _Lat.Size[0], Ly = _Lat.Size[1];
    ASSERT_ALLWAYS(Lx > 1 && Ly > 1, "System size should be bigger than 1!");
    int spin_up = SpinIndex(UP,  //InOfW/InOfVertex
                            UP,  //InOfW/OutOfVertex
                            UP,  //OutOfW/InOfVertex
                            UP); //OutOfW/OutOfVertex

    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (!_Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _Lat.Vec2Index({0, 0});
        MeasureWeight[spin_up][sub][coor] = Complex(1.0, 0.0);
        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = Complex(1.0, 0.0);
            SmoothWeight[spin_up][sub][coor][tau] = weight;
        }
    }

    DeltaTWeight = 0.0;
    BareWeight = 0.0;
}

void W::_InitialBareJ1J2()
{
    BareWeight = 0.0;
    int Lx = _Lat.Size[0], Ly = _Lat.Size[1];
    ASSERT_ALLWAYS(Lx > 1 && Ly > 1, "System size should be bigger than 1!");
    ASSERT_ALLWAYS(_Lat.LatticeType == CHECKBOARD, "J1J2 lattice should be checkboard!");
    ASSERT_ALLWAYS(_Model == J1J2, ToString(int(_Model)) + " is not J1J2 model!");
    int spinindex = SpinIndex(UP,  //InOfW/InOfVertex
                              UP,  //InOfW/OutOfVertex
                              UP,  //OutOfW/InOfVertex
                              UP); //OutOfW/OutOfVertex

    vector<initializer_list<int>> coord;
    int sublatA2B = _Lat.Sublat2Index(0, 1);
    int sublatA2A = _Lat.Sublat2Index(0, 0);
    {
        //J1 interaction A-->B
        coord.push_back({0, 0});
        coord.push_back({0, Ly - 1});
        coord.push_back({Lx - 1, 0});
        coord.push_back({Lx - 1, Ly - 1});
        for (auto e : coord)
            BareWeight[spinindex][sublatA2B][_Lat.Vec2Index(e)] = _Interaction[0];
        //J2 interaction A-->A
        coord.clear();
        coord.push_back({0, 1});
        coord.push_back({1, 0});
        coord.push_back({0, Ly - 1});
        coord.push_back({Lx - 1, 0});
        for (auto e : coord)
            BareWeight[spinindex][sublatA2A][_Lat.Vec2Index(e)] = _Interaction[1];
    }
    int sublatB2A = _Lat.Sublat2Index(1, 0);
    int sublatB2B = _Lat.Sublat2Index(1, 1);
    {
        //J1 interaction B-->A
        coord.clear();
        coord.push_back({0, 0});
        coord.push_back({0, 1});
        coord.push_back({1, 0});
        coord.push_back({1, 1});
        for (auto e : coord)
            BareWeight[spinindex][sublatB2A][_Lat.Vec2Index(e)] = _Interaction[0];
        //J2 interaction B-->B
        coord.clear();
        coord.push_back({0, 1});
        coord.push_back({1, 0});
        coord.push_back({0, Ly - 1});
        coord.push_back({Lx - 1, 0});
        for (auto e : coord)
            BareWeight[spinindex][sublatB2B][_Lat.Vec2Index(e)] = _Interaction[1];
    }
    //Generate other non-zero spin configuration
    {
        for (auto i : GetSpinIndexVector_FourSpinsFileter(UpUp2UpUp))
            BareWeight[i] = BareWeight[spinindex];

        for (auto i : GetSpinIndexVector_FourSpinsFileter(UpDown2UpDown)) {
            BareWeight[i] = BareWeight[spinindex];
            BareWeight[i] *= -1.0;
        }
        for (auto i : GetSpinIndexVector_FourSpinsFileter(UpDown2DownUp)) {
            BareWeight[i] = BareWeight[spinindex];
            BareWeight[i] *= 2.0;
        }
    }
    DeltaTWeight = BareWeight;
    SmoothWeight = 0.0;
    MeasureWeight = Complex(1.0, 0.0);
}

void W::_InitialBareHubbard()
{
    BareWeight = 0.0;
    int Lx = _Lat.Size[0], Ly = _Lat.Size[1];
    ASSERT_ALLWAYS(Lx > 1 && Ly > 1, "System size should be bigger than 1!");
    ASSERT_ALLWAYS(_Lat.LatticeType == CHECKBOARD, "J1J2 lattice should be checkboard!");
    ASSERT_ALLWAYS(_Model == J1J2, ToString(int(_Model)) + " is not J1J2 model!");
    DeltaTWeight = BareWeight;
    SmoothWeight = 0.0;
    MeasureWeight = Complex(1.0, 0.0);
}
