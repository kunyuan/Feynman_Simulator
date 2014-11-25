//
//  weight_builder.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/22/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_initializer.h"
using namespace weight;

string InitialErrorMessage(model _Model, const Lattice &_Lat)
{
    return "Model" + ToString((int)_Model) + " on Lattice " +
           ToString((int)_Lat.LatticeType) +
           " has not yet been implemented!";
}

void GInitializer::BuildNew()
{
    switch (_G._Model) {
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
            if (_G._Lat.LatticeType == SQUARE)
                _InitialBareHubbardSquare();
            else
                ABORT(InitialErrorMessage(_G._Model, _G._Lat));
            break;
        default:
            ABORT(InitialErrorMessage(_G._Model, _G._Lat));
    }
}
// interaction

void WInitializer::BuildNew()
{
    switch (_W._Model) {
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
        default:
            ABORT(InitialErrorMessage(_W._Model, _W._Lat));
    }
}

void GInitializer::_InitialTest()
{
    //Spin independent; dr==0; exp(i*tau)
    _G._SmoothTWeight = 0.0;
    int spin_down = _Map.SpinIndex(DOWN, DOWN);
    int spin_up = _Map.SpinIndex(UP, UP);
    for (int sub = 0; sub < _G.GetShape()[SUB]; sub++) {
        if (!_G._Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _G._Lat.Vec2Index({0, 0});

        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = exp(Complex(0.0, _Map.IndexToTau(tau)));
            _G._SmoothTWeight[spin_down][sub][coor][tau] = weight;
            _G._SmoothTWeight[spin_up][sub][coor][tau] = weight;
        }
    }
}

void GInitializer::_InitialDiagCounter()
{
    _G._SmoothTWeight = 0.0;
    int spin_down = _Map.SpinIndex(DOWN, DOWN);
    int spin_up = _Map.SpinIndex(UP, UP);

    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (!_G._Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _G._Lat.Vec2Index({0, 0});
        //        MeasureWeight[spin_up][sub][coor] = Complex(1.0, 0.0);
        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = Complex(1.0, 0.0);
            _G._SmoothTWeight[spin_down][sub][coor][tau] = weight;
            _G._SmoothTWeight[spin_up][sub][coor][tau] = weight;
        }
    }
}

void GInitializer::_InitialBareSpin()
{
    _G._SmoothTWeight = 0.0;
    Complex mu = Complex(0.0, PI / 2.0 / _G._Beta);
    int spin_down = _Map.SpinIndex(DOWN, DOWN);
    int spin_up = _Map.SpinIndex(UP, UP);
    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (_G._Lat.IsOnSameSubLat(sub)) {
            int coor = 0;
            for (int tau = 0; tau < _Shape[TAU]; tau++) {
                Complex weight = exp(mu * _Map.IndexToTau(tau)) / Complex(1.0, 1.0);
                _G._SmoothTWeight[spin_down][sub][coor][tau] = weight;
                _G._SmoothTWeight[spin_up][sub][coor][tau] = weight;
            }
        }
    }
}

void GInitializer::_InitialBareHubbardSquare()
{
    _G._SmoothTWeight = 0.0;
    int Lx = _G._Lat.Size[0], Ly = _G._Lat.Size[1];
    ASSERT_ALLWAYS(Lx > 1 && Ly > 1, "System size should be bigger than 1!");
    ASSERT_ALLWAYS(_G._Lat.LatticeType == SQUARE, "lattice should be square!");
    ASSERT_ALLWAYS(_G._Model == HUBBARD, int(_G._Model) << " is not Hubbard model!");
    //Initialize in momentum space first
    //    int spin_down = SpinIndex(DOWN, DOWN);
    //    int spin_up = SpinIndex(UP, UP);
    //    int sub = _Lat.Sublat2Index(0, 0);

    //            _W._BareWeight[spin_down][sub][coor][tau] = weight;
    //            _W._BareWeight[spin_up][sub][coor][tau] = weight;
}

void WInitializer::_InitialTest()
{
    //spin conserved; dr==0; exp(-i*tau)
    _W._DeltaTWeight = 0.0;
    _W._SmoothTWeight = 0.0;
    int Lx = _W._Lat.Size[0], Ly = _W._Lat.Size[1];
    ASSERT_ALLWAYS(Lx > 1 && Ly > 1, "System size should be bigger than 1!");
    int spin_up = _Map.SpinIndex(UP,  //InOfW/InOfVertex
                                 UP,  //InOfW/OutOfVertex
                                 UP,  //OutOfW/InOfVertex
                                 UP); //OutOfW/OutOfVertex

    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (!_W._Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _W._Lat.Vec2Index({0, 0});

        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = exp(Complex(0.0, -_Map.IndexToTau(tau)));
            _W._SmoothTWeight[spin_up][sub][coor][tau] = weight;
        }
    }

    for (auto i : _Map.GetSpinIndexVector(UpUp2UpUp))
        _W._SmoothTWeight[i] = _W._SmoothTWeight[spin_up];

    for (auto i : _Map.GetSpinIndexVector(UpDown2UpDown)) {
        _W._SmoothTWeight[i] = _W._SmoothTWeight[spin_up];
        _W._SmoothTWeight[i] *= -1.0;
    }
    for (auto i : _Map.GetSpinIndexVector(UpDown2DownUp)) {
        _W._SmoothTWeight[i] = _W._SmoothTWeight[spin_up];
        _W._SmoothTWeight[i] *= 2.0;
    }
}

void WInitializer::_InitialDiagCounter()
{
    //spin==UP,UP,UP,UP; dr==0; independent of tau
    _W._DeltaTWeight = 0.0;
    _W._SmoothTWeight = 0.0;
    int Lx = _W._Lat.Size[0], Ly = _W._Lat.Size[1];
    ASSERT_ALLWAYS(Lx > 1 && Ly > 1, "System size should be bigger than 1!");
    int spin_up = _Map.SpinIndex(UP,  //InOfW/InOfVertex
                                 UP,  //InOfW/OutOfVertex
                                 UP,  //OutOfW/InOfVertex
                                 UP); //OutOfW/OutOfVertex

    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (!_W._Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _W._Lat.Vec2Index({0, 0});
        for (int tau = 0; tau < _Shape[TAU]; tau++) {
            Complex weight = Complex(1.0, 0.0);
            _W._SmoothTWeight[spin_up][sub][coor][tau] = weight;
        }
    }
}

void WInitializer::_InitialBareJ1J2()
{
    _W._DeltaTWeight = 0.0;
    _W._SmoothTWeight = 0.0;
    int Lx = _W._Lat.Size[0], Ly = _W._Lat.Size[1];
    ASSERT_ALLWAYS(Lx > 1 && Ly > 1, "System size should be bigger than 1!");
    ASSERT_ALLWAYS(_W._Lat.LatticeType == CHECKBOARD, "J1J2 lattice should be checkboard!");
    ASSERT_ALLWAYS(_W._Model == J1J2, ToString(int(_W._Model)) + " is not J1J2 model!");
    int spinindex = _Map.SpinIndex(UP,  //InOfW/InOfVertex
                                   UP,  //InOfW/OutOfVertex
                                   UP,  //OutOfW/InOfVertex
                                   UP); //OutOfW/OutOfVertex

    vector<initializer_list<int>> coord;
    int sublatA2B = _W._Lat.Sublat2Index(0, 1);
    int sublatA2A = _W._Lat.Sublat2Index(0, 0);
    {
        //J1 interaction A-->B
        coord.push_back({0, 0});
        coord.push_back({0, Ly - 1});
        coord.push_back({Lx - 1, 0});
        coord.push_back({Lx - 1, Ly - 1});
        for (auto e : coord)
            _W._DeltaTWeight[spinindex][sublatA2B][_W._Lat.Vec2Index(e)] = _W._Interaction[0];
        //J2 interaction A-->A
        coord.clear();
        coord.push_back({0, 1});
        coord.push_back({1, 0});
        coord.push_back({0, Ly - 1});
        coord.push_back({Lx - 1, 0});
        for (auto e : coord)
            _W._DeltaTWeight[spinindex][sublatA2A][_W._Lat.Vec2Index(e)] = _W._Interaction[1];
    }
    int sublatB2A = _W._Lat.Sublat2Index(1, 0);
    int sublatB2B = _W._Lat.Sublat2Index(1, 1);
    {
        //J1 interaction B-->A
        coord.clear();
        coord.push_back({0, 0});
        coord.push_back({0, 1});
        coord.push_back({1, 0});
        coord.push_back({1, 1});
        for (auto e : coord)
            _W._DeltaTWeight[spinindex][sublatB2A][_W._Lat.Vec2Index(e)] = _W._Interaction[0];
        //J2 interaction B-->B
        coord.clear();
        coord.push_back({0, 1});
        coord.push_back({1, 0});
        coord.push_back({0, Ly - 1});
        coord.push_back({Lx - 1, 0});
        for (auto e : coord)
            _W._DeltaTWeight[spinindex][sublatB2B][_W._Lat.Vec2Index(e)] = _W._Interaction[1];
    }
    //Generate other non-zero spin configuration
    {
        for (auto i : _Map.GetSpinIndexVector(UpUp2UpUp))
            _W._DeltaTWeight[i] = _W._DeltaTWeight[spinindex];

        for (auto i : _Map.GetSpinIndexVector(UpDown2UpDown)) {
            _W._DeltaTWeight[i] = _W._DeltaTWeight[spinindex];
            _W._DeltaTWeight[i] *= -1.0;
        }
        for (auto i : _Map.GetSpinIndexVector(UpDown2DownUp)) {
            _W._DeltaTWeight[i] = _W._DeltaTWeight[spinindex];
            _W._DeltaTWeight[i] *= 2.0;
        }
    }
}

void WInitializer::_InitialBareHubbard()
{
    int Lx = _W._Lat.Size[0], Ly = _W._Lat.Size[1];
    ASSERT_ALLWAYS(_W._Model == HUBBARD, int(_W._Model) << " is not Hubbard model!");
    _W._DeltaTWeight = 0.0;
    _W._SmoothTWeight = 0.0;
}