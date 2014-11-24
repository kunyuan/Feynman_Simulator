//
//  weight_builder.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/22/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight_initializer.h"
#include "module/weight/component.h"
using namespace weight0;

string InitialErrorMessage(model _Model, const Lattice &_Lat)
{
    return "Model" + ToString((int)_Model) + " on Lattice " +
           ToString((int)_Lat.LatticeType) +
           " has not yet been implemented!";
}

void GInitializer::BuildTest()
{
    _InitialTest();
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

void GInitializer::_InitialTest()
{
    //Spin independent; dr==0; exp(i*tau)
    _G._SmoothTWeight = 0.0;
    int spin_down = _G.SpinIndex(DOWN, DOWN);
    int spin_up = _G.SpinIndex(UP, UP);
    for (int sub = 0; sub < _G.GetShape()[SUB]; sub++) {
        if (!_G._Lat.IsOnSameSubLat(sub))
            continue;
        int coor = _G._Lat.Vec2Index({0, 0});

        for (int tau = 0; tau < _G.GetShape()[TAU]; tau++) {
            Complex weight = exp(Complex(0.0, _G.IndexToTau(tau)));
            _G._SmoothTWeight[spin_down][sub][coor][tau] = weight;
            _G._SmoothTWeight[spin_up][sub][coor][tau] = weight;
        }
    }
    _G._BareWeight = 0.0;
}

void GInitializer::_InitialDiagCounter()
{
    uint *_Shape = _G.GetShape().data();
    _G._SmoothTWeight = 0.0;
    int spin_down = _G.SpinIndex(DOWN, DOWN);
    int spin_up = _G.SpinIndex(UP, UP);

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
    _G._BareWeight = 0.0;
}

void GInitializer::_InitialBareSpin()
{
    uint *_Shape = _G.GetShape().data();
    _G._BareWeight = 0.0;
    Complex mu = Complex(0.0, PI / 2.0 / _G._Beta);
    int spin_down = _G.SpinIndex(DOWN, DOWN);
    int spin_up = _G.SpinIndex(UP, UP);
    for (int sub = 0; sub < _Shape[SUB]; sub++) {
        if (_G._Lat.IsOnSameSubLat(sub)) {
            int coor = 0;
            for (int tau = 0; tau < _Shape[TAU]; tau++) {
                Complex weight = exp(mu * _G.IndexToTau(tau)) / Complex(1.0, 1.0);
                _G._BareWeight[spin_down][sub][coor][tau] = weight;
                _G._BareWeight[spin_up][sub][coor][tau] = weight;
            }
        }
    }
    _G._SmoothTWeight = _G._BareWeight;
}

void GInitializer::_InitialBareHubbardSquare()
{
    _G._BareWeight = 0.0;
    int Lx = _G._Lat.Size[0], Ly = _G._Lat.Size[1];
    ASSERT_ALLWAYS(Lx > 1 && Ly > 1, "System size should be bigger than 1!");
    ASSERT_ALLWAYS(_G._Lat.LatticeType == SQUARE, "lattice should be square!");
    ASSERT_ALLWAYS(_G._Model == HUBBARD, int(_G._Model) << " is not Hubbard model!");
    //Initialize in momentum space first
    //    int spin_down = SpinIndex(DOWN, DOWN);
    //    int spin_up = SpinIndex(UP, UP);
    //    int sub = _Lat.Sublat2Index(0, 0);

    //            BareWeight[spin_down][sub][coor][tau] = weight;
    //            BareWeight[spin_up][sub][coor][tau] = weight;

    _G._SmoothTWeight = _G._BareWeight;
}