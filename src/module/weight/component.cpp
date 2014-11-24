//
//  component.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/22/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "component.h"
#include "weight_initializer.h"
#include "utility/cnpy.h"

using namespace weight0;
using namespace std;

template <typename T>
bool LoadMatrix(T &matrix, const string &FileName, const string &Name)
{
    cnpy::NpyArray weight = cnpy::npz_load(FileName, Name);
    if (weight.data == nullptr) {
        ABORT("Can't find " << Name << ".Weight in .npz data file!");
        return false;
    }
    //assignment here will copy data in weight.data into *this
    matrix = (reinterpret_cast<Complex *>(weight.data));
    return true;
}

template <typename T>
void SaveMatrix(T &matrix, const string &FileName, const std::string Mode,
                const string &Name, const vector<uint> &Shape, int Dim)
{
    cnpy::npz_save(FileName, Name, matrix(), Shape.data(), Dim, Mode);
}

G::G(model Model, const Lattice &lat, real beta,
     const std::vector<real> &Hopping,
     const std::vector<real> &RealChemicalPotential,
     real ExternalField, TauSymmetry Symmetry)
    : weight0::BasicWithTwoSpins(lat, beta, Model, Symmetry, "G")
{
    _Hopping = Hopping;
    _ExternalField = ExternalField;
    _RealChemicalPotential = RealChemicalPotential;
    _SmoothTWeight.Allocate(GetShape().data());
    _BareWeight.Allocate(GetShape().data());
    //use _Shape[SP] to _Shape[TAU] to construct array3
}

void G::BuildNew()
{
    GInitializer(*this).BuildNew();
}

void G::BuildTest()
{
    GInitializer(*this).BuildTest();
}

bool G::Load(const std::string &FileName)
{
    return LoadMatrix(_SmoothTWeight, FileName, _Name + ".Smooth") &&
           LoadMatrix(_BareWeight, FileName, _Name + ".Bare");
}
void G::Save(const std::string &FileName, const std::string Mode)
{
    SaveMatrix(_SmoothTWeight, FileName, Mode, _Name + ".Smooth", GetShape(), 4);
    SaveMatrix(_BareWeight, FileName, "a", _Name + ".Bare", GetShape(), 3);
}

W::W(model Model, const Lattice &lat, real Beta,
     const vector<real> &Interaction_, real ExternalField)
    : weight0::BasicWithFourSpins(lat, Beta, Model, TauSymmetric, "W")
{
    _Interaction = Interaction_;
    _ExternalField = ExternalField;
    _BareWeight.Allocate(GetShape().data());
    _SmoothTWeight.Allocate(GetShape().data());
    //use _Shape[SP] to _Shape[VOL] to construct array3
}
