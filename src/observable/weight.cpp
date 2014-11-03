//
//  weight.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight.h"

using namespace std;

weight::Weight::Weight()
{
    Sigma = nullptr;
    Polar = nullptr;
    G = nullptr;
    W = nullptr;
    WormWeight = nullptr;
}

weight::Weight::~Weight()
{
    delete Sigma;
    delete Polar;
    delete G;
    delete W;
    delete WormWeight;
}

bool weight::Weight::BuildNew(const Lattice &Lat, real Beta, int order)
{
    _AllocateResources(Lat, Beta, order);
    G->InitialWithBare();
    W->InitialWithBare();
    Sigma->ClearStatistics();
    Polar->ClearStatistics();
    return true;
}

bool weight::Weight::Load(const std::string &InputFile, const Lattice &Lat, real Beta, int order)
{
    _AllocateResources(Lat, Beta, order);
    G->Load(InputFile);
    W->Load(InputFile);
    Sigma->Load(InputFile);
    Polar->Load(InputFile);
    return true;
}

void weight::Weight::Save(const std::string &InputFile, const std::string &Mode)
{
    G->Save(InputFile, Mode);
    W->Save(InputFile, "a");
    Sigma->Save(InputFile, "a");
    Polar->Save(InputFile, "a");
}

void weight::Weight::_AllocateResources(const Lattice &Lat, real Beta, int order)
{
    //make sure old Sigma/Polar/G/W are released before assigning new memory
    delete Sigma;
    Sigma = new weight::Sigma(Lat, Beta, order);
    delete Polar;
    Polar = new weight::Polar(Lat, Beta, order);
    delete G;
    G = new weight::G(Lat, Beta, order);
    delete W;
    W = new weight::W(Lat, Beta, order);
    delete WormWeight;
    WormWeight = new weight::Worm();
}
