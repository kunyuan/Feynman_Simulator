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
}

weight::Weight::~Weight()
{
    delete Sigma;
    delete Polar;
    delete G;
    delete W;
}
/**
*  Build G, W, Sigma, Polar from file, you may use flag weight::GW and weight::SigmaPolar to control which group to load. Notice those in unflaged group will remain the same.
*
*  @param _flag     weight::GW to allocate memory for GW, weight::SigmaPolar to allocate memory for Sigma and Polar
*  @param Lat       Lattice
*  @param Beta      Beta
*  @param order     order
*/

bool weight::Weight::BuildNew(flag _flag, const Parameter &para)
{
    if (_flag & weight::GW) {
        _AllocateGW(para.Lat, para.Beta, para.Order);
        G->InitialWithBare();
        W->InitialWithBare();
    }
    if (_flag & weight::SigmaPolar) {
        _AllocateSigmaPolar(para.Lat, para.Beta, para.Order);
        Sigma->ClearStatistics();
        Polar->ClearStatistics();
    }
    return true;
}

void weight::Weight::ReWeight(flag _flag, const Parameter &para)
{
    if (_flag & weight::GW) {
        //TODO: reweight G, W
    }
    if (_flag & weight::SigmaPolar) {
        //TODO: reweight SigmaPolar
    }
}
/**
*  Load G, W, Sigma, Polar from file, you may use flag weight::GW and weight::SigmaPolar to control which group to load. Notice those in unflaged group will remain the same.
*
*  @param InputFile Input File Name
*  @param _flag     weight::GW|weight::SigmaPolar
*  @param Lat       Lattice
*  @param Beta      Beta
*  @param order     order
*
*  @return alway true for now
*/
bool weight::Weight::Load(const std::string &InputFile, flag _flag, const Parameter &para)
{
    if (_flag & weight::GW) {
        _AllocateGW(para.Lat, para.Beta, para.Order);
        G->Load(InputFile);
        W->Load(InputFile);
    }
    if (_flag & weight::SigmaPolar) {
        _AllocateSigmaPolar(para.Lat, para.Beta, para.Order);
        Sigma->Load(InputFile);
        Polar->Load(InputFile);
    }
    return true;
}
/**
*  Save G, W, Sigma, Polar to file, you may use flag weight::GW and weight::SigmaPolar to control which group to save
*
*  @param InputFile output file name
*  @param _flag     weight::GW|weight::SigmaPolar
*  @param Mode      "a" or "w"
*/
void weight::Weight::Save(const string &FileName, flag _flag, string Mode)
{
    if (_flag & weight::GW) {
        G->Save(FileName, Mode);
        W->Save(FileName, "a");
        Mode = "a";
    }
    if (_flag & weight::SigmaPolar) {
        Sigma->Save(FileName, Mode);
        Polar->Save(FileName, "a");
        Mode = "a";
    }
}

void weight::Weight::SetTest(const Parameter &para)
{
    _AllocateGW(para.Lat, para.Beta, para.Order);
    _AllocateSigmaPolar(para.Lat, para.Beta, para.Order);
    G->SetTest();
    W->SetTest();
    Sigma->SetTest();
    Polar->SetTest();
}

void weight::Weight::_AllocateGW(const Lattice &Lat, real Beta, int order)
{
    //make sure old Sigma/Polar/G/W are released before assigning new memory
    delete G;
    G = new weight::G(Lat, Beta, order);
    delete W;
    W = new weight::W(Lat, Beta, order);
}

void weight::Weight::_AllocateSigmaPolar(const Lattice &Lat, real Beta, int order)
{
    delete Sigma;
    Sigma = new weight::Sigma(Lat, Beta, order);
    delete Polar;
    Polar = new weight::Polar(Lat, Beta, order);
}
