//
//  weight.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight.h"
#include "component.h"
#include "module/parameter/parameter.h"

using namespace std;
using namespace para;

weight::Weight::Weight(bool IsAllSymmetric)
{
    _IsAllSymmetric = IsAllSymmetric;
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

bool weight::Weight::BuildNew(flag _flag, const ParaMC& para)
{
    if (para.Order == 0)
        ABORT("Order can not be zero!!!");
    //GW can only be loaded
    if (_flag & weight::SigmaPolar) {
        _AllocateSigmaPolar(para);
        Sigma->BuildNew();
        Polar->BuildNew();
    }
    return true;
}

void weight::Weight::ReWeight(flag _flag, const ParaMC& para)
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
bool weight::Weight::Load(const std::string& InputFile, flag _flag, const ParaMC& para)
{
    if (_flag & weight::GW) {
        _AllocateGW(para);
        G->Load(InputFile);
        W->Load(InputFile);
    }
    if (_flag & weight::SigmaPolar) {
        _AllocateSigmaPolar(para);
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
void weight::Weight::Save(const string& FileName, flag _flag, string Mode)
{
    //GW never saved here
    if (_flag & weight::SigmaPolar) {
        Sigma->Save(FileName, Mode);
        Polar->Save(FileName, "a");
        Mode = "a";
    }
}

int weight::Weight::UpdateSigmaPolarWeight(int OrderAccepted, real ErrorThreshold)
{
    int SigmaOrder = Sigma->Estimator.OrderAcceptable(OrderAccepted, ErrorThreshold);
    int PolarOrder = Polar->Estimator.OrderAcceptable(OrderAccepted, ErrorThreshold);
    int NewOrderAccepted = (SigmaOrder < PolarOrder ? SigmaOrder : PolarOrder);
    Sigma->UpdateWeight(NewOrderAccepted);
    Polar->UpdateWeight(NewOrderAccepted);
    return NewOrderAccepted;
}

void weight::Weight::SetTest(const ParaMC& para)
{
    _AllocateGW(para);
    _AllocateSigmaPolar(para);
    G->BuildTest(model::Trivial);
    W->BuildTest(model::Trivial);
}

void weight::Weight::SetDiagCounter(const ParaMC& para)
{
    _AllocateGW(para);
    _AllocateSigmaPolar(para);
    G->BuildTest(model::DiagCount);
    W->BuildTest(model::DiagCount);
}

void weight::Weight::_AllocateGW(const ParaMC& para)
{
    //make sure old Sigma/Polar/G/W are released before assigning new memory
    delete G;
    auto symmetry = _IsAllSymmetric ? TauSymmetric : TauAntiSymmetric;
    G = new weight::G(para.Lat, para.Beta, para.MaxTauBin, symmetry);
    delete W;
    W = new weight::W(para.Lat, para.Beta, para.MaxTauBin);
}

void weight::Weight::_AllocateSigmaPolar(const ParaMC& para)
{
    auto symmetry = _IsAllSymmetric ? TauSymmetric : TauAntiSymmetric;
    delete Sigma;
    Sigma = new weight::Sigma(para.Lat, para.Beta, para.MaxTauBin, para.Order, symmetry);
    delete Polar;
    Polar = new weight::Polar(para.Lat, para.Beta, para.MaxTauBin, para.Order);
}
