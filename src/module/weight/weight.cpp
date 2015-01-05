//
//  weight.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight.h"
#include "component.h"
#include "utility/dictionary.h"
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

bool weight::Weight::FromDict(const Dictionary& dict, flag _flag, const para::ParaMC& para)
{
    if (_flag & weight::GW) {
        _AllocateGW(para);
        G->FromDict(dict.Get<Dictionary>("G"));
        W->FromDict(dict.Get<Dictionary>("W"));
    }
    if (_flag & weight::SigmaPolar) {
        _AllocateSigmaPolar(para);
        Sigma->FromDict(dict.Get<Dictionary>("Sigma"));
        Polar->FromDict(dict.Get<Dictionary>("Polar"));
    }
    return true;
}
Dictionary weight::Weight::ToDict(flag _flag)
{
    Dictionary dict;
    if (_flag & weight::GW) {
        dict["G"] = G->ToDict();
        dict["W"] = W->ToDict();
    }
    if (_flag & weight::SigmaPolar) {
        dict["Sigma"] = Sigma->ToDict();
        dict["Polar"] = Polar->ToDict();
    }
    return dict;
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
