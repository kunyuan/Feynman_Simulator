//
//  weight.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight.h"
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

bool weight::Weight::BuildNew(flag _flag, const Parameter &para)
{
    if (para.Order == 0)
        ABORT("Order can not be zero!!!");
    if (_flag & weight::GW) {
        _AllocateGW(para);
        G->Initial(MODEL);
        W->Initial(MODEL);
    }
    if (_flag & weight::SigmaPolar) {
        _AllocateSigmaPolar(para);
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

int weight::Weight::UpdateSigmaPolarWeight(int OrderAccepted, real ErrorThreshold)
{
    int SigmaOrder = Sigma->OrderAcceptable(OrderAccepted, ErrorThreshold);
    int PolarOrder = Polar->OrderAcceptable(OrderAccepted, ErrorThreshold);
    int NewOrderAccepted = (SigmaOrder < PolarOrder ? SigmaOrder : PolarOrder);
    Sigma->UpdateWeight(NewOrderAccepted);
    Polar->UpdateWeight(NewOrderAccepted);
    return NewOrderAccepted;
}

void weight::Weight::SetTest(const Parameter &para)
{
    _AllocateGW(para);
    _AllocateSigmaPolar(para);
    G->Initial(model::TEST);
    W->Initial(model::TEST);
}

void weight::Weight::SetDiagCounter(const Parameter &para)
{
    _AllocateGW(para);
    _AllocateSigmaPolar(para);
    G->Initial(model::DIAGRAMCOUNTER);
    W->Initial(model::DIAGRAMCOUNTER);
}

void weight::Weight::_AllocateGW(const Parameter &para)
{
    //make sure old Sigma/Polar/G/W are released before assigning new memory
    delete G;
    G = new weight::G(para.Lat, para.Beta, para.Order,
                      para.Hopping, para.ExternalField,
                      para.RealChemicalPotential, _IsAllSymmetric);
    delete W;
    W = new weight::W(para.Lat, para.Beta, para.Order,
                      para.Interaction, para.ExternalField);
}

void weight::Weight::_AllocateSigmaPolar(const Parameter &para)
{
    delete Sigma;
    Sigma = new weight::Sigma(para.Lat, para.Beta, para.Order, _IsAllSymmetric);
    delete Polar;
    Polar = new weight::Polar(para.Lat, para.Beta, para.Order);
}
