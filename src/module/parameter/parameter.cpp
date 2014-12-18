//
//  state.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "parameter.h"
#include "utility/utility.h"
using namespace para;

Message Parameter::GenerateMessage()
{
    Message Message_;
    Message_.Version = Version;
    Message_.Beta = Beta;
    return Message_;
}

void Parameter::UpdateWithMessage(const Message& Message_)
{
    Version = Message_.Version;
    Beta = Message_.Beta;
    T = 1.0 / Beta;
}

bool Parameter::_BuildNew(const std::string& InputFile)
{
    _para.ParseFile(InputFile);
    GetPara(_para, L);
    GetPara(_para, InitialBeta);
    GetPara(_para, DeltaBeta);
    GetPara(_para, FinalBeta);
    GetPara(_para, Order);
    GetPara(_para, NSublat)

        Lat.Initialize(L, NSublat);
    Version = 0;
    Beta = InitialBeta;
    T = 1.0 / Beta;
    if (Order >= MAX_ORDER)
        ABORT("Order can not be bigger than " << MAX_ORDER);
    return true;
}

bool Parameter::_Load(const std::string& InputFile)
{
    _para.ParseFile(InputFile);
    GetPara(_para, Version);
    GetPara(_para, L);
    GetPara(_para, InitialBeta);
    GetPara(_para, DeltaBeta);
    GetPara(_para, FinalBeta);
    //!!!Beta should be a part of state, so it will be stored
    GetPara(_para, Beta);
    GetPara(_para, Order);
    GetPara(_para, NSublat)

        Lat.Initialize(L, NSublat);
    T = 1.0 / Beta;
    if (Order >= MAX_ORDER)
        ABORT("Order can not be bigger than " << MAX_ORDER);
    return true;
}

/**
*  subclass may have more parameters to write to _para, so the function does not really save _para to file here
*/
void Parameter::_SavePreparation()
{
    _para.clear();
    SetPara(_para, Version);
    SetPara(_para, L);
    SetPara(_para, InitialBeta);
    SetPara(_para, DeltaBeta);
    SetPara(_para, FinalBeta);
    SetPara(_para, Beta);
    SetPara(_para, Order);
    SetPara(_para, NSublat);
}

bool ParaMC::BuildNew(const std::string& InputFile)
{
    Parameter::_BuildNew(InputFile);
    GetPara(_para, Toss);
    GetPara(_para, Sample);
    GetPara(_para, Sweep);
    GetPara(_para, Seed);
    GetPara(_para, WormSpaceReweight);
    GetPara(_para, OrderReWeight);
    GetPara(_para, MaxTauBin);

    ASSERT_ALLWAYS(OrderReWeight.size() == Order + 1, "OrderReWeight should have Order+1 elementes!");
    Counter = 0;
    this->RNG.Reset(Seed);
    return true;
}

bool ParaMC::Load(const std::string& InputFile)
{
    Parameter::_Load(InputFile);
    GetPara(_para, Counter);
    GetPara(_para, Toss);
    GetPara(_para, Sample);
    GetPara(_para, Sweep);
    GetPara(_para, WormSpaceReweight);
    GetPara(_para, OrderReWeight);
    GetPara(_para, RNG);
    GetPara(_para, MaxTauBin);
    return true;
}

void ParaMC::Save(const std::string& OutputFile, string Mode)
{
    Parameter::_SavePreparation();
    SetPara(_para, Counter);
    SetPara(_para, Toss);
    SetPara(_para, Sample);
    SetPara(_para, Sweep);
    SetPara(_para, WormSpaceReweight);
    SetPara(_para, OrderReWeight);
    SetPara(_para, RNG);
    SetPara(_para, MaxTauBin);
    ASSERT_ALLWAYS(OrderReWeight.size() == Order + 1, "OrderReWeight should have Order+1 elementes!");
    _para.SaveToFile(OutputFile, Mode);
    //save with append mode, so that it will not overwrite stuff wroten by Parameter:SaveParameter
}

void ParaMC::SetTest()
{
    Version = 0;
    int size[2] = { 8, 8 };
    NSublat = 2;
    L = Vec<int>(size);
    Lat = Lattice(L, NSublat);
    InitialBeta = 1.0;
    DeltaBeta = 0.0;
    FinalBeta = 1.0;
    Beta = 1.0;
    Order = 3;
    Toss = 10000;
    Sample = 5000000;
    Seed = 519180543;
    WormSpaceReweight = 0.1;
    OrderReWeight = { 1, 1, 1, 1 };
    T = 1.0 / Beta;
    Counter = 0;
    MaxTauBin = 32;
}
