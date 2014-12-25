//
//  state.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "parameter.h"
#include "utility/utility.h"
#include "utility/dictionary.h"
using namespace para;

const std::string KEYNAME = "Para";

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
    Dictionary _para;
    _para.Load(InputFile, KEYNAME);
    GET(_para, L);
    GET(_para, InitialBeta);
    GET(_para, DeltaBeta);
    GET(_para, FinalBeta);
    GET(_para, Order);
    GET(_para, NSublat);

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
    Dictionary _para;
    _para.Load(InputFile, KEYNAME);
    GET(_para, Version);
    GET(_para, L);
    GET(_para, InitialBeta);
    GET(_para, DeltaBeta);
    GET(_para, FinalBeta);
    //!!!Beta should be a part of state, so it will be stored
    GET(_para, Beta);
    GET(_para, Order);
    GET(_para, NSublat)

    Lat.Initialize(L, NSublat);
    T = 1.0 / Beta;
    if (Order >= MAX_ORDER)
        ABORT("Order can not be bigger than " << MAX_ORDER);
    return true;
}

bool ParaMC::BuildNew(const std::string& InputFile)
{
    Parameter::_BuildNew(InputFile);
    Dictionary _para;
    _para.Load(InputFile, KEYNAME);
    GET(_para, Toss);
    GET(_para, Sample);
    GET(_para, Sweep);
    GET(_para, Seed);
    GET(_para, WormSpaceReweight);
    GET(_para, OrderReWeight);
    GET(_para, MaxTauBin);

    ASSERT_ALLWAYS(OrderReWeight.size() == Order + 1, "OrderReWeight should have Order+1 elementes!");
    Counter = 0;
    this->RNG.Reset(Seed);
    return true;
}

bool ParaMC::Load(const std::string& InputFile)
{
    Parameter::_Load(InputFile);
    Dictionary _para;
    _para.Load(InputFile, KEYNAME);
    GET(_para, Counter);
    GET(_para, Toss);
    GET(_para, Sample);
    GET(_para, Sweep);
    GET(_para, WormSpaceReweight);
    GET(_para, OrderReWeight);
    GET(_para, MaxTauBin);
    GET(_para, RNG);
    return true;
}

void ParaMC::Save(const std::string& OutputFile, string Mode)
{
    Dictionary _para;
    SET(_para, Version);
    SET(_para, L);
    SET(_para, InitialBeta);
    SET(_para, DeltaBeta);
    SET(_para, FinalBeta);
    SET(_para, Beta);
    SET(_para, Order);
    SET(_para, NSublat);
    SET(_para, Counter);
    SET(_para, Toss);
    SET(_para, Sample);
    SET(_para, Sweep);
    SET(_para, WormSpaceReweight);
    SET(_para, OrderReWeight);
    SET(_para, MaxTauBin);
    SET(_para, RNG);
    ASSERT_ALLWAYS(OrderReWeight.size() == Order + 1, "OrderReWeight should have Order+1 elementes!");
    _para.Save(OutputFile, Mode, KEYNAME);
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
