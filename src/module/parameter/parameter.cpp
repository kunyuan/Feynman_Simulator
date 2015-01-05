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

bool Parameter::_FromDict(const Dictionary& Para)
{
    GET_WITH_DEFAULT(Para, Version, 0);
    auto _para = Para.Get<Dictionary>("Tau");
    GET(_para, Beta);
    GET(_para, DeltaBeta);
    GET(_para, FinalBeta);
    GET(_para, MaxTauBin);
    _para = Para.Get<Dictionary>("Lattice");
    GET(_para, L);
    GET(_para, NSublat);

    Lat.Initialize(L, NSublat);
    T = 1.0 / Beta;

    ASSERT_ALLWAYS(Beta <= FinalBeta, "Beta should be <= Final Beta");
    return true;
}
Dictionary Parameter::_ToDict()
{
    Dictionary Para, _para;
    SET(_para, L);
    SET(_para, NSublat);
    Para["Lattice"] = _para;
    _para.Clear();
    SET(_para, Beta);
    SET(_para, DeltaBeta);
    SET(_para, FinalBeta);
    SET(_para, MaxTauBin);
    Para["Tau"] = _para;
    SET(Para, Version);
    return Para;
}

bool ParaMC::BuildNew(const std::string& InputFile)
{
    Dictionary _para;
    _para.Load(InputFile);
    FromDict(_para);
    return true;
}

bool ParaMC::FromDict(const Dictionary& Para)
{
    Parameter::_FromDict(Para);
    auto _para = Para.Get<Dictionary>("Markov");
    GET(_para, Toss);
    GET(_para, Sample);
    GET(_para, Sweep);
    GET(_para, WormSpaceReweight);
    GET(_para, OrderReWeight);
    GET(_para, Order);
    GET_WITH_DEFAULT(_para, Counter, 0);
    GET_WITH_DEFAULT(_para, Seed, 0);
    if (_para.HasKey("RNG"))
        GET(_para, RNG);
    else
        RNG.Reset(Seed);
    ASSERT_ALLWAYS(Order < MAX_ORDER, "Order can not be bigger than " << MAX_ORDER);
    ASSERT_ALLWAYS(OrderReWeight.size() == Order + 1, "OrderReWeight should have Order+1 elementes!");
    return true;
}
Dictionary ParaMC::ToDict()
{
    Dictionary _para;
    SET(_para, Toss);
    SET(_para, Sample);
    SET(_para, Sweep);
    SET(_para, WormSpaceReweight);
    SET(_para, OrderReWeight);
    SET(_para, Counter);
    SET(_para, RNG);
    SET(_para, Order);
    Dictionary Para;
    Para["Markov"] = _para;
    Para.Update(Parameter::_ToDict());
    return Para;
}

void ParaMC::SetTest()
{
    Version = 0;
    int size[2] = { 8, 8 };
    NSublat = 2;
    L = Vec<int>(size);
    Lat = Lattice(L, NSublat);
    Beta = 1.0;
    DeltaBeta = 0.0;
    FinalBeta = 1.0;
    Beta = 1.0;
    Order = 1;
    OrderReWeight = { 1, 1 };
    Toss = 10000;
    Sample = 5000000;
    Seed = 519180543;
    WormSpaceReweight = 0.1;
    T = 1.0 / Beta;
    Counter = 0;
    MaxTauBin = 32;
}
