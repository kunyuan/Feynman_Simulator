//
//  state.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "parameter.h"
using namespace para;

status Parameter::GetStatus()
{
    status Status;
    Status.Version = Version;
    Status.Jcp = Jcp;
    Status.Beta = Beta;
    return Status;
}

void Parameter::SetStatus(const status &Status)
{
    Version = Status.Version;
    Jcp = Status.Jcp;
    Beta = Status.Beta;
    T = 1.0 / Beta;
}

bool Parameter::_BuildNew(const std::string &InputFile)
{
    _para.ParseFile(InputFile);
    GetPara(_para, L);
    GetPara(_para, Jcp);
    GetPara(_para, InitialBeta);
    GetPara(_para, DeltaBeta);
    GetPara(_para, FinalBeta);
    GetPara(_para, Order);

    Lat.Reset(L);
    Version = 0;
    Beta = InitialBeta;
    T = 1.0 / Beta;
    if (Order >= MAX_ORDER)
        ABORT("Order can not be bigger than " << MAX_ORDER);
    return true;
}

bool Parameter::_Load(const std::string &InputFile)
{
    _para.ParseFile(InputFile);
    GetPara(_para, Version);
    GetPara(_para, L);
    GetPara(_para, Jcp);
    GetPara(_para, InitialBeta);
    GetPara(_para, DeltaBeta);
    GetPara(_para, FinalBeta);
    //!!!Beta should be a part of state, so it will be stored
    GetPara(_para, Beta);
    GetPara(_para, Order);

    Lat.Reset(L);
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
    SetPara(_para, Jcp);
    SetPara(_para, InitialBeta);
    SetPara(_para, DeltaBeta);
    SetPara(_para, FinalBeta);
    //!!!Beta should be a part of state, so it will be stored
    SetPara(_para, Beta);
    SetPara(_para, Order);
}

bool ParaMC::BuildNew(const std::string &InputFile)
{
    Parameter::_BuildNew(InputFile);
    GetPara(_para, Toss);
    GetPara(_para, Sample);
    GetPara(_para, Sweep);
    GetPara(_para, Seed);
    GetPara(_para, WormSpaceReweight);
    GetParaArray(_para, OrderReWeight, Order);

    Counter = 0;
    this->RNG.Reset(Seed);
    return true;
}

bool ParaMC::Load(const std::string &InputFile)
{
    Parameter::_Load(InputFile);
    GetPara(_para, Counter);
    GetPara(_para, Toss);
    GetPara(_para, Sample);
    GetPara(_para, Sweep);
    GetPara(_para, WormSpaceReweight);
    GetParaArray(_para, OrderReWeight, Order);
    GetPara(_para, RNG);
    return true;
}

void ParaMC::Save(const std::string &OutputFile, string Mode)
{
    Parameter::_SavePreparation();
    SetPara(_para, Counter);
    SetPara(_para, Toss);
    SetPara(_para, Sample);
    SetPara(_para, Sweep);
    SetPara(_para, WormSpaceReweight);
    SetParaArray(_para, OrderReWeight, Order);
    SetPara(_para, RNG);
    _para.SaveToFile(OutputFile, Mode);
    //save with append mode, so that it will not overwrite stuff wroten by Parameter:SaveParameter
}

void ParaMC::SetTest()
{
    Version = 0;
    int size[2] = {8, 8};
    L = Vec<int>(size);
    Lat.Reset(L);
    Jcp = 1.0;
    InitialBeta = 1.0;
    DeltaBeta = 0.0;
    FinalBeta = 1.0;
    Beta = 1.0;
    Order = 4;
    Toss = 10000;
    Sample = 5000000;
    Seed = 519180543;
    WormSpaceReweight = 0.1;
    OrderReWeight[0] = 1.5;
    OrderReWeight[1] = 1;
    OrderReWeight[2] = 3.0;
    OrderReWeight[3] = 4.0;
    T = 1.0 / Beta;
}

bool ParaDyson::BuildNew(const string &InputFile)
{
    Parameter::_BuildNew(InputFile);
    GetPara(_para, OrderAccepted);
    GetPara(_para, ErrorThreshold);
    GetPara(_para, SleepTime);
    return true;
}

bool ParaDyson::Load(const string &InputFile)
{
    Parameter::_Load(InputFile);
    GetPara(_para, OrderAccepted);
    GetPara(_para, ErrorThreshold);
    GetPara(_para, SleepTime);
    return true;
}

void ParaDyson::Save(const std::string &OutputFile, string Mode)
{
    Parameter::_SavePreparation();
    SetPara(_para, OrderAccepted);
    SetPara(_para, ErrorThreshold);
    SetPara(_para, SleepTime);
    _para.SaveToFile(OutputFile, Mode);
}