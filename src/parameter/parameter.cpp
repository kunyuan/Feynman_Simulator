//
//  state.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "parameter.h"

bool Parameter::BuildNew(const std::string &InputFile)
{
    _para.LoadFromFile(InputFile);
    GetPara(_para, L);
    GetPara(_para, Jcp);
    GetPara(_para, InitialBeta);
    GetPara(_para, DeltaBeta);
    GetPara(_para, FinalBeta);
    GetPara(_para, Order);

    Lat.Reset(L);
    Beta = InitialBeta;
    T = 1.0 / Beta;
    if (Order >= MAX_ORDER)
        ABORT("Order can not be bigger than " << MAX_ORDER);
    return true;
}

bool Parameter::Load(const std::string &InputFile)
{
    _para.LoadFromFile(InputFile);
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

void Parameter::Save(const std::string &OutputFile, string Mode)
{
    _para.clear();
    SetPara(_para, L);
    SetPara(_para, Jcp);
    SetPara(_para, InitialBeta);
    SetPara(_para, DeltaBeta);
    SetPara(_para, FinalBeta);
    //!!!Beta should be a part of state, so it will be stored
    SetPara(_para, Beta);
    SetPara(_para, Order);
    _para.SaveToFile(OutputFile, Mode);
}

bool ParameterMC::BuildNew(const std::string &InputFile)
{
    Parameter::BuildNew(InputFile);
    GetPara(_para, Toss);
    GetPara(_para, Sample);
    GetPara(_para, Sweep);
    GetPara(_para, Seed);
    GetPara(_para, WormSpaceReweight);
    GetParaArray(_para, OrderReWeight, Order);

    Version = 0;
    Counter = 0;
    this->RNG.Reset(Seed);
    return true;
}

bool ParameterMC::Load(const std::string &InputFile)
{
    Parameter::Load(InputFile);
    GetPara(_para, Version);
    GetPara(_para, Counter);
    GetPara(_para, Toss);
    GetPara(_para, Sample);
    GetPara(_para, Sweep);
    GetPara(_para, WormSpaceReweight);
    GetParaArray(_para, OrderReWeight, Order);
    GetPara(_para, RNG);
    return true;
}

void ParameterMC::Save(const std::string &OutputFile, string Mode)
{
    Parameter::Save(OutputFile, Mode);
    _para.clear();
    SetPara(_para, Version);
    SetPara(_para, Counter);
    SetPara(_para, Toss);
    SetPara(_para, Sample);
    SetPara(_para, Sweep);
    SetPara(_para, WormSpaceReweight);
    SetParaArray(_para, OrderReWeight, Order);
    SetPara(_para, RNG);
    _para.SaveToFile(OutputFile, "a");
    //save with append mode, so that it will not overwrite stuff wroten by Parameter:SaveParameter
}