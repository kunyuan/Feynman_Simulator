//
//  state.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "parameter.h"
#include "utility.h"
using namespace para;

Message Parameter::GenerateMessage()
{
    Message Message_;
    Message_.Version = Version;
    Message_.Interaction = Interaction;
    Message_.ExternalField = ExternalField;
    Message_.Beta = Beta;
    return Message_;
}

void Parameter::UpdateWithMessage(const Message &Message_)
{
    Version = Message_.Version;
    Interaction = Message_.Interaction;
    ExternalField = Message_.ExternalField;
    Beta = Message_.Beta;
    T = 1.0 / Beta;
}

bool Parameter::_BuildNew(const std::string &InputFile)
{
    _para.ParseFile(InputFile);
    GetPara(_para, Hopping);
    GetPara(_para, Interaction);
    GetPara(_para, ChemicalPotential);
    GetPara(_para, ExternalField);
    GetPara(_para, L);
    GetPara(_para, InitialBeta);
    GetPara(_para, DeltaBeta);
    GetPara(_para, FinalBeta);
    GetPara(_para, Order);

    Lat.Initialize(L, LATTICE);
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
    GetPara(_para, Hopping);
    GetPara(_para, Interaction);
    GetPara(_para, ChemicalPotential);
    GetPara(_para, ExternalField);
    GetPara(_para, L);
    GetPara(_para, InitialBeta);
    GetPara(_para, DeltaBeta);
    GetPara(_para, FinalBeta);
    //!!!Beta should be a part of state, so it will be stored
    GetPara(_para, Beta);
    GetPara(_para, Order);

    Lat.Initialize(L, LATTICE);
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
    SetPara(_para, Hopping);
    SetPara(_para, Interaction);
    SetPara(_para, ChemicalPotential);
    SetPara(_para, ExternalField);
    SetPara(_para, L);
    SetPara(_para, InitialBeta);
    SetPara(_para, DeltaBeta);
    SetPara(_para, FinalBeta);
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
    GetPara(_para, OrderReWeight);

    ASSERT_ALLWAYS(OrderReWeight.size() == Order + 1, "OrderReWeight should have Order+1 elementes!");
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
    GetPara(_para, OrderReWeight);
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
    SetPara(_para, OrderReWeight);
    SetPara(_para, RNG);
    ASSERT_ALLWAYS(OrderReWeight.size() == Order + 1, "OrderReWeight should have Order+1 elementes!");
    _para.SaveToFile(OutputFile, Mode);
    //save with append mode, so that it will not overwrite stuff wroten by Parameter:SaveParameter
}

void ParaMC::SetTest()
{
    Version = 0;
    int size[2] = {8, 8};
    L = Vec<int>(size);
    Lat = Lattice(L, CHECKBOARD);
    Hopping.push_back(0.0);
    Interaction.push_back(1.0);
    ChemicalPotential.push_back(0.0);
    ChemicalPotential.push_back(0.0);
    ExternalField = 0.0;
    InitialBeta = 1.0;
    DeltaBeta = 0.0;
    FinalBeta = 1.0;
    Beta = 1.0;
    Order = 3;
    Toss = 10000;
    Sample = 5000000;
    Seed = 519180543;
    WormSpaceReweight = 0.1;
    OrderReWeight = {1, 1, 1, 1};
    T = 1.0 / Beta;
    Counter = 0;
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

void ParaDyson::SetTest()
{
    Version = 0;
    int size[2] = {8, 8};
    L = Vec<int>(size);
    Lat = Lattice(L, CHECKBOARD);
    Hopping.push_back(0.0);
    Interaction.push_back(1.0);
    ChemicalPotential.push_back(0.0);
    ChemicalPotential.push_back(0.0);
    ExternalField = 0.0;
    InitialBeta = 1.0;
    DeltaBeta = 0.0;
    FinalBeta = 1.0;
    Beta = 1.0;
    Order = 3;
    T = 1.0 / Beta;
    OrderAccepted = 3;
    ErrorThreshold = 0.05;
    SleepTime = 300;
}