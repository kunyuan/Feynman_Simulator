//
//  measure.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/19/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "markov_monitor.h"
#include "module/diagram/diagram.h"
#include "module/parameter/parameter.h"
#include "module/weight/weight.h"
#include "module/weight/component.h"
#include "utility/dictionary.h"

using namespace std;
using namespace diag;
using namespace para;
using namespace mc;

MarkovMonitor::MarkovMonitor()
{
}

bool MarkovMonitor::BuildNew(ParaMC& para, Diagram& diag, weight::Weight& weight)
{
    Para = &para;
    Diag = &diag;
    Weight = &weight;
    for (int i = 0; i <= Para->Order; i++) {
        WormEstimator.AddEstimator("Order" + ToString(i));
        PhyEstimator.AddEstimator("Order" + ToString(i));
    }
    WormEstimator.ClearStatistics();
    PhyEstimator.ClearStatistics();
    SigmaEstimator.ClearStatistics();
    PolarEstimator.ClearStatistics();
    return true;
}

void MarkovMonitor::Reset(ParaMC& para, Diagram& diag, weight::Weight& weight)
{
    Para = &para;
    Diag = &diag;
    Weight = &weight;
}

bool MarkovMonitor::FromDict(const Dictionary& dict, ParaMC& para, Diagram& diag, weight::Weight& weight)
{
    Para = &para;
    Diag = &diag;
    Weight = &weight;
    for (int i = 0; i <= Para->Order; i++) {
        WormEstimator.AddEstimator("Order" + ToString(i));
        PhyEstimator.AddEstimator("Order" + ToString(i));
    }
    return WormEstimator.FromDict(dict.Get<Dictionary>("WormEstimator")) 
        || PhyEstimator.FromDict(dict.Get<Dictionary>("PhyEstimator"))
        ||SigmaEstimator.FromDict(dict.Get<Dictionary>("SigmaEstimator"))
        ||PolarEstimator.FromDict(dict.Get<Dictionary>("PolarEstimator"));
}
Dictionary MarkovMonitor::ToDict()
{
    Dictionary dict;
    dict["WormEstimator"] = WormEstimator.ToDict();
    dict["PhyEstimator"] = PhyEstimator.ToDict();
    dict["SigmaEstimator"] = SigmaEstimator.ToDict();
    dict["PolarEstimator"] = PolarEstimator.ToDict();
    return dict;
}

void MarkovMonitor::SqueezeStatistics(real factor)
{
    Weight->Sigma->Estimator.SqueezeStatistics(factor);
    Weight->Polar->Estimator.SqueezeStatistics(factor);
    WormEstimator.SqueezeStatistics(factor);
    PhyEstimator.SqueezeStatistics(factor);
    SigmaEstimator.SqueezeStatistics(factor);
    PolarEstimator.SqueezeStatistics(factor);
}

bool MarkovMonitor::AdjustOrderReWeight()
{
    for (int i = 0; i <= Para->Order; i++) {
        if (PhyEstimator[i].Norm() < 1000.0)
            return false;
    }
    if(PolarEstimator.Norm() <1000.0)
        return false;
    Para->OrderReWeight[0] = 1.0;
    real weight[Para->Order + 1];
    real wormweight = 0.0, phyweight = 0.0;
    weight[0] = WormEstimator[0].Value() + PhyEstimator[0].Value();
    for (int i = 1; i <= Para->Order; i++) {
        weight[i] = WormEstimator[i].Value() + PhyEstimator[i].Value();
        wormweight += WormEstimator[i].Value();
        phyweight += PhyEstimator[i].Value();

        Para->OrderReWeight[i] = Para->OrderTimeRatio[i] * weight[0] / weight[i];
    }
    Para->WormSpaceReweight = phyweight / wormweight;
    Para->PolarReweight = SigmaEstimator.Value()/PolarEstimator.Value();
    return true;
}

void MarkovMonitor::Measure()
{
    if (Diag->Worm.Exist) {
        WormEstimator[Diag->Order].Measure(1.0 / Para->OrderReWeight[Diag->Order] / Para->WormSpaceReweight);
        if(Diag->MeasureGLine)
            SigmaEstimator.Measure(1.0 / Para->OrderReWeight[Diag->Order] / Para->WormSpaceReweight);
        else
            PolarEstimator.Measure(1.0 / Para->OrderReWeight[Diag->Order] / Para->WormSpaceReweight);
    }
    else {
        PhyEstimator[Diag->Order].Measure(1.0 / Para->OrderReWeight[Diag->Order]);
        if (Diag->MeasureGLine) {
            SigmaEstimator.Measure(1.0/Para->OrderReWeight[Diag->Order]);
            if (Diag->Order == 0) {
                Weight->Sigma->Estimator.MeasureNorm();
            }
            else {
                gLine g = Diag->GMeasure;
                vertex vin = g->NeighVer(OUT);
                vertex vout = g->NeighVer(IN);
                Weight->Sigma->Measure(vin->R, vout->R, vin->Tau, vout->Tau, g->Spin(OUT), g->Spin(IN), Diag->Order, Diag->Phase / Para->OrderReWeight[Diag->Order]);
            }
        }
        else {
            PolarEstimator.Measure(1.0/Para->OrderReWeight[Diag->Order]);
            if (Diag->Order == 0) {
                Weight->Polar->Estimator.MeasureNorm();
            }
            else {
                wLine w = Diag->WMeasure;
                vertex vin = w->NeighVer(OUT);
                vertex vout = w->NeighVer(IN);
                Weight->Polar->Measure(vin->R, vout->R, vin->Tau, vout->Tau, vin->Spin(), vout->Spin(), Diag->Order, -Diag->Phase / Para->OrderReWeight[Diag->Order]);
            }
        }
    }
}

void MarkovMonitor::AddStatistics()
{
}
