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
#include "module/weight/gamma3.h"
#include "utility/dictionary.h"

using namespace std;
using namespace diag;
using namespace para;
using namespace mc;
using namespace weight;

MarkovMonitor::MarkovMonitor()
{
}

bool MarkovMonitor::BuildNew(ParaMC &para, Diagram &diag, weight::Weight &weight)
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

void MarkovMonitor::Reset(ParaMC &para, Diagram &diag, weight::Weight &weight)
{
    Para = &para;
    Diag = &diag;
    Weight = &weight;
}

bool MarkovMonitor::FromDict(const Dictionary &dict, ParaMC &para, Diagram &diag, weight::Weight &weight)
{
    Para = &para;
    Diag = &diag;
    Weight = &weight;
    for (int i = 0; i <= Para->Order; i++) {
        WormEstimator.AddEstimator("Order" + ToString(i));
        PhyEstimator.AddEstimator("Order" + ToString(i));
    }
    bool flag = true;
    flag &= WormEstimator.FromDict(dict.Get<Dictionary>("WormEstimator"),
                                   true //allow failure
                                   );
    flag &= PhyEstimator.FromDict(dict.Get<Dictionary>("PhyEstimator"),
                                  true //allow failure
                                  );
    flag &= SigmaEstimator.FromDict(dict.Get<Dictionary>("SigmaEstimator"));
    flag &= PolarEstimator.FromDict(dict.Get<Dictionary>("PolarEstimator"));
    return flag;
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
    Weight->GammaG->SqueezeStatistics(factor);
    Weight->GammaW->SqueezeStatistics(factor);
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
    if (PolarEstimator.Norm() < 1000.0)
        return false;
    if (Zero(PolarEstimator.Value()) || Zero(SigmaEstimator.Value()))
        return false;

    real weight[Para->Order + 1];
    real wormweight = 0.0, phyweight = 0.0;
    weight[0] = WormEstimator[0].Value() + PhyEstimator[0].Value();

    LOG_INFO("Number of confs in different Order => Reweight Ratio:");
    stringstream output;
    for (int i = 1; i <= Para->Order; i++) {
        weight[i] = WormEstimator[i].Value() + PhyEstimator[i].Value();
        wormweight += WormEstimator[i].Value();
        phyweight += PhyEstimator[i].Value();
        if (Zero(weight[i]))
            continue;
        Para->OrderReWeight[i] = Para->OrderTimeRatio[i] * weight[0] / weight[i];
        output << "Order0 :" << weight[0] << " Order" << i << " :" << weight[i] << " => " << Para->OrderReWeight[i] << endl;
    }
    Para->OrderReWeight[0] = 1.0;

    if (Zero(wormweight))
        return false;

    Para->WormSpaceReweight = phyweight / wormweight;
    output << "Worm confs:" << wormweight << " Physics confs: " << phyweight << " => " << Para->WormSpaceReweight << endl;
    LOG_INFO(output.str());

    Para->PolarReweight = SigmaEstimator.Value() / PolarEstimator.Value();
    return true;
}

void MarkovMonitor::Measure()
{
    real OrderReWeight = Para->OrderReWeight[Diag->Order];
    if (Diag->Worm.Exist) {
        real WormWeight = 1.0 / OrderReWeight / Para->WormSpaceReweight;
        WormEstimator[Diag->Order].Measure(WormWeight);
        if (Diag->HasGammaGW==0 ) {
            if (Diag->MeasureGLine)
                SigmaEstimator.Measure(WormWeight);
            else if (!Diag->MeasureGLine)
                PolarEstimator.Measure(WormWeight);
        }
    }
    else {
        real OrderWeight = 1.0 / OrderReWeight;

        PhyEstimator[Diag->Order].Measure(OrderWeight);

        if (Diag->Order == 0){
            if (Diag->MeasureGLine) {
                Weight->Sigma->Estimator.MeasureNorm(OrderWeight);
                Weight->GammaG->MeasureNorm(OrderWeight);
            }else{
                Weight->Polar->Estimator.MeasureNorm(OrderWeight);
                Weight->GammaW->MeasureNorm(OrderWeight);
            }
        }

        if (Diag->HasGammaGW==0 ){
            if (Diag->MeasureGLine) {
                SigmaEstimator.Measure(OrderWeight);
                if (Diag->Order != 0) {
                    gLine g = Diag->GMeasure;
                    vertex vin = g->NeighVer(OUT);
                    vertex vout = g->NeighVer(IN);
                    Weight->Sigma->Measure(vin->R, vout->R, vin->Tau, vout->Tau, g->Spin(OUT), g->Spin(IN), Diag->Order, Diag->Phase * OrderWeight);
                }
            }
            else {
                PolarEstimator.Measure(OrderWeight);
                if (Diag->Order != 0) {
                    wLine w = Diag->WMeasure;
                    vertex vin = w->NeighVer(OUT);
                    vertex vout = w->NeighVer(IN);
                    Weight->Polar->Measure(vin->R, vout->R, vin->Tau, vout->Tau, vin->Spin(), vout->Spin(),
                                           Diag->Order, -1.0* Diag->Phase * OrderWeight);
                }
            }
//            if (Diag->MeasureGLine && Diag->Order != 0) {
//                gLine g = Diag->GMeasure;
//                vertex vin = g->NeighVer(OUT);
//                vertex vout = g->NeighVer(IN);
//                ExtPoint& Ext = Diag->UExt;
//
//                Weight->GammaG->Measure(vin->R, vout->R, vin->R, vin->Tau, vout->Tau, 0, g->Spin(OUT), g->Spin(IN), UP,
//                                        Diag->Phase * OrderWeight);
//
////                Weight->GammaG->Measure(vin->R, vout->R, vin->R, vin->Tau, vout->Tau, 0, g->Spin(OUT), g->Spin(IN), UP,
////                                        Diag->Phase * OrderWeight);
//            }
//        }else if(Diag->MeasureGLine && Diag->MeasureGammaGW==1) {
        }else if(Diag->MeasureGLine && Diag->MeasureGammaGW==1 && Diag->HasGammaGW== 2) {
            if (Diag->Order != 0) {
                gLine g = Diag->GMeasure;
                vertex vin = g->NeighVer(OUT);
                vertex vout = g->NeighVer(IN);
                ExtPoint& Ext = Diag->UExt;
                Weight->GammaG->Measure(vin->R, vout->R, Ext.R, vin->Tau, vout->Tau, Ext.Tau, g->Spin(OUT), g->Spin(IN), Ext.Spin,
                                       Diag->Phase * OrderWeight);
            }
        }else if(!Diag->MeasureGLine && Diag->MeasureGammaGW==2 && Diag->HasGammaGW==1) {
//        }else if(!Diag->MeasureGLine && Diag->MeasureGammaGW==0 && Diag->HasGammaGW==1) {
            if (Diag->Order != 0) {
                wLine w = Diag->WMeasure;
                vertex vin = w->NeighVer(OUT);
                vertex vout = w->NeighVer(IN);
                ExtPoint& Ext = Diag->UExt;
                Weight->GammaW->Measure(vin->R, vout->R, Ext.R, vin->Tau, vout->Tau, Ext.Tau, vin->Spin(), vout->Spin(), Ext.Spin,
                                        -1.0* Diag->Phase * OrderWeight);
            }
        }
    }
}

void MarkovMonitor::AddStatistics()
{
}
