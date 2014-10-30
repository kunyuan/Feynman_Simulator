//
//  measure.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/19/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "markov_monitor.h"

MarkovMonitor::MarkovMonitor(EnvMonteCarlo *env)
{
    _Env = env;
}

void MarkovMonitor::Annealing()
{
    SqueezeStatistics();
}

void MarkovMonitor::SqueezeStatistics()
{
}

void MarkovMonitor::ReWeightEachOrder()
{
}

void MarkovMonitor::Measure()
{
    //    cEstimator[0].Measure(<#const Complex &#>);
    //    _Env->cEstimator["1"].Measure(<#const Complex &#>);
}

void MarkovMonitor::AddStatistics()
{
    _Env->cEstimator.AddStatistics();
    _Env->rEstimator.AddStatistics();
}