//
//  initialization.h
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Fermion_Simulator__initialization__
#define __Fermion_Simulator__initialization__

#include "definition.h"
#include "array.h"
#include "logger.h"
#include "rng.h"
#include "cnpy.h"
#include "timer.h"
#include "test.h"
#include "environment.h"
#include "lattice.h"
#include "job.h"
#include "diagram.h"
#include "weight.h"

void InitGlobalUtility();

void InitEveryOneNeedsIt();
void InitMonteCarlo();
void InitDyson();
// define all objects to start the simulation here

#endif /* defined(__Fermion_Simulator__initialization__) */
