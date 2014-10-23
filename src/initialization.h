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
#include "test.h"
#include "utility/array.h"
#include "utility/logger.h"
#include "utility/rng.h"
#include "utility/cnpy.h"
#include "utility/timer.h"
#include "environment/environment.h"
#include "lattice/lattice.h"
#include "diagram/diagram.h"
#include "observable/weight.h"
#include "job/job.h"

void InitGlobalUtility();

void InitEveryOneNeedsIt();
void InitMonteCarlo();
void InitDyson();
// define all objects to start the simulation here

#endif /* defined(__Fermion_Simulator__initialization__) */
