//
//  main.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/2/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

/********************** include files *****************************************/
#include <iostream>
#include <unistd.h>
#include "test.h"
#include "environment/environment.h"
#include "job/job.h"
#include <Python/Python.h>
using namespace std;
using namespace para;

void MonteCarlo(const Job&);
int main(int argc, const char* argv[])
{
    Py_Initialize();
    RunTest();
    return 0;
    string InputFile = "infile/_in_MC_2";
    para::Job Job(InputFile);

    if (Job.Type == "MC") {
        MonteCarlo(Job);
    }
    else if (Job.Type == "DYSON") {
        cout << "Not Defined" << endl;
    }
    Py_Finalize();
    return 0;
}

void MonteCarlo(const para::Job& Job)
{
    EnvMonteCarlo PaddyField(Job);
    if (Job.DoesLoad)
        PaddyField.Load();
    else
        PaddyField.BuildNew();

    auto& Para = PaddyField.Para;
    auto& Grasshopper = PaddyField.Grasshopper;
    auto& Scarecrow = PaddyField.Scarecrow;

    while (Para.Counter < 10) {
        Para.Counter++;
        Grasshopper.Hop(0);

        Scarecrow.Measure();

        if (Para.Counter % 10 == 0) {
            //            Env.AddStatistics();
            Scarecrow.ReWeightEachOrder();
        }
        if (Para.Counter % 100000 == 0) {
            PaddyField.Save();
        }
        if (Para.Counter % 20 == 0) {
            PaddyField.ListenToMessage();
        }
    }
}
