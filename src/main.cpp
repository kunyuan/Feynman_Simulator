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
#include "utility/pyglue/pywrapper.h"

using namespace std;
using namespace para;

const string HelpStr = "Usage:"
                       "-p N / --PID N   use N to construct input file path."
                       "or -f / --file PATH   use PATH as the input file path.";

void MonteCarlo(const Job&);
int main(int argc, const char* argv[])
{
    Python::Initialize();
    Python::ArrayInitialize();
    RunTest();
    ASSERT_ALLWAYS(argc == 3, HelpStr);
    string InputFile;
    if (strcmp(argv[1], "-p") == 0 || strcmp(argv[1], "--PID") == 0)
        InputFile = string("infile/_in_MC_") + argv[2];
    else if (strcmp(argv[1], "-f") == 0 || strcmp(argv[1], "--file") == 0)
        InputFile = argv[2];
    else
        ABORT("Unable to parse arguments!\n" + HelpStr);

    para::Job Job(InputFile);

    if (Job.Type == "MC")
        MonteCarlo(Job);
    else
        cout << "Not Defined" << endl;
    Python::Finalize();
    return 0;
}

void MonteCarlo(const para::Job& Job)
{
    EnvMonteCarlo PaddyField(Job);
    if (Job.DoesLoad)
        PaddyField.Load();
    else
        PaddyField.BuildNew();

    auto& Grasshopper = PaddyField.Grasshopper;
    auto& Scarecrow = PaddyField.Scarecrow;
    auto& Para = PaddyField.Para;

    int icount = 0;
    //Don't use Para.Counter as counter

    while (icount < Para.Sample) {
        icount++;
        Grasshopper.Hop(Para.Sweep);

        Scarecrow.Measure();
        
        if(icount %10==0)
            Scarecrow.AddStatistics();

        if (icount % 1000 == 0) {
            PaddyField.Diag.CheckDiagram();
            if (!PaddyField.Diag.Worm.Exist)
                PaddyField.Diag.WriteDiagram2gv("diagram/" + ToString(Para.Counter) + ".gv");

            Scarecrow.ReWeightEachOrder();
        }
        if (icount % 10000 == 0) {
            PaddyField.Save();
            Grasshopper.PrintDetailBalanceInfo();
        }
        if (icount % 200 == 0) {
            PaddyField.ListenToMessage();
        }
    }
}
