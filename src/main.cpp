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
#include "utility/pyglue/pywrapper.h"
#include "job/job.h"
#include "utility/timer.h"

using namespace std;
using namespace para;

const string HelpStr = "Usage:"
                       "-p N / --PID N   use N to construct input file path."
                       "or -f / --file PATH   use PATH as the input file path.";

void MonteCarlo(const Job &);
int main(int argc, const char *argv[])
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

void MonteCarlo(const para::Job &Job)
{
    EnvMonteCarlo Env(Job);
    if (Job.DoesLoad)
        Env.Load();
    else
        Env.BuildNew();

    auto &Markov = Env.Markov;
    auto &MarkovMonitor = Env.MarkovMonitor;
    auto &Para = Env.Para;

    LOG_INFO("Markov is started!");
    timer PrinterTimer, DiskWriterTimer, MessageTimer;
    PrinterTimer.start();
    DiskWriterTimer.start();
    MessageTimer.start();

    int sigma[MAX_ORDER] = {0};
    int polar[MAX_ORDER] = {0};

    Env.ListenToMessage();
    Env.Diag.WriteDiagram2gv("diagram/" + ToString(Para.Counter) + ".gv");
    Env.Diag.CheckDiagram();

    for (uint Step = 0; Step < Job.Sample; Step++) {
        //Don't use Para.Counter as counter
        Markov.Hop(Para.Sweep);
        MarkovMonitor.Measure();
        if (!Markov.Diag->Worm.Exist) {
            if (Markov.Diag->MeasureGLine)
                sigma[Markov.Diag->Order]++;
            else
                polar[Markov.Diag->Order]++;
        }

        if (Step % 1000 == 0) {
            MarkovMonitor.AddStatistics();
            //            if (!Env.Diag.Worm.Exist && Env.Diag.Order == 3)
            //                Env.Diag.WriteDiagram2gv("diagram/" + ToString(Para.Counter) + ".gv");
            if (PrinterTimer.check(10)) {
                Env.Diag.CheckDiagram();
                Markov.PrintDetailBalanceInfo();
            }
            if (DiskWriterTimer.check(60)) {
                Env.AdjustOrderReWeight();
                Env.Save();
            }
            if (MessageTimer.check(10))
                Env.ListenToMessage();
        }
    }

    //    Markov.PrintDetailBalanceInfo();
    Env.Save();
    LOG_INFO("Markov is ended!");
}
