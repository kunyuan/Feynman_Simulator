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
    InterruptHandler Interrupt;
    EnvMonteCarlo Env(Job);
    if (Job.DoesLoad)
        Env.Load();
    else
        Env.BuildNew();

    auto& Markov = Env.Markov;
    auto& MarkovMonitor = Env.MarkovMonitor;
    auto& Para = Env.Para;

    LOG_INFO("Markov is started!");
    timer ReweightTimer, PrinterTimer, DiskWriterTimer, MessageTimer;
    PrinterTimer.start();
    DiskWriterTimer.start();
    MessageTimer.start();
    ReweightTimer.start();

    int sigma[MAX_ORDER] = { 0 };
    int polar[MAX_ORDER] = { 0 };

    Env.ListenToMessage();

    for (uint Step = 0; Step < Para.Toss; Step++) {
        Markov.Hop(Para.Sweep);
    }

    //    for (uint i = 0; i < 1000; i++) {
    //        for (uint Step = 0; Step < Job.Sample; Step++) {
    uint Step = 0;
    while (true) {
        //Don't use Para.Counter as counter
        Step++;
        Markov.Hop(Para.Sweep);
        MarkovMonitor.Measure();
        if (!Markov.Diag->Worm.Exist) {
            if (Markov.Diag->MeasureGLine)
                sigma[Markov.Diag->Order]++;
            else
                polar[Markov.Diag->Order]++;
        }

        if (Step % 100 == 0) {
            MarkovMonitor.AddStatistics();
	    //if (Markov.Diag->Order>2) {
		//Markov.Diag->WriteDiagram2gv("diagram/" + ToString(Para.Counter) + ".gv");
	    //}

            if (PrinterTimer.check(Para.PrinterTimer)) {
                Env.Diag.CheckDiagram();
                Markov.PrintDetailBalanceInfo();
            }

            if (DiskWriterTimer.check(Para.DiskWriterTimer)) {
                Interrupt.Delay();
                Env.Save();
                Interrupt.Resume();
            }

            if (MessageTimer.check(Para.MessageTimer))
                Env.ListenToMessage();

            if (ReweightTimer.check(Para.ReweightTimer))
                Env.AdjustOrderReWeight();
        }
    }
    LOG_INFO("Markov is ended!");
}
