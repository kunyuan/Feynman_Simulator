//
//  component_test.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/7/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "diagram.h"
#include "../utility/sput.h"
using namespace std;

void Test_Diagram_Component();
void Test_Diagram_IO();

int TestDiagram()
{
    sput_start_testing();
    sput_enter_suite("Test Underlying Diagram Data Structure");
    sput_run_test(Test_Diagram_Component);
    sput_run_test(Test_Diagram_IO);
    sput_finish_testing();
    return sput_get_return_value();
}

void Test_Diagram_Component()
{
    Bundle<GLine> gbundle("GLine");
    GLine &New1 = gbundle.Add();
    GLine &New2 = gbundle.Add();
    GLine &New3 = gbundle.Add();
    gbundle.Remove(New2);
    GLine &New4 = gbundle.Add();
    sput_fail_unless(&gbundle[0] == &New1 && New1.Name == 0, "Check Bundle[0]");
    sput_fail_unless(&gbundle[1] == &New3 && New3.Name == 1, "Check Bundle[1]");
    sput_fail_unless(&gbundle[2] == &New4 && New4.Name == 2, "Check Bundle[2]");
    sput_fail_unless(gbundle.HowMany() == 3, "Check Bundle counter");
    gbundle.Remove(New3);
    gbundle.Remove(New4);
    gbundle.Recover(2);
    sput_fail_unless(&gbundle[1] == &New4 && New4.Name == 1, "Check Bundle[1] after reversing");
    sput_fail_unless(&gbundle[2] == &New3 && New3.Name == 2, "Check Bundle[2] after reversing");
    bool flag = true;
    for (int i = 0; i < gbundle.HowMany(); i++) {
        flag = (gbundle[i].Name == i);
        if (!flag)
            break;
    }
    sput_fail_unless(flag, "Check all objects have right names");
    sput_fail_unless(!gbundle.Exist(gbundle[3]), "Check existence of first useless object");
    sput_fail_unless(gbundle.Exist(gbundle[2]), "Check existence of last useful object");
    GLine Outside;
    Outside.Name = 0;
    sput_fail_unless(!gbundle.Exist(Outside), "Check the existence of an outside G");
    return;
}

void Test_Diagram_IO()
{
    Diagram Diag;
    stringstream strtemp;
    GLine NewG;
    NewG.Name = 2;
    NewG.Vertex[IN] = 0;
    NewG.Vertex[OUT] = 1;
    NewG.K = 1;
    NewG.Weight = Complex(1.0, 3.0e9);
    NewG.SaveConfig(strtemp);
    NewG.LoadConfig(strtemp);
    sput_fail_unless(Equal(NewG.Vertex[OUT], 1), "Check GLine reading");

    WLine NewW;
    NewW.Vertex[IN] = 2;
    NewW.Vertex[OUT] = 3;
    NewW.K = 2;
    NewW.SaveConfig(strtemp);
    NewW.LoadConfig(strtemp);
    sput_fail_unless(Equal(NewW.Vertex[OUT], 3), "Check WLine reading");

    Vertex NewV;
    NewV.Name = 0;
    NewV.R.Sublattice = 1;
    NewV.R.Coordinate[0] = 12;
    NewV.R.Coordinate[1] = 2;
    NewV.Tau = 0.51;
    NewV.Spin[IN] = UP;
    NewV.Spin[OUT] = DOWN;
    NewV.SaveConfig(strtemp);
    NewV.LoadConfig(strtemp);
    sput_fail_unless(Equal(NewV.Tau, 0.51), "Check Vertex reading");

    Lattice lat;
    Weight::G G(lat, 1.0, 1);
    Weight::W W(lat, 1.0, 1);
    G.InitializeState();
    W.InitializeState();

    Diag.Build(&lat, &G, &W);

    Diag.LoadConfig("../src/diagram/diagram_template.config");
    LOG_INFO(Diag.PrettyString(Diag.Ver[0]));
    sput_fail_unless(Diag.CheckDiagram(), "Check diagram reading");
    sput_fail_unless(Equal(Diag.Weight, Complex(64.0, 0.0)), "Check diagram reading");
    Diag.SaveConfig("diagram_template.config", "w");
    Diag.WriteDiagram2gv("./test.gv");
}
