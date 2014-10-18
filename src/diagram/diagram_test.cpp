//
//  component_test.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/7/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include <sstream>
#include "sput.h"
#include "component_bundle.h"
#include "utility.h"
#include "diagram.h"
#include "logger.h"
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
    NewG.Weight = Complex(1.0, 3.0e9);
    LOG_INFO(NewG);
    strtemp << NewG;

    WLine NewW;
    NewW.Vertex[IN] = 2;
    NewW.Vertex[OUT] = 3;
    LOG_INFO(NewW);
    strtemp << NewW;

    Vertex NewV;
    NewV.Name = 0;
    NewV.R.Sublattice = 1;
    NewV.R.Coordinate[0] = 12;
    NewV.R.Coordinate[1] = 2;
    NewV.tau = 0.51;
    NewV.Spin[IN]=UP;
    NewV.Spin[OUT]=DOWN;
    LOG_INFO(NewV);
    strtemp << NewV;

    strtemp >> NewG;
    sput_fail_unless(Equal(Diag.Spin(NewG,OUT), DOWN), "Check GLine reading");
    strtemp >> NewW;
    sput_fail_unless(Equal(NewW.Vertex[OUT], 3), "Check WLine reading");
    strtemp >> NewV;
    sput_fail_unless(Equal(NewV.tau, 0.51), "Check Vertex reading");

    Diag.ReadDiagram("../src/diagram/diagram_template.config");
    sput_fail_unless(Diag.CheckDiagram(), "Check diagram reading");
    Diag.WriteDiagram(cout);
    Diag.WriteDiagram2gv("./test.gv");
}