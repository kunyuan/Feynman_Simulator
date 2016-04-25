//
//  component_test.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/7/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "diagram.h"
#include "utility/sput.h"
#include "module/weight/component.h"
using namespace std;
using namespace diag;

void Test_Diagram_Component();
void Test_Diagram_Component_Bundle();
void Test_Diagram_IO();

int diag::TestDiagram()
{
    sput_start_testing();
    sput_enter_suite("Test Underlying Diagram Data Structure");
    sput_run_test(Test_Diagram_Component);
    sput_run_test(Test_Diagram_Component_Bundle);
    sput_run_test(Test_Diagram_IO);
    sput_finish_testing();
    return sput_get_return_value();
}

bool CheckNameInBundle(Bundle<GLine>& b)
{
    bool flag = true;
    for (int i = 0; i < b.HowMany(); i++) {
        flag = (b[i].Name == i);
        if (!flag)
            break;
    }
    return flag;
}

void Test_Diagram_Component()
{
    Bundle<GLine> gbundle("GLine");
    gLine New1 = gbundle.Add();
    gLine New2 = gbundle.Add();
    gLine New3 = gbundle.Add();
    //New1, New2, New3,...
    gbundle.Remove(New2);
    //New1, New3,...
    gLine New4 = gbundle.Add();
    //New1, New3, New4...
    sput_fail_unless(gbundle(0) == New1 && New1->Name == 0, "Check Bundle[0]");
    sput_fail_unless(gbundle(1) == New3 && New3->Name == 1, "Check Bundle[1]");
    sput_fail_unless(gbundle(2) == New4 && New4->Name == 2, "Check Bundle[2]");
    sput_fail_unless(gbundle.HowMany() == 3, "Check Bundle counter");
    gbundle.Remove(New3);
    gbundle.Remove(New4);
    //New1...
    gbundle.Recover(2);
    //New1, New3, New4...
    sput_fail_unless(gbundle(1) == New4 && New4->Name == 1, "Check Bundle[1] after reversing");
    sput_fail_unless(gbundle(2) == New3 && New3->Name == 2, "Check Bundle[2] after reversing");
    sput_fail_unless(CheckNameInBundle(gbundle), "Check all objects have right names");

    sput_fail_unless(gbundle.Exist(gbundle(2)), "Check existence of last useful object");
    return;
}

void Test_Diagram_Component_Bundle()
{
    Bundle<GLine> gbundle("GLine");
    for (int i = 0; i < MAX_BUNDLE; i++)
        gbundle.Add();
    sput_fail_unless(gbundle.HowMany() == MAX_BUNDLE, "Check maximum capacity");
    for (int i = 0; i < MAX_BUNDLE; i++)
        gbundle.Remove(0);
    sput_fail_unless(gbundle.HowMany() == 0, "Check maximum capacity");
    sput_fail_unless(CheckNameInBundle(gbundle), "Check all objects have right names");
}

void Test_Diagram_IO()
{
    Lattice lat(Vec<int>(8));
    weight::G G(lat, 1.0, 32);
    weight::W W(lat, 1.0, 32);
    G.BuildTest();
    W.BuildTest();
    Diagram Diag;

    Diag.BuildNew(lat, G, W);
    LOG_INFO(Diag.Ver(0)->PrettyString());
    LOG_INFO(Diag.Ver(1)->PrettyString());
    sput_fail_unless(Diag.CheckDiagram(), "Check diagram G,W,Ver and Weight");
    Diag.WriteDiagram2gv("./test.gv");
    system("rm ./test.gv");
}
