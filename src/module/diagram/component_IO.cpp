//
//  component.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/11/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

//TODO: G, W IsDelta, IsMeasure
#include "diagram.h"
#include "utility/abort.h"

#define SEP ' '
using namespace std;
using namespace diag;

#define READ(is, thing)                              \
    {                                                \
        is >> (thing);                               \
        if (is.fail())                               \
            ABORT("Fail to read " << #thing << "!"); \
    }

/*******************   Read/write component to dat file  ********************************/
//bool Diagram::LoadConfig(istream &is, WormClass &worm)
//{
//    //format: i/Ira/Masha/dSpin/K
//    name ira, masha;
//    READ(is, ira);
//    worm.Ira = Ver(ira);
//    READ(is, masha);
//    worm.Masha = Ver(masha);
//    READ(is, worm.dSpin);
//    READ(is, worm.K);
//    return true;
//}
//
//void Diagram::SaveConfig(ostream &os, WormClass &worm)
//{
//    os << worm.Ira->Name << SEP << worm.Masha->Name << SEP << worm.dSpin << SEP << worm.K << endl;
//}
//
//bool Diagram::LoadConfig(istream &is, gLine g)
//{
//    //format: g/start/end/K
//    name g_in, g_out;
//    READ(is, g_in);
//    g->nVer[IN] = Ver(g_in);
//    READ(is, g_out);
//    g->nVer[OUT] = Ver(g_out);
//    READ(is, g->K);
//    AddGHash(g->K);
//    READ(is, g->IsMeasure);
//    if(g->IsMeasure)
//    {
//        MeasureGLine = true;
//        GMeasure = g;
//        WMeasure = nullptr;
//    }
//    return true;
//}
//
//void Diagram::SaveConfig(ostream &os, gLine g)
//{
//    os << g->nVer[IN]->Name << SEP << g->nVer[OUT]->Name << SEP << g->K << SEP << g->IsMeasure << endl;
//}
//
//bool Diagram::LoadConfig(istream &is, wLine w)
//{
//    //format: w/start/end
//    name w_in, w_out;
//    READ(is, w_in);
//    w->nVer[IN] = Ver(w_in);
//    READ(is, w_out);
//    w->nVer[OUT] = Ver(w_out);
//    READ(is, w->K);
//    AddWHash(w->K);
//    READ(is, w->IsDelta);
//    READ(is, w->IsMeasure);
//    if(w->IsMeasure)
//    {
//        MeasureGLine = false;
//        GMeasure = nullptr;
//        WMeasure = w;
//    }
//    return true;
//}
//
//void Diagram::SaveConfig(ostream &os, wLine w)
//{
//    os << w->nVer[IN]->Name << SEP << w->nVer[OUT]->Name << SEP << w->K << SEP << w->IsDelta << SEP << w->IsMeasure << endl;
//}
//
//bool Diagram::LoadConfig(istream &is, vertex v)
//{
//    int spinin, spinout;
//    //format: name sublattice r tau
//    READ(is, v->Name);
//    READ(is, v->R.Sublattice);
//    READ(is, v->R.Coordinate);
//    READ(is, v->Tau);
//    READ(is, spinin);
//    READ(is, spinout);
//    if (is.good()) {
//        v->_spin[IN] = spin(spinin);
//        v->_spin[OUT] = spin(spinout);
//    }
//    return true;
//}
//
//void Diagram::SaveConfig(ostream &os, vertex v)
//{
//    os << v->Name << SEP << v->R.Sublattice << SEP << v->R.Coordinate << SEP << v->Tau << SEP << int(v->Spin(IN)) << SEP << int(v->Spin(OUT)) << endl;
//}
