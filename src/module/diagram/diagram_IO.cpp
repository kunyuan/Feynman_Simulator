//
//  diagram_IO_global.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/9/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//
//TODO: G, W IsDelta, IsMeasure

#include <vector>
#include <iostream>
#include "diagram.h"
#include "utility/dictionary.h"
#include "utility/abort.h"
#include "utility/scopeguard.h"
using namespace std;
using namespace diag;

/*******************  Read/write diagram to dat file ****************/
Dictionary Diagram::_ToDict(WormClass worm)
{
    Dictionary WormDict;
    WormDict["Ira"] = worm.Ira->Name;
    WormDict["Masha"] = worm.Masha->Name;
    WormDict["dSpin"] = worm.dSpin;
    WormDict["K"] = worm.K;
    return WormDict;
}
void Diagram::_FromDict(const Dictionary& WormDict, WormClass& worm)
{
    name ira, masha;
    WormDict.Get("Ira", ira);
    worm.Ira = Ver(ira);
    WormDict.Get("Masha", masha);
    worm.Masha = Ver(masha);
    WormDict.Get("dSpin", worm.dSpin);
    WormDict.Get("K", worm.K);
}
Dictionary Diagram::_ToDict(gLine g)
{
    Dictionary GDict;
    GDict["IN"] = g->nVer[IN]->Name;
    GDict["OUT"] = g->nVer[OUT]->Name;
    GDict["K"] = g->K;
    GDict["IsMeasure"] = g->IsMeasure;
    return GDict;
}
void Diagram::_FromDict(const Dictionary& GDict, gLine g)
{
    name g_in, g_out;
    GDict.Get("IN", g_in);
    g->nVer[IN] = Ver(g_in);
    GDict.Get("OUT", g_out);
    g->nVer[OUT] = Ver(g_out);
    GDict.Get("K", g->K);
    AddGHash(g->K);
    GDict.Get("IsMeasure", g->IsMeasure);
    if (g->IsMeasure) {
        MeasureGLine = true;
        GMeasure = g;
        WMeasure = nullptr;
    }
}

Dictionary Diagram::_ToDict(wLine w)
{
    Dictionary WDict;
    WDict["IN"] = w->nVer[IN]->Name;
    WDict["OUT"] = w->nVer[OUT]->Name;
    WDict["K"] = w->K;
    WDict["IsDelta"] = w->IsDelta;
    WDict["IsMeasure"] = w->IsMeasure;
    return WDict;
}
void Diagram::_FromDict(const Dictionary& WDict, wLine w)
{
    name w_in, w_out;
    WDict.Get("IN", w_in);
    w->nVer[IN] = Ver(w_in);
    WDict.Get("OUT", w_out);
    w->nVer[OUT] = Ver(w_out);
    WDict.Get("K", w->K);
    AddWHash(w->K);
    WDict.Print();
    WDict.Get("IsDelta", w->IsDelta);
    WDict.Get("IsMeasure", w->IsMeasure);
    if (w->IsMeasure) {
        MeasureGLine = false;
        GMeasure = nullptr;
        WMeasure = w;
    }
}
Dictionary Diagram::_ToDict(vertex v)
{
    Dictionary VerDict;
    VerDict["Name"] = v->Name;
    VerDict["Sublat"] = v->R.Sublattice;
    VerDict["Coordi"] = v->R.Coordinate;
    VerDict["Tau"] = v->Tau;
    VerDict["SpinIn"] = (int)v->Spin(IN);
    VerDict["SpinOut"] = (int)v->Spin(OUT);
    return VerDict;
}
void Diagram::_FromDict(const Dictionary& VerDict, vertex v)
{
    VerDict.Get("Name", v->Name);
    VerDict.Get("Sublat", v->R.Sublattice);
    VerDict.Get("Coordi", v->R.Coordinate);
    VerDict.Get("Tau", v->Tau);
    int spinin, spinout;
    VerDict.Get("SpinIn", spinin);
    VerDict.Get("SpinOut", spinout);
    v->_spin[IN] = spin(spinin);
    v->_spin[OUT] = spin(spinout);
}

Dictionary Diagram::ToDict()
{
    Dictionary Config;
    vector<Dictionary> VerList, GList, WList;
    for (int index = 0; index < Ver.HowMany(); index++)
        VerList.push_back(_ToDict(Ver(index)));
    Config["Ver"] = VerList;
    for (int index = 0; index < G.HowMany(); index++)
        GList.push_back(_ToDict(G(index)));
    Config["G"] = GList;
    for (int index = 0; index < W.HowMany(); index++)
        WList.push_back(_ToDict(W(index)));
    Config["W"] = WList;
    if (Worm.Exist) {
        Config["Worm"] = _ToDict(Worm);
    }
    Config["SignFermiLoop"] = SignFermiLoop;
    return Config;
}

bool Diagram::FromDict(const Dictionary& dict, Lattice& lat, weight::G& g, weight::W& w)
{
    Reset(lat, g, w);
    return FromDict(dict);
}

bool Diagram::FromDict(const Dictionary& Config)
{
    ClearDiagram();
    for (auto& dict : Config.Get<vector<Dictionary> >("Ver"))
        _FromDict(dict, Ver.Add());
    for (auto& dict : Config.Get<vector<Dictionary> >("W"))
        _FromDict(dict, W.Add());
    for (auto& dict : Config.Get<vector<Dictionary> >("G"))
        _FromDict(dict, G.Add());
    Worm.Exist = Config.HasKey("Worm");
    if (Config.HasKey("Worm")){
        _FromDict(Config.Get<Dictionary>("Worm"), Worm);
    }
    SignFermiLoop = Config.Get<real>("SignFermiLoop");
    FixDiagram();
    return true;
}

void Diagram::BuildNew(Lattice& lat, weight::G& g, weight::W& w)
{
    Reset(lat, g, w);
    Dictionary Config;
    string coord;
    if (D == 2)
        coord = "[1,0]";
    else if (D == 3)
        coord = "[1,0,0]";
    else
        ABORT("Not implemented!");
    Config.LoadFromString(
        "{'SignFermiLoop': 1.0,"
        "'Ver': "
        "[{'Name': 0, 'Sublat': 0, 'Coordi':"+coord+", 'Tau': 0.0, 'SpinIn': 1, 'SpinOut' :1},"
        "{'Name': 1, 'Sublat': 0, 'Coordi':"+coord+", 'Tau': 0.0, 'SpinIn': 1, 'SpinOut' :1}],"
        "'G':"
        "[{'IN': 0, 'OUT': 1, 'K': 1, 'IsMeasure': True},"
        "{'IN': 1, 'OUT': 0, 'K': 2, 'IsMeasure': False}],"
        "'W':"
        "[{'IN': 0, 'OUT': 1, 'K': 1, 'IsDelta': False, 'IsMeasure': False}]}");
    if (!FromDict(Config))
        ABORT("Faile to construct diagram!");
}

void Diagram::SetTest(Lattice& lat, weight::G& g, weight::W& w)
{
    Reset(lat, g, w);
    Dictionary Config;
    string coord;
    if (D == 2)
        coord = "[1,0]";
    else if (D == 3)
        coord = "[1,0,0]";
    else
        ABORT("Not implemented!");
    Config.LoadFromString(
        "{'SignFermiLoop': 1.0,"
        "'Ver': "
        "[{'Name': 0, 'Sublat': 0, 'Coordi':"+coord+", 'Tau': 0.0, 'SpinIn': 1, 'SpinOut' :1},"
        "{'Name': 1, 'Sublat': 0, 'Coordi':"+coord+", 'Tau': 0.0, 'SpinIn': 1, 'SpinOut' :1}],"
        "'G':"
        "[{'IN': 0, 'OUT': 1, 'K': 1, 'IsMeasure': True},"
        "{'IN': 1, 'OUT': 0, 'K': 2, 'IsMeasure': False}],"
        "'W':"
        "[{'IN': 0, 'OUT': 1, 'K': 1, 'IsDelta': False, 'IsMeasure': False}]}");
    if (!FromDict(Config))
        ABORT("Faile to construct diagram!");
}

/************************   write component to gv ****************************************/
string GLineStyle(bool IsMeasure, spin in, spin out)
{
    string color = "";
    if (IsMeasure)
        color = "color=\"green\"";
    else {
        string str[2] = { "blue", "red" };
        color = "color=\"" + str[in] + ":" + str[out] + ";0.5\"";
    }
    return "[" + color + "]";
}

string WLineStyle(bool IsMeasure)
{
    string color = "";
    if (IsMeasure)
        color = "color=green,";
    return "[" + color + "style=dashed arrowhead=none]";
}

string VertexStyle(bool IsWorm, int sublattice)
{
    string shape = "";
    if (IsWorm)
        shape = "shape=square,";
    string colorstr;
    if (sublattice == 0)
        colorstr = "grey";
    else if (sublattice == 1)
        colorstr = "palegreen";
    else
        colorstr = "palegreen";
    return "[" + shape + "fillcolor=" + colorstr + "]";
}

/**
*  write diagram object to .gv file so that it can be visualized by Graphviz
*
*  @param os ostream& as the target stream to output
*/
void Diagram::WriteDiagram2gv(string path)
{
    ofstream os(path, ios::out);
    if (!os) {
        ABORT("Cannot open " + path);
    }
    string head = "digraph Feynman{\n";
    string tail = "}\n";
    string node_attribute = "    node [margin=0.1 fillcolor=grey fontcolor=black fontsize=10 width=0.2 shape=circle style=filled fixedsize=true]\n";

    //Prettystring
    os << "//Order=" << Order << ", Weight=" << Weight
       << ", FermiLoop=" << SignFermiLoop << ", WormExist=" << Worm.Exist << endl;
    os << "//" << Worm.PrettyString() << endl;
    for (int i = 0; i < Ver.HowMany(); i++)
        os << "//" << Ver(i)->PrettyString() << endl;
    for (int i = 0; i < G.HowMany(); i++)
        os << "//" << G(i)->PrettyString() << endl;
    for (int i = 0; i < W.HowMany(); i++)
        os << "//" << W(i)->PrettyString() << endl;
    os << endl;

    //gv file
    os << head << node_attribute << endl;
    os << "    //" << Ver.BundleName() << endl;
    for (int index = 0; index < Ver.HowMany(); index++) {
        os << "    " << Ver(index)->Name << " "
           << VertexStyle(IsWorm(Ver(index)), Ver(index)->R.Sublattice) << ";" << endl;
    }
    os << "    //" << G.BundleName() << endl;
    for (int i = 0; i < G.HowMany(); i++) {
        os << "    " << G(i)->nVer[IN]->Name << "->" << G(i)->nVer[OUT]->Name;
        os << " " << GLineStyle(G(i)->IsMeasure, G(i)->Spin(IN), G(i)->Spin(OUT))
           << ";" << endl;
    }
    os << "    //" << W.BundleName() << endl;
    for (int i = 0; i < W.HowMany(); i++) {
        os << "    " << W(i)->nVer[IN]->Name << "->" << W(i)->nVer[OUT]->Name;
        os << " " << WLineStyle(W(i)->IsMeasure) << ";" << endl;
    }
    os << tail;
}
