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

const char SEP = ' ';
const char COMMENT = '#';
const string SEP_LINE = "#######################################################";
const string SEP_LINE_SHORT = "####";

/*******************  Read/write diagram to dat file ****************/
Dictionary GetDict(WormClass worm)
{
    Dictionary WormDict("Worm");
    WormDict.Set("Ira", worm.Ira->Name);
    WormDict.Set("Masha", worm.Masha->Name);
    WormDict.Set("dSpin", (int)worm.dSpin);
    WormDict.Set("K", worm.K.K);
    return WormDict;
}
Dictionary GetDict(gLine g)
{
    Dictionary GDict("G");
    GDict.Set("IN", g->nVer[IN]->Name);
    GDict.Set("OUT", g->nVer[OUT]->Name);
    GDict.Set("K", g->K.K);
    GDict.Set("IsMeasure", g->IsMeasure);
    return GDict;
}
Dictionary GetDict(wLine w)
{
    Dictionary WDict("W");
    WDict.Set("IN", w->nVer[IN]->Name);
    WDict.Set("OUT", w->nVer[OUT]->Name);
    WDict.Set("K", w->K.K);
    WDict.Set("IsDelta", w->IsDelta);
    WDict.Set("IsMeasure", w->IsMeasure);
    return WDict;
}
Dictionary GetDict(vertex v)
{
    Dictionary VerDict("Ver");
    VerDict.Set("Name", v->Name);
    VerDict.Set("Sublat", v->R.Sublattice);
    VerDict.Set("Coordi", v->R.Coordinate);
    VerDict.Set("Tau", v->Tau);
    VerDict.Set("SpinIn", (int)v->Spin(IN));
    VerDict.Set("SpinOut", (int)v->Spin(OUT));
    return VerDict;
}
void Diagram::Save(const std::string& FileName, string Mode)
{
    Dictionary Config("Config");
    for (int index = 0; index < Ver.HowMany(); index++) {
        Config.Set("Ver", GetDict(Ver(index)));
    }
    for (int index = 0; index < G.HowMany(); index++) {
        Config.Set("G", GetDict(G(index)));
    }
    for (int index = 0; index < W.HowMany(); index++) {
        Config.Set("W", GetDict(W(index)));
    }
    if (Worm.Exist) {
        Config.Set("Worm", GetDict(Worm));
    }
}

bool Diagram::_Load(istream& ifs)
{
    string line;
    //locate the last configration block
    streampos lastBlockPos = 0;
    while (getline(ifs, line)) {
        if (line.compare(0, SEP_LINE_SHORT.size(), SEP_LINE_SHORT) == 0)
            lastBlockPos = ifs.tellg();
    }
    ifs.clear();
    ifs.seekg(lastBlockPos);

    ClearDiagram();
    char head;
    string temp;
    while (!ifs.eof()) {
        ifs >> head;
        if (head == COMMENT) {
            getline(ifs, temp);
            continue;
        }
        else if (head == 'g')
            LoadConfig(ifs, G.Add());
        else if (head == 'w')
            LoadConfig(ifs, W.Add());
        else if (head == 'v')
            LoadConfig(ifs, Ver.Add());
        else if (head == 'i')
            LoadConfig(ifs, Worm);
        else if (head == 's')
            ifs >> SignFermiLoop;
        else
            ABORT("Error in reading diagram! Get " + ToString(head) + " as the head!");
        head = COMMENT;
    }
    FixDiagram();
    return true;
}

bool Diagram::Load(const std::string& FileName)
{

    ifstream ifs(FileName, ios::in);
    ON_SCOPE_EXIT([&] {ifs.close(); });
    if (!ifs.is_open()) {
        ABORT("Cannot open " + FileName);
        return false;
    }
    return _Load(ifs);
}

bool Diagram::Load(const std::string& FileName, Lattice& lat,
                   weight::G& g, weight::W& w)
{
    Reset(lat, g, w);
    return Load(FileName);
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
