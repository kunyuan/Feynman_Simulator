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
#include "utility/abort.h"
#include "utility/scopeguard.h"
using namespace std;
using namespace diag;

const char SEP = ' ';
const char COMMENT = '#';
const string SEP_LINE = "#######################################################";
const string SEP_LINE_SHORT = "####";

/*******************  Read/write diagram to dat file ****************/
void Diagram::Save(const std::string &FileName, string Mode)
{
    ofstream os;
    if (Mode == "w")
        os.open(FileName, ios::out);
    else if (Mode == "a")
        os.open(FileName, ios::app);
    else
        ABORT("What is Mode=" << Mode << "?");
    ON_SCOPE_EXIT([&] {os.close(); });
    if (!os.is_open()) {
        ABORT("Cannot open " + FileName);
    }
    else {
        os << SEP_LINE << endl;
        os << COMMENT << Ver.BundleName() << endl;
        for (int index = 0; index < Ver.HowMany(); index++)
            SaveConfig(os << 'v' << SEP, Ver(index));

        os << COMMENT << G.BundleName() << endl;
        for (int index = 0; index < G.HowMany(); index++)
            SaveConfig(os << 'g' << SEP, G(index));

        os << COMMENT << W.BundleName() << endl;
        for (int index = 0; index < W.HowMany(); index++)
            SaveConfig(os << 'w' << SEP, W(index));

        if (Worm.Exist) {
            os << COMMENT << "Worm" << endl;
            SaveConfig(os << 'i' << SEP, Worm);
        }
    }
}

bool Diagram::_Load(istream &ifs)
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

bool Diagram::Load(const std::string &FileName)
{

    ifstream ifs(FileName, ios::in);
    ON_SCOPE_EXIT([&] {ifs.close(); });
    if (!ifs.is_open()) {
        ABORT("Cannot open " + FileName);
        return false;
    }
    return _Load(ifs);
}

bool Diagram::Load(const std::string &FileName, Lattice &lat, RandomFactory &rng, weight::G *g, weight::W *w)
{
    Reset(lat, rng, g, w);
    return Load(FileName);
}

/************************   write component to gv ****************************************/
string GLineStyle(bool IsMeasure, spin in, spin out)
{
    string color = "";
    if (IsMeasure)
        color = "color=\"green\"";
    else {
        string str[2] = {"blue", "red"};
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
