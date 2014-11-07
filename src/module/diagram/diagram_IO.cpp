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
            //TODO read from Worm
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

    ifstream ifs;
    ifs.open(FileName, ios::in);
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
string EdgeColor(spin in, spin out)
{
    string str[2] = {"blue", "red"};
    return "[color=\"" + str[in] + ":" + str[out] + ";0.5\"]";
}

string VertexColor(int sublattice)
{
    string colorstr;
    if (sublattice == 0)
        colorstr = "grey";
    else if (sublattice == 1)
        colorstr = "palegreen";
    else
        colorstr = "palegreen";
    return "[fillcolor=" + colorstr + "]";
}

ostream &Diagram::Component2gv(ostream &os, gLine r)
{
    os << r->nVer[IN]->Name << "->" << r->nVer[OUT]->Name;
    os << " " << EdgeColor(r->Spin(IN), r->Spin(OUT)) << ";";
    os << "  //" << r->PrettyString() << endl;
    return os;
}

ostream &Diagram::Component2gv(ostream &os, wLine r)
{
    os << r->nVer[IN]->Name << "->" << r->nVer[OUT]->Name;
    os << " [style=dashed arrowhead=none]; ";
    os << "  //" << r->PrettyString() << endl;
    return os;
}

ostream &Diagram::Component2gv(ostream &os, vertex r)
{
    os << r->Name << " " << VertexColor(r->R.Sublattice) << ";";
    os << "  //" << r->PrettyString() << endl;
    return os;
}

/********************  write diagram to gv ********************/

template <typename T>
ostream &Diagram::Bundle2gv(ostream &os, Bundle<T> &r)
{
    os << "    //" << r.BundleName() << endl;
    for (int index = 0; index < r.HowMany(); index++) {
        Component2gv(os << "    ", r(index));
    }
    os << endl;
    return os;
}

template ostream &Diagram::Bundle2gv(ostream &os, Bundle<GLine> &r);
template ostream &Diagram::Bundle2gv(ostream &os, Bundle<WLine> &r);
template ostream &Diagram::Bundle2gv(ostream &os, Bundle<Vertex> &r);

/**
*  write diagram object to .gv file so that it can be visualized by Graphviz
*
*  @param os ostream& as the target stream to output
*/
void Diagram::WriteDiagram2gv(ostream &os)
{
    string head = "digraph Feynman{\n";
    string tail = "}\n";
    string dpi = "graph[dpi=200];\n";
    string node_attribute = "    node [margin=0.1 fillcolor=grey fontcolor=black fontsize=10 width=0.2 shape=circle style=filled fixedsize=true]\n";
    os << head << dpi << node_attribute << endl;
    Bundle2gv(os, Ver);
    Bundle2gv(os, G);
    Bundle2gv(os, W);
    os << tail;
}

void Diagram::WriteDiagram2gv(string path)
{
    ofstream os;
    os.open(path, ios::out);
    if (!os) {
        ABORT("Cannot open " + path);
    }
    else {
        WriteDiagram2gv(os);
    }
}
