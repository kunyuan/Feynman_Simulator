//
//  diagram_IO_global.cpp
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/9/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include <vector>
#include <sstream>
#include "component.h"
#include "component_bundle.h"
#include "diagram.h"
#include "pystring.h"
#include "abort.h"
using namespace std;

#define SEP ' '
#define COMMENT '#'
#define SEP_LINE "#######################################################"
#define SEP_LINE_SHORT "####"

/*******************  Read/write diagram to dat file ****************/

template <typename T>
ostream &operator<<(ostream &os, Bundle<T> &r)
{
    os << SEP_LINE << endl;
    os << COMMENT << r.Name() << endl;
    for (int index = 0; index < r.HowMany(); index++) {
        if (r.Name() == "GLine")
            os << 'g' << SEP;
        else if (r.Name() == "WLine")
            os << 'w' << SEP;
        else if (r.Name() == "Vertex")
            os << 'v' << SEP;
        else
            ABORT("What the hell is " + ToString(r.Name()) + "?");
        os << r[index];
    }
    return os;
}

template ostream &operator<<(ostream &os, Bundle<GLine> &r);
template ostream &operator<<(ostream &os, Bundle<WLine> &r);
template ostream &operator<<(ostream &os, Bundle<Vertex> &r);

void Diagram::WriteDiagram(ostream &os)
{
    os << Ver << G << W;
}

void Diagram::WriteDiagram(string path)
{
    ofstream os;
    os.open(path, ios::out);
    if (!os) {
        ABORT("Cannot open " + path);
    }
    else {
        WriteDiagram(os);
    }
}

bool Diagram::ReadDiagram(istream &is)
{
    char head;
    string temp;
    while (!is.eof()) {
        is >> head;
        if (head == COMMENT) {
            getline(is, temp);
            continue;
        }
        else if (head == 'g')
            is >> G.Add();
        else if (head == 'w')
            is >> W.Add();
        else if (head == 'v')
            is >> Ver.Add();
        else
            ABORT("Error in reading diagram! Get " + ToString(head) + " as the head!");
    }
    FixDiagram();
    return true;
}

bool Diagram::ReadDiagram(string path)
{
    ifstream is;
    is.open(path, ios::in);
    if (!is) {
        ABORT("Cannot find " + path);
        return false;
    }
    else {
        return ReadDiagram(is);
    }
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

ostream &Component2gv(ostream &os, GLine &r)
{
    os << r.Vertex[IN] << "->" << r.Vertex[OUT];
    os << " " << EdgeColor(r.Spin[IN], r.Spin[OUT]) << ";";
    os << "  //Name=" << r.Name << "; Weight=" << r.Weight << endl;
    return os;
}

ostream &Component2gv(ostream &os, WLine &r)
{
    os << r.Vertex[IN] << "->" << r.Vertex[OUT];
    os << " [style=dashed arrowhead=none]; ";
    os << "  //Name=" << r.Name << "; Weight=" << r.Weight << endl;
    return os;
}

ostream &Component2gv(ostream &os, Vertex &r)
{
    os << r.Name << " " << VertexColor(r.Sublattice) << ";";
    os << "  //InG=" << r.G[IN] << "; OutG=" << r.G[OUT] << "; W=" << r.W << endl;
    return os;
}

/********************  write diagram to gv ********************/
template <typename T>
ostream &Bundle2gv(ostream &os, Bundle<T> &r)
{
    os << "    //" << r.Name() << endl;
    for (int index = 0; index < r.HowMany(); index++) {
        Component2gv(os << "    ", r[index]);
    }
    os << endl;
    return os;
}

template ostream &Bundle2gv(ostream &os, Bundle<GLine> &r);
template ostream &Bundle2gv(ostream &os, Bundle<WLine> &r);
template ostream &Bundle2gv(ostream &os, Bundle<Vertex> &r);

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