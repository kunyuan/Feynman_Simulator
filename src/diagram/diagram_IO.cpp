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
void Diagram::WriteDiagram(ostream &os)
{
    os << SEP_LINE << endl;
    os << COMMENT << Ver.Name() << endl;
    for (int index = 0; index < Ver.HowMany(); index++)
        os << 'v' << SEP<<Ver[index];
    
    os << SEP_LINE << endl;
    os << COMMENT << G.Name() << endl;
    for (int index = 0; index < G.HowMany(); index++)
        os << 'g' << SEP<<G[index];
    
    os << SEP_LINE << endl;
    os << COMMENT << W.Name() << endl;
    for (int index = 0; index < W.HowMany(); index++)
        os << 'w' << SEP<<W[index];
    
    if(Worm.Exist)
    {
        os << SEP_LINE << endl;
        os << COMMENT << "Worm" << endl;
        os << "worm" << SEP<< Worm.Ira <<SEP <<Worm.Masha <<SEP<<Worm.dSpin
            <<SEP<<Worm.K<<endl;
    }
    
    os << SEP_LINE << endl;
    os << COMMENT << "Weight" << Weight << endl;
    os << endl;
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
    ClearDiagram();
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

ostream& Diagram::Component2gv(ostream &os, GLine &r)
{
    os << r.Vertex[IN] << "->" << r.Vertex[OUT];
    os << " " << EdgeColor(Spin(r,IN), Spin(r,OUT)) << ";";
    os << "  //Name=" << r.Name << "; K="<<r.K<< "; Weight=" << r.Weight << endl;
    return os;
}

ostream& Diagram::Component2gv(ostream &os, WLine &r)
{
    os << r.Vertex[IN] << "->" << r.Vertex[OUT];
    os << " [style=dashed arrowhead=none]; ";
    os << "  //Name=" << r.Name << "; K="<<r.K<<"; Weight=" << r.Weight << endl;
    return os;
}

ostream& Diagram::Component2gv(ostream &os, Vertex &r)
{
    os << r.Name << " " << VertexColor(r.R.Sublattice) << ";";
    os << "  //InG=" << r.G[IN] << "; OutG=" << r.G[OUT] << "; W=" << r.W << endl;
    return os;
}

/********************  write diagram to gv ********************/

template <typename T>
ostream& Diagram::Bundle2gv(ostream &os, Bundle<T> &r)
{
    os << "    //" << r.Name() << endl;
    for (int index = 0; index < r.HowMany(); index++) {
        Component2gv(os << "    ", r[index]);
    }
    os << endl;
    return os;
}

template ostream& Diagram::Bundle2gv(ostream &os, Bundle<GLine> &r);
template ostream& Diagram::Bundle2gv(ostream &os, Bundle<WLine> &r);
template ostream& Diagram::Bundle2gv(ostream &os, Bundle<Vertex> &r);

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