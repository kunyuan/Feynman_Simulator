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
#include "abort.h"
using namespace std;

const char SEP = ' ';
const char COMMENT = '#';
const string SEP_LINE = "#######################################################";
const string SEP_LINE_SHORT = "####";

/*******************  Read/write diagram to dat file ****************/
void Diagram::SaveConfig(const std::string &FileName, string Mode)
{
    ofstream os;
    if (Mode == "w")
        os.open(FileName, ios::out);
    else if (Mode == "a")
        os.open(FileName, ios::app);
    else
        ABORT("What is Mode=" << Mode << "?");
    if (!os.is_open()) {
        ABORT("Cannot open " + FileName);
    }
    else {
        os << SEP_LINE << endl;
        os << COMMENT << Ver.Name() << endl;
        for (int index = 0; index < Ver.HowMany(); index++)
            Ver[index].SaveConfig(os << 'v' << SEP);

        os << COMMENT << G.Name() << endl;
        for (int index = 0; index < G.HowMany(); index++)
            G[index].SaveConfig(os << 'g' << SEP);

        os << COMMENT << W.Name() << endl;
        for (int index = 0; index < W.HowMany(); index++)
            W[index].SaveConfig(os << 'w' << SEP);

        if(Worm.Exist)
        {
            os << COMMENT << "Worm" << endl;
            Worm.SaveConfig(os << 'i' <<SEP);
        }
    }
}

bool Diagram::LoadConfig(const std::string &FileName)
{
    ifstream ifs;
    ifs.open(FileName, ios::in);
    if (!ifs) {
        ABORT("Cannot find " + FileName);
        return false;
    }
    else {
        string line;
        int i = 0;
        //locate the last configration block
        streampos lastBlockPos = 0;
        while (std::getline(ifs, line)) {
            i++;
            if (line.compare(0, SEP_LINE_SHORT.size(), SEP_LINE_SHORT) == 0)
                lastBlockPos = ifs.tellg();
        }
        ifs.clear();
        ifs.seekg(lastBlockPos);
        ClearDiagram();
        char head;
        string temp;
        i=0;
        while (!ifs.eof()) {
            ifs >> head;
            if (head == COMMENT) {
                getline(ifs, temp);
                continue;
            }
            else if (head == 'g')
                G.Add().LoadConfig(ifs);
            else if (head == 'w')
                W.Add().LoadConfig(ifs);
            else if (head == 'v')
                Ver.Add().LoadConfig(ifs);
            else if (head == 'i')
                //TODO read from Worm
                Worm.LoadConfig(ifs);
            else
                ABORT("Error in reading diagram! Get " + ToString(head) + " as the head!");
            head=COMMENT;
        }
        FixDiagram();
        return true;
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

ostream &Diagram::Component2gv(ostream &os, GLine &r)
{
    os << r.Vertex[IN] << "->" << r.Vertex[OUT];
    os << " " << EdgeColor(Spin(r, IN), Spin(r, OUT)) << ";";
    os << "  //" << PrettyString(r);
    return os;
}

ostream &Diagram::Component2gv(ostream &os, WLine &r)
{
    os << r.Vertex[IN] << "->" << r.Vertex[OUT];
    os << " [style=dashed arrowhead=none]; ";
    os << "  //" << PrettyString(r);
    return os;
}

ostream &Diagram::Component2gv(ostream &os, Vertex &r)
{
    os << r.Name << " " << VertexColor(r.R.Sublattice) << ";";
    os << "  //" << PrettyString(r);
    return os;
}

/********************  write diagram to gv ********************/

template <typename T>
ostream &Diagram::Bundle2gv(ostream &os, Bundle<T> &r)
{
    os << "    //" << r.Name() << endl;
    for (int index = 0; index < r.HowMany(); index++) {
        Component2gv(os << "    ", r[index]);
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
