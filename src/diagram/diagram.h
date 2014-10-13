//
//  diagram_global.h
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/9/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Fermion_Simulator__diagram_global__
#define __Fermion_Simulator__diagram_global__

#include <iostream>
#include "logger.h"
#include "component_bundle.h"

class Diagram {
  public:
    Diagram();
    Bundle<GLine> G;
    Bundle<WLine> W;
    Bundle<Vertex> Ver;

    //GLine
    spin G_Spin(GLine &g, int dir);
    int G_Sublattice(GLine &g, int dir);

    //Vertex
    spin VerSpin(Vertex &, int dir);

    //Diagram
    Vertex &NeighVer(GLine &, int dir);
    Vertex &NeighVer(WLine &, int dir);
    GLine &NeighG(Vertex &, int dir);
    WLine &NeighW(Vertex &);
    bool FixDiagram();

    //Diagram Check
    bool CheckG();
    bool CheckW();
    bool CheckVer();
    bool CheckDiagram();

    //Diagram IO
    bool ReadDiagram(std::string);
    bool ReadDiagram(std::istream &);
    void WriteDiagram(std::string);
    void WriteDiagram(std::ostream &);
    void WriteDiagram2gv(std::string);
    void WriteDiagram2gv(std::ostream &);
};
#endif /* defined(__Fermion_Simulator__diagram_global__) */
