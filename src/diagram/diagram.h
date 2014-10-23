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
#include "component_bundle.h"
#include "../utility/rng.h"
#include "../observable/weight.h"

class Diagram {
  public:
    Diagram();

    void SetGWWeight(Weight::G *, Weight::W *);
    Weight::G *GWeight;
    Weight::W *WWeight;

    void SetLat(Lattice *);
    Lattice *Lat;

    int Order;
    Complex Phase, Weight;

    Bundle<GLine> G;
    Bundle<WLine> W;
    Bundle<Vertex> Ver;

    Worm Worm;
    bool IsWorm(const Vertex &);

    spin Spin(GLine &, const int &);
    spin Spin(WLine &, const int &, const int &);
    spin Spin(Vertex &, const int &);

    int Sublattice(GLine &, const int &);
    int Sublattice(WLine &, const int &);
    int Sublattice(Vertex &);

    string PrettyString(GLine &);
    string PrettyString(WLine &);
    string PrettyString(Vertex &);

    //Randomly Pick
    GLine &RandomPickG();
    WLine &RandomPickW();
    Vertex &RandomPickVer();

    //Diagram
    Vertex &NeighVer(GLine &, int dir);
    Vertex &NeighVer(WLine &, int dir);
    GLine &NeighG(Vertex &, int dir);
    WLine &NeighW(Vertex &);
    bool FixDiagram();
    void ClearDiagram();

    //Diagram Check
    bool CheckG();
    bool CheckW();
    bool CheckVer();
    bool CheckWeight();
    bool CheckDiagram();

    //Diagram IO
    bool LoadConfig(const std::string &FileName);
    void SaveConfig(const std::string &FileName, std::string Mode = "a");

    void WriteDiagram2gv(std::string);
    void WriteDiagram2gv(std::ostream &);

  private:
    std::ostream &Component2gv(std::ostream &, GLine &);
    std::ostream &Component2gv(std::ostream &, WLine &);
    std::ostream &Component2gv(std::ostream &, Vertex &);

    template <typename T>
    std::ostream &Bundle2gv(std::ostream &, Bundle<T> &);
};
#endif /* defined(__Fermion_Simulator__diagram_global__) */
