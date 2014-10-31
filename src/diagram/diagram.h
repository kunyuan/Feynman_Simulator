//
//  diagram_global.h
//  Fermion_Simulator
//
//  Created by Kun Chen on 10/9/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Fermion_Simulator__diagram_global__
#define __Fermion_Simulator__diagram_global__

#include <iosfwd>
#include "component_bundle.h"
namespace Weight {
class G;
class W;
}

class Diagram {
  public:
    Diagram();

    void Build(Lattice *, Weight::G *, Weight::W *);

    Weight::G *GWeight;
    Weight::W *WWeight;
    Lattice *Lat;

    int Order;
    Complex Phase, Weight;

    Bundle<GLine> G;
    Bundle<WLine> W;
    Bundle<Vertex> Ver;

    WormClass Worm;
    bool IsWorm(vertex);

    spin Spin(gLine, int);
    spin Spin(wLine, int, int);
    spin Spin(vertex, int);

    int Sublattice(gLine, int);
    int Sublattice(wLine, int);
    int Sublattice(vertex);

    string PrettyString(gLine);
    string PrettyString(wLine);
    string PrettyString(vertex);

    //Randomly Pick
    gLine RandomPickG();
    wLine RandomPickW();
    vertex RandomPickVer();

    //Diagram
    vertex NeighVer(gLine, int);
    vertex NeighVer(wLine, int);
    gLine NeighG(vertex, int);
    wLine NeighW(vertex);
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
