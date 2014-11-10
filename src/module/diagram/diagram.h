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
#include "utility/rng.h"
namespace weight {
class G;
class W;
}

namespace diag {
class Diagram {
  public:
    Diagram();

    void BuildNew(Lattice &, RandomFactory &, weight::G *, weight::W *);
    bool Load(const std::string &FileName, Lattice &, RandomFactory &, weight::G *, weight::W *);
    bool Load(const std::string &FileName);
    void Save(const std::string &FileName, std::string Mode = "a");
    void Reset(Lattice &, RandomFactory &, weight::G *, weight::W *);
    void SetTest(Lattice &, RandomFactory &, weight::G *, weight::W *);

    weight::G *GWeight;
    weight::W *WWeight;
    Lattice *Lat;
    RandomFactory *RNG;

    int Order;
    Complex Phase, Weight;
    real SignFermiLoop;

    Bundle<GLine> G;
    Bundle<WLine> W;
    Bundle<Vertex> Ver;

    WormClass Worm;
    bool IsWorm(vertex);
    
    bool MeasureGLine;
    gLine GMeasure;
    wLine WMeasure;

    spin Spin(gLine, int); ///For gline with different spins on two ends

    spin Spin(gLine); ///For gline with same spin on two ends
    void FlipGSpin(gLine);
    spin Spin(wLine, int, int);
    spin Spin(vertex, int);

    int Sublattice(gLine, int);
    int Sublattice(wLine, int);
    int Sublattice(vertex);

    //Randomly Pick
    gLine RandomPickG();
    wLine RandomPickW();
    vertex RandomPickVer();

    //Diagram
    bool FixDiagram();
    void ClearDiagram();

    //Diagram Check
    bool CheckG();
    bool CheckW();
    bool CheckVer();
    bool CheckWeight();
    bool CheckDiagram();

    //Diagram IO
    std::string PrettyString(gLine);
    std::string PrettyString(wLine);
    std::string PrettyString(vertex);

    void WriteDiagram2gv(std::string);
    void WriteDiagram2gv(std::ostream &);

  private:
    bool _Load(std::istream &);

    std::ostream &Component2gv(std::ostream &, gLine);
    std::ostream &Component2gv(std::ostream &, wLine);
    std::ostream &Component2gv(std::ostream &, vertex);

    bool LoadConfig(std::istream &is, WormClass &);
    bool LoadConfig(std::istream &is, gLine);
    bool LoadConfig(std::istream &is, wLine);
    bool LoadConfig(std::istream &is, vertex);

    void SaveConfig(std::ostream &is, WormClass &);
    void SaveConfig(std::ostream &is, gLine);
    void SaveConfig(std::ostream &is, wLine);
    void SaveConfig(std::ostream &is, vertex);

    template <typename T>
    std::ostream &Bundle2gv(std::ostream &, Bundle<T> &);
};

int TestDiagram();
}
#endif /* defined(__Fermion_Simulator__diagram_global__) */
