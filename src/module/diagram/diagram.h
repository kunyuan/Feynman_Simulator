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

    void BuildNew(Lattice &, weight::G &, weight::W &);
    bool Load(const std::string &FileName, Lattice &,
              weight::G &, weight::W &);
    bool Load(const std::string &FileName);
    void Save(const std::string &FileName, std::string Mode = "a");
    void Reset(Lattice &, weight::G &, weight::W &);
    void SetTest(Lattice &, weight::G &, weight::W &);
    bool CheckDiagram();
    bool FixDiagram();

    Lattice *Lat;
    weight::G *GWeight;
    weight::W *WWeight;

    int Order;
    Complex Phase, Weight;
    real SignFermiLoop;

    Bundle<GLine> G;
    Bundle<WLine> W;
    Bundle<Vertex> Ver;
    
    bool GHash[2*MAX_K+1];
    bool WHash[MAX_K+1];

    WormClass Worm;
    bool IsWorm(vertex);

    bool MeasureGLine;
    gLine GMeasure;
    wLine WMeasure;

    //Diagram
    void ClearDiagram();

    //Diagram Hash Table Check
    bool GHashCheck(Momentum);
    bool WHashCheck(Momentum);
    
    void AddGHash(Momentum);
    void AddWHash(Momentum);
    void RemoveGHash(Momentum);
    void RemoveWHash(Momentum);
    void ReplaceGHash(Momentum, Momentum);
    void ReplaceWHash(Momentum, Momentum);

    void WriteDiagram2gv(std::string);

  private:
    bool _CheckTopo();
    bool _CheckStatus();
    bool _CheckK();
    bool _CheckSpin();
    bool _CheckWeight();

    bool _Load(std::istream &);
    bool LoadConfig(std::istream &is, WormClass &);
    bool LoadConfig(std::istream &is, gLine);
    bool LoadConfig(std::istream &is, wLine);
    bool LoadConfig(std::istream &is, vertex);

    void SaveConfig(std::ostream &is, WormClass &);
    void SaveConfig(std::ostream &is, gLine);
    void SaveConfig(std::ostream &is, wLine);
    void SaveConfig(std::ostream &is, vertex);
};

int TestDiagram();
}
#endif /* defined(__Fermion_Simulator__diagram_global__) */
