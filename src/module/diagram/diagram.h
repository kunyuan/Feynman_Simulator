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
class Dictionary;

namespace diag {
class Diagram {
public:
    Diagram();

    void BuildNew(Lattice&, weight::G&, weight::W&);
    bool Load(const std::string& FileName, Lattice&,
              weight::G&, weight::W&);
    bool Load(const std::string& FileName);
    void Save(const std::string& FileName, std::string Mode = "a");
    void Reset(Lattice&, weight::G&, weight::W&);
    void SetTest(Lattice&, weight::G&, weight::W&);
    bool CheckDiagram();
    bool FixDiagram();

    Lattice* Lat;
    weight::G* GWeight;
    weight::W* WWeight;

    int Order;
    Complex Phase, Weight;
    real SignFermiLoop;

    Bundle<GLine> G;
    Bundle<WLine> W;
    Bundle<Vertex> Ver;

    bool GHash[2 * MAX_K + 1];
    bool WHash[MAX_K + 1];

    WormClass Worm;
    bool IsWorm(vertex);

    bool MeasureGLine;
    gLine GMeasure;
    wLine WMeasure;
    Complex ConstSigma=Complex(1.0, 0.0);
    Complex ConstPolar=Complex(1.0, 0.0);

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

    bool _Load(const Dictionary&);
    void _FromDict(const Dictionary&, wLine);
    void _FromDict(const Dictionary&, gLine);
    void _FromDict(const Dictionary&, vertex);
    void _FromDict(const Dictionary& WormDict, WormClass& worm);
    Dictionary _ToDict(WormClass);
    Dictionary _ToDict(wLine);
    Dictionary _ToDict(gLine);
    Dictionary _ToDict(vertex);
};

int TestDiagram();
}
#endif /* defined(__Fermion_Simulator__diagram_global__) */
