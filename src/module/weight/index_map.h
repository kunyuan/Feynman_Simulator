//
//  index_map.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/24/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__index_map__
#define __Feynman_Simulator__index_map__

#include "utility/convention.h"
#include "lattice/lattice.h"
#include <vector>

namespace weight {
const uint MAX_TAU_BIN = 32;

enum SPIN4Filter { UpUp2UpUp,
                   UpDown2UpDown,
                   UpDown2DownUp };

class IndexMap {
  public:
    IndexMap(real Beta, const Lattice &Lat);
    int TauIndex(real tau);
    int TauIndex(real t_in, real t_out);
    real IndexToTau(int TauIndex);

    int SublatIndex(const Distance &dist);
    int CoordiIndex(const Distance &dist);

  protected:
    real _Beta;
    real _dBeta;
    real _dBetaInverse;
    Lattice _Lat;
};

class IndexMapSPIN2 : public IndexMap {
  public:
    using IndexMap::IndexMap;
    static int SpinIndex(spin SpinIn, spin SpinOut);
    static bool IsSameSpin(int spindex);
    void Map(uint *result, spin in, spin out,
             const Site &rin, const Site &rout, real tin, real tout);
    void MapDeltaT(uint *result, spin in, spin out, const Site &rin, const Site &rout);
};

class IndexMapSPIN4 : public IndexMap {
  public:
    using IndexMap::IndexMap;
    //First In/Out: direction of WLine; Second In/Out: direction of Vertex
    static int SpinIndex(spin SpinInIn, spin SpinInOut, spin SpinOutIn, spin SpinOutOut);
    static int SpinIndex(const spin *TwoSpinIn, const spin *TwoSpinOut);
    static std::vector<int> GetSpinIndexVector(SPIN4Filter filter);
    void Map(uint *result, const spin *in, const spin *out,
             const Site &rin, const Site &rout, real tin, real tout);
    void MapDeltaT(uint *result, const spin *in, const spin *out,
                   const Site &rin, const Site &rout);
};
}

#endif /* defined(__Feynman_Simulator__index_map__) */
