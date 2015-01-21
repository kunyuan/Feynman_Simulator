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

namespace weight {

enum Dim {
    SP1 = 0,
    SUB1,
    SP2,
    SUB2,
    VOL,
    TAU,
};
enum TauSymmetry {
    TauSymmetric = 1,
    TauAntiSymmetric = -1
};

const uint DELTA_T_SIZE = 5;
const uint SMOOTH_T_SIZE = 6;

enum SPIN4Filter { UpUp2UpUp,
                   UpDown2UpDown,
                   UpDown2DownUp };

class IndexMap {
public:
    IndexMap(real Beta, uint MaxTauBin, const Lattice& Lat, TauSymmetry Symmetry);
    int GetTauSymmetryFactor(real t_in, real t_out) const;
    const uint* GetShape() const; //the shape of internal weight array
    real Beta;
    Lattice Lat;
    uint MaxTauBin;
    TauSymmetry Symmetry;
    int TauIndex(real tau) const;
    int TauIndex(real t_in, real t_out) const;
    real IndexToTau(int TauIndex) const;

protected:
    void _UpdateCache();
    uint _Shape[SMOOTH_T_SIZE];
    uint _CacheDeltaT[DELTA_T_SIZE];
    uint _CacheSmoothT[SMOOTH_T_SIZE];
    uint _SizeDeltaT;
    uint _SizeSmoothT;
    int _TauSymmetryFactor;
    real _dBeta;
    real _dBetaInverse;
};

class IndexMapSPIN2 : public IndexMap {
public:
    IndexMapSPIN2(real Beta, uint MaxTauBin, const Lattice& Lat, TauSymmetry Symmetry);
    static bool IsSameSpin(int spindex);
    uint GetIndex(spin in, spin out,
                  const Site& rin, const Site& rout,
                  real tin, real tout) const;
    uint GetIndex(spin in, spin out,
                  const Site& rin, const Site& rout) const;

private:
    static int SpinIndex(spin SpinIn, spin SpinOut);
};

class IndexMapSPIN4 : public IndexMap {
public:
    IndexMapSPIN4(real Beta, uint MaxTauBin, const Lattice& Lat, TauSymmetry Symmetry);
    //First In/Out: direction of WLine; Second In/Out: direction of Vertex
    uint GetIndex(const spin* in, const spin* out,
                  const Site& rin, const Site& rout,
                  real tin, real tout) const;
    uint GetIndex(const spin* in, const spin* out,
                  const Site& rin, const Site& rout) const;

private:
    static int SpinIndex(const spin* Spin);
    static int SpinIndex(spin SpinInIn, spin SpinInOut, spin SpinOutIn, spin SpinOutOut);
    static int SpinIndex(const spin* TwoSpinIn, const spin* TwoSpinOut);
};
}

#endif /* defined(__Feynman_Simulator__index_map__) */
