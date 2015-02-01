//
//  index_map.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/24/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "index_map.h"
#include "utility/logger.h"
#include "utility/abort.h"
#include <math.h>

using namespace weight;

IndexMap::IndexMap(real Beta_, uint MaxTauBin_, const Lattice& lat, TauSymmetry Symmetry_)
{
    MaxTauBin = MaxTauBin_;
    Beta = Beta_;
    _dBeta = Beta / MaxTauBin;
    _dBetaInverse = 1.0 / _dBeta;
    Lat = lat;
    Symmetry = Symmetry_;
    _TauSymmetryFactor = int(Symmetry);
    _Shape[SUB1] = (uint)Lat.SublatVol;
    _Shape[SUB2] = (uint)Lat.SublatVol;
    _Shape[VOL] = (uint)Lat.Vol;
    _Shape[TAU] = MaxTauBin;
}

int IndexMap::GetTauSymmetryFactor(real t_in, real t_out) const
{
    return (t_out > t_in) ? 1 : _TauSymmetryFactor;
}

int IndexMap::TauIndex(real tau) const
{
    if (DEBUGMODE && (tau < -Beta || tau >= Beta))
        LOG_INFO("tau=" << tau << " is out of the range ["
                        << -Beta << "," << Beta << ")");
    //TODO: mapping between tau and bin

    int bin = tau < 0 ? floor(tau * _dBetaInverse) + MaxTauBin
                      : floor(tau * _dBetaInverse);
    if (DEBUGMODE && (bin < 0 || bin >= MaxTauBin)) {
        LOG_INFO("tau=" << tau << " is out of the range ["
                        << -Beta << "," << Beta << ")");
        LOG_INFO("bin=" << bin << " is out of the range ["
                        << 0 << "," << MaxTauBin << "]");
    }
    return bin;
}

int IndexMap::TauIndex(real t_in, real t_out) const
{
    return TauIndex(t_out - t_in);
}

real IndexMap::IndexToTau(int Bin) const
{
    //TODO: mapping between tau and bin
    return Bin * _dBeta + _dBeta / 2;
}

const uint* IndexMap::GetShape() const
{
    return _Shape;
}

void IndexMap::_UpdateCache()
{
    _SizeDeltaT = 1;
    for (auto i = 0; i < DELTA_T_SIZE; i++) {
        _CacheDeltaT[DELTA_T_SIZE - 1 - i] = _SizeDeltaT;
        _SizeDeltaT *= _Shape[DELTA_T_SIZE - 1 - i];
    }
    _SizeSmoothT = 1;
    for (auto i = 0; i < SMOOTH_T_SIZE; i++) {
        _CacheSmoothT[SMOOTH_T_SIZE - 1 - i] = _SizeSmoothT;
        _SizeSmoothT *= _Shape[SMOOTH_T_SIZE - 1 - i];
    }
}

IndexMapSPIN2::IndexMapSPIN2(real Beta, uint MaxTauBin, const Lattice& Lat, TauSymmetry Symmetry)
    : IndexMap(Beta, MaxTauBin, Lat, Symmetry)
{
    _Shape[SP1] = 2;
    _Shape[SP2] = 2;
    _UpdateCache();
}

int IndexMapSPIN2::SpinIndex(spin SpinIn, spin SpinOut)
{
    return SpinIn * SPIN + SpinOut;
}

bool IndexMapSPIN2::IsSameSpin(int spindex)
{
    return (spindex == 0 || spindex == 2);
}

uint IndexMapSPIN2::GetIndex(spin in, spin out, const Site& rin, const Site& rout,
                             real tin, real tout) const
{
    auto coord = Lat.CoordiIndex(rin, rout);
    uint Index = in * _CacheSmoothT[SP1] + rin.Sublattice * _CacheSmoothT[SUB1]
                 + out * _CacheSmoothT[SP2] + rout.Sublattice * _CacheSmoothT[SUB2]
                 + coord * _CacheSmoothT[VOL] + TauIndex(tin, tout);
    if (DEBUGMODE && Index >= _SizeSmoothT)
        THROW_ERROR(IndexInvalid, "exceed array bound!");
    return Index;
}

uint IndexMapSPIN2::GetIndex(spin in, spin out,
                             const Site& rin, const Site& rout) const
{
    auto coord = Lat.CoordiIndex(rin, rout);
    uint Index = in * _CacheDeltaT[SP1] + rin.Sublattice * _CacheDeltaT[SUB1]
                 + out * _CacheDeltaT[SP2] + rout.Sublattice * _CacheDeltaT[SUB2]
                 + coord;
    if (DEBUGMODE && Index >= _SizeDeltaT)
        THROW_ERROR(IndexInvalid, "exceed array bound!");
    return Index;
}

IndexMapSPIN4::IndexMapSPIN4(real Beta, uint MaxTauBin, const Lattice& Lat, TauSymmetry Symmetry)
    : IndexMap(Beta, MaxTauBin, Lat, Symmetry)
{
    _Shape[SP1] = 4;
    _Shape[SP2] = 4;
    _UpdateCache();
}

//First In/Out: direction of WLine; Second In/Out: direction of Vertex
int IndexMapSPIN4::SpinIndex(spin SpinInIn, spin SpinInOut, spin SpinOutIn, spin SpinOutOut)
{
    return SpinInIn * SPIN3 + SpinInOut * SPIN2 + SpinOutIn * SPIN + SpinOutOut;
}
int IndexMapSPIN4::SpinIndex(const spin* TwoSpinIn, const spin* TwoSpinOut)
{
    return SpinIndex(TwoSpinIn[0], TwoSpinIn[1],
                     TwoSpinOut[0], TwoSpinOut[1]);
}

int IndexMapSPIN4::SpinIndex(const spin* Spin)
{
    return Spin[IN] * SPIN + Spin[OUT];
}

uint IndexMapSPIN4::GetIndex(const spin* SpinIn, const spin* SpinOut, const Site& rin, const Site& rout, real tin, real tout) const
{
    auto coord = Lat.CoordiIndex(rin, rout);
    uint Index = SpinIndex(SpinIn) * _CacheSmoothT[SP1] + rin.Sublattice * _CacheSmoothT[SUB1]
                 + SpinIndex(SpinOut) * _CacheSmoothT[SP2] + rout.Sublattice * _CacheSmoothT[SUB2] + coord * _CacheSmoothT[VOL] + TauIndex(tin, tout);
    if (DEBUGMODE && Index >= _SizeSmoothT)
        THROW_ERROR(IndexInvalid, "exceed array bound!");
    return Index;
}

uint IndexMapSPIN4::GetIndex(const spin* SpinIn, const spin* SpinOut, const Site& rin, const Site& rout) const
{
    auto coord = Lat.CoordiIndex(rin, rout);
    uint Index = SpinIndex(SpinIn) * _CacheDeltaT[SP1] + rin.Sublattice * _CacheDeltaT[SUB1]
                 + SpinIndex(SpinOut) * _CacheDeltaT[SP2] + rout.Sublattice * _CacheDeltaT[SUB2]
                 + coord;
    if (DEBUGMODE && Index >= _SizeDeltaT)
        THROW_ERROR(IndexInvalid, "exceed array bound!");
    return Index;
}