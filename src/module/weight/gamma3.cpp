//
// Created by kun on 10/18/17.
//

#include "component.h"
#include "gamma3.h"
#include "utility/dictionary.h"

using namespace weight;
using namespace std;


GammaGClass::GammaGClass(const Lattice &lat, real beta, uint MaxTauBin, real Norm)
    : _Map(IndexMapSPIN2(beta, MaxTauBin, lat, TauAntiSymmetric))
{
    int Vol = _Map.Lat.Vol;
    _Beta = _Map.Beta;
    _Norm = Norm * (_Map.MaxTauBin / _Beta) / _Beta / Vol;
    uint SubVol=(uint)lat.SublatVol;
    uint MeaShape[6]={2,SubVol,SubVol,(uint)Vol,MaxTauBin,MaxTauBin};
    //Gspin, Gsub, Usub, G_r1, Gtau1,dtau2
    _WeightAccu.Allocate(MeaShape, SMOOTH);
    _WeightSize = _WeightAccu.GetSize();
    ClearStatistics();
}

Complex GammaGClass::Weight(const Site &Gr1, const Site &Gr2, const Site &Ur,
                            real Gt1, real Gt2, real Ut, spin Gspin1, spin Gspin2, spin Uspin) const
{
    return Complex(1.0,0.0);
}

void GammaGClass::MeasureNorm(real weight)
{
    _NormAccu += weight;
}

void GammaGClass::Measure(const Site &Gr1, const Site &Gr2, const Site &Ur,
                          real Gt1, real Gt2, real Ut1, spin Gspin1, spin Gspin2, spin Uspin,
              const Complex &weight)
{
    //TODO: calculate the Index
//    _WeightAccu[Index] += weight;
}


void GammaGClass::ClearStatistics()
{
    _NormAccu = 0.0;
    _WeightAccu.Assign(0.0);
}
//TODO: you may have to replace int with size_t here

void GammaGClass::SqueezeStatistics(real factor)
{
    ASSERT_ALLWAYS(factor > 0, "factor=" << factor << "<=0!");
    _NormAccu /= factor;
    _WeightAccu *= 1.0 / factor;
}

/**********************   Weight IO ****************************************/

bool GammaGClass::FromDict(const Dictionary& dict)
{
    _Norm = dict.Get<real>("Norm");
    _NormAccu = dict.Get<real>("NormAccu");
    auto arr = dict.Get<Python::ArrayObject>("WeightAccu");
    //assert estimator shape except order dimension
    ASSERT_ALLWAYS(Equal(arr.Shape().data() + 1, _WeightAccu.GetShape() + 1, _WeightAccu.GetDim() - 1), "Shape should match!");
    _WeightAccu.Assign(0.0);
    _WeightAccu.Assign(arr.Data<Complex>(), arr.Size());
    return true;
}

Dictionary GammaGClass::ToDict()
{
    Dictionary dict;
    dict["Norm"] = _Norm;
    dict["NormAccu"] = _NormAccu;
    dict["WeightAccu"] = Python::ArrayObject(_WeightAccu.Data(), _WeightAccu.GetShape(), _WeightAccu.GetDim());
    return dict;
}


GammaWClass::GammaWClass(const Lattice &lat, real beta, uint MaxTauBin, real Norm)
        : _Map(IndexMapSPIN4(beta, MaxTauBin, lat, TauSymmetric))
{
    int Vol = _Map.Lat.Vol;
    _Beta = _Map.Beta;
    _Norm = Norm * (_Map.MaxTauBin / _Beta) / _Beta / Vol;
    uint SubVol=(uint)lat.SublatVol;
    uint MeaShape[8]={6,SubVol,SubVol,(uint)Vol,(uint)Vol, MaxTauBin,MaxTauBin};
    //Wspin12, Wsub1, Wsub2, Usub, W_r1, dW_r2, Wtau1,dtau2
    _WeightAccu.Allocate(MeaShape, SMOOTH);
    _WeightSize = _WeightAccu.GetSize();
//    ClearStatistics();
}

Complex GammaWClass::Weight(const Site &Wr1, const Site &Wr2, const Site &Ur,
                            real Wt1, real Wt2, real Ut, spin* Wspin1, spin* Wspin2, spin Uspin) const
{
    return Complex(1.0,0.0);
}

void GammaWClass::Measure(const Site &Wr1, const Site &Wr2, const Site &Ur, real Wt1, real Wt2,
                          real Ut, spin* Wspin1, spin* Wspin2, spin Uspin, const Complex &Weight)
{
    //TODO: calculate the Index
//    _WeightAccu[Index] += weight;
}

void GammaWClass::MeasureNorm(real weight)
{
    _NormAccu += weight;
}

void GammaWClass::ClearStatistics()
{
    _NormAccu = 0.0;
    _WeightAccu.Assign(0.0);
}
//TODO: you may have to replace int with size_t here

void GammaWClass::SqueezeStatistics(real factor)
{
    ASSERT_ALLWAYS(factor > 0, "factor=" << factor << "<=0!");
    _NormAccu /= factor;
    _WeightAccu *= 1.0 / factor;
}

bool GammaWClass::FromDict(const Dictionary& dict)
{
    _Norm = dict.Get<real>("Norm");
    _NormAccu = dict.Get<real>("NormAccu");
    auto arr = dict.Get<Python::ArrayObject>("WeightAccu");
    //assert estimator shape except order dimension
    ASSERT_ALLWAYS(Equal(arr.Shape().data() + 1, _WeightAccu.GetShape() + 1, _WeightAccu.GetDim() - 1), "Shape should match!");
    _WeightAccu.Assign(0.0);
    _WeightAccu.Assign(arr.Data<Complex>(), arr.Size());
    return true;
}

Dictionary GammaWClass::ToDict()
{
    Dictionary dict;
    dict["Norm"] = _Norm;
    dict["NormAccu"] = _NormAccu;
    dict["WeightAccu"] = Python::ArrayObject(_WeightAccu.Data(), _WeightAccu.GetShape(), _WeightAccu.GetDim());
    return dict;
}

