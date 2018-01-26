//
// Created by kun on 10/18/17.
//

#include "gamma3.h"
#include "utility/dictionary.h"
#include "utility/complex.h"
#include "iostream"

using namespace weight;
using namespace std;


GammaGClass::GammaGClass(const Lattice &lat, real beta, uint MaxTauBin, real Norm)
    : _Map(IndexMapSPIN2(beta, MaxTauBin, lat, TauAntiSymmetric))
{
    int Vol = _Map.Lat.Vol;
    _Beta = _Map.Beta;
    _MaxTauBin = MaxTauBin;
    _dBetaInverse = _MaxTauBin/_Beta;

    _Norm = Norm * pow(_Map.MaxTauBin / _Beta, 2.0) / _Beta / Vol;
//    _Norm = pow(_Map.MaxTauBin / _Beta, 1.0) / _Beta / pow(Vol, 1.0)*Norm;

//    uint SubVol=(uint)lat.SublatVol;
    uint MeaShape[4] = {2, (uint)Vol, MaxTauBin, MaxTauBin};
    //Gspin, Gsub, Usub, G_r1, Gtau1,dtau2
    _Weight.Allocate(MeaShape, SMOOTH);
    _WeightAccu.Allocate(MeaShape, SMOOTH);
    _WeightSize = _WeightAccu.GetSize();
    _CacheIndex[0] = 1;
    _CacheIndex[1] = MaxTauBin;
    _CacheIndex[2] = MaxTauBin*MaxTauBin;
    _CacheIndex[3] = MaxTauBin*MaxTauBin*Vol;
    ClearStatistics();
}

Complex GammaGClass::Weight(const Site &Gr_in, const Site &Gr_out, const Site &Ur,
                            real Gt_in, real Gt_out, real Ut, spin Gspin_in, spin Gspin_out, spin Uspin) const
{
    auto coord = _Map.Lat.CoordiIndex(Gr_in, Ur);

    real sign=1;
    auto t1 = Gt_out - Ut;
    if(t1<0){
        sign = -sign;
        t1 += _Beta;
    }
    int t1Index = floor(t1*_dBetaInverse);

    auto t2 = Gt_in - Ut;
    if(t2 < 0){
        sign = -sign;
        t2 += _Beta;
    }
    int t2Index = floor(t2*_dBetaInverse);

    //if(t1Index==t2Index){
        //if(t1>=t2) 
            //t1Index++;
        //else
            //t2Index++;
        //if(t1Index==_MaxTauBin) t1Index=0;
        //if(t2Index==_MaxTauBin) t2Index=0;
    //}

    uint Index = Gspin_in * _CacheIndex[3] + coord * _CacheIndex[2] + t1Index*_CacheIndex[1] + t2Index;
    return _Weight(Index) * sign;
}

void GammaGClass::MeasureNorm(real weight)
{
    _NormAccu += weight;
}

void GammaGClass::Measure(const Site &Gr_in, const Site &Gr_out, const Site &Ur,
                          real Gt_in, real Gt_out, real Ut, spin Gspin_in, spin Gspin_out, spin Uspin,
              const Complex &weight)
{
    auto coord = _Map.Lat.CoordiIndex(Gr_in, Ur);

//    int sign=1;
//    auto t1=Gt_out-Ut;
//    if(t1<0){
//        sign=-sign;
//        t1+=_Beta;
//    }
//    int t1Index=floor(t1*_dBetaInverse);

//    //TODO: Test for Sigma
//    real sign = 1.0;
//    int t1Index = 0;
//
//    auto dt = Gt_out - Gt_in;
//    if(dt < 0){
//        sign = -sign;
//        dt += _Beta;
//    }
//    int dtIndex = floor(dt * _dBetaInverse);
//
//    uint Index = Gspin_in * _CacheIndex[3] + coord * _CacheIndex[2] + t1Index*_CacheIndex[1] + dtIndex;

    real sign=1;
    auto t1 = Gt_out - Ut;
    if(t1<0){
        sign = -sign;
        t1 += _Beta;
    }
    int t1Index = floor(t1*_dBetaInverse);

    auto t2 = Gt_in - Ut;
    if(t2 < 0){
        sign = -sign;
        t2 += _Beta;
    }
    int t2Index = floor(t2*_dBetaInverse);

    uint Index = Gspin_in * _CacheIndex[3] + coord * _CacheIndex[2] + t1Index*_CacheIndex[1] + t2Index;
    _WeightAccu[Index] += weight * sign;
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

bool GammaGClass::WeightFromDict(const Dictionary& dict) {
    _Weight.FromDict(dict);
    return true;
}

bool GammaGClass::StatisFromDict(const Dictionary &dict)
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

Dictionary GammaGClass::StatisToDict()
{
    Dictionary dict;
    dict["Norm"] = _Norm;
    dict["NormAccu"] = _NormAccu;
    dict["WeightAccu"] = Python::ArrayObject(_WeightAccu.Data(), _WeightAccu.GetShape(), _WeightAccu.GetDim());
    return dict;
}

Dictionary GammaGClass::WeightToDict()
{
    return _Weight.ToDict();
}


GammaWClass::GammaWClass(const Lattice &lat, real beta, uint MaxTauBinTiny,  std::vector<basis> & BoseBasis, real Norm)
        : _Map(IndexMapSPIN4(beta, MaxTauBinTiny, lat, TauSymmetric))
{
    int Vol = _Map.Lat.Vol;
    _Vol=Vol;
    _Beta = _Map.Beta;
    _Norm = Norm * pow(_Map.MaxTauBin / _Beta, 2.0) / _Beta / Vol;

    _MaxTauBinTiny = MaxTauBinTiny;
    _dBetaInverse = _MaxTauBinTiny/_Beta;
//    uint SubVol=(uint)lat.SublatVol;
    _SpaceTimeSize=(MaxTauBinTiny*Vol+1)*(MaxTauBinTiny*Vol)/2;

//    uint WeightShape[5] = {2, (uint)Vol, (uint)Vol, MaxTauBin, MaxTauBin};
    //uint WeightShape[5] = {2, _SpaceTimeSize};
    //Wspin12, W_r1, dW_r2, Wtau1,dtau2
    //_Weight.Allocate(WeightShape, SMOOTH);

    _BasisVec.swap(BoseBasis);
    _BasisNum=_BasisVec.size();
    _BasisMaxTauBin=_BasisVec[0].size();
    _dBetaInverse = _BasisMaxTauBin/_Beta;

//    uint MeaShape[5] = {2, (uint)Vol, (uint)Vol, _MaxTauBin, _MaxTauBin};
//    uint MeaShape[5] = {2, _SpaceTimeSize};
    //_WeightAccu.Allocate(WeightShape, SMOOTH);
    //_WeightSize = _WeightAccu.GetSize();

    _CacheIndex[0] = 1;
    _CacheIndex[1] = MaxTauBinTiny;
    _CacheIndex[2] = MaxTauBinTiny*MaxTauBinTiny;
    _CacheIndex[3] = MaxTauBinTiny*MaxTauBinTiny*Vol;
    _CacheIndex[4] = MaxTauBinTiny*MaxTauBinTiny*Vol*Vol;

    _CacheIndexBasis[0] = 1;
    _CacheIndexBasis[1] = _BasisNum;
    _CacheIndexBasis[2] = _BasisNum*_BasisNum;
    _CacheIndexBasis[3] = _BasisNum*_BasisNum*Vol;
    _CacheIndexBasis[4] = _BasisNum*_BasisNum*Vol*Vol;


    ClearStatistics();
}

uint GammaWClass::_SpinIndex(spin * L, spin * R) const {
    if(L[0] == UP && L[1] == UP && R[0] == UP && R[1] == UP)
        return 0;
    else if (L[0] == DOWN && L[1] == DOWN && R[0] == DOWN && R[1] == DOWN)
        return 1;
    else if (L[0] == UP && L[1] == UP && R[0] == DOWN && R[1] ==DOWN)
        return 2;
    else if (L[0] == DOWN && L[1] == DOWN && R[0] == UP && R[1] ==UP)
        return 3;
    else if (L[0] == UP && L[1] == DOWN && R[0] == DOWN && R[1] ==UP)
        return 4;
    else if (L[0] == DOWN && L[1] == UP && R[0] == UP && R[1] ==DOWN)
        return 5;
    else
        return -1;
}

uint GammaWClass::_GetIndex(int r1, int r2, real t1, real t2, int& SymmetryFactor, bool& DoesMirrored) const {
    if(t1<0){
        t1+=_Beta;
    }
    int t1Index=floor(t1*_dBetaInverse);

    if(t2<0){
        t2+=_Beta;
    }
    int t2Index=floor(t2*_dBetaInverse);

    auto Row=r1*_MaxTauBinTiny+t1Index;
    auto Col=r2*_MaxTauBinTiny+t2Index;

    if(Row>Col){
        auto temp=Row;
        Row=Col;
        Col=temp;
        DoesMirrored=true;
    } else{
        DoesMirrored=false;
    }
    if(Row==Col)
        SymmetryFactor=1.0;
    else
        SymmetryFactor=2.0;
//    uint Index = coord_r1 * _CacheIndex[3] + coord_r2 * _CacheIndex[2]+ t1Index * _CacheIndex[1] + t2Index;
    uint Index=(_MaxTauBinTiny*_Vol*Row)+Col-(Row*(Row+1))/2;
    return Index;
}

//Complex GammaWClass::Weight(const Site &Wr_in, const Site &Wr_out, const Site &Ur,
//                            real Wt_in, real Wt_out, real Ut, spin* Wspin_in, spin* Wspin_out, spin Uspin) const
//{
//    auto coord_r1 = _Map.Lat.CoordiIndex(Wr_out, Ur);
//    auto coord_r2 = _Map.Lat.CoordiIndex(Wr_in,  Ur);
//
//    auto t1 = Wt_out - Ut;
//    if(t1<0){
//        t1+=_Beta;
//    }
//    int t1Index=floor(t1*_dBetaInverse);
//
//    auto t2=Wt_in-Ut;
//    if(t2<0){
//        t2+=_Beta;
//    }
//    int t2Index=floor(t2*_dBetaInverse);
//    auto SpinIndex=_SpinIndex(Wspin_out, Wspin_in);
//    uint Index = coord_r1 * _CacheIndex[3] + coord_r2 * _CacheIndex[2]+ t1Index * _CacheIndex[1] + t2Index;
//    if(SpinIndex==0||SpinIndex==1)
//        return _Weight(Index);
//    else if(SpinIndex==2||SpinIndex==3)
//        return -_Weight(Index);
//    else if(SpinIndex==4){
//        Index+= _CacheIndex[4];
//        return -conjugate(_Weight(Index));
//    }else{
//        //SpinIndex==5
//        Index+= _CacheIndex[4];
//        return _Weight(Index);
//    }
//}

Complex GammaWClass::Weight(const Site &Wr_in, const Site &Wr_out, const Site &Ur,
                            real Wt_in, real Wt_out, real Ut, spin* Wspin_in, spin* Wspin_out, spin Uspin) const
{
    int coord_r1 = _Map.Lat.CoordiIndex(Wr_out, Ur);
    int coord_r2 = _Map.Lat.CoordiIndex(Wr_in,  Ur);

    auto t1 = Wt_out - Ut;
    auto t2=Wt_in-Ut;
    auto SpinIndex=_SpinIndex(Wspin_out, Wspin_in);

    int SymmetryFactor;
    bool DoesMirrored;
    uint Index=_GetIndex(coord_r1, coord_r2, t1, t2, SymmetryFactor, DoesMirrored);
    if(SpinIndex>=4)
        Index+= _SpaceTimeSize;
    auto NewWeight=_Weight(Index);
    if(DoesMirrored)
        NewWeight=-conjugate(NewWeight);

    if(SpinIndex==0||SpinIndex==1) {
        return NewWeight;
    }
    else if(SpinIndex==2||SpinIndex==3) {
        return -NewWeight;
    }
    else if(SpinIndex==4){
        return -conjugate(NewWeight);
    }else{
        //SpinIndex==5
        return NewWeight;

    }
}

//void GammaWClass::Measure(const Site &Wr_in, const Site &Wr_out, const Site &Ur, real Wt_in, real Wt_out,
//                          real Ut, spin* Wspin_in, spin* Wspin_out, spin Uspin, const Complex &Weight)
//{
//    auto coord_r1 = _Map.Lat.CoordiIndex(Wr_out, Ur);
//    auto coord_r2 = _Map.Lat.CoordiIndex(Wr_in,  Ur);
//
//    auto t1=Wt_out-Ut;
//    if(t1<0){
//        t1+=_Beta;
//    }
//    int t1Index=floor(t1*_dBetaInverse);
//
//    auto t2=Wt_in-Ut;
//    if(t2<0){
//        t2+=_Beta;
//    }
//    int t2Index=floor(t2*_dBetaInverse);
//    auto SpinIndex=_SpinIndex(Wspin_out, Wspin_in);
//    uint Index = coord_r1 * _CacheIndex[3] +coord_r2*_CacheIndex[2]+ t1Index*_CacheIndex[1]+t2Index;
////    uint Index=_GetIndex(coord_r1, coord_r2, t1, t2, SymmetryFactor);
//    if(SpinIndex==0||SpinIndex==1)
//        _WeightAccu[Index]+=Weight;
//    else if(SpinIndex==2||SpinIndex==3)
//        _WeightAccu[Index]+=-Weight;
//    else if(SpinIndex==4){
//        Index+= _CacheIndex[4];
//        _WeightAccu[Index]+=Weight/2.0;
//    }
//    else if(SpinIndex==5){
//        Index+= _CacheIndex[4];
//        _WeightAccu[Index]+=-conjugate(Weight)/2.0;
//    }
//}
//
void GammaWClass::Measure(const Site &Wr_in, const Site &Wr_out, const Site &Ur, real Wt_in, real Wt_out,
                          real Ut, spin* Wspin_in, spin* Wspin_out, spin Uspin, const Complex &Weight)
{
    auto coord_r1 = _Map.Lat.CoordiIndex(Wr_out, Ur);
    auto coord_r2 = _Map.Lat.CoordiIndex(Wr_in,  Ur);

    auto t1=Wt_out-Ut;
    auto t2=Wt_in-Ut;
    auto SpinIndex=_SpinIndex(Wspin_out, Wspin_in);
//    uint Index = coord_r1 * _CacheIndex[3] +coord_r2*_CacheIndex[2]+ t1Index*_CacheIndex[1]+t2Index;
    int SymmetryFactor;
    bool DoesMirrored;
    uint Index=_GetIndex(coord_r1, coord_r2, t1, t2, SymmetryFactor, DoesMirrored);
    auto NewWeight=Weight/SymmetryFactor;
    if(DoesMirrored)
        NewWeight=-conjugate(NewWeight);
    if(SpinIndex==0||SpinIndex==1)
        _WeightAccu[Index]+=NewWeight;
    else if(SpinIndex==2||SpinIndex==3)
        _WeightAccu[Index]+=-NewWeight;
    else if(SpinIndex==4){
//        Index+= _CacheIndex[4];
        Index+=_SpaceTimeSize;
        _WeightAccu[Index]+=NewWeight/2.0;
    }
    else if(SpinIndex==5){
//        Index+= _CacheIndex[4];
        Index+=_SpaceTimeSize;
        _WeightAccu[Index]+=-conjugate(NewWeight)/2.0;
    }
}

//void GammaWClass::Measure(const Site &Wr_in, const Site &Wr_out, const Site &Ur, real Wt_in, real Wt_out,
//                          real Ut, spin* Wspin_in, spin* Wspin_out, spin Uspin, const Complex &Weight)
//{
//    auto coord_r1 = _Map.Lat.CoordiIndex(Wr_out, Ur);
//    auto coord_r2 = _Map.Lat.CoordiIndex(Wr_in,  Ur);
//
//    auto t1=Wt_out-Ut;
//    if(t1<0){
//        t1+=_Beta;
//    }
//    int t1Index=floor(t1*_dBetaInverse);
//
//    auto t2=Wt_in-Ut;
//    if(t2<0){
//        t2+=_Beta;
//    }
//    int t2Index=floor(t2*_dBetaInverse);
//    auto SpinIndex=_SpinIndex(Wspin_out, Wspin_in);
//    uint Index0 = coord_r1 * _CacheIndex[3] +coord_r2*_CacheIndex[2];
//    for(int b1=0;b1<_BasisNum/2;b1++){
//        for(int b2=0;b2<_BasisNum/2;b2++){
//            uint Index = Index0 + b1*_CacheIndex[1]+b2;
//            auto Factor=_BasisVec[b1][t1Index]*_BasisVec[b2][t2Index];
//            if(SpinIndex==0||SpinIndex==1)
//                _WeightAccu[Index]+=Weight*Factor;
//            else if(SpinIndex==2||SpinIndex==3)
//                _WeightAccu[Index]+=-Weight*Factor;
//            else if(SpinIndex==4){
//                Index+= _CacheIndex[4];
//                _WeightAccu[Index]+=Weight*Factor/2.0;
//            }
//            else if(SpinIndex==5){
//                Index+= _CacheIndex[4];
//                _WeightAccu[Index]+=-conjugate(Weight*Factor)/2.0;
//            }
//        }
//    }
//}

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

/**********************   Weight IO ****************************************/

bool GammaWClass::WeightFromDict(const Dictionary& dict) {
    _RSqueeze=dict.Get<vector<vector<int>>>("RSqueeze");
    _TauSqueeze=dict.Get<vector<vector<int>>>("TauSqueeze");
    _RSymFactor=dict.Get<vector<vector<int>>>("RSymFactor");
    _TauSymFactor=dict.Get<vector<vector<int>>>("TauSymFactor");
    _RSize=dict.Get<uint>("RSize");
    _TauSize=dict.Get<uint>("TauSize");

    uint WeightShape[5] = {2, _RSize, _TauSize};
    //Wspin12, W_r1, dW_r2, Wtau1,dtau2
    _Weight.Allocate(WeightShape, SMOOTH);
//    uint MeaShape[5] = {2, (uint)Vol, (uint)Vol, _MaxTauBin, _MaxTauBin};
//    uint MeaShape[5] = {2, _SpaceTimeSize};
    _WeightAccu.Allocate(WeightShape, SMOOTH);
    _WeightSize = _WeightAccu.GetSize();

    _Weight.FromDict(dict);
    return true;
}

bool GammaWClass::StatisFromDict(const Dictionary &dict)
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

Dictionary GammaWClass::StatisToDict()
{
    Dictionary dict;
    dict["Norm"] = _Norm;
    dict["NormAccu"] = _NormAccu;
    dict["WeightAccu"] = Python::ArrayObject(_WeightAccu.Data(), _WeightAccu.GetShape(), _WeightAccu.GetDim());
    return dict;
}

Dictionary GammaWClass::WeightToDict()
{
    return _Weight.ToDict();
}
