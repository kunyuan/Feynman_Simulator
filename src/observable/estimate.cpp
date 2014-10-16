//
//  Estimate.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "estimate.h"
#include "utility.h"
#include "cnpy.h"
#include "abort.h"

using namespace std;

template <typename T>
Estimate<T>::Estimate():Mean(),Error(){}

/**
*  Estimate constructor of type T
*
*  @param mean __reference__ of a T number as the mean value
*  @param error __reference__ of a T number as the error
*
*/
template <typename T>
Estimate<T>::Estimate(const T& mean, const T& error):Mean(mean),Error(error){}

/**
*  \brief Pretty output of Complex Estimate
*/
ostream& operator<<(ostream &os, Estimate<Complex> & e)
{
    os.setf(ios::showpoint);
    os<<"("<<e.Mean.Re<<"+/-"<<e.Error.Re<<","<<e.Mean.Im<<"+/-"<<e.Error.Im<<")"<<endl;
    return os;
}

ostream& operator<<(ostream &os, Estimate<real> & e)
{
    os.setf(ios::showpoint);
    os<<e.Mean<<"+/-"<<e.Error<<endl;
    return os;
}

template class Estimate<Complex>;
template class Estimate<real>;

/**
*  Normalization factor is initialized as 1.0 when Estimator is constructed.
*
*/
template<typename T>
Estimator<T>::Estimator(string Name)
{
    _name=Name;
    ClearStatistics();
}

template<typename T>
void Estimator<T>::ClearStatistics()
{
    _history.clear();
    _ratio=1.0;
    _norm=1.0;
}
/**
*  \brief Using statistics from ThrowRatio*100% to 100% to estimate the error bar
*/
const real ThrowRatio=1.0/3;
template<>
void Estimator<real>::_update()
{
    int size=(int)_history.size();
    if(size==0) return;
    real Min=MaxReal,Max=MinReal;
    int MinIndex=0, MaxIndex=0;
    for(int i=size*ThrowRatio;i<size;i++)
    {
        if(Min>_history[i])
        {
            Min=_history[i];
            MinIndex=i;
        }
        if(Max<_history[i])
        {
            Max=_history[i];
            MaxIndex=i;
        }
    }
    _value.Error=abs(Max-Min)/2.0;
    _value.Mean=_accumulator/_norm;
    _ratio=(MaxIndex-MinIndex)/(real)size*(1.0-ThrowRatio);
}

template<>
void Estimator<Complex>::_update()
{
    int size=(int)_history.size();
    if(size==0) return;
    Complex Min(MaxReal,MaxReal),Max(MinReal,MinReal);
    int MinIndexRe=0, MaxIndexRe=0;
    int MinIndexIm=0, MaxIndexIm=0;
    for(int i=size*ThrowRatio;i<size;i++)
    {
        if(Min.Re>_history[i].Re)
        {
            Min.Re=_history[i].Re;
            MinIndexRe=i;
        }
        if(Max.Re<_history[i].Re)
        {
            Max.Re=_history[i].Re;
            MaxIndexRe=i;
        }
        if(Min.Im>_history[i].Im)
        {
            Min.Im=_history[i].Im;
            MinIndexIm=i;
        }
        if(Max.Im<_history[i].Im)
        {
            Max.Im=_history[i].Im;
            MaxIndexIm=i;
        }
    }
    _value.Mean=_accumulator/_norm;
    _value.Error.Re=abs(Max.Re-Min.Re)/2.0;
    _value.Error.Im=abs(Max.Im-Min.Im)/2.0;
    if(MaxIndexRe-MinIndexRe<MaxIndexIm-MinIndexIm)
        _ratio=(MaxIndexIm-MinIndexIm)/(real)size*(1.0-ThrowRatio);
    else
        _ratio=(MaxIndexRe-MinIndexRe)/(real)size*(1.0-ThrowRatio);
}

template <typename T>
void Estimator<T>::AddStatistics(const T& t)
{
    _accumulator+=t;
    _norm+=1.0;
    _history.push_back(_accumulator/_norm);
}

template <typename T>
Estimate<T> Estimator<T>::Estimate()
{
    _update();
    return _value;
}

template <typename T>
real Estimator<T>::Ratio()
{
    return _ratio;
}

template <typename T>
bool Estimator<T>::ReadState(cnpy::npz_t NpzMap)
{
    cnpy::NpyArray history=NpzMap[_name];
    T* start = reinterpret_cast<T*>(history.data);
    if(start==NULL) ABORT("Can't find estimator "<<_name<<" in .npz data file!"<<endl);
    ClearStatistics();
    size_t size=history.shape[0];
    _history.assign(start, start+size);
    _norm=real(size+1);
    _accumulator=_history[size-1]*_norm;
    _update();
    return true;
}

template <typename T>
void Estimator<T>::SaveState(const string FileName, string Mode)
{
    unsigned int shape[1];
    shape[0]=(unsigned int)_history.size();
    //!!!Assume _norm==_history.size()+1 here, so don't have to store _norm
    cnpy::npz_save(cnpy::npz_name(FileName),_name,_history.data(),shape,1,Mode);
}

template class Estimator<real>;
template class Estimator<Complex>;

template <typename T>
bool EstimatorVector<T>::ReadState(const string FileName)
{
    cnpy::npz_t NpzMap=cnpy::npz_load(cnpy::npz_name(FileName));
    for(unsigned int i=0;i<vector<Estimator<T>>::size();i++)
    {
        vector<Estimator<T>>::at(i).ReadState(NpzMap);
    }
    return true;
}

template <typename T>
void EstimatorVector<T>::SaveState(const string FileName, string Mode)
{
    string Mod=Mode;
    for(unsigned int i=0;i<vector<Estimator<T>>::size();i++)
    {
        vector<Estimator<T>>::at(i).SaveState(FileName, Mod);
        if(i==0&&Mod=="w") Mod="a"; //the second and the rest elements will be wrote as appended
    }
}
/**
*  \brief clear all statistics of the elements in the EstimatorVector.
*   __memory__ of the vector will not be freed!
*
*/
template <typename T>
void EstimatorVector<T>::ClearStatistics()
{
    for(unsigned int i=0;i<vector<Estimator<T>>::size();i++)
    {
        vector<Estimator<T>>::at(i).ClearStatistics();
    }
}

template class EstimatorVector<Complex>;
template class EstimatorVector<real>;