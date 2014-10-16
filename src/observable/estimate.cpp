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
    _value.Error=fabs(Max-Min)/2.0;
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
    _value.Error.Re=fabs(Max.Re-Min.Re)/2.0;
    _value.Error.Im=fabs(Max.Im-Min.Im)/2.0;
    if(MaxIndexRe-MinIndexRe<MaxIndexIm-MinIndexIm)
        _ratio=(MaxIndexIm-MinIndexIm)/(real)size*(1.0-ThrowRatio);
    else
        _ratio=(MaxIndexRe-MinIndexRe)/(real)size*(1.0-ThrowRatio);
}

template <typename T>
void Estimator<T>::Measure(const T& t)
{
    _accumulator+=t;
    _norm+=1.0;
}

template <typename T>
void Estimator<T>::AddStatistics()
{
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
    _history.assign(start,start+history.shape[0]);
    //read normalization factor
    cnpy::NpyArray norm=NpzMap[_name+"_Norm"];
    real* start_Norm = reinterpret_cast<real*>(norm.data);
    if(start_Norm==NULL) ABORT("Can't find estimator "<<_name<<"_Norm in .npz data file!"<<endl);
    _norm=*start_Norm;
    cnpy::NpyArray accu=NpzMap[_name+"_Accu"];
    T* start_accu = reinterpret_cast<T*>(accu.data);
    if(start_accu==NULL) ABORT("Can't find estimator "<<_name<<"_Accu in .npz data file!"<<endl);
    _accumulator=*start_accu;
    _update();
    return true;
}

template <typename T>
void Estimator<T>::SaveState(const string FileName, string Mode)
{
    unsigned int shape[1];
    shape[0]=(unsigned int)_history.size();
    cnpy::npz_save(cnpy::npz_name(FileName),_name,_history.data(),shape,1,Mode);
    shape[0]=1;
    cnpy::npz_save(cnpy::npz_name(FileName),_name+"_Norm",&_norm,shape,1,"a");
    cnpy::npz_save(cnpy::npz_name(FileName),_name+"_Accu",&_accumulator,shape,1,"a");
}

template class Estimator<real>;
template class Estimator<Complex>;

template <typename T>
void EstimatorVector<T>::AddEstimator(string name)
{
    vector<Estimator<T>>::push_back(Estimator<T>(name));
}

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

/**
*  \brief This insertion is perferred than other ways provided by unorder_map class, because it
*   takes care of Estimator<T>._name automatically.
*  @param name The name of the estimator
*/
template <typename T>
void EstimatorMap<T>::AddEstimator(string name)
{
//    Estimator<T> NewEstimator(name);
//    unordered_map<string, Estimator<T>>::insert(make_pair(name, NewEstimator));
}

template <>
void EstimatorMap<Complex>::AddEstimator(string name)
{
    Estimator<Complex> NewEstimator(name);
    unordered_map<string, Estimator<Complex>>::insert(make_pair(name, NewEstimator));
}
template <typename T>
bool EstimatorMap<T>::ReadState(const string FileName)
{
    cnpy::npz_t NpzMap=cnpy::npz_load(cnpy::npz_name(FileName));
    auto end=unordered_map<string,Estimator<T>>::end();
    for(auto it=unordered_map<string,Estimator<T>>::begin();it!=end;it++)
    {
        it->second.ReadState(NpzMap);
    }
    return true;
}
template class EstimatorMap<Complex>;
template class EstimatorMap<real>;