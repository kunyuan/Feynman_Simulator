//
//  Estimate.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__Estimate__
#define __Feynman_Simulator__Estimate__

#include <iostream>
#include <vector>
#include "component.h"
#include "cnpy.h"
/**
*  \brief estimate with mean value and standard error
*/
template <typename T>
class Estimate
{
public:
    Estimate();
    Estimate(const T &m, const T &e);
    T Mean;
    T Error;
    friend std::ostream& operator<<(std::ostream &, Estimate &);
};

/**
*  \brief maintain a history of observables, where an Estimate object can be calculated
*/

template <typename T>
class Estimator
{
private:
    std::vector<T> _history;
    std::string _name;
    T _accumulator;
    real _norm;
    real _ratio;
    Estimate<T> _value;
    void _update();
public:
    Estimator(std::string);
    void AddStatistics(const T&);
    Estimate<T> Estimate();
    double Ratio();
    bool ReadState(cnpy::npz_t);
    void SaveState(const std::string FileName, const std::string Mode="a");
    void ClearStatistics();
};

/**
*  \brief an extension to the std::vector<Estimator<T>> to contain Estimator<T>
*/

template <typename T>
class EstimatorVector: public std::vector<Estimator<T>>
{
public:
    bool ReadState(const std::string FileName);
    void SaveState(const std::string FileName, const std::string Mode="a");
    void ClearStatistics();
};

#endif /* defined(__Feynman_Simulator__Estimate__) */
