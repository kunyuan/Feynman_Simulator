//
//  Estimate.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__Estimate__
#define __Feynman_Simulator__Estimate__

#include <iosfwd>
#include <vector>
#include <unordered_map>
#include "utility/complex.h"

#define SIZE 10000

class Dictionary;
/**
*  \brief estimate with mean value and standard error
*/
template <typename T>
class EstimateClass {
public:
    EstimateClass();
    EstimateClass(const T& m, const T& e);
    T Mean;
    T Error;
    T RelativeError();
    friend std::ostream& operator<<(std::ostream&, const EstimateClass&);
};

/**
*  \brief maintain a history of observables, where an Estimate object can be calculated
*/

template <typename T>
class Estimator {
private:
    uint _interval;
    uint _counter;
    uint _end;
    T _accumulator;
    real _norm;
    real _ratio;
    EstimateClass<T> _value;
    void _update();
    T _history[SIZE];

public:
    Estimator();
    Estimator(std::string);
    std::string Name;
    void Measure(const T&);
    void AddStatistics();
    T Value();
    real Norm();
    double Ratio();
    EstimateClass<T> Estimate();
    bool FromDict(const Dictionary&);
    Dictionary ToDict();
    void ClearStatistics();
    void SqueezeStatistics(real factor);
};

/**
*  \brief an extension to the std::vector<Estimator<T>> to contain Estimator<T>
*/

template <typename T>
class EstimatorBundle {
private:
    typedef Estimator<T> EstimatorT;
    std::vector<EstimatorT> _EstimatorVector;
    std::unordered_map<std::string, unsigned long> _EstimatorMap;
    bool _MakeSureKeyNotExists(std::string key);

public:
    void AddEstimator(const std::string);
    void AddEstimator(const EstimatorT&);
    void RemoveEstimator(const std::string);
    void AddStatistics();
    int HowMany();
    //if AllowFailure is set to true, then FromDict will not throw exception even if some of the estimators fail to read
    bool FromDict(const Dictionary&, bool AllowFailure = false);
    Dictionary ToDict();
    EstimatorT& operator[](int);
    EstimatorT& operator[](std::string);
    void ClearStatistics();
    void SqueezeStatistics(real factor);
};

int TestEstimator();

#endif /* defined(__Feynman_Simulator__Estimate__) */
