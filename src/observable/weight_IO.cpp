//
//  weight_IO.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/20/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "weight.h"
#include "abort.h"

using namespace std;
using namespace Array;
using namespace Weight;

/**********************   Weight IO ****************************************/
void Sigma::SaveState(const std::string &FileName, const std::string &Mode)
{
    unsigned int shape[1] = {1};
    cnpy::npz_save(cnpy::npz_name(FileName), _Name + "_Norm", &_Norm, shape, 1, Mode);
    cnpy::npz_save(cnpy::npz_name(FileName), _Name + "_Accu", _WeightAccu->Data(), _Shape, 5, "a");
    unsigned int shape2[2];
    shape2[0] = _Shape[ORDER];
    shape2[1] = _Shape[TAU];
    cnpy::npz_save(cnpy::npz_name(FileName), _Name + "_SquareAccu", _WeightSquareAccu->Data(), shape2, 2, "a");
}

bool Sigma::LoadState(const std::string &FileName)
{
    cnpy::npz_t NpzMap = cnpy::npz_load(cnpy::npz_name(FileName));

    cnpy::NpyArray sigma_accu = NpzMap[_Name + "_Accu"];
    Complex *start = reinterpret_cast<Complex *>(sigma_accu.data);
    if (start == NULL)
        ABORT("Can't find estimator " << _Name << " _Accu in .npz data file!" << endl);
    _WeightAccu->Set(start);

    cnpy::NpyArray sigma_squareaccu = NpzMap[_Name + "_SquareAccu"];
    start = reinterpret_cast<Complex *>(sigma_squareaccu.data);
    if (start == NULL)
        ABORT("Can't find estimator " << _Name << " _Accu in .npz data file!" << endl);
    _WeightSquareAccu->Set(start);

    //read normalization factor
    cnpy::NpyArray norm = NpzMap[_Name + "_Norm"];
    real *start_Norm = reinterpret_cast<real *>(norm.data);
    if (start_Norm == NULL)
        ABORT("Can't find estimator " << _Name << "_Norm in .npz data file!" << endl);
    _Norm = *start_Norm;

    NpzMap.destruct();
    return true;
}