//
//  dyson.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/4/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "dyson.h"
#include "module/parameter/parameter.h"
#include "module/weight/weight.h"
#include "module/weight/weight_inherit.h"
#include "utility/utility.h"

using namespace dyson;
using namespace weight;

bool Dyson::BuildNew(para::ParaDyson &para, Weight &weight)
{
    Beta = para.Beta;

    G = weight.G;
    W = weight.W;
    Sigma = weight.Sigma;
    Polar = weight.Polar;
    return true;
}

void Dyson::DeriveG()
{
    unsigned int *GShape = G->Shape();
    Array::array4<Complex> &GSmooth = G->SmoothWeight;
    GSmooth = G->BareWeight; ///TODO: need test
    Array::array1<Complex> DiscretCorrection(GShape[TAU]);
    for (int omega = 0; omega < GShape[TAU]; omega++)
        DiscretCorrection[omega] = cos((omega + 0.5) * PI / GShape[TAU]);

    for (int sp = 0; sp < GShape[SP]; sp++)
        if (G->IsSameSpin(sp)) {
            MatrixInverse(GSmooth[sp](), GShape[VOL] * GShape[TAU]);
            GSmooth[sp] += Sigma->SmoothWeight[sp];
            for (int sub = 0; sub < GShape[SUB]; sub++)
                for (int k = 0; k < GShape[VOL]; k++)
                    GSmooth[sp][sub][k] += DiscretCorrection * Sigma->DeltaTWeight[sp][sub][k];
            MatrixInverse(GSmooth[sp](), GShape[VOL] * GShape[TAU]);
        }
        else
            GSmooth[sp] = 0.0; ///TODO:need test
}

/**
*  Polar's value is untouched!!!
*/
void Dyson::DeriveW()
{
    unsigned int *WShape = W->Shape();
    Array::array4<Complex> &WSmooth = W->SmoothWeight;
    WSmooth = Polar->SmoothWeight;
    Array::array1<Complex> DiscretCorrection(WShape[TAU]);
    for (int omega = 0; omega < WShape[TAU]; omega++)
        DiscretCorrection[omega] = cos((omega)*PI / WShape[TAU]);

    int SpaceTimeVol = WShape[VOL] * WShape[TAU];
    for (int sp = 0; sp < WShape[SP]; sp++) {
        MatrixMultiply(WSmooth[sp](), W->BareWeight[sp](), SpaceTimeVol, WShape[VOL]);
        ///TODO: need test
        MatrixInverse(WSmooth[sp](), WShape[VOL] * WShape[TAU]);
        for (int sub = 0; sub < WShape[SUB]; sub++)
            for (int k = 0; k < WShape[VOL]; k++)
                WSmooth[sp][sub][k] -= DiscretCorrection;
        MatrixInverse(WSmooth[sp](), WShape[VOL] * WShape[TAU]);
        MatrixMultiply(WSmooth[sp](), W->BareWeight[sp](), SpaceTimeVol, WShape[VOL]);
    }
}

/**
 *  2*2*INTERVAL matrix inverse ( 0  1
 *                                2  3)
 *
 *  @param Complex* original matrix
 *
 *  @return  the inversed matrix
 */
#define POS0(x) (0)
#define POS1(x) (x)
#define POS2(x) (x + x)
#define POS3(x) (x + x + x)
void dyson::MatrixInverse(Complex *matrix, int interval)
{
    for (int i = 0; i < interval; i++) {
        Complex inverse_denominator = 1.0 / (matrix[i + POS0(interval)] * matrix[i + POS3(interval)] - matrix[i + POS1(interval)] * matrix[i + POS2(interval)]);
        Complex tmp = matrix[i + POS3(interval)] * inverse_denominator;
        matrix[i + POS3(interval)] = matrix[i + POS0(interval)] * inverse_denominator;
        matrix[i + POS0(interval)] = tmp;
        matrix[i + POS1(interval)] *= -inverse_denominator;
        matrix[i + POS2(interval)] *= -inverse_denominator;
    }
}

void dyson::MatrixMultiply(Complex *matrix1, Complex *matrix2, int interval)
{
    Complex tmp[4];
    for (int i = 0; i < interval; i++) {
        tmp[0] = matrix1[i + POS0(interval)];
        tmp[1] = matrix1[i + POS1(interval)];
        tmp[2] = matrix1[i + POS2(interval)];
        tmp[3] = matrix1[i + POS3(interval)];

        matrix1[i + POS0(interval)] = tmp[0] * matrix2[i + POS0(interval)] + tmp[1] * matrix2[i + POS2(interval)];
        matrix1[i + POS1(interval)] = tmp[0] * matrix2[i + POS1(interval)] + tmp[1] * matrix2[i + POS3(interval)];
        matrix1[i + POS2(interval)] = tmp[2] * matrix2[i + POS0(interval)] + tmp[3] * matrix2[i + POS2(interval)];
        matrix1[i + POS3(interval)] = tmp[2] * matrix2[i + POS1(interval)] + tmp[3] * matrix2[i + POS3(interval)];
    }
}

//interval2 is smaller than interval1
void dyson::MatrixMultiply(Complex *matrix1, Complex *matrix2, int interval1, int interval2)
{
    if (interval1 < interval2 || interval1 % interval2 != 0)
        return;
    int step = interval1 / interval2;
    int j;

    Complex tmp[4];
    for (int i = 0; i < interval1; i++) {
        j = i / step;
        tmp[0] = matrix1[i + POS0(interval1)];
        tmp[1] = matrix1[i + POS1(interval1)];
        tmp[2] = matrix1[i + POS2(interval1)];
        tmp[3] = matrix1[i + POS3(interval1)];

        matrix1[i + POS0(interval1)] = tmp[0] * matrix2[j + POS0(interval2)] + tmp[1] * matrix2[j + POS2(interval2)];
        matrix1[i + POS1(interval1)] = tmp[0] * matrix2[j + POS1(interval2)] + tmp[1] * matrix2[j + POS3(interval2)];
        matrix1[i + POS2(interval1)] = tmp[2] * matrix2[j + POS0(interval2)] + tmp[3] * matrix2[j + POS2(interval2)];
        matrix1[i + POS3(interval1)] = tmp[2] * matrix2[j + POS1(interval2)] + tmp[3] * matrix2[j + POS3(interval2)];
    }
}