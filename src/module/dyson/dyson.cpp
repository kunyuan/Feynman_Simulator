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

/**
 *  G = 1/(G_0^-1-Sigma)
 */
void Dyson::DeriveG()
{
    unsigned int *GShape = G->Shape();
    Array::array4<Complex> &GSmooth = G->SmoothWeight;
    
    Array::array1<Complex> DiscretCorrection(GShape[TAU]);
    for (int omega = 0; omega < GShape[TAU]; omega++)
        DiscretCorrection[omega] = cos((omega + 0.5) * PI / GShape[TAU]);
    
    int SpaceTimeVol = GShape[VOL] * GShape[TAU];

    GSmooth = G->BareWeight;
    for (int sp = 0; sp < GShape[SP]; sp++)
        if (G->IsSameSpin(sp)) {
            MatrixInverseSUB(GSmooth[sp](), SpaceTimeVol);
            GSmooth[sp] += Sigma->SmoothWeight[sp];
            for (int sub = 0; sub < GShape[SUB]; sub++)
                for (int k = 0; k < GShape[VOL]; k++)
                    GSmooth[sp][sub][k] += DiscretCorrection * Sigma->DeltaTWeight[sp][sub][k];
            MatrixInverseSUB(GSmooth[sp](), SpaceTimeVol);
        }
        else
            GSmooth[sp] = 0.0;
}

/**
 *  need refactoring
 *  W = J*Polar/(J^-1 - Polar)
 */
void Dyson::DeriveW()
{
    unsigned int *WShape = W->Shape();
    Array::array4<Complex> &WSmooth = W->SmoothWeight;
    Array::array3<Complex> WBareInverse = W->BareWeight;
    
    Array::array1<Complex> DiscretCorrection(WShape[TAU]);
    for (int omega = 0; omega < WShape[TAU]; omega++)
        DiscretCorrection[omega] = -cos((omega)*PI/WShape[TAU]);

    int SpaceTimeVol = WShape[VOL] * WShape[TAU];
    WSmooth = Polar->SmoothWeight;
    
    for(int sp=0; sp<WShape[SP]; sp++)
        MatrixInverseSUB(WBareInverse[sp](), WShape[VOL]);
    MatrixInverseSP(WBareInverse(), WShape[SUB]*WShape[VOL]);
    
    for(int sp=0; sp<WShape[SP]; sp++)
        for(int sub = 0; sub < WShape[SUB]; sub++)
            for(int k = 0; k < WShape[VOL]; k++)
                for(int omega = 0; omega < WShape[TAU]; omega++)
                {
                    WSmooth[sp][sub][k][omega] *= DiscretCorrection[omega];
                    WSmooth[sp][sub][k][omega] += WBareInverse[sp][sub][k];
                }
    
    for(int sp=0; sp<WShape[SP]; sp++)
        MatrixInverseSUB(WSmooth[sp](), SpaceTimeVol);
    MatrixInverseSP(WSmooth(), WShape[SUB]*SpaceTimeVol);
    
    for(int sp=0; sp<WShape[SP]; sp++)
        MatrixMultiplySUB(WSmooth[sp](), Polar->SmoothWeight[sp](), SpaceTimeVol, SpaceTimeVol);
    MatrixMultiplySP(WSmooth(), W->BareWeight(), WShape[SUB]*SpaceTimeVol, WShape[SUB]*SpaceTimeVol);
    
    for(int sp=0; sp<WShape[SP]; sp++)
        MatrixMultiplySUB(WSmooth[sp](), W->BareWeight[sp](), SpaceTimeVol, WShape[VOL]);
    MatrixMultiplySP(WSmooth(), W->BareWeight(), WShape[SUB]*SpaceTimeVol, WShape[SUB]*WShape[VOL]);
}

void Dyson::DeriveChi()
{
//    unsigned int *WShape = W->Shape();
//    Array::array4<Complex> WSmooth = ;
    
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
void dyson::MatrixInverseSUB(Complex *matrix, int interval)
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

void dyson::MatrixMultiplySUB(Complex *matrix1, Complex *matrix2, int interval)
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
void dyson::MatrixMultiplySUB(Complex *matrix1, Complex *matrix2, int interval1, int interval2)
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

/**
 *  4*4 matrix inverse
 *
 *  @param matrix   4*4 matrix
 *  @param interval the inner dimensions
 */
void dyson::MatrixInverseSP(Complex *matrix, int interval)
{
    Complex denom1, tmp1, denom2, tmp2;
    for(int i=0; i<interval; i++)
    {
        denom1 = matrix[i]*matrix[i+5*interval]-matrix[i+interval]*matrix[i+4*interval];
        tmp1 = matrix[i];
        matrix[i] = matrix[i+5*interval]/denom1;
        matrix[i+5*interval] = tmp1/denom1;
        matrix[i+interval] /= (-denom1);
        matrix[i+4*interval] /= (-denom1);
        
        denom2 = matrix[i+10*interval]*matrix[i+15*interval]-matrix[i+11*interval]*matrix[i+14*interval];
        tmp2 = matrix[i+10*interval];
        matrix[i+10*interval] = matrix[i+15*interval]/denom2;
        matrix[i+15*interval] = tmp2/denom2;
        matrix[i+11*interval] /= (-denom2);
        matrix[i+14*interval] /= (-denom2);
    }
    
}

void dyson::MatrixMultiplySP(Complex *matrix1, Complex *matrix2, int interval)
{
    Complex tmp[16]={0.0};
    for(int i=0; i<interval; i++)
    {
        for(int k=0; k<16; k++)
        {
            tmp[k] = matrix1[i+k*interval];
            matrix1[i+k*interval] = 0.0;
        }
        
        matrix1[i] = tmp[0]*matrix2[i]
                    + tmp[1]*matrix2[i+4*interval];
        matrix1[i+interval]   = tmp[0]*matrix2[i+interval]
                              + tmp[1]*matrix2[i+5*interval];
        matrix1[i+4*interval] = tmp[4]*matrix2[i]
                              + tmp[5]*matrix2[i+4*interval];
        matrix1[i+5*interval] = tmp[4]*matrix2[i+interval]
                              + tmp[5]*matrix2[i+5*interval];
        
        matrix1[i+10*interval]= tmp[10]*matrix2[i+10*interval]
                               +tmp[11]*matrix2[i+14*interval];
        matrix1[i+11*interval]= tmp[10]*matrix2[i+11*interval]
                               +tmp[11]*matrix2[i+15*interval];
        matrix1[i+14*interval]= tmp[14]*matrix2[i+10*interval]
                               +tmp[15]*matrix2[i+14*interval];
        matrix1[i+15*interval]= tmp[14]*matrix2[i+11*interval]
                               +tmp[15]*matrix2[i+15*interval];
    }
}

void dyson::MatrixMultiplySP(Complex *matrix1, Complex *matrix2, int interval1, int interval2)
{
    if (interval1 < interval2 || interval1 % interval2 != 0)
        return;
    Complex tmp[16]={0.0};
    int step = interval1/interval2;
    int j;
    
    for(int i=0; i<interval1; i++)
    {
        j = i / step;
        for(int k=0; k<16; k++)
        {
            tmp[k] = matrix1[i+k*interval1];
            matrix1[i+k*interval1] = 0.0;
        }
        
        matrix1[i] = tmp[0]*matrix2[j]
                    + tmp[1]*matrix2[j+4*interval2];
        matrix1[i+interval1] =   tmp[0]*matrix2[j+interval2]
                                +tmp[1]*matrix2[j+5*interval2];
        matrix1[i+4*interval1] = tmp[4]*matrix2[j]
                                +tmp[5]*matrix2[j+4*interval2];
        matrix1[i+5*interval1] = tmp[4]*matrix2[j+interval2]
                                +tmp[5]*matrix2[j+5*interval2];
        
        matrix1[i+10*interval1]= tmp[10]*matrix2[j+10*interval2]
                                +tmp[11]*matrix2[j+14*interval2];
        matrix1[i+11*interval1]= tmp[10]*matrix2[j+11*interval2]
                                +tmp[11]*matrix2[j+15*interval2];
        matrix1[i+14*interval1]= tmp[14]*matrix2[j+10*interval2]
                                +tmp[15]*matrix2[j+14*interval2];
        matrix1[i+15*interval1]= tmp[14]*matrix2[j+11*interval2]
                                +tmp[15]*matrix2[j+15*interval2];
    }
}