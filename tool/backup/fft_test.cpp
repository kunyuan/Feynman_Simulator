//
//  fft_test.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/23/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "fft.h"
#include "sput.h"
#include "complex.h"
#include <iostream>
using namespace fft;

void Test_fft1D();
void Test_fft2D();
void Test_fft3D();
void Test_fftnD();
void NaiveFFT3D(Complex *array, int X, int Y, int Z, bool DoItX, bool DoItY, bool DotItZ);

int TestFFT()
{
    sput_start_testing();
    sput_enter_suite("Test fft...");
    sput_run_test(Test_fft1D);

    //    sput_enter_suite("Test fft2D...");
    //    sput_run_test(Test_fft2D);
    //
    sput_run_test(Test_fft3D);
    sput_run_test(Test_fftnD);
    sput_finish_testing();
    return sput_get_return_value();
    return 0;
}

bool CheckArray(Complex *source, Complex *target, int size)
{
    for (int i = 0; i < size; i++) {
        if (!Equal(source[i], target[i], 1e-5)) {
            std::cout << source[i] << "!=" << target[i] << std::endl;
            return false;
        }
    }
    return true;
}

void OutputAll(Complex *data, int size)
{
    for (int i = 0; i < size; i++) {
        std::cout.precision(15);
        std::cout << std::fixed << data[i] << std::endl;
    }
}

void Test_fft1D()
{
    Complex in[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    Complex out[8] =
        {{28.000000000000000, 0.000000000000000},
         {-4.000000000000000, 9.656854249492380},
         {-4.000000000000000, 4.000000000000000},
         {-4.000000000000000, 1.656854249492380},
         {-4.000000000000000, 0.000000000000000},
         {-4.000000000000000, -1.656854249492380},
         {-4.000000000000000, -4.000000000000000},
         {-3.999999999999999, -9.656854249492380}};
    fft::fft(in, 8, FORTH);
    sput_fail_unless(CheckArray(in, out, 8), "Check fft 1d");
}

void Test_fft3D()
{
    const int Nx = 4;
    const int Ny = 4;
    const int Nz = 4;
    unsigned int shape[] = {Nx, Ny, Nz};
    Complex in[Nx][Ny][Nz];
    Complex out3D[Nx][Ny][Nz];
    for (int i = 0;
         i < shape[0]; i++)
        for (int j = 0; j < shape[1]; j++)
            for (int k = 0; k < shape[2]; k++) {
                in[i][j][k] = (i * i + j * j) * cos(k);
                out3D[i][j][k] = in[i][j][k];
            }
    NaiveFFT3D((Complex *)out3D, Nx, Ny, Nz, true, true, true);

    fft3D((Complex *)in, Nx, Ny, Nz, FORTH);
    sput_fail_unless(CheckArray((Complex *)in, (Complex *)out3D, Nx * Ny * Nz), "Check fft 3d");
    fft3D((Complex *)in, Nx, Ny, Nz, BACK);
    fft3D((Complex *)in, Nx, Ny, Nz, FORTH);
    sput_fail_unless(CheckArray((Complex *)in, (Complex *)out3D, Nx * Ny * Nz), "Check fft 3d backforth");
}

void Test_fftnD()
{
    const int Nx = 4;
    const int Ny = 8;
    const int Nz = 16;
    bool DoX = true;
    bool DoY = false;
    bool DoZ = true;
    unsigned int shape[] = {Nx, Ny, Nz};
    bool mask[3] = {DoX, DoY, DoZ};
    Complex in[Nx][Ny][Nz];
    Complex out3D[Nx][Ny][Nz];
    for (int i = 0;
         i < shape[0]; i++)
        for (int j = 0; j < shape[1]; j++)
            for (int k = 0; k < shape[2]; k++) {
                in[i][j][k] = (i * i + j * j) * cos(k);
                out3D[i][j][k] = in[i][j][k];
            }
    NaiveFFT3D((Complex *)out3D, Nx, Ny, Nz, DoX, DoY, DoZ);

    fftnD((Complex *)in, (int *)shape, 3, FORTH, (bool *)mask);
    sput_fail_unless(CheckArray((Complex *)in, (Complex *)out3D, Nx * Ny * Nz), "Check fft nd");
    fftnD((Complex *)in, (int *)shape, 3, BACK, (bool *)mask);
    fftnD((Complex *)in, (int *)shape, 3, FORTH, (bool *)mask);
    sput_fail_unless(CheckArray((Complex *)in, (Complex *)out3D, Nx * Ny * Nz), "Check fft nd backforth");
}

void NaiveFFT3D(Complex *array, int X, int Y, int Z, bool DoItX, bool DoItY, bool DoItZ)
{
    Complex *temp = new Complex[X * Y * Z];
    int windex;
    int index;

    for (int i = 0; i < X * Y * Z; i++)
        temp[i] = array[i];
    for (int iw = 0; iw < X; iw++)
        for (int jw = 0; jw < Y; jw++)
            for (int kw = 0; kw < Z; kw++) {
                windex = iw * Y * Z + jw * Z + kw;
                array[windex] = 0.0;
                Complex factorX;
                Complex factorY;
                Complex factorZ;
                for (int i = 0; i < X; i++) {
                    if (DoItX)
                        factorX = exp(-2 * PI * Complex(0.0, iw * i) / X);
                    else
                        factorX = (i == iw ? 1.0 : 0.0);
                    for (int j = 0; j < Y; j++) {
                        if (DoItY)
                            factorY = exp(-2 * PI * Complex(0.0, jw * j) / Y);
                        else
                            factorY = (j == jw ? 1.0 : 0.0);
                        for (int k = 0; k < Z; k++) {
                            if (DoItZ)
                                factorZ = exp(-2 * PI * Complex(0.0, kw * k) / Z);
                            else
                                factorZ = (k == kw ? 1.0 : 0.0);
                            index = i * Y * Z + j * Z + k;
                            array[windex] += temp[index] * factorX * factorY * factorZ;
                        }
                    }
                }
            }
    delete[] temp;
}
