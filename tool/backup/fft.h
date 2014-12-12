//
//  fft.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/23/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef Feynman_Simulator_fft_h
#define Feynman_Simulator_fft_h

class Complex;
namespace fft {
enum Dir {
    FORTH = 1,
    BACK = -1
};

void fft(Complex *x, int n, Dir);
void fft2D(Complex *x, int n1, int n2, Dir);
void fft3D(Complex *x, int n1, int n2, int n3, Dir);
void fftnD(Complex *x, int *size, int dim, Dir, bool *Mask = nullptr);
}

void stockham(Complex x[], int n, int flag, int n2, Complex y[]);
void cooley_tukey(Complex x[], int n, int flag, int n2);
int TestFFT();
#endif
