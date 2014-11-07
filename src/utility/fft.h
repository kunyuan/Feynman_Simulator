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
    FORTH,
    BACK
};

void fft(Complex *x, int n, Dir, bool *Mask=nullptr);
void fft2D(Complex *x, int n1, int n2, Dir, bool *Mask = nullptr);
void fft3D(Complex *x, int n1, int n2, int n3, Dir, bool *Mask = nullptr);
void fft4D(Complex *x, int n1, int n2, int n3, int n4, Dir, bool *Mask = nullptr);

void fft(Complex *x, int *size, int dim, Dir, bool *Mask = nullptr);
}

int TestFFT();
#endif
