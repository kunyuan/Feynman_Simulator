//
//  fft.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 10/23/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef Feynman_Simulator_fft_h
#define Feynman_Simulator_fft_h

#include "complex.h"
#include "convention.h"

void fft(Complex x[], int n, int flag);
void fft2D(Complex x[], int n1, int n2, int flag);
void fft3D(Complex x[], int n1, int n2, int n3, int flag);

int TestFFT();
#endif
