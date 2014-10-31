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
void fft(Complex x[], int n, int flag);
void fft2D(Complex x[], int n1, int n2, int flag);
void fft3D(Complex x[], int n1, int n2, int n3, int flag);
void fft4D(Complex x[], int n1, int n2, int n3, int n4, int flag);

void fft(Complex x[], int *size, int dim, int flag);

int TestFFT();
#endif
