#include "complex.h"
#include <complex>
#include "timer.h"
#include <iostream>
#include "fft.h"

void complex_benchmark();
void fft_benchmark();

int main(int argc, const char *argv[])
{
    //complex_benchmark();
    fft_benchmark();
    return 0;
}

void fft_benchmark()
{
    const int X=1024;
    const int Y=32;
    const int num=1000;

    Complex a[X][Y];
    Complex b[X][Y];
    Complex c[X][Y];
    timer T;
    T.start();
    for(int t=0;t<num;t++)
    {
        for(int i=0;i<X;i++)
            for(int j=0;j<Y;j++)
                a[i][j]=Complex(i,j);

        cooley_tukey((Complex *)a, X*Y, 1, Y);
    }
    T.stop();
    std::cout<<T<<std::endl;

    T.restart();
    for(int t=0;t<num;t++)
    {
        for(int i=0;i<X;i++)
            for(int j=0;j<Y;j++)
                b[i][j]=Complex(i,j);

        stockham((Complex *)b, X*Y, 1, Y, (Complex*) c);
    }
    T.stop();
    std::cout<<T<<std::endl;


    //for(int i=0;i<X;i++)
        //for(int j=0;j<Y;j++)
            //if(!Equal(a[i][j],b[i][j]/X, 1e-6))
            //{
                //std::cout<<"wrong calling fft!"<<std::endl;
                //std::cout<<a[i][j]<<", "<<b[i][j]<<std::endl;
                //return;
            //}
}

void complex_benchmark()
{
    timer T;
    T.start();
    const int NUM=100000000;
    Complex c;
    for(int i=0;i<NUM;i++)
    {
        Complex a=Complex(1.0,2*i);
        Complex a0=Complex(3.0,2*i);
        c+=arg(a/a0);
    }
    T.stop();
    std::cout<<T<<", "<<c<<std::endl;

    T.restart();
    std::complex<double> d;
    for(int i=0;i<NUM;i++)
    {
        std::complex<double> b=std::complex<double>(1.0,2*i);
        std::complex<double> b0=std::complex<double>(3.0,2*i);
        d+=std::arg(b/b0);
    }
    T.stop();
    std::cout<<T<<",  "<<d<<std::endl;
}


