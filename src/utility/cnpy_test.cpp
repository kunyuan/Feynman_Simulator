#include "cnpy.h"
#include "complex.h"
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "rng.h"

const int Nx = 128;
const int Ny = 64;
const int Nz = 32;

void Testcnpy()
{
    RandomFactory RNG;
    //create random data
    Complex *data = new Complex[Nx * Ny * Nz];
    for (int i = 0; i < Nx * Ny * Nz; i++)
        data[i] = Complex(RNG.urn(), RNG.urn());

    //save it to file
    const unsigned int shape[] = {Nx, Ny, Nz};
    cnpy::npy_save("arr1.npy", data, shape, 3, "w");

    //load it into a new array
    cnpy::NpyArray arr = cnpy::npy_load("arr1.npy");
    Complex *loaded_data = reinterpret_cast<Complex *>(arr.data);

    //make sure the loaded data matches the saved data
    assert(arr.word_size == sizeof(Complex));
    assert(arr.shape.size() == 3 && arr.shape[0] == Nx && arr.shape[1] == Ny && arr.shape[2] == Nz);
    for (int i = 0; i < Nx * Ny * Nz; i++)
        assert(Equal(data[i], loaded_data[i]));
    //append the same data to file
    //npy array on file now has shape (Nz+Nz,Ny,Nx)
    cnpy::npy_save("arr1.npy", data, shape, 3, "a");

    //now write to an npz file
    //non-array variables are treated as 1D arrays with 1 element
    double myVar1 = 1.2;
    char myVar2 = 'a';
    unsigned int shape2[] = {1};
    cnpy::npz_save("out.npz", "myVar1", &myVar1, shape2, 1, "w"); //"w" overwrites any existing file
    cnpy::npz_save("out.npz", "myVar2", &myVar2, shape2, 1, "a"); //"a" appends to the file we created above
    cnpy::npz_save("out.npz", "arr1", data, shape, 3, "a");       //"a" appends to the file we created above

    //load a single var from the npz file
    cnpy::NpyArray arr2 = cnpy::npz_load("out.npz", "arr1");
    //check that the loaded arr1 matches myVar1
    assert(arr2.word_size == sizeof(Complex));
    assert(arr2.shape.size() == 3 && arr2.shape[0] == Nx && arr2.shape[1] == Ny && arr2.shape[2] == Nz);

    //load the entire npz file
    cnpy::npz_t my_npz = cnpy::npz_load("out.npz");

    //check that the loaded myVar1 matches myVar1
    cnpy::NpyArray arr_mv1 = my_npz["myVar1"];
    double *mv1 = reinterpret_cast<double *>(arr_mv1.data);
    assert(arr_mv1.shape.size() == 1 && arr_mv1.shape[0] == 1);
    assert(mv1[0] == myVar1);

    //cleanup: note that we are responsible for deleting all loaded data
    delete[] data;
    //Better to use the destruct in NpyArray
    arr.destruct();
    arr2.destruct();
    my_npz.destruct();
}
