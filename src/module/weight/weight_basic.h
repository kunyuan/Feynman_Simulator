//
//  weight_basic.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/21/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weight_basic__
#define __Feynman_Simulator__weight_basic__

#include "utility/convention.h"
#include "lattice/lattice.h"
#include <vector>
class Lattice;

namespace weight {

class Basic {
  public:
    static int SpinIndex(spin SpinIn, spin SpinOut);
    static bool IsSameSpin(int spindex)
    {
        return (spindex == 0 || spindex == 2);
    }

    //First In/Out: direction of WLine; Second In/Out: direction of Vertex
    static int SpinIndex(spin SpinInIn, spin SpinInOut, spin SpinOutIn, spin SpinOutOut);
    static int SpinIndex(spin *TwoSpinIn, spin *TwoSpinOut);

    enum SpinFilter { UpUp2UpUp,
                      UpDown2UpDown,
                      UpDown2DownUp };
    static std::vector<int> GetSpinIndexVector_FourSpinsFileter(SpinFilter);
    /**
    *  @return return 1 if weight is symmtric in tau, otherwise, return -1
    */
    int TauToBin(real tau);
    int TauToBin(real t_in, real t_out);
    real BinToTau(int Bin);

    void Reset(real beta);

  protected:
    Basic(const Lattice &lat, real Beta, model Model,
          bool IsTauSymmetric, std::string);
    model _Model;
    bool _IsTauSymmetric;
    std::string _Name;
    real _Beta;
    real _dBeta; //_Beta/MAX_TAU
    real _dBetaInverse;
    Lattice _Lat;
};
}

#endif /* defined(__Feynman_Simulator__weight_basic__) */
