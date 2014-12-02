//
//  weight_builder.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/22/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weight_builder__
#define __Feynman_Simulator__weight_builder__

#include "module/weight/component.h"

namespace weight {
class G;
class W;
class GInitializer {
  public:
    GInitializer(G &G_)
        : _G(G_), _Map(G_._Map), _Shape(G_.GetShape())
    {
    }
    void BuildNew();

  private:
    G &_G;
    IndexMapSPIN2 &_Map;
    uint *_Shape;
    void _InitialTest();
    void _InitialDiagCounter();
    void _InitialBareSpin();
    void _InitialBareHubbardSquare();
};

class WInitializer {
  public:
    WInitializer(W &W_)
        : _W(W_), _Map(_W._Map), _Shape(W_.GetShape())
    {
    }
    void BuildNew();

  private:
    W &_W;
    IndexMapSPIN4 &_Map;
    uint *_Shape;
    void _InitialTest();
    void _InitialDiagCounter();
    void _InitialBareJ1J2();
    void _InitialBareHubbard();
};
}

#endif /* defined(__Feynman_Simulator__weight_builder__) */
