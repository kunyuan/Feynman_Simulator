//
//  weight_builder.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/22/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weight_builder__
#define __Feynman_Simulator__weight_builder__

namespace weight0 {
class G;
class GInitializer {
  public:
    GInitializer(G &G_)
        : _G(G_)
    {
    }
    void BuildNew();
    void BuildTest();

  private:
    G &_G;
    void _InitialTest();
    void _InitialDiagCounter();
    void _InitialBareSpin();
    void _InitialBareHubbardSquare();
};
}

#endif /* defined(__Feynman_Simulator__weight_builder__) */
