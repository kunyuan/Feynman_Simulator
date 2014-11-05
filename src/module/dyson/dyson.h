//
//  dyson.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/4/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__dyson__
#define __Feynman_Simulator__dyson__

#include "utility/convention.h"
namespace weight {
class Weight;
class G;
class W;
class Sigma;
class Polar;
}
namespace para {
class ParaDyson;
}

namespace dyson {
class Dyson {
  public:
    real Beta;
    weight::G *G;
    weight::W *W;
    weight::Sigma *Sigma;
    weight::Polar *Polar;
    bool BuildNew(para::ParaDyson &, weight::Weight &);
};
}

#endif /* defined(__Feynman_Simulator__dyson__) */
