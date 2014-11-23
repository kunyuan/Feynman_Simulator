//
//  weight_builder.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/22/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__weight_builder__
#define __Feynman_Simulator__weight_builder__

namespace weight {
class G;
class GInitializer {
  public:
    GInitializer(const G &G_);
    void BuildNew();

  private:
    G &_G;
};
}

#endif /* defined(__Feynman_Simulator__weight_builder__) */
