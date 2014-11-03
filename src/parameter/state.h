//
//  control.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__control__
#define __Feynman_Simulator__control__

#include "parameter_map.h"
#include <string>

class State {
  public:
    int Version;
    real Jcp;
    real Beta;
    int OrderAccepted;

    bool Load();
    void Save();

  private:
    const std::string _StateFile = "state.txt";
    ParameterMap _Para;
};

#endif /* defined(__Feynman_Simulator__control__) */
