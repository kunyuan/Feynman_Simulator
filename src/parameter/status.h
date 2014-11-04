//
//  control.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__status__
#define __Feynman_Simulator__status__

#include "parser.h"
#include <string>

class Status {
  public:
    int Version;
    real Jcp;
    real Beta;
    int OrderAccepted;

    bool Load();
    void Save();

  private:
    const std::string _StatusFile = "status.txt";
    SimpleParser _Para;
};

#endif /* defined(__Feynman_Simulator__control__) */
