//
//  control.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 11/3/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__status__
#define __Feynman_Simulator__status__

#include <string>
#include "utility/convention.h"

namespace para {

class Message {
public:
    int Version;
    real Beta;
    real SqueezeFactor;

    bool Load(const std::string& FileName);
    void Save(const std::string& FileName);
    std::string PrettyString();
};
}

#endif /* defined(__Feynman_Simulator__control__) */
