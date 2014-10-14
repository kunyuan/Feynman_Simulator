//
//  vector.h
//  Feynman_Simulator
//
//  Created by yuan on 10/14/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__vector__
#define __Feynman_Simulator__vector__

#include <stdio.h>
#include "convention.h"

template <typename T>
class Vec{
private:
    T _Arrary[D];
public:
    Vec(T t)
    {
        for(int i=0; i<D; i++)
            _Arrary[i]=t;
    }
    
    T& operator[](int index)
    {
        return _Arrary[index];
    }
};
#endif /* defined(__Feynman_Simulator__vector__) */
