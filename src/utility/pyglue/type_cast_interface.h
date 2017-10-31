//
//  type_cast_interface.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/23/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef Feynman_Simulator_type_cast_interface_h
#define Feynman_Simulator_type_cast_interface_h
#include "object.h"
namespace Python {
class ITypeCast {
public:
    virtual Object ToPy() const = 0;
    virtual bool FromPy(const Object&) = 0;

protected:
    ~ITypeCast() {}
};
}

//Tempate:
/*
class MyType : public ITypeCast {
public:
    virtual Object CastToPy() const
    {
        //convert MyType to Object
    }
    virtual bool Convert(Object)
    {
        //convert object to MyType
    }
};
     */
#endif
