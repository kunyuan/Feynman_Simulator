//
//  pyarraywrapper.h
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/25/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#ifndef __Feynman_Simulator__pyarraywrapper__
#define __Feynman_Simulator__pyarraywrapper__

#include "utility/abort.h"
#include <vector>
#include "object.h"

namespace Python {
void ArrayInitialize();
class ArrayObject;
bool Convert(Object, ArrayObject&);
class ArrayObject : public Object {
public:
    ArrayObject()
        : Object()
    {
    }
    ArrayObject(const Object& obj);
    ArrayObject(PyObject* obj, OwnerShip ownership = NewRef);
    template <typename T>
    ArrayObject(T* data, uint* Shape, int Dim);
    void* Data();
    std::vector<uint> Shape();
    int Dim();
    ArrayObject& operator=(const ArrayObject& obj)
    {
        Object::operator=(obj);
        return *this;
    }
};
}

#endif /* defined(__Feynman_Simulator__pyarraywrapper__) */
