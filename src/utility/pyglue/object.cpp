//
//  object.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/21/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "object.h"
#include "utility/abort.h"
#include <Python/Python.h>
namespace Python {
void IncreaseRef(PyObject* obj, Ownership os)
{
    if (obj != nullptr)
        if (os == Received)
            Py_INCREF(obj);
}
void DecreaseRef(PyObject* obj, Ownership os)
{
    if (os == Received)
        Py_XDECREF(obj);
}

void Initialize()
{
    Py_Initialize();
}

void Finalize()
{
    Py_Finalize();
}

void ClearError()
{
    PyErr_Clear();
}

void PrintError()
{
    PyErr_Print();
}

void MakeSureNoPyError(ERRORCODE e)
{
    if (PyErr_Occurred()) {
        PrintError();
        ERRORCODEABORT(e, "Python runtime error!");
    }
}

Object::Object(PyObject* ptr, Ownership os)
{
    _PyPtr = ptr;
    _OwnerShip = os;
    IncreaseRef(_PyPtr, _OwnerShip);
}
Object::Object(const Object& obj, Ownership os)
{
    _PyPtr = obj.GetPtr();
    _OwnerShip = obj.GetOwnerShip();
    IncreaseRef(_PyPtr, _OwnerShip);
}
void Object::Destroy()
{
    DecreaseRef(_PyPtr, _OwnerShip);
    _PyPtr = nullptr;
}

void Object::Print()
{
    PyObject_Print(_PyPtr, stdout, 0);
}

std::string Object::PrettyString()
{
    Object result = PyObject_Repr(_PyPtr);
    MakeSureNoPyError(ERR_VALUE_INVALID);
    return PyString_AsString(result.GetPtr());
}

void Object::MakeSureNotNull()
{
    if (_PyPtr == nullptr)
        ERRORCODEABORT(ERR_VALUE_INVALID, "PyObject* is null!");
}
}