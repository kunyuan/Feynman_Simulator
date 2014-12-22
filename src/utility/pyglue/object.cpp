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

void PrintPyObject(PyObject* obj)
{
    PyObject_Print(obj, stdout, 0);
}
void MakeSureNoPyError(ERRORCODE e)
{
    if (PyErr_Occurred()) {
        PrintError();
        ERRORCODEABORT(e, "Python runtime error!");
    }
}

void IncreaseRef(PyObject* obj)
{
    if (obj != nullptr)
        Py_INCREF(obj);
}
void DecreaseRef(PyObject* obj, bool IsOwner)
{
    if (IsOwner)
        Py_XDECREF(obj);
}

Object::Object(const Object& obj)
{
    _PyPtr = obj.Borrow();
    _IsOwner = true;
    IncreaseRef(_PyPtr);
}

Object& Object::operator=(const Object& obj)
{
    _PyPtr = obj.Borrow();
    _IsOwner = true;
    IncreaseRef(_PyPtr);
    return *this;
}

Object::Object(PyObject* ptr, bool IsOwner)
{
    _PyPtr = ptr;
    _IsOwner = IsOwner;
}

Object Object::Borrow(PyObject* ptr)
{
    return Object(ptr, false);
}
Object Object::Steal(PyObject* ptr)
{
    return Object(ptr, true);
}

PyObject* Object::Borrow() const
{
    return _PyPtr;
}

PyObject* Object::Steal()
{
    if (IsOwner() == false)
        ERRORCODEABORT(ERR_VALUE_INVALID, "Object does not have ownership!");
    return _PyPtr;
}

void Object::Destroy()
{
    DecreaseRef(_PyPtr, IsOwner());
    _PyPtr = nullptr;
    _IsOwner = false;
}

void Object::Print()
{
    PyObject_Print(_PyPtr, stdout, 0);
}

void Object::_PrintDebug() const
{
    LOG_INFO("PyObject ref=" << _PyPtr->ob_refcnt);
}
std::string Object::PrettyString()
{
    Object result = Object::Steal(PyObject_Repr(_PyPtr));
    MakeSureNoPyError(ERR_VALUE_INVALID);
    return PyString_AsString(result.Borrow());
}

void Object::MakeSureNotNull()
{
    if (_PyPtr == nullptr)
        ERRORCODEABORT(ERR_VALUE_INVALID, "PyObject* is null!");
}
}