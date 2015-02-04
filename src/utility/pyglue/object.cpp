//
//  object.cpp
//  Feynman_Simulator
//
//  Created by Kun Chen on 12/21/14.
//  Copyright (c) 2014 Kun Chen. All rights reserved.
//

#include "object.h"
#include <iostream>
#include "utility/abort.h"
#include "utility/scopeguard.h"
//#include <Python/Python.h>
#include <Python.h>

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

std::string GetPyErrorStr()
{
    PyObject* ptype;
    PyObject* pvalue;
    PyObject* ptraceback;
    PyErr_Fetch(&ptype, &pvalue, &ptraceback);
    PyErr_NormalizeException(&ptype, &pvalue, &ptraceback);
    Object type = PyObject_Str(ptype);
    Object value = PyObject_Str(pvalue);
    auto result = std::string(PyString_AsString(type.Get())) + "\n"
                  + std::string(PyString_AsString(value.Get()));
    Py_XDECREF(ptype);
    Py_XDECREF(pvalue);
    Py_XDECREF(ptraceback);
    return result;
}

void PropagatePyError()
{
    ON_SCOPE_EXIT([&] {PyErr_Clear(); });
    //clear Python error state once the error is propagated
    if (PyErr_Occurred()) {
        if (PyErr_ExceptionMatches(PyExc_IOError))
            THROW(IOInvalid, GetPyErrorStr(), WARNING);
        else if (PyErr_ExceptionMatches(PyExc_KeyError))
            THROW(KeyInvalid, GetPyErrorStr(), WARNING);
        else if (PyErr_ExceptionMatches(PyExc_IndexError))
            THROW(IndexInvalid, GetPyErrorStr(), WARNING);
        else if (PyErr_ExceptionMatches(PyExc_MemoryError))
            THROW(MemoryException, GetPyErrorStr(), WARNING);
        else
            THROW(RunTimeException, GetPyErrorStr(), WARNING);
    }
}
void IncreaseRef(PyObject* obj)
{
    if (obj != nullptr)
        Py_INCREF(obj);
}
void DecreaseRef(PyObject* obj)
{
    Py_XDECREF(obj);
}
Object Object::Copy() const
{
    IncreaseRef(_PyPtr);
    return Object(_PyPtr, NewRef);
}

Object::Object(const Object& obj)
{
    _PyPtr = obj.Get(NewRef);
}

Object& Object::operator=(const Object& obj)
{
    if (this == &obj)
        return *this;
    // The PyObject behind obj._PyPtr may be deleted once DecreaseRef(this->_PyPtr) is called. Thus DecreaseRef() should be called after the assigment is done. Otherwise, it can be dangerous when the old this->_PyPtr==obj._PyPtr
    PyObject* temp = _PyPtr;
    _PyPtr = obj.Get(NewRef);
    DecreaseRef(temp);
    return *this;
}

Object::Object(PyObject* obj, OwnerShip os)
{
    _PyPtr = obj;
    if (os != NewRef)
        IncreaseRef(_PyPtr);
}

PyObject* Object::Get(OwnerShip os) const
{
    if (os == NewRef)
        IncreaseRef(_PyPtr);
    return _PyPtr;
}

void Object::Destroy()
{
    DecreaseRef(_PyPtr);
    _PyPtr = nullptr;
}

long Object::RefCount() const
{
    if (_PyPtr != nullptr)
        return _PyPtr->ob_refcnt;
    else
        return -1;
}

void Object::Print() const
{
    if (_PyPtr != nullptr) {
        PyObject_Print(_PyPtr, stdout, 0);
        std::cout << std::endl;
    }
    else
        std::cout << "nullptr" << std::endl;
    //    PyRun_SimpleString("\n");
}

void Object::_PrintDebug() const
{
    if (_PyPtr != nullptr)
        LOG_INFO("PyObject ref=" << _PyPtr->ob_refcnt);
    Print();
}
std::string Object::PrettyString()
{
    Object result = PyObject_Repr(_PyPtr);
    PropagatePyError();
    return PyString_AsString(result.Get());
}

void Object::MakeSureNotNull()
{
    if (_PyPtr == nullptr)
        ABORT("PyObject* is null!");
}
}