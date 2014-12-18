# distutils: language = c++
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
cimport cython
from cython.operator cimport dereference, preincrement
cimport cython.operator as co
ctypedef public double real

cdef extern from "../complex.h":
    cppclass Complex:
        real Re
        real Im

ctypedef fused generic:
    bool
    int
    real
    string
    Complex

cdef public object Cast(generic value):
    if generic is Complex:
        return value.Re+1j*value.Im
    else:
        return value
cdef public void As(object value, generic* result):
    if generic is Complex:
        result[0].Re=value.real
        result[0].Im=value.imag
    elif generic is string:
        result[0].append(<char *>value)
    else:
        result[0]=value

cdef public object CastBool(bool value):
    return value
cdef public object CastInt(int value):
    return value
cdef public object CastReal(real value):
    return value
cdef public object CastStr(string value):
    return value
cdef public object CastComplex(Complex value):
    return value.Re+1j*value.Im
cdef public bool AsBool(value):
    return value
cdef public void AsInt(object value, int& result):
    (&result)[0]=value
cdef public real AsReal(object value, real& result):
    (&result)[0]=value
cdef public real AsComplex(object value, Complex& result):
    result.Re=value.real
    result.Im=value.imag
cdef public void AsStr(value, string& result):
    result.append(<char *>value)

cdef public NewDict():
    node={}
    return node

cdef public void DictSet(node, string key, value):
    node[key]=value

cdef public DictGet(node, string key) except +:
    return node[key]

cdef public void DictDel(node, string key, value):
    node.pop(key, None)

cdef public void DictPrint(node):
    print node

def hello():
    print "hello"

cdef void world():
    print "world"

