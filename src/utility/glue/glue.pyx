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

cdef inline object Cast(generic value):
    if generic is Complex:
        return value.Re+1j*value.Im
    else:
        return value
cdef inline void GetAs(object value, generic* result):
    if generic is Complex:
        result[0].Re=value.real
        result[0].Im=value.imag
    elif generic is string:
        result[0].append(<char *>value)
    else:
        result[0]=value

cdef public NewDict():
    node={}
    return node

cdef public void DictSet_int(node, string key, int value):
    node[key]=Cast(value)

cdef public void DictGet_int(node, string key, int* result) except +:
    GetAs(node[key], result)

cdef public void DictDel(node, string key, value):
    node.pop(key, None)

cdef public void DictPrint(node):
    print node

def hello():
    print "hello"

cdef void world():
    print "world"

