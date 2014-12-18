#ifndef __PYX_HAVE__glue
#define __PYX_HAVE__glue

/* "glue.pyx":8
 * from cython.operator cimport dereference, preincrement
 * cimport cython.operator as co
 * ctypedef public double real             # <<<<<<<<<<<<<<
 * 
 * cdef extern from "../complex.h":
 */

typedef double real;

#ifndef __PYX_HAVE_API__glue

#ifndef __PYX_EXTERN_C
#ifdef __cplusplus
#define __PYX_EXTERN_C extern "C"
#else
#define __PYX_EXTERN_C extern
#endif
#endif

__PYX_EXTERN_C DL_IMPORT(PyObject) * CastBool(bool);
__PYX_EXTERN_C DL_IMPORT(PyObject) * CastInt(int);
__PYX_EXTERN_C DL_IMPORT(PyObject) * CastReal(real);
__PYX_EXTERN_C DL_IMPORT(PyObject) * CastStr(std::string);
__PYX_EXTERN_C DL_IMPORT(PyObject) * CastComplex(Complex);
__PYX_EXTERN_C DL_IMPORT(bool) AsBool(PyObject*);
__PYX_EXTERN_C DL_IMPORT(void) AsInt(PyObject*, int&);
__PYX_EXTERN_C DL_IMPORT(real) AsReal(PyObject*, real&);
__PYX_EXTERN_C DL_IMPORT(real) AsComplex(PyObject*, Complex&);
__PYX_EXTERN_C DL_IMPORT(void) AsStr(PyObject*, std::string&);
__PYX_EXTERN_C DL_IMPORT(PyObject) * NewDict(void);
__PYX_EXTERN_C DL_IMPORT(void) DictSet(PyObject*, std::string, PyObject*);
__PYX_EXTERN_C DL_IMPORT(PyObject) * DictGet(PyObject*, std::string);
__PYX_EXTERN_C DL_IMPORT(void) DictDel(PyObject*, std::string, PyObject*);
__PYX_EXTERN_C DL_IMPORT(void) DictPrint(PyObject*);
__PYX_EXTERN_C DL_IMPORT(PyObject) * __pyx_fuse_0Cast(bool);
__PYX_EXTERN_C DL_IMPORT(PyObject) * __pyx_fuse_1Cast(int);
__PYX_EXTERN_C DL_IMPORT(PyObject) * __pyx_fuse_2Cast(real);
__PYX_EXTERN_C DL_IMPORT(PyObject) * __pyx_fuse_3Cast(std::string);
__PYX_EXTERN_C DL_IMPORT(PyObject) * __pyx_fuse_4Cast(Complex);
__PYX_EXTERN_C DL_IMPORT(void) __pyx_fuse_0As(PyObject*, bool*);
__PYX_EXTERN_C DL_IMPORT(void) __pyx_fuse_1As(PyObject*, int*);
__PYX_EXTERN_C DL_IMPORT(void) __pyx_fuse_2As(PyObject*, real*);
__PYX_EXTERN_C DL_IMPORT(void) __pyx_fuse_3As(PyObject*, std::string*);
__PYX_EXTERN_C DL_IMPORT(void) __pyx_fuse_4As(PyObject*, Complex*);

#endif /* !__PYX_HAVE_API__glue */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initglue(void);
#else
PyMODINIT_FUNC PyInit_glue(void);
#endif

#endif /* !__PYX_HAVE__glue */
