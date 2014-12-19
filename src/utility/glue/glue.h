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

__PYX_EXTERN_C DL_IMPORT(PyObject) *NewDict(void);
__PYX_EXTERN_C DL_IMPORT(void) DictSet_int(PyObject *, std::string, int);
__PYX_EXTERN_C DL_IMPORT(void) DictGet_int(PyObject *, std::string, int *);
__PYX_EXTERN_C DL_IMPORT(void) DictDel(PyObject *, std::string, PyObject *);
__PYX_EXTERN_C DL_IMPORT(void) DictPrint(PyObject *);

#endif /* !__PYX_HAVE_API__glue */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initglue(void);
#else
PyMODINIT_FUNC PyInit_glue(void);
#endif

#endif /* !__PYX_HAVE__glue */
