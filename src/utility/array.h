/* Array.h:  A high-performance multi-dimensional C++ array class
Copyright (C) 1997-2010 John C. Bowman

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

#ifndef __Array_h__
#define __Array_h__ 1

#define __ARRAY_H_VERSION__ 1.47

// Defining NDEBUG improves optimization but disables argument checking.
// Defining __NOARRAY2OPT inhibits special optimization of Array2[].
// Constructor (including copy constructor, as array3 A=new array3())
// always has reference sematic
// Assignment always has deep copy sematic

#include "abort.h"
#include <sstream>
#include <climits>
#include <cstdlib>
#include <cerrno>
#include <iostream>

#ifdef NDEBUG
#define __check(i, n, dim, m)
#define __checkSize()
#define __checkEqual(a, b, dim, m)
#define __checkActivate(i, align) this->Activate(align)
#else
#define __check(i, n, dim, m) this->Check(i, n, dim, m)
#define __checkSize() this->CheckSize()
#define __checkEqual(a, b, dim, m) this->CheckEqual(a, b, dim, m)
#define __checkActivate(i, align) this->CheckActivate(i, align)
#ifndef __NOARRAY2OPT
#define __NOARRAY2OPT
#endif
#endif

#ifndef HAVE_POSIX_MEMALIGN

#ifdef __GLIBC_PREREQ
#if __GLIBC_PREREQ(2, 3)
#define HAVE_POSIX_MEMALIGN
#endif
#else
#ifdef _POSIX_SOURCE
#define HAVE_POSIX_MEMALIGN
#endif
#endif

#else

#ifdef _AIX
extern "C" int posix_memalign(void **memptr, size_t alignment, size_t size);
#endif

#endif

int TestArray();

namespace Array {
inline std::ostream &_newl(std::ostream &s)
{
    s << '\n';
    return s;
}

inline void ArrayExit(const char *x);

#ifndef __ExternalArrayExit
inline void ArrayExit(const char *x)
{
    ABORT(_newl << "ERROR: " << x << "." << std::endl);
}
#endif

#ifndef __fftwpp_h__

// Adapted from FFTW aligned malloc/free.  Assumes that malloc is at least
// sizeof(void*)-aligned. Allocated memory must be freed with free0.
inline int posix_memalign0(void **memptr, size_t alignment, size_t size)
{
    if (alignment % sizeof(void *) != 0 || (alignment & (alignment - 1)) != 0)
        return EINVAL;
    void *p0 = malloc(size + alignment);
    if (!p0)
        return ENOMEM;
    void *p = (void *)(((size_t)p0 + alignment) & ~(alignment - 1));
    *((void **)p - 1) = p0;
    *memptr = p;
    return 0;
}

inline void free0(void *p)
{
    if (p)
        free(*((void **)p - 1));
}

template <class T>
inline void newAlign(T *&v, size_t len, size_t align)
{
    void *mem = NULL;
    const char *invalid = "Invalid alignment requested";
    const char *nomem = "Memory limits exceeded";
#ifdef HAVE_POSIX_MEMALIGN
    int rc = posix_memalign(&mem, align, len * sizeof(T));
#else
    int rc = posix_memalign0(&mem, align, len * sizeof(T));
#endif
    if (rc == EINVAL)
        Array::ArrayExit(invalid);
    if (rc == ENOMEM)
        Array::ArrayExit(nomem);
    v = (T *)mem;
    for (size_t i = 0; i < len; i++)
        new (v + i) T;
}

template <class T>
inline void deleteAlign(T *v, size_t len)
{
    for (size_t i = len - 1; i > 0; i--)
        v[i]. ~T();
    v[0]. ~T();
#ifdef HAVE_POSIX_MEMALIGN
    free(v);
#else
    free0(v);
#endif
}

#endif

template <class T>
class array1 {
  protected:
    T *v;
    unsigned int size;
    mutable int state;

  public:
    enum alloc_state { unallocated = 0,
                       allocated = 1,
                       temporary = 2,
                       aligned = 4 };
    virtual unsigned int Size() const
    {
        return size;
    }
    void CheckSize() const
    {
        if (!test(allocated) && size == 0)
            ArrayExit("Operation attempted on unallocated array");
    }
    void CheckEqual(int a, int b, unsigned int dim, unsigned int m) const
    {
        if (a != b) {
            std::ostringstream buf;
            buf << "Array" << dim << " index ";
            if (m)
                buf << m << " ";
            buf << "is incompatible in assignment (" << a << " != " << b << ")";
            ArrayExit(buf.str().c_str());
        }
    }

    int test(int flag) const
    {
        return state & flag;
    }
    void clear(int flag) const
    {
        //clear flag from the state
        state &= ~flag;
    }
    void set(int flag) const
    {
        state |= flag;
    }
    void Activate(size_t align = 0)
    {
        if (align) {
            newAlign(v, size, align);
            set(allocated | aligned);
        }
        else {
            v = new T[size];
            set(allocated);
        }
    }
    void CheckActivate(int dim, size_t align = 0)
    {
        if (test(allocated)) {
            std::ostringstream buf;
            buf << "Reallocation of Array" << dim
                << " attempted (must Deallocate first)";
            ArrayExit(buf.str().c_str());
        }
        Activate(align);
    }
    void Deallocate() const
    {
        if (test(allocated)) {
            if (test(aligned))
                deleteAlign(v, size);
            else
                delete[] v;
            state = unallocated;
        }
    }
    virtual void Dimension(unsigned int nx0)
    {
        size = nx0;
    }
    void Dimension(unsigned int nx0, T *v0)
    {
        //if the data pointer is passed from outside, then this array has not been allocated by itself
        Dimension(nx0);
        v = v0;
        clear(allocated);
    }
    void Dimension(const array1<T> &A)
    {
        Dimension(A.size, A.v);
        state = A.test(temporary);
    }

    void CheckActivate(size_t align = 0)
    {
        __checkActivate(1, align);
    }

    void Allocate(unsigned int nx0, size_t align = 0)
    {
        Dimension(nx0);
        CheckActivate(align);
    }

    void Reallocate(unsigned int nx0, size_t align = 0)
    {
        Deallocate();
        Allocate(nx0, align);
    }

    array1()
        : size(0), state(unallocated)
    {
    }
    array1(const void *)
        : size(0), state(unallocated)
    {
    }
    array1(unsigned int nx0, size_t align = 0)
        : state(unallocated)
    {
        Allocate(nx0, align);
    }
    array1(unsigned int nx0, T *v0)
        : state(unallocated)
    {
        Dimension(nx0, v0);
    }
    array1(T *v0)
        : state(unallocated)
    {
        Dimension(INT_MAX, v0);
    }
    array1(const array1<T> &A)
        : v(A.v), size(A.size),
          state(A.test(temporary))
    {
        //if A is temporary, then this array is also temporary
        //if A is activte or aligned, then this array is nothing but a reference to A, thus it is unallocated
    }

    virtual ~array1()
    {
        Deallocate();
    }

    void Freeze()
    {
        state = unallocated;
    }
    void Hold()
    {
        if (test(allocated)) {
            state = temporary;
        }
    }
    void Purge() const
    {
        if (test(temporary)) {
            Deallocate();
            state = unallocated;
        }
    }

    virtual void Check(int i, int n, unsigned int dim, unsigned int m,
                       int o = 0) const
    {
        if (i < 0 || i >= n) {
            std::ostringstream buf;
            buf << "Array" << dim << " index ";
            if (m)
                buf << m << " ";
            buf << "is out of bounds (" << i + o;
            if (n == 0)
                buf << " index given to empty array";
            else {
                if (i < 0)
                    buf << " < " << o;
                else
                    buf << " > " << n + o - 1;
            }
            buf << ")";
            ArrayExit(buf.str().c_str());
        }
    }

    unsigned int Nx() const
    {
        return size;
    }

#ifdef NDEBUG
    typedef T *opt;
#else
    typedef array1<T> opt;
#endif

    T &operator[](int ix) const
    {
        __check(ix, size, 1, 1);
        return v[ix];
    }
    T &operator()(int ix) const
    {
        __check(ix, size, 1, 1);
        return v[ix];
    }
    T *operator()() const
    {
        return v;
    }
    operator T *() const
    {
        return v;
    }

    array1<T> operator+(int i) const
    {
        return array1<T>(size - i, v + i);
    }

    T *Data() const
    {
        return v;
    }

    void Load(T a) const
    {
        __checkSize();
        for (unsigned int i = 0; i < size; i++)
            v[i] = a;
    }
    void Load(const T *a) const
    {
        for (unsigned int i = 0; i < size; i++)
            v[i] = a[i];
    }
    void Store(T *a) const
    {
        for (unsigned int i = 0; i < size; i++)
            a[i] = v[i];
    }
    void Set(T *a)
    {
        v = a;
        clear(allocated);
    }
    T Min()
    {
        if (size == 0)
            ArrayExit("Cannot take minimum of empty array");
        T min = v[0];
        for (unsigned int i = 1; i < size; i++)
            if (v[i] < min)
                min = v[i];
        return min;
    }
    T Max()
    {
        if (size == 0)
            ArrayExit("Cannot take maximum of empty array");
        T max = v[0];
        for (unsigned int i = 1; i < size; i++)
            if (v[i] > max)
                max = v[i];
        return max;
    }

    std::istream &Input(std::istream &s) const
    {
        __checkSize();
        for (unsigned int i = 0; i < size; i++)
            s >> v[i];
        return s;
    }

    array1<T> &operator=(T a)
    {
        Load(a);
        return *this;
    }
    array1<T> &operator=(const T *a)
    {
        Load(a);
        return *this;
    }
    array1<T> &operator=(const array1<T> &A)
    {
        __checkEqual(size, A.Size(), 1, 1);
        Load(A());
        A.Purge();
        return *this;
    }

    array1<T> &operator+=(const array1<T> &A)
    {
        __checkSize();
        for (unsigned int i = 0; i < size; i++)
            v[i] += A(i);
        return *this;
    }
    array1<T> &operator-=(const array1<T> &A)
    {
        __checkSize();
        for (unsigned int i = 0; i < size; i++)
            v[i] -= A(i);
        return *this;
    }
    array1<T> &operator*=(const array1<T> &A)
    {
        __checkSize();
        for (unsigned int i = 0; i < size; i++)
            v[i] *= A(i);
        return *this;
    }
    array1<T> &operator/=(const array1<T> &A)
    {
        __checkSize();
        for (unsigned int i = 0; i < size; i++)
            v[i] /= A(i);
        return *this;
    }

    array1<T> &operator+=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < size; i++)
            v[i] += a;
        return *this;
    }
    array1<T> &operator-=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < size; i++)
            v[i] -= a;
        return *this;
    }
    array1<T> &operator*=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < size; i++)
            v[i] *= a;
        return *this;
    }
    array1<T> &operator/=(T a)
    {
        __checkSize();
        T ainv = 1.0 / a;
        for (unsigned int i = 0; i < size; i++)
            v[i] *= ainv;
        return *this;
    }

    double L1() const
    {
        __checkSize();
        double norm = 0.0;
        for (unsigned int i = 0; i < size; i++)
            norm += abs(v[i]);
        return norm;
    }
#ifdef __ArrayExtensions
    double Abs2() const
    {
        __checkSize();
        double norm = 0.0;
        for (unsigned int i = 0; i < size; i++)
            norm += abs2(v[i]);
        return norm;
    }
    double L2() const
    {
        return sqrt(Abs2());
    }
    double LInfinity() const
    {
        __checkSize();
        double norm = 0.0;
        for (unsigned int i = 0; i < size; i++) {
            T a = abs(v[i]);
            if (a > norm)
                norm = a;
        }
        return norm;
    }
    double LMinusInfinity() const
    {
        __checkSize();
        double norm = DBL_MAX;
        for (unsigned int i = 0; i < size; i++) {
            T a = abs(v[i]);
            if (a < norm)
                norm = a;
        }
        return norm;
    }
#endif
};

template <class T>
array1<T> operator+(const array1<T> &t, const T &a)
{
    auto s = t;
    s += a;
    return s;
}

template <class T>
array1<T> operator+(const T &a, const array1<T> &t)
{
    auto s = t;
    s += a;
    return s;
}

template <class T>
array1<T> operator-(const array1<T> &t, const T &a)
{
    auto s = t;
    s -= a;
    return s;
}

template <class T>
array1<T> operator*(const array1<T> &t, const T &a)
{
    auto s = t;
    s *= a;
    return s;
}
template <class T>
array1<T> operator*(const T &a, const array1<T> &t)
{
    auto s = t;
    s *= a;
    return s;
}
template <class T>
array1<T> operator/(const array1<T> &t, const T &a)
{
    auto s = t;
    s /= a;
    return s;
}

template <class T>
void swaparray(T &A, T &B)
{
    T C;
    C.Dimension(A);
    A.Dimension(B);
    B.Dimension(C);
}

template <class T>
void leftshiftarray(T &A, T &B, T &C)
{
    T D;
    D.Dimension(A);
    A.Dimension(B);
    B.Dimension(C);
    C.Dimension(D);
}

template <class T>
void rightshiftarray(T &A, T &B, T &C)
{
    T D;
    D.Dimension(C);
    C.Dimension(B);
    B.Dimension(A);
    A.Dimension(D);
}

template <class T>
std::ostream &operator<<(std::ostream &s, const array1<T> &A)
{
    T *p = A();
    for (unsigned int i = 0; i < A.Nx(); i++) {
        s << *(p++) << " ";
    }
    return s;
}

template <class T>
std::istream &operator>>(std::istream &s, const array1<T> &A)
{
    return A.Input(s);
}

template <class T>
class array2 : public array1<T> {
  protected:
    unsigned int nx;
    unsigned int ny;

  public:
    using array1<T>::Dimension;

    void Dimension(unsigned int nx0, unsigned int ny0)
    {
        nx = nx0;
        ny = ny0;
        this->size = nx * ny;
    }
    void Dimension(unsigned int nx0, unsigned int ny0, T *v0)
    {
        Dimension(nx0, ny0);
        this->v = v0;
        this->clear(this->allocated);
    }
    void Dimension(const array1<T> &A)
    {
        ArrayExit("Operation not implemented");
    }

    void Allocate(unsigned int nx0, unsigned int ny0, size_t align = 0)
    {
        Dimension(nx0, ny0);
        __checkActivate(2, align);
    }

    void Allocate(unsigned int *n)
    {
        Allocate(n[0], n[1]);
    }

    array2()
        : nx(0), ny(0)
    {
    }
    array2(unsigned int nx0, unsigned int ny0, size_t align = 0)
    {
        Allocate(nx0, ny0, align);
    }
    array2(unsigned int nx0, unsigned int ny0, T *v0)
    {
        Dimension(nx0, ny0, v0);
    }
    array2(unsigned int *n)
    {
        Allocate(n[0], n[1]);
    }

    unsigned int Nx() const
    {
        return nx;
    }
    unsigned int Ny() const
    {
        return ny;
    }

#ifndef __NOARRAY2OPT
    T *operator[](int ix) const
    {
        return this->v + ix * ny;
    }
#else
    array1<T> operator[](int ix) const
    {
        __check(ix, nx, 2, 1);
        return array1<T>(ny, this->v + ix * ny);
    }
#endif
    T &operator()(int ix, int iy) const
    {
        __check(ix, nx, 2, 1);
        __check(iy, ny, 2, 2);
        return this->v[ix * ny + iy];
    }
    T &operator()(int i) const
    {
        __check(i, this->size, 2, 0);
        return this->v[i];
    }
    T *operator()() const
    {
        return this->v;
    }

    array2<T> &operator=(T a)
    {
        this->Load(a);
        return *this;
    }
    array2<T> &operator=(T *a)
    {
        this->Load(a);
        return *this;
    }
    array2<T> &operator=(const array2<T> &A)
    {
        __checkEqual(nx, A.Nx(), 2, 1);
        __checkEqual(ny, A.Ny(), 2, 2);
        this->Load(A());
        A.Purge();
        return *this;
    }

    array2<T> &operator+=(const array2<T> &A)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] += A(i);
        return *this;
    }
    array2<T> &operator-=(const array2<T> &A)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] -= A(i);
        return *this;
    }
    array2<T> &operator*=(const array2<T> &A);

    array2<T> &operator+=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] += a;
        return *this;
    }
    array2<T> &operator-=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] -= a;
        return *this;
    }
    array2<T> &operator*=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] *= a;
        return *this;
    }
    array2<T> &operator/=(T a)
    {
        __checkSize();
        T ainv = 1.0 / a;
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] *= ainv;
        return *this;
    }

    void Identity()
    {
        Load((T)0);
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] = (T)1;
    }
};

template <class T>
std::ostream &operator<<(std::ostream &s, const array2<T> &A)
{
    T *p = A();
    for (unsigned int i = 0; i < A.Nx(); i++) {
        for (unsigned int j = 0; j < A.Ny(); j++) {
            s << *(p++) << " ";
        }
        s << _newl;
    }
    s << std::flush;
    return s;
}

template <class T>
std::istream &operator>>(std::istream &s, const array2<T> &A)
{
    return A.Input(s);
}

template <class T>
class array3 : public array1<T> {
  protected:
    unsigned int nx;
    unsigned int ny;
    unsigned int nz;
    unsigned int nyz;

  public:
    using array1<T>::Dimension;

    void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0)
    {
        nx = nx0;
        ny = ny0;
        nz = nz0;
        nyz = ny * nz;
        this->size = nx * nyz;
    }
    void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0, T *v0)
    {
        Dimension(nx0, ny0, nz0);
        this->v = v0;
        this->clear(this->allocated);
    }

    void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                  size_t align = 0)
    {
        Dimension(nx0, ny0, nz0);
        __checkActivate(3, align);
    }

    void Allocate(unsigned int *n)
    {
        Allocate(n[0], n[1], n[2]);
    }

    array3()
        : nx(0), ny(0), nz(0), nyz(0)
    {
    }
    array3(unsigned int nx0, unsigned int ny0, unsigned int nz0,
           size_t align = 0)
    {
        Allocate(nx0, ny0, nz0, align);
    }
    array3(unsigned int nx0, unsigned int ny0, unsigned int nz0, T *v0)
    {
        Dimension(nx0, ny0, nz0, v0);
    }
    array3(unsigned int *n)
    {
        Allocate(n[0], n[1], n[2]);
    }

    unsigned int Nx() const
    {
        return nx;
    }
    unsigned int Ny() const
    {
        return ny;
    }
    unsigned int Nz() const
    {
        return nz;
    }

    array2<T> operator[](int ix) const
    {
        __check(ix, nx, 3, 1);
        return array2<T>(ny, nz, this->v + ix * nyz);
    }
    T &operator()(int ix, int iy, int iz) const
    {
        __check(ix, nx, 3, 1);
        __check(iy, ny, 3, 2);
        __check(iz, nz, 3, 3);
        return this->v[ix * nyz + iy * nz + iz];
    }
    T &operator()(int i) const
    {
        __check(i, this->size, 3, 0);
        return this->v[i];
    }
    T *operator()() const
    {
        return this->v;
    }

    array3<T> &operator=(T a)
    {
        this->Load(a);
        return *this;
    }
    array3<T> &operator=(T *a)
    {
        this->Load(a);
        return *this;
    }
    array3<T> &operator=(const array3<T> &A)
    {
        __checkEqual(nx, A.Nx(), 3, 1);
        __checkEqual(ny, A.Ny(), 3, 2);
        __checkEqual(nz, A.Nz(), 3, 3);
        this->Load(A());
        A.Purge();
        return *this;
    }

    array3<T> &operator+=(const array3<T> &A)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] += A(i);
        return *this;
    }
    array3<T> &operator-=(const array3<T> &A)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] -= A(i);
        return *this;
    }

    array3<T> &operator+=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] += a;
        return *this;
    }
    array3<T> &operator-=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] -= a;
        return *this;
    }
    array3<T> &operator*=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] *= a;
        return *this;
    }
    array3<T> &operator/=(T a)
    {
        __checkSize();
        T ainv = 1.0 / a;
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] *= ainv;
        return *this;
    }
};

template <class T>
std::ostream &operator<<(std::ostream &s, const array3<T> &A)
{
    T *p = A();
    for (unsigned int i = 0; i < A.Nx(); i++) {
        for (unsigned int j = 0; j < A.Ny(); j++) {
            for (unsigned int k = 0; k < A.Nz(); k++) {
                s << *(p++) << " ";
            }
            s << _newl;
        }
        s << _newl;
    }
    s << std::flush;
    return s;
}

template <class T>
std::istream &operator>>(std::istream &s, const array3<T> &A)
{
    return A.Input(s);
}

template <class T>
class array4 : public array1<T> {
  protected:
    unsigned int nx;
    unsigned int ny;
    unsigned int nz;
    unsigned int nw;
    unsigned int nyz;
    unsigned int nzw;
    unsigned int nyzw;

  public:
    using array1<T>::Dimension;

    void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                   unsigned int nw0)
    {
        nx = nx0;
        ny = ny0;
        nz = nz0;
        nw = nw0;
        nzw = nz * nw;
        nyzw = ny * nzw;
        this->size = nx * nyzw;
    }
    void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                   unsigned int nw0, T *v0)
    {
        Dimension(nx0, ny0, nz0, nw0);
        this->v = v0;
        this->clear(this->allocated);
    }

    void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                  unsigned int nw0, size_t align = 0)
    {
        Dimension(nx0, ny0, nz0, nw0);
        __checkActivate(4, align);
    }

    void Allocate(unsigned int *n)
    {
        Allocate(n[0], n[1], n[2], n[3]);
    }

    array4()
        : nx(0), ny(0), nz(0), nw(0), nyz(0), nzw(0), nyzw(0)
    {
    }
    array4(unsigned int nx0, unsigned int ny0, unsigned int nz0,
           unsigned int nw0, size_t align = 0)
    {
        Allocate(nx0, ny0, nz0, nw0, align);
    }
    array4(unsigned int nx0, unsigned int ny0, unsigned int nz0,
           unsigned int nw0, T *v0)
    {
        Dimension(nx0, ny0, nz0, nw0, v0);
    }
    array4(unsigned int *n)
    {
        Allocate(n[0], n[1], n[2], n[3]);
    }

    unsigned int Nx() const
    {
        return nx;
    }
    unsigned int Ny() const
    {
        return ny;
    }
    unsigned int Nz() const
    {
        return nz;
    }
    unsigned int N4() const
    {
        return nw;
    }

    array3<T> operator[](int ix) const
    {
        __check(ix, nx, 3, 1);
        return array3<T>(ny, nz, nw, this->v + ix * nyzw);
    }
    T &operator()(int ix, int iy, int iz, int iw) const
    {
        __check(ix, nx, 4, 1);
        __check(iy, ny, 4, 2);
        __check(iz, nz, 4, 3);
        __check(iw, nw, 4, 4);
        return this->v[ix * nyzw + iy * nzw + iz * nw + iw];
    }
    T &operator()(int i) const
    {
        __check(i, this->size, 4, 0);
        return this->v[i];
    }
    T *operator()() const
    {
        return this->v;
    }

    array4<T> &operator=(T a)
    {
        this->Load(a);
        return *this;
    }
    array4<T> &operator=(T *a)
    {
        this->Load(a);
        return *this;
    }

    array4<T> &operator=(const array4<T> &A)
    {
        __checkEqual(nx, A.Nx(), 4, 1);
        __checkEqual(ny, A.Ny(), 4, 2);
        __checkEqual(nz, A.Nz(), 4, 3);
        __checkEqual(nw, A.N4(), 4, 4);
        this->Load(A());
        A.Purge();
        return *this;
    }

    array4<T> &operator+=(const array4<T> &A)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] += A(i);
        return *this;
    }

    array4<T> &operator-=(const array4<T> &A)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] -= A(i);
        return *this;
    }

    array4<T> &operator+=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] += a;
        return *this;
    }
    array4<T> &operator-=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] -= a;
        return *this;
    }
    array4<T> &operator*=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] *= a;
        return *this;
    }
    array4<T> &operator/=(T a)
    {
        __checkSize();
        T ainv = 1.0 / a;
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] *= ainv;
        return *this;
    }
};

template <class T>
std::ostream &operator<<(std::ostream &s, const array4<T> &A)
{
    T *p = A;
    for (unsigned int i = 0; i < A.Nx(); i++) {
        for (unsigned int j = 0; j < A.Ny(); j++) {
            for (unsigned int k = 0; k < A.Nz(); k++) {
                for (unsigned int l = 0; l < A.N4(); l++) {
                    s << *(p++) << " ";
                }
                s << _newl;
            }
            s << _newl;
        }
        s << _newl;
    }
    s << std::flush;
    return s;
}

template <class T>
std::istream &operator>>(std::istream &s, const array4<T> &A)
{
    return A.Input(s);
}

template <class T>
class array5 : public array1<T> {
  protected:
    unsigned int nx;
    unsigned int ny;
    unsigned int nz;
    unsigned int nw;
    unsigned int nv;
    unsigned int nwv;
    unsigned int nzwv;
    unsigned int nyzwv;

  public:
    using array1<T>::Dimension;

    void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                   unsigned int nw0, unsigned int nv0)
    {
        nx = nx0;
        ny = ny0;
        nz = nz0;
        nw = nw0;
        nv = nv0;
        nwv = nw * nv;
        nzwv = nz * nwv;
        nyzwv = ny * nzwv;
        this->size = nx * nyzwv;
    }
    void Dimension(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                   unsigned int nw0, unsigned int nv0, T *v0)
    {
        Dimension(nx0, ny0, nz0, nw0, nv0);
        this->v = v0;
        this->clear(this->allocated);
    }

    void Allocate(unsigned int nx0, unsigned int ny0, unsigned int nz0,
                  unsigned int nw0, unsigned int nv0, size_t align = 0)
    {
        Dimension(nx0, ny0, nz0, nw0, nv0);
        __checkActivate(5, align);
    }

    void Allocate(unsigned int *n)
    {
        Allocate(n[0], n[1], n[2], n[3], n[4]);
    }

    array5()
        : nx(0), ny(0), nz(0), nw(0), nv(0), nwv(0), nzwv(0), nyzwv(0)
    {
    }
    array5(unsigned int nx0, unsigned int ny0, unsigned int nz0,
           unsigned int nw0, unsigned int nv0, size_t align = 0)
    {
        Allocate(nx0, ny0, nz0, nw0, nv0, align);
    }
    array5(unsigned int nx0, unsigned int ny0, unsigned int nz0,
           unsigned int nw0, unsigned int nv0, T *v0)
    {
        Dimension(nx0, ny0, nz0, nw0, nv0, v0);
    }
    array5(unsigned int *n)
    {
        Allocate(n[0], n[1], n[2], n[3], n[4]);
    }

    unsigned int Nx() const
    {
        return nx;
    }
    unsigned int Ny() const
    {
        return ny;
    }
    unsigned int Nz() const
    {
        return ny;
    }
    unsigned int N4() const
    {
        return nw;
    }
    unsigned int N5() const
    {
        return nv;
    }

    array4<T> operator[](int ix) const
    {
        __check(ix, nx, 4, 1);
        return array4<T>(ny, nz, nw, nv, this->v + ix * nyzwv);
    }
    T &operator()(int ix, int iy, int iz, int iw, int iv) const
    {
        __check(ix, nx, 5, 1);
        __check(iy, ny, 5, 2);
        __check(iz, nz, 5, 3);
        __check(iw, nw, 5, 4);
        __check(iv, nv, 5, 5);
        return this->v[ix * nyzwv + iy * nzwv + iz * nwv + iw * nv + iv];
    }
    T &operator()(int i) const
    {
        __check(i, this->size, 5, 0);
        return this->v[i];
    }
    T *operator()() const
    {
        return this->v;
    }

    array5<T> &operator=(T a)
    {
        this->Load(a);
        return *this;
    }
    array5<T> &operator=(T *a)
    {
        this->Load(a);
        return *this;
    }

    array5<T> &operator=(const array5<T> &A)
    {
        __checkEqual(nx, A.Nx(), 5, 1);
        __checkEqual(ny, A.Ny(), 5, 2);
        __checkEqual(nz, A.Nz(), 5, 3);
        __checkEqual(nw, A.N4(), 5, 4);
        __checkEqual(nv, A.N5(), 5, 5);
        this->Load(A());
        A.Purge();
        return *this;
    }

    array5<T> &operator+=(const array5<T> &A)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] += A(i);
        return *this;
    }
    array5<T> &operator-=(const array5<T> &A)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] -= A(i);
        return *this;
    }

    array5<T> &operator+=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] += a;
        return *this;
    }
    array5<T> &operator-=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] -= a;
        return *this;
    }
    array5<T> &operator*=(T a)
    {
        __checkSize();
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] *= a;
        return *this;
    }
    array5<T> &operator/=(T a)
    {
        __checkSize();
        T ainv = 1.0 / a;
        for (unsigned int i = 0; i < this->size; i++)
            this->v[i] *= ainv;
        return *this;
    }
};

template <class T>
std::ostream &operator<<(std::ostream &s, const array5<T> &A)
{
    T *p = A;
    for (unsigned int i = 0; i < A.Nx(); i++) {
        for (unsigned int j = 0; j < A.Ny(); j++) {
            for (unsigned int k = 0; k < A.Nz(); k++) {
                for (unsigned int l = 0; l < A.N4(); l++) {
                    for (unsigned int l = 0; l < A.N5(); l++) {
                        s << *(p++) << " ";
                    }
                    s << _newl;
                }
                s << _newl;
            }
            s << _newl;
        }
        s << _newl;
    }
    s << std::flush;
    return s;
}

template <class T>
std::istream &operator>>(std::istream &s, const array5<T> &A)
{
    return A.Input(s);
}
}

#undef __check
#undef __checkSize
#undef __checkActivate

#endif
