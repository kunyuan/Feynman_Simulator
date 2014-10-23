/* zlib.h -- interface of the 'zlib' general purpose compression library
  version 1.2.8, April 28th, 2013

  Copyright (C) 1995-2013 Jean-loup Gailly and Mark Adler

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Jean-loup Gailly        Mark Adler
  jloup@gzip.org          madler@alumni.caltech.edu


  The data format used by the zlib library is described by RFCs (Request for
  Comments) 1950 to 1952 in the files http://tools.ietf.org/html/rfc1950
  (zlib format), rfc1951 (deflate format) and rfc1952 (gzip format).
*/

#ifndef ZLIB_H
#define ZLIB_H

#if defined(__MSDOS__) && !defined(MSDOS)
#define MSDOS
#endif
#if (defined(OS_2) || defined(__OS2__)) && !defined(OS2)
#define OS2
#endif
#if defined(_WINDOWS) && !defined(WINDOWS)
#define WINDOWS
#endif
#if defined(_WIN32) || defined(_WIN32_WCE) || defined(__WIN32__)
#ifndef WIN32
#define WIN32
#endif
#endif
#if (defined(MSDOS) || defined(OS2) || defined(WINDOWS)) && !defined(WIN32)
#if !defined(__GNUC__) && !defined(__FLAT__) && !defined(__386__)
#ifndef SYS16BIT
#define SYS16BIT
#endif
#endif
#endif

/*
 * Compile with -DMAXSEG_64K if the alloc function cannot allocate more
 * than 64k bytes at a time (needed on systems with 16-bit int).
 */
#ifdef SYS16BIT
#define MAXSEG_64K
#endif
#ifdef MSDOS
#define UNALIGNED_OK
#endif

#ifdef __STDC_VERSION__
#ifndef STDC
#define STDC
#endif
#if __STDC_VERSION__ >= 199901L
#ifndef STDC99
#define STDC99
#endif
#endif
#endif
#if !defined(STDC) && (defined(__STDC__) || defined(__cplusplus))
#define STDC
#endif
#if !defined(STDC) && (defined(__GNUC__) || defined(__BORLANDC__))
#define STDC
#endif
#if !defined(STDC) && (defined(MSDOS) || defined(WINDOWS) || defined(WIN32))
#define STDC
#endif
#if !defined(STDC) && (defined(OS2) || defined(__HOS_AIX__))
#define STDC
#endif

#if defined(__OS400__) && !defined(STDC) /* iSeries (formerly AS/400). */
#define STDC
#endif

#ifndef STDC
#ifndef const /* cannot use !defined(STDC) && !defined(const) on Mac */
#define const /* note: need a more gentle solution here */
#endif
#endif

#if defined(ZLIB_CONST) && !defined(z_const)
#define z_const const
#else
#define z_const
#endif

/* Some Mac compilers merge all .h files incorrectly: */
#if defined(__MWERKS__) || defined(applec) || defined(THINK_C) || defined(__SC__)
#define NO_DUMMY_DECL
#endif
#if !defined(__MACTYPES__)
typedef unsigned char Byte; /* 8 bits */
#endif
typedef unsigned int uInt;   /* 16 bits or more */
typedef unsigned long uLong; /* 32 bits or more */

typedef Byte Bytef;
typedef char charf;
typedef int intf;
typedef uInt uIntf;
typedef uLong uLongf;

#if !defined(Z_U4) && !defined(STDC)
#include <limits.h>
#if (UINT_MAX == 0xffffffffUL)
#define Z_U4 unsigned
#elif(ULONG_MAX == 0xffffffffUL)
#define Z_U4 unsigned long
#elif(USHRT_MAX == 0xffffffffUL)
#define Z_U4 unsigned short
#endif
#endif

#ifdef Z_U4
typedef Z_U4 z_crc_t;
#else
typedef unsigned long z_crc_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* checksum functions */

extern uLong crc32(uLong crc, const Bytef *buf, uInt len);
/*
     Update a running CRC-32 with the bytes buf[0..len-1] and return the
   updated CRC-32.  If buf is Z_NULL, this function returns the required
   initial value for the crc.  Pre- and post-conditioning (one's complement) is
   performed within this function so it shouldn't be done by the application.

   Usage example:

     uLong crc = crc32(0L, Z_NULL, 0);

     while (read_buffer(buffer, length) != EOF) {
       crc = crc32(crc, buffer, length);
     }
     if (crc != original_crc) error();
*/

int TestCRC32();

#ifdef __cplusplus
}
#endif

#endif /* ZLIB_H */
