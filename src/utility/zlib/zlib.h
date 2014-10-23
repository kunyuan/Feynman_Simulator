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

#if !defined(Z_U4) && !defined(Z_SOLO) && defined(STDC)
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
