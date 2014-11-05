#include "crc32.h"

/* Standard CRC32 checksum: fast public domain implementation for
 * little-endian architectures.  Written for compilation with an
 * optimizer set to perform loop unwinding.  Outputs the checksum for
 * each file given as a command line argument.  Invalid file names and
 * files that cause errors are silently skipped.  The program reads
 * from stdin if it is called with no arguments. */

#include <stdint.h>

uint32_t crc32_for_byte(uint32_t r)
{
    for (int j = 0; j < 8; ++j)
        r = (r & 1 ? 0 : (uint32_t)0xEDB88320L) ^ r >> 1;
    return r ^ (uint32_t)0xFF000000L;
}

/* Any unsigned integer type with at least 32 bits may be used as
 * accumulator type for fast crc32-calulation, but unsigned long is
 * probably the optimal choice for most systems. */
typedef unsigned long accum_t;

void init_tables(uint32_t *table, uint32_t *wtable)
{
    for (size_t i = 0; i < 0x100; ++i)
        table[i] = crc32_for_byte((uint32_t)i);
    for (size_t k = 0; k < sizeof(accum_t); ++k)
        for (size_t w, i = 0; i < 0x100; ++i) {
            for (size_t j = w = 0; j < sizeof(accum_t); ++j)
                w = table[(uint8_t)(j == k ? w ^ i : w)] ^ w >> 8;
            wtable[(k << 8) + i] = (uint32_t)w ^ (k ? wtable[0] : 0);
        }
}

uint32_t crc32(uint32_t crc, unsigned char *data, size_t n_bytes)
{
    if (data == NULL)
        return 0;
    static uint32_t table[0x100], wtable[0x100 * sizeof(accum_t)];
    size_t n_accum = n_bytes / sizeof(accum_t);
    if (!*table)
        init_tables(table, wtable);
    for (size_t i = 0; i < n_accum; ++i) {
        accum_t a = crc ^ ((accum_t *)data)[i];
        for (size_t j = crc = 0; j < sizeof(accum_t); ++j)
            crc ^= wtable[(j << 8) + (uint8_t)(a >> 8 * j)];
    }
    for (size_t i = n_accum * sizeof(accum_t); i < n_bytes; ++i)
        crc = table[(uint8_t)crc ^ ((uint8_t *)data)[i]] ^ crc >> 8;
    return crc;
}