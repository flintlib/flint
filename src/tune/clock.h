/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FLINT_CLOCK_H
#define FLINT_CLOCK_H

#include <stdlib.h>
#include <time.h>
#include "flint.h"

#define FLINT_CPU_FREQUENCY_DEFAULT 3.2e9

#if FLINT64 && defined(__amd64__)
typedef ulong flint_time_t[1];

FLINT_FORCE_INLINE
double flint_time_nsec_diff(flint_time_t t1, flint_time_t t0)
{
    char * str = getenv("FLINT_CPU_FREQUENCY");
    double freq;
    double seconds;

    if (str == NULL)
        freq = FLINT_CPU_FREQUENCY_DEFAULT;
    else
        freq = strtod(str, NULL);

    seconds = (double) (*t1 - *t0) / freq;

    return seconds * 10e9;
}

FLINT_FORCE_INLINE
void flint_time_get(flint_time_t t0)
{
    __asm__ volatile (
            "rdtsc\n\t"
            "shl $32, %%rdx\n\t"
            "or %%rdx, %0"
            : "=a" (*t0) : : "rdx");
}
#else
typedef struct timespec flint_time_t[1];

FLINT_FORCE_INLINE
double flint_time_nsec_diff(flint_time_t t1, flint_time_t t0)
{
    return 1000000000.0 * (t1->tv_sec - t0->tv_sec)
        + (double) (t1->tv_nsec - t0->tv_nsec);
}

FLINT_FORCE_INLINE
void flint_time_get(flint_time_t t0)
{
    return clock_gettime(CLOCK_PROCESS_CPUTIME_ID, t0);
}
#endif

#endif /* FLINT_CLOCK_H */
