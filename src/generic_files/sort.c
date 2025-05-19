/*
   Copyright (C) 2025 Marc Mezzarobba

   This file is part of FLINT.

   FLINT is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License (LGPL) as published
   by the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>

#include "flint.h"

typedef int (* cmp_t) (const void *, const void *, void *);

static void
merge_sort(unsigned char * buf0, unsigned char * alt0, slong len, slong size,
           cmp_t cmp, void * data)
{
    if (len <= 1)
        return;

    slong l0 = len / 2;
    slong l1 = len - l0;

    unsigned char * buf1 = buf0 + l0 * size;
    unsigned char * alt1 = alt0 + l0 * size;

    merge_sort(alt0, buf0, l0, size, cmp, data);
    merge_sort(alt1, buf1, l1, size, cmp, data);

    slong i0 = 0, i1 = 0;

    for (slong k = 0; k < len; k++)
    {
        unsigned char * elt;

        if (i1 >= l1 || (i0 < l0 && cmp(alt0 + i0 * size,
                                        alt1 + i1 * size, data) <= 0))
            elt = alt0 + i0++ * size;
        else
            elt = alt1 + i1++ * size;

        memcpy(buf0 + k * size, elt, size);
    }
}

void
flint_merge_sort(void * buf, slong len, slong size,
                 cmp_t cmp, void * data)
{
    unsigned char * alt = flint_malloc(len * size);
    memcpy(alt, buf, len * size);
    merge_sort(buf, alt, len, size, cmp, data);
    flint_free(alt);
}

void
flint_sort(void * buf, slong len, slong size,
           cmp_t cmp, void * data)
{
#ifdef __GLIBC__
    qsort_r(buf, len, size, cmp, data);
#else
    flint_merge_sort(buf, len, size, cmp, data);
#endif
}
