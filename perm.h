/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#ifndef PERM_H
#define PERM_H

#undef ulong /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long

#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Memory management *********************************************************/

static __inline__ long * _perm_init(long n)
{
    long i, *vec;

    vec = (long *) flint_malloc(n * sizeof(long));

    if (!vec)
    {
        printf("ERROR (_perm_init).\n\n");
        abort();
    }

    for (i = 0; i < n; i++)
        vec[i] = i;

    return vec;
}

static __inline__ void _perm_clear(long * vec)
{
    flint_free(vec);
}

/* Assignment ****************************************************************/

static __inline__ long _perm_equal(const long *vec1, const long *vec2, long n)
{
    long i;

    for (i = 0; i < n; i++)
        if (vec1[i] != vec2[i])
            return 0;

    return 1;
}

static __inline__ void _perm_set(long *res, const long *vec, long n)
{
    long i;

    for (i = 0; i < n; i++)
        res[i] = vec[i];
}

static __inline__ void _perm_set_one(long *vec, long n)
{
    long i;

    for (i = 0; i < n; i++)
        vec[i] = i;
}

static __inline__ void
 _perm_inv(long *res, const long *vec, long n)
{
    long i;

    if (res == vec)
    {
        long *t = (long *) flint_malloc(n * sizeof(long));

        if (!t)
        {
            printf("ERROR (_perm_inv).\n\n");
            abort();
        }

        for (i = 0; i < n; i++)
            t[i] = vec[i];
        for (i = 0; i < n; i++)
            res[t[i]] = i;

        flint_free(t);
    }
    else
    {
        for (i = 0; i < n; i++)
            res[vec[i]] = i;
    }
}

/* Composition ***************************************************************/

static __inline__ void
_perm_compose(long *res, const long *vec1, const long *vec2, long n)
{
    long i;

    if (res == vec1)
    {
        long *t = (long *) flint_malloc(n * sizeof(long));

        for (i = 0; i < n; i++)
            t[i] = vec1[i];
        for (i = 0; i < n; i++)
            res[i] = t[vec2[i]];

        flint_free(t);
    }
    else
    {
        for (i = 0; i < n; i++)
            res[i] = vec1[vec2[i]];
    }
}

/* Randomisation *************************************************************/

int _perm_randtest(long * vec, long n, flint_rand_t state);

/* Parity ********************************************************************/

int _perm_parity(const long * vec, long n);

/* Input and output **********************************************************/

static __inline__  int _long_vec_print(const long * vec, long len)
{
    long i;

    printf("%ld", len);
    if (len > 0)
    {
        printf(" ");
        for (i = 0; i < len; i++)
            printf(" %ld", vec[i]);
    }

    return 1;
}

static __inline__ int _perm_print(const long * vec, long n)
{
    long i;

    printf("%ld", n);
    if (n > 0)
    {
        printf(" ");
        for (i = 0; i < n; i++)
            printf(" %ld", vec[i]);
    }

    return 1;
}

#ifdef __cplusplus
}
#endif

#endif
