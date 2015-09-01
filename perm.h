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

#ifdef PERM_INLINES_C
#define PERM_INLINE FLINT_DLL
#else
#define PERM_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Memory management *********************************************************/

PERM_INLINE slong * _perm_init(slong n)
{
    slong i, *vec;

    vec = (slong *) flint_malloc(n * sizeof(slong));

    if (!vec)
    {
        flint_printf("ERROR (_perm_init).\n\n");
        abort();
    }

    for (i = 0; i < n; i++)
        vec[i] = i;

    return vec;
}

PERM_INLINE void _perm_clear(slong * vec)
{
    flint_free(vec);
}

/* Assignment ****************************************************************/

PERM_INLINE slong _perm_equal(const slong *vec1, const slong *vec2, slong n)
{
    slong i;

    for (i = 0; i < n; i++)
        if (vec1[i] != vec2[i])
            return 0;

    return 1;
}

PERM_INLINE void _perm_set(slong *res, const slong *vec, slong n)
{
    slong i;

    for (i = 0; i < n; i++)
        res[i] = vec[i];
}

PERM_INLINE void _perm_set_one(slong *vec, slong n)
{
    slong i;

    for (i = 0; i < n; i++)
        vec[i] = i;
}

PERM_INLINE void
 _perm_inv(slong *res, const slong *vec, slong n)
{
    slong i;

    if (res == vec)
    {
        slong *t = (slong *) flint_malloc(n * sizeof(slong));

        if (!t)
        {
            flint_printf("ERROR (_perm_inv).\n\n");
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

PERM_INLINE void
_perm_compose(slong *res, const slong *vec1, const slong *vec2, slong n)
{
    slong i;

    if (res == vec1)
    {
        slong *t = (slong *) flint_malloc(n * sizeof(slong));

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

FLINT_DLL int _perm_randtest(slong * vec, slong n, flint_rand_t state);

/* Parity ********************************************************************/

FLINT_DLL int _perm_parity(const slong * vec, slong n);

/* Input and output **********************************************************/

PERM_INLINE  int _long_vec_print(const slong * vec, slong len)
{
    slong i;

    flint_printf("%wd", len);
    if (len > 0)
    {
        flint_printf(" ");
        for (i = 0; i < len; i++)
            flint_printf(" %wd", vec[i]);
    }

    return 1;
}

PERM_INLINE int _perm_print(const slong * vec, slong n)
{
    slong i;

    flint_printf("%wd", n);
    if (n > 0)
    {
        flint_printf(" ");
        for (i = 0; i < n; i++)
            flint_printf(" %wd", vec[i]);
    }

    return 1;
}

#ifdef __cplusplus
}
#endif

#endif
