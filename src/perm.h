/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef PERM_H
#define PERM_H

#ifdef PERM_INLINES_C
#define PERM_INLINE
#else
#define PERM_INLINE static inline
#endif


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
        flint_throw(FLINT_ERROR, "ERROR (_perm_init).\n\n");
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

PERM_INLINE void _perm_one(slong *vec, slong n)
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
            flint_throw(FLINT_ERROR, "ERROR (_perm_inv).\n\n");
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

int _perm_randtest(slong * vec, slong n, flint_rand_t state);

/* Parity ********************************************************************/

int _perm_parity(const slong * vec, slong n);

#define _long_vec_print _Pragma("GCC error \"'_long_vec_print(vec, len)' is deprecated. Use 'flint_printf(\"%{slong*}\", vec, len)' instead.\"")
#define _perm_print _Pragma("GCC error \"'_perm_print(vec, len)' is deprecated. Use 'flint_printf(\"%{slong*}\", vec, len)' instead.\"")
#define _perm_set_one _Pragma("GCC error \"'_perm_set_one' is deprecated. Use '_perm_one' instead.\"")

#ifdef __cplusplus
}
#endif

#endif
