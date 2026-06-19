/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2026 Ricardo Buring
    Copyright (C) 2026 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
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

PERM_INLINE slong _perm_is_one(const slong *vec, slong n)
{
    slong i;

    for (i = 0; i < n; i++)
        if (vec[i] != i)
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

// if p is encoded by vec1 and q is encoded by vec2, then p*q is encoded by res
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

// if a is encoded by vec1 and b is encoded by vec2, then (a^-1)*b is encoded by res
PERM_INLINE void
_perm_compose_inv1(slong *res, const slong *vec1, const slong *vec2, slong n)
{
    slong i;

    // FIXME: do this without allocating memory, if non-aliasing
    if (vec1 == vec2)
    {
        _perm_one(res, n);
    }
    if (res == vec2)
    {
        slong *t = (slong *) flint_malloc(n * sizeof(slong));

        for (i = 0; i < n; i++)
            t[vec1[i]] = i;
        for (i = 0; i < n; i++)
            res[i] = t[vec2[i]];

        flint_free(t);
    }
    else
    {
        _perm_inv(res, vec1, n);
        _perm_compose(res, res, vec2, n);
    }
}

// if a is encoded by vec1 and b is encoded by vec2, then a*(b^-1) is encoded by res
PERM_INLINE void
_perm_compose_inv2(slong *res, const slong *vec1, const slong *vec2, slong n)
{
    slong i;

    if (vec1 == vec2)
    {
        _perm_one(res, n);
    }
    else if (res == vec1)
    {
        slong *t = (slong *) flint_malloc(n * sizeof(slong));

        for (i = 0; i < n; i++)
            t[i] = vec1[i];
        for (i = 0; i < n; i++)
            res[vec2[i]] = t[i];

        flint_free(t);
    }
    else if (res == vec2)
    {
        slong *t = (slong *) flint_malloc(n * sizeof(slong));

        for (i = 0; i < n; i++)
            t[vec2[i]] = i;
        for (i = 0; i < n; i++)
            res[i] = vec1[t[i]];

        flint_free(t);
    }
    else
    {
        for (i = 0; i < n; i++)
            res[vec2[i]] = vec1[i];
    }
}

/* Randomisation *************************************************************/

int _perm_randtest(slong * vec, slong n, flint_rand_t state);

/* Parity ********************************************************************/

int _perm_parity(const slong * vec, slong n);

/* Iteration *****************************************************************/

PERM_INLINE int _perm_next_lex(slong *perm, slong n)
{
    slong k = n-2;
    while (k >= 0 && perm[k] >= perm[k+1])
        k--;
    if (k < 0)
        return 0;
    slong l = n-1;
    while (l >= k && perm[k] >= perm[l])
        l--;
    FLINT_SWAP(slong, perm[k], perm[l]);
    slong i = k+1, j = n-1;
    while (i < j) {
        FLINT_SWAP(slong, perm[i], perm[j]);
        i++;
        j--;
    }
    return 1;
}

PERM_INLINE int _perm_next_heap(slong *perm, slong n, slong *stack)
{
    slong sp = 1;
    while (sp < n && stack[sp] >= sp)
        stack[sp++] = 0;
    if (sp == n) return 0;
    FLINT_SWAP(slong, perm[sp % 2 ? stack[sp] : 0], perm[sp]);
    stack[sp]++;
    return 1;
}

#define _long_vec_print _Pragma("GCC error \"'_long_vec_print(vec, len)' is deprecated. Use 'flint_printf(\"%{slong*}\", vec, len)' instead.\"")
#define _perm_print _Pragma("GCC error \"'_perm_print(vec, len)' is deprecated. Use 'flint_printf(\"%{slong*}\", vec, len)' instead.\"")
#define _perm_set_one _Pragma("GCC error \"'_perm_set_one' is deprecated. Use '_perm_one' instead.\"")

#ifdef __cplusplus
}
#endif

#endif
