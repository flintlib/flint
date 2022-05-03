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
#define PERM_INLINE FLINT_DLL
#else
#define PERM_INLINE static __inline__
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
        flint_throw(FLINT_ALLOC, "Could not allocate in _perm_init\n");

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
            flint_throw(FLINT_ALLOC, "Could not allocate in _perm_inv\n");

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

#if defined (FILE)                  \
  || defined (H_STDIO)              \
  || defined (_H_STDIO)             \
  || defined (_STDIO_H)             \
  || defined (_STDIO_H_)            \
  || defined (__STDIO_H)            \
  || defined (__STDIO_H__)          \
  || defined (_STDIO_INCLUDED)      \
  || defined (__dj_include_stdio_h_)\
  || defined (_FILE_DEFINED)        \
  || defined (__STDIO__)            \
  || defined (_MSL_STDIO_H)         \
  || defined (_STDIO_H_INCLUDED)    \
  || defined (_ISO_STDIO_ISO_H)     \
  || defined (__STDIO_LOADED)       \
  || defined (_STDIO)               \
  || defined (__DEFINED_FILE)

#include "flint-impl.h"

PERM_INLINE  int _long_vec_print(const slong * vec, slong len)
{
    slong i;

    printf(WORD_FMT "d", len);
    if (len > 0)
    {
        printf(" ");
        for (i = 0; i < len; i++)
            printf(" " WORD_FMT "d", vec[i]);
    }

    return 1;
}

PERM_INLINE int _perm_print(const slong * vec, slong n)
{
    slong i;

    printf(WORD_FMT "d", n);
    if (n > 0)
    {
        printf(" ");
        for (i = 0; i < n; i++)
            printf(" " WORD_FMT "d", vec[i]);
    }

    return 1;
}
#endif

#ifdef __cplusplus
}
#endif

#endif
