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
#define PERM_INLINE static __inline__
#endif

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Memory management *********************************************************/

slong * _perm_init(slong n);
void _perm_clear(slong * vec);

/* Assignment ****************************************************************/

PERM_INLINE
slong _perm_equal(const slong * vec1, const slong * vec2, slong n)
{
    slong i;

    for (i = 0; i < n; i++)
        if (vec1[i] != vec2[i])
            return 0;

    return 1;
}

PERM_INLINE
void _perm_set(slong * res, const slong * vec, slong n)
{
    slong i;

    for (i = 0; i < n; i++)
        res[i] = vec[i];
}

PERM_INLINE
void _perm_set_one(slong * vec, slong n)
{
    slong i;

    for (i = 0; i < n; i++)
        vec[i] = i;
}

void _perm_inv(slong * res, const slong * vec, slong n);

/* Composition ***************************************************************/

void _perm_compose(slong * res, const slong * vec1, const slong * vec2, slong n);

/* Randomisation *************************************************************/

int _perm_randtest(slong * vec, slong n, flint_rand_t state);

/* Parity ********************************************************************/

int _perm_parity(const slong * vec, slong n);

/* Input and output **********************************************************/

int _perm_print(const slong * vec, slong n);

#ifdef __cplusplus
}
#endif

#endif
