/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "perm.h"

void
_perm_compose(slong * res, const slong * vec1, const slong * vec2, slong n)
{
    slong i;

    if (res == vec1)
    {
        slong * t = (slong *) flint_malloc(n * sizeof(slong));

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

slong * _perm_init(slong n)
{
    slong i, *vec;

    vec = (slong *) flint_malloc(n * sizeof(slong));

    for (i = 0; i < n; i++)
        vec[i] = i;

    return vec;
}

void _perm_clear(slong * vec)
{
    flint_free(vec);
}

void
 _perm_inv(slong * res, const slong * vec, slong n)
{
    slong i;

    if (res == vec)
    {
        slong * t = (slong *) flint_malloc(n * sizeof(slong));

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

int
_perm_parity(const slong *vec, slong n)
{
    slong i, k;
    int * encountered;
    int parity;
    TMP_INIT;

    if (n <= 1)
        return 0;

    TMP_START;

    parity = 0;
    encountered = (int *) TMP_ALLOC(n*sizeof(int));

    for (i = 0; i < n; i++)
       encountered[i] = 0;

    for (i = 0; i < n; i++)
    {
        if (encountered[i] != 0)
        {
            parity ^= 1;
        }
        else
        {
            k = i;
            do
            {
                k = vec[k];
                encountered[k] = 1;
            }
            while (k != i);
        }
    }

    TMP_END;

    return parity;
}

int
_perm_print(const slong * vec, slong n)
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

int _perm_randtest(slong *vec, slong n, flint_rand_t state)
{
    slong i, j, t;

    int parity = 0;

    for (i = 0; i < n; i++)
        vec[i] = i;

    for (i = n - 1; i > 0; i--)
    {
        j = n_randint(state, i + 1);
        parity ^= (i != j);
        t = vec[i];
        vec[i] = vec[j];
        vec[j] = t;
    }

    return parity;
}
