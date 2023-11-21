/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "bernoulli.h"

FLINT_TLS_PREFIX slong bernoulli_cache_num = 0;

FLINT_TLS_PREFIX fmpq * bernoulli_cache = NULL;

void
bernoulli_cleanup(void)
{
    slong i;

    for (i = 0; i < bernoulli_cache_num; i++)
        fmpq_clear(bernoulli_cache + i);

    flint_free(bernoulli_cache);
    bernoulli_cache = NULL;
    bernoulli_cache_num = 0;
}

void
bernoulli_cache_compute(slong n)
{
    slong old_num = bernoulli_cache_num;

    if (old_num < n)
    {
        slong i, new_num;

        if (old_num == 0)
        {
            flint_register_cleanup_function(bernoulli_cleanup);
        }

        if (n <= 128)
            new_num = FLINT_MAX(old_num + 32, n);
        else
            new_num = FLINT_MAX(old_num + 128, n);

        bernoulli_cache = flint_realloc(bernoulli_cache, new_num * sizeof(fmpq));
        for (i = old_num; i < new_num; i++)
            fmpq_init(bernoulli_cache + i);

        if (new_num <= 128)
        {
            /* todo: use recursion, but only compute new entries */
            arith_bernoulli_number_vec(bernoulli_cache, new_num);
        }
        else
        {
            bernoulli_fmpq_vec_no_cache(bernoulli_cache + old_num, old_num, new_num - old_num);
        }

        bernoulli_cache_num = new_num;
    }
}

