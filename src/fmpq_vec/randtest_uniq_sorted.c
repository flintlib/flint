/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpq.h"
#include "fmpq_vec.h"

void _fmpq_vec_randtest_uniq_sorted(fmpq * vec, flint_rand_t state, slong len, flint_bitcnt_t bits)
{
    slong i;
    int do_again;

    /* if 2^bits < len we are too likely to have collision */
    if (4 * n_sizeinbase(len, 2) > bits)
        flint_throw(FLINT_ERROR, "bits too small in %s\n", __func__);

    _fmpq_vec_randtest(vec, state, len, bits);
    if (len <= 1) return;

    do
    {
        do_again = 0;
        _fmpq_vec_sort(vec, len);
        for (i = 0; i < len - 1; i++)
        {
            if (fmpq_equal(vec + i, vec + i + 1))
            {
                fmpq_randtest(vec + i, state, bits);
                do_again = 1;
            }
        }
    } while (do_again);
}
