/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

/* Small enough that a numerical computation is safe */
#define DOUBLE_CUTOFF (UWORD(1) << 21)

void
fmpq_dedekind_sum_coprime(fmpq_t s, const fmpz_t h, const fmpz_t k)
{
    if (fmpz_cmp_ui(k, DOUBLE_CUTOFF) < 0)
    {
        double t;

        t = fmpq_dedekind_sum_coprime_d(*h, *k) * (6 * (*k));

        /* Round to nearest after truncation */
        if (t > 0)
            t += 0.5;
        else
            t -= 0.5;

        fmpz_set_d(fmpq_numref(s), t);
        fmpz_set_ui(fmpq_denref(s), UWORD(6) * (*k));
        fmpq_canonicalise(s);
    }
    else
    {
        fmpq_dedekind_sum_coprime_large(s, h, k);
    }
}


