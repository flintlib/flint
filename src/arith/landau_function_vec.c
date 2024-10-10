/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "arith.h"

#ifdef __GNUC__
# define log __builtin_log
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

void arith_landau_function_vec(fmpz * res, slong len)
{
    ulong p, pmax;
    ulong pk, pkhi;
    fmpz_t a;
    slong k;
    ulong n;

    if (len < 1)
        return;

    for (k = 0; k < len; k++)
        fmpz_one(res + k);

    pmax = 1.328 * sqrt(len*log(len) + 1);

    fmpz_init(a);

    for (p = UWORD(2); p <= pmax; p = n_nextprime(p, 0))
    {
        for (n = len - 1; n >= p; n--)
        {
            pk = p;
            pkhi = UWORD(0);

            for (k = 1; k <= len; k++)
            {
                if (pk > n || pkhi)
                    break;

                fmpz_mul_ui(a, res + n - pk, pk);
                if (fmpz_cmp(res + n, a) < 0)
                    fmpz_set(res + n, a);

                umul_ppmm(pkhi, pk, pk, p);
            }
        }
    }

    fmpz_clear(a);
}
