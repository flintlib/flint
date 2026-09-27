/*
    Copyright (C) 2015, Elias Tsigaridas
    Copyright (C) 2016, Vincent Delecroix
    Copyright (C) 2026, Fredrik Johansson

    The implementation was inspired from the SLV library version 0.5 by Elias
    Tsigaridas (namely the function Descartes_test in the file vca_solver_1.c
    lines 67-125)

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/* We compute the Taylor shift of the reversal of p row by row (as in
   Horner's scheme), checking the sign variations of the finished
   coefficients as we go so that we can stop early. The coefficients are
   stored as d-limb two's complement integers: all intermediate values
   are sums of coefficients of p times binomial coefficients, so they
   have at most bits + len + 1 bits. */

#define SGN(j) \
    ((d == 1) ? FLINT_SGN((slong) q[j]) : \
    (((slong) q[d * (j) + d - 1] < 0) ? -1 : (flint_mpn_zero_p(q + d * (j), d) ? 0 : 1)))

/* q[j + 1] += q[j] */
#define ADD_NEXT(j) \
    do { \
        if (d == 1) \
            q[(j) + 1] += q[j]; \
        else \
            mpn_add_n(q + d * ((j) + 1), q + d * ((j) + 1), q + d * (j), d); \
    } while (0)

static slong
_fmpz_poly_descartes_bound_0_1_mpn(const fmpz * p, slong len, slong bound, slong d)
{
    slong V = 0;
    slong i, j;
    int s, t;
    slong deg = len - 1;
    nn_ptr q;
    TMP_INIT;

    TMP_START;
    q = TMP_ALLOC(d * len * sizeof(ulong));

    for (j = 0; j < len; j++)
        fmpz_get_signed_ui_array(q + d * j, d, p + j);

    for (j = 0; j <= deg - 1; j++)
        ADD_NEXT(j);

    s = SGN(deg);  /* = sign(p(1)) */

    for (i = 1; i <= deg - 1; i++)
    {
        j = deg - i;
        t = s;
        while ((j >= 0) && (t == 0))
        {
            t = SGN(j);
            j--;
        }

        while ((j >= 0) && ((SGN(j) == t) || (SGN(j) == 0)))
            j--;

        if (j < 0)
        {
            /* all coefficients of q have the same sign */
            TMP_END;
            return V;
        }

        for (j = 0; j <= deg - i - 1; j++)
            ADD_NEXT(j);

        if (s == 0)
            s = SGN(deg - i);
        else if (s == -SGN(deg - i))
        {
            if (V == bound)
            {
                TMP_END;
                return WORD_MAX;
            }
            V++;
            s = -s;
        }
    }

    if (s == -SGN(0))
    {
        if (V == bound)
        {
            TMP_END;
            return WORD_MAX;
        }
        V++;
    }

    TMP_END;
    return V;
}

slong _fmpz_poly_descartes_bound_0_1(const fmpz * p, slong len, slong bound)
{
    slong j, bits;
    int t;

    if (len <= 1)
        return 0;

    /* quick exit: all coefficients have the same sign */
    j = len - 1;
    t = fmpz_sgn(p + j);
    while ((j >= 0) && ((fmpz_sgn(p + j) == t) || fmpz_sgn(p + j) == 0))
        j--;
    if (j < 0)
        return 0;

    bits = FLINT_ABS(_fmpz_vec_max_bits(p, len));
    bits = bits + len + 1;

    return _fmpz_poly_descartes_bound_0_1_mpn(p, len, bound,
        (bits + FLINT_BITS - 1) / FLINT_BITS);
}
