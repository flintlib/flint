/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "mpn_extras.h"
#include "fixed.h"
#include "profiler.h"

#define TIME(expr, res) \
    do { timeit_t __t; slong __r; expr; \
         TIMEIT_REPEAT(__t, __r) { expr; } TIMEIT_END_REPEAT(__t, __r); \
         res = (double) __t->wall * 0.001 / __r; } while (0)

int main(int argc, char * argv[])
{
    flint_rand_t state;
    mp_size_t an, bn, ratio;
    double tg, tn, tu, tp, ti;
    flint_rand_init(state);

    flint_printf("Euclidean division: an x bn, times in seconds\n");
    flint_printf("%8s %8s %10s %10s %10s %10s %10s\n", "bn", "an", "gmp", "newton", "unbal", "preinvn", "inv");

    for (ratio = 1; ratio <= 16; ratio *= 2)
    {
        for (bn = 32; bn <= 32768; bn *= 2)
        {
            mp_ptr a, b, q, r, dinv, a2, q2;
            an = bn + ratio * bn;
            if (ratio == 1) an = 2 * bn;
            if (an > 400000) continue;

            a = flint_malloc(an * sizeof(mp_limb_t));
            a2 = flint_malloc((an + 1) * sizeof(mp_limb_t));
            b = flint_malloc(bn * sizeof(mp_limb_t));
            q = flint_malloc((an - bn + 2) * sizeof(mp_limb_t));
            q2 = flint_malloc((an + 2) * sizeof(mp_limb_t));
            r = flint_malloc(bn * sizeof(mp_limb_t));
            dinv = flint_malloc(bn * sizeof(mp_limb_t));
            flint_mpn_rrandom(a, state, an);
            flint_mpn_rrandom(b, state, bn);
            b[bn - 1] |= (UWORD(1) << (FLINT_BITS - 1));

            TIME(mpn_tdiv_qr(q, r, 0, a, an, b, bn), tg);
            TIME(_flint_mpn_tdiv_qr_newton(q, r, a, an, b, bn), tn);
            if (an > 2 * bn)
                TIME(_flint_mpn_tdiv_qr_unbalanced(q, r, a, an, b, bn), tu);
            else
                tu = 0;
            /* existing FLINT preinvn division (normalised divisor) */
            TIME((flint_mpn_preinvn(dinv, b, bn), flint_mpn_divrem_preinvn(q, a2, a, an, b, bn, dinv)), tp);
            /* reciprocal floor(B^an / b) via inv vs tdiv */
            TIME(flint_mpn_inv(q2, b, bn, an), ti);

            flint_printf("%8wd %8wd %10.3e %10.3e %10.3e %10.3e %10.3e\n", bn, an, tg, tn, tu, tp, ti);
            flint_free(a); flint_free(a2); flint_free(b); flint_free(q); flint_free(q2); flint_free(r); flint_free(dinv);
        }
        flint_printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
