/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"
#include "ball_helpers.h"

/* _mp_real_submul_bounded against mp_real_mul + mp_real_sub: a = b c + d with
   |d| < B^E, at various precisions, exact and inexact operands */
static void
test_submul(flint_rand_t state, slong iters)
{
    slong iter;
    for (iter = 0; iter < iters; iter++)
    {
        mp_real_t a, b, c, d, p, r1, r2;
        arb_t x, y, xa;
        slong n = 1 + n_randint(state, (iter % 10 == 0) ? 300 : 20);
        slong sb = 1 + n_randint(state, 2 * n + 1), sc = 1 + n_randint(state, 2 * n + 1);
        slong E, sd;
        nn_ptr tmp;

        mp_real_init(a); mp_real_init(b); mp_real_init(c); mp_real_init(d);
        mp_real_init(p); mp_real_init(r1); mp_real_init(r2);
        arb_init(x); arb_init(y); arb_init(xa);

        tmp = flint_malloc((2 * n + 2) * sizeof(ulong));
        flint_mpn_rrandom(tmp, state, sb);
        if (tmp[sb - 1] == 0) tmp[sb - 1] = 1;
        _mp_real_set_mpn_2exp(b, tmp, sb, (slong) n_randint(state, 300) - 150);
        flint_mpn_rrandom(tmp, state, sc);
        if (tmp[sc - 1] == 0) tmp[sc - 1] = 1;
        _mp_real_set_mpn_2exp(c, tmp, sc, (slong) n_randint(state, 300) - 150);
        if (n_randint(state, 2)) mp_real_neg(b, b);
        if (n_randint(state, 2)) mp_real_neg(c, c);
        if (n_randint(state, 3) == 0) _mp_real_add_error_ulps(b, 1.0 + n_randint(state, 5));
        if (n_randint(state, 3) == 0) _mp_real_add_error_ulps(c, 1.0 + n_randint(state, 5));

        /* the exact product of the midpoints (with inexact b or c, a
           product of the balls, and a sum with it, would round the
           midpoint to the accuracy of the radius, above d), and a
           small d */
        {
            mp_real_t bm, cm;
            mp_real_init(bm);
            mp_real_init(cm);
            mp_real_set(bm, b);
            mp_real_set(cm, c);
            bm->err = 0;
            cm->err = 0;
            mp_real_mul(p, bm, cm, sb + sc + 2);
            mp_real_clear(bm);
            mp_real_clear(cm);
        }
        sd = 1 + n_randint(state, n + 1);
        flint_mpn_rrandom(tmp, state, sd);
        if (tmp[sd - 1] == 0) tmp[sd - 1] = 1;
        /* d = (tmp, sd) B^(p.exp - n - 1 - sd + t) < B^(p.exp - n - 1 + t) */
        {
            slong t = n_randint(state, 3);
            _mp_real_set_mpn_2exp(d, tmp, sd, FLINT_BITS * (p->exp - n - 1 - sd + t));
            E = p->exp - n - 1 + t;
        }
        if (n_randint(state, 2)) mp_real_neg(d, d);
        if (n_randint(state, 5) == 0) mp_real_zero(d);
        /* a = p + d exactly (p spans sb + sc limbs, d reaches n + 1 + sd
           <= 2 n + 2 limbs below p's top), the radius added after */
        mp_real_add(a, p, d, sb + sc + 2 * n + 8);
        if (n_randint(state, 3) == 0) _mp_real_add_error_ulps(a, 1.0 + n_randint(state, 5));
        /* the bound |a - b c| < B^E on the midpoints: d's magnitude */
        E = E + 1;

        _mp_real_submul_bounded(r1, a, b, c, E, n);
        mp_real_mul(r2, b, c, sb + sc + 2);
        mp_real_sub(r2, a, r2, sb + sc + n + 4);

        mp_real_get_arb(x, r1);
        mp_real_get_arb(y, r2);
        if (!arb_overlaps(x, y))
        {
            flint_printf("FAIL: submul_bounded (iter %wd), n = %wd, E = %wd\n",
                iter, n, E);
            mp_real_print(a); mp_real_print(b); mp_real_print(c);
            mp_real_print(r1); mp_real_print(r2);
            flint_abort();
        }
        /* accuracy: the limbs of the residual above B^(E - n), its
           weight window, less a partial top limb and the 3 units */
        mp_real_get_arb(xa, a);
        if (a->err == 0 && b->err == 0 && c->err == 0
            && r1->size > 0
            && arb_rel_accuracy_bits(x) < FLINT_BITS * (r1->exp - (E - n) - 1) - 8)
        {
            flint_printf("FAIL: submul_bounded accuracy (iter %wd), n = %wd: "
                "%wd bits\n", iter, n, arb_rel_accuracy_bits(x));
            mp_real_print(r1); mp_real_print(r2);
            flint_abort();
        }

        flint_free(tmp);
        mp_real_clear(a); mp_real_clear(b); mp_real_clear(c); mp_real_clear(d);
        mp_real_clear(p); mp_real_clear(r1); mp_real_clear(r2);
        arb_clear(x); arb_clear(y); arb_clear(xa);
    }
}

TEST_FUNCTION_START(mp_real_submul_bounded, state)
{
    test_submul(state, 300 + 300 * flint_test_multiplier());

    TEST_FUNCTION_END(state);
}
