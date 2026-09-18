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
#include "ulong_extras.h"
#include "fmpq.h"
#include "arb.h"
#include "fixed.h"

/* every relation of every precomputed table holds for the actual
   logarithms and angles (the row sums match the epsilons to a
   relative 1e-12), the Machin-type sets reproduce the logarithms and arguments, and
   the cache/generation interface behaves */

TEST_FUNCTION_START(fixed_rel_tab, state)
{
    slong i, g, j, k;
    slong prec = 800;

    for (g = 0; g < 2; g++)
    {
        for (i = 0; i < _fixed_rel_static_num; i++)
        {
            const fixed_rel_static_struct * st = _fixed_rel_static + i;
            const fixed_rel_struct * t;
            arb_ptr alpha;

            if (st->gaussian != g)
                continue;

            if (!fixed_rel_table_is_cached(g, st->num))
                TEST_FUNCTION_FAIL("precomputed table not reported cached: g = %wd, num = %wd\n", g, st->num);
            t = fixed_rel_table(g, st->num);
            if (t->rows != st->rows || t->d != st->d || !t->is_static)
                TEST_FUNCTION_FAIL("table mismatch: g = %wd, num = %wd\n", g, st->num);

            alpha = _arb_vec_init(st->num);
            if (g)
                _fixed_atan_gauss_vec(alpha, st->num, prec);
            else
                arb_log_primes_vec_bsplit(alpha, st->num, prec);

            for (j = 0; j < t->rows; j++)
            {
                arb_t s;
                double e;
                arb_init(s);
                for (k = 0; k < t->num; k++)
                    arb_addmul_si(s, alpha + k, t->d[j * t->num + k], prec);
                e = arf_get_d(arb_midref(s), ARF_RND_NEAR);
                /* (the epsilons are not strictly monotone: the
                   generator keeps a row whenever it improves on the
                   rows before it, which lets |epsilon| step up
                   occasionally; the descent does not depend on
                   monotonicity) */
                if (!(fabs(e - t->epsilon[j]) <= 1e-12 * fabs(t->epsilon[j]))
                    || fabs(t->epsilon_inv[j] * t->epsilon[j] - 1.0) > 1e-15)
                    TEST_FUNCTION_FAIL("relation: g = %wd, num = %wd, row %wd: sum %g, epsilon %g\n",
                        g, st->num, j, e, t->epsilon[j]);
                arb_clear(s);
            }
            if (t->weights[0] != 0.0f || (g == 0 && t->primes[0] != 2))
                TEST_FUNCTION_FAIL("weights/primes: g = %wd, num = %wd\n", g, st->num);
            _arb_vec_clear(alpha, st->num);
        }

        /* a generated table */
        {
            slong num = 5 + n_randint(state, 10);
            const fixed_rel_struct * t;
            arb_ptr alpha;
            int was = fixed_rel_table_is_cached(g, num);

            t = fixed_rel_table(g, num);
            if (!fixed_rel_table_is_cached(g, num) || t->num != num)
                TEST_FUNCTION_FAIL("generated table: g = %wd, num = %wd\n", g, num);
            if (t->is_static != was && was)
                TEST_FUNCTION_FAIL("static flag: g = %wd, num = %wd\n", g, num);
            alpha = _arb_vec_init(num);
            if (g)
                _fixed_atan_gauss_vec(alpha, num, prec);
            else
                arb_log_primes_vec_bsplit(alpha, num, prec);
            for (j = 0; j < t->rows; j += 7)
            {
                arb_t s;
                double e;
                arb_init(s);
                for (k = 0; k < num; k++)
                    arb_addmul_si(s, alpha + k, t->d[j * num + k], prec);
                e = arf_get_d(arb_midref(s), ARF_RND_NEAR);
                if (!(fabs(e - t->epsilon[j]) <= 1e-12 * fabs(t->epsilon[j])))
                    TEST_FUNCTION_FAIL("generated relation: g = %wd, num = %wd, row %wd\n", g, num, j);
                arb_clear(s);
            }
            _arb_vec_clear(alpha, num);
        }
    }

    /* the Machin-type sets */
    for (g = 0; g < 2; g++)
    {
        slong num;
        for (num = 2; num <= fixed_machin_table_max(g); num++)
        {
            const fixed_machin_struct * mt = fixed_machin_table(g, num);
            arb_ptr y;
            fmpz_t p, q;

            /* the chosen set has exactly num terms whenever one is
               tabulated for that size; otherwise it is the largest
               below (few followups) or the next one up */
            if (mt == NULL || (num <= 32 && mt->num != FLINT_MAX(num, g ? 3 : 4)))
                TEST_FUNCTION_FAIL("machin lookup: g = %wd, num = %wd, got %wd\n",
                    g, num, mt == NULL ? WORD(-1) : mt->num);
            if (num != mt->num)
                continue;   /* checked at its own size */

            y = _arb_vec_init(mt->num);
            fmpz_init(p); fmpz_init(q);
            for (j = 0; j < mt->num; j++)
            {
                fmpz_one(p);
                fixed_machin_get_x(q, mt, j);
                arb_atan_frac_bsplit(y + j, p, q, g == 0, prec);
            }
            for (i = 0; i < mt->num; i++)
            {
                arb_t s, r;
                arb_init(s); arb_init(r);
                {
                    fmpz_t cc;
                    fmpz_init(cc);
                    arb_zero(s);
                    for (j = 0; j < mt->num; j++)
                    {
                        fixed_machin_get_c(cc, mt, i, j);
                        arb_addmul_fmpz(s, y + j, cc, prec);
                    }
                    arb_div_ui(s, s, mt->den, prec);
                    fmpz_clear(cc);
                }
                if (g == 0)
                    arb_log_ui(r, n_nth_prime(i + 1), prec);
                else
                {
                    fmpq_t fr;
                    fmpq_init(fr);
                    fmpq_set_si(fr, _fixed_gaussian_primes[2 * i + 1], _fixed_gaussian_primes[2 * i]);
                    arb_set_fmpq(r, fr, prec);
                    arb_atan(r, r, prec);
                    fmpq_clear(fr);
                }
                if (!arb_overlaps(s, r) || arb_rel_accuracy_bits(s) < prec - 100)
                    TEST_FUNCTION_FAIL("machin set: g = %wd, num = %wd, entry %wd\n", g, mt->num, i);
                arb_clear(s); arb_clear(r);
            }
            _arb_vec_clear(y, mt->num);
            fmpz_clear(p); fmpz_clear(q);
        }
    }

    TEST_FUNCTION_END(state);
}
