/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <float.h>
#include "test_helpers.h"
#include "arb.h"
#include "mp_real.h"
#include "dfloat.h"
#include "ball_helpers.h"

/* mp_real_set_dfloat: the ball contains the dfloat ball (its exact
   endpoints), with the exact sum of the components as the midpoint when
   the radius is zero (compared with _dfloat_get_arb);
   mp_real_get_dfloat: the dfloat ball contains the mp_real ball, each
   component is the double nearest to the remainder (|x_{i+1}| <=
   ulp(x_i)/2, no overlap), and for an exact input of a modest size the
   radius is below 2^-(53 n - 3) relative; out of range gives the whole
   line; round trips enclose. */

/* b contains the exact endpoints of x (mp_real_get_arb would round the
   radius up to 30 bits) */
static int
_dfc_contains_exact(const arb_t b, const mp_real_t x)
{
    arf_t m, r, t;
    int ok;
    mp_real_t y;
    arb_t a;

    arf_init(m); arf_init(r); arf_init(t);
    mp_real_init(y); arb_init(a);
    mp_real_set(y, x);
    y->err = 0;
    mp_real_get_arb(a, y);          /* the midpoint, exactly */
    arf_set(m, arb_midref(a));
    arf_set_ui(r, x->err);
    arf_mul_2exp_si(r, r, FLINT_BITS * (x->exp - x->size));
    arf_add(t, m, r, ARF_PREC_EXACT, ARF_RND_DOWN);
    ok = arb_contains_arf(b, t);
    arf_sub(t, m, r, ARF_PREC_EXACT, ARF_RND_DOWN);
    ok = ok && arb_contains_arf(b, t);
    arf_clear(m); arf_clear(r); arf_clear(t);
    mp_real_clear(y); arb_clear(a);
    return ok;
}

TEST_FUNCTION_START(mp_real_dfloat_conv, state)
{
    slong iter;

    if (!dfloat_is_supported())
    {
        FLINT_TEST_CLEAR(state);
        printf("%.*s(" _YELLOW_B "SKIPPED" _RESET ")\n", 54, _test_io_string_);
        return 0;
    }

    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        slong n = 1 + n_randint(state, 4), i;
        double d[5], r, e[5], er;
        mp_real_t x;
        arb_t a, b;

        mp_real_init(x);
        arb_init(a);
        arb_init(b);

        /* a random dfloat ball (finite) */
        switch (n)
        {
            case 1: { d1b_t t; d1b_randtest(t, state); for (i = 0; i < 1; i++) d[i] = t->d[i]; r = t->rad; break; }
            case 2: { d2b_t t; d2b_randtest(t, state); for (i = 0; i < 2; i++) d[i] = t->d[i]; r = t->rad; break; }
            case 3: { d3b_t t; d3b_randtest(t, state); for (i = 0; i < 3; i++) d[i] = t->d[i]; r = t->rad; break; }
            default: { d4b_t t; d4b_randtest(t, state); for (i = 0; i < 4; i++) d[i] = t->d[i]; r = t->rad; break; }
        }
        if (n_randint(state, 3) == 0)
            r = 0.0;
        if (!(r <= DBL_MAX))
            r = 1.0;
        for (i = 0; i < n; i++)
            if (!(fabs(d[i]) <= DBL_MAX))
                d[i] = 0.0;

        mp_real_set_dfloat(x, d, n, r);
        mp_real_get_arb(a, x);
        _dfloat_get_arb(b, d, (int) n, 0.0);
        {
            /* the exact endpoints mid +- r */
            arf_t lo, hi, rr;
            int ok;
            arf_init(lo); arf_init(hi); arf_init(rr);
            arf_set_d(rr, r);
            arf_sub(lo, arb_midref(b), rr, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_add(hi, arb_midref(b), rr, ARF_PREC_EXACT, ARF_RND_DOWN);
            ok = arb_contains_arf(a, lo) && arb_contains_arf(a, hi)
                && (r != 0.0 || arf_equal(arb_midref(a), arb_midref(b)));
            arf_clear(lo); arf_clear(hi); arf_clear(rr);
            if (!ok)
                arb_zero(b);    /* force the failure below */
            else
                arb_set(b, a);
        }
        if (!arb_contains(a, b) || (r == 0.0 && !arb_is_exact(a)))
        {
            flint_printf("FAIL: set_dfloat, n = %wd\n", n);
            arb_printd(a, 40); flint_printf("\n");
            arb_printd(b, 40); flint_printf("\n");
            flint_abort();
        }

        /* the round trip */
        mp_real_get_dfloat(e, &er, n, x);
        _dfloat_get_arb(b, e, (int) n, er);
        if (er < HUGE_VAL && !_dfc_contains_exact(b, x))
        {
            flint_printf("FAIL: round trip, n = %wd\n", n);
            arb_printd(a, 40); flint_printf("\n");
            arb_printd(b, 40); flint_printf("\n");
            flint_abort();
        }

        mp_real_clear(x);
        arb_clear(a);
        arb_clear(b);
    }

    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        slong n = 1 + n_randint(state, 4), i;
        double e[4], er;
        mp_real_t x;
        arb_t a, b;
        int exact;

        mp_real_init(x);
        arb_init(a);
        arb_init(b);

        mp_real_randtest(x, state, 1 + n_randint(state, 8), n_randint(state, 2));
        if (n_randint(state, 4) == 0)
            x->exp += (slong) n_randint(state, 40) - 20;       /* out of range */
        exact = (x->err == 0);
        mp_real_get_arb(a, x);

        mp_real_get_dfloat(e, &er, n, x);

        if (er == HUGE_VAL)
        {
            /* only beyond the range, and as the whole line */
            for (i = 0; i < n; i++)
                if (e[i] != 0.0)
                    TEST_FUNCTION_FAIL("whole line with nonzero components\n");
            if (x->size != 0 && FLINT_BITS * (x->exp - 1) < 1000
                && mp_real_abs_bound_lt_2exp_si(x) < 1000)
                TEST_FUNCTION_FAIL("spurious whole line\n");
        }
        else
        {
            _dfloat_get_arb(b, e, (int) n, er);
            if (!_dfc_contains_exact(b, x))
            {
                flint_printf("FAIL: get_dfloat containment, n = %wd\n", n);
                arb_printd(a, 40); flint_printf("\n");
                arb_printd(b, 40); flint_printf("\n");
                flint_abort();
            }
            for (i = 0; i + 1 < n; i++)
            {
                int ex;
                if (e[i] == 0.0)
                {
                    if (e[i + 1] != 0.0)
                        TEST_FUNCTION_FAIL("zero component before a nonzero one\n");
                    continue;
                }
                frexp(e[i], &ex);
                if (fabs(e[i + 1]) > ldexp(1.0, ex - 54))
                    TEST_FUNCTION_FAIL("overlapping components\n");
            }
            /* accuracy for exact inputs well inside the normal range */
            if (exact && x->size != 0 && FLINT_BITS * (x->exp - 1) < 700
                && FLINT_BITS * (x->exp - 1) > -700 && !arb_is_exact(b))
            {
                arb_t c;
                arb_init(c);
                arb_set_arf(c, arb_midref(b));
                mag_set_d(arb_radref(c), er);
                if (arb_rel_accuracy_bits(b) < 53 * n - 3)
                {
                    flint_printf("FAIL: get_dfloat accuracy, n = %wd, acc %wd\n", n, arb_rel_accuracy_bits(b));
                    arb_printd(a, 60); flint_printf("\n");
                    arb_printd(b, 60); flint_printf("\n");
                    flint_abort();
                }
                arb_clear(c);
            }
        }

        mp_real_clear(x);
        arb_clear(a);
        arb_clear(b);
    }

    /* special values */
    {
        double d[2] = { 0.0, 0.0 }, e[2], er;
        mp_real_t x;
        mp_real_init(x);
        mp_real_set_dfloat(x, d, 2, 0.0);
        if (x->size != 0 || x->err != 0)
            TEST_FUNCTION_FAIL("zero\n");
        mp_real_get_dfloat(e, &er, 2, x);
        if (e[0] != 0.0 || e[1] != 0.0 || er != 0.0)
            TEST_FUNCTION_FAIL("zero back\n");
        d[0] = 0x1p-1074; d[1] = 0.0;
        mp_real_set_dfloat(x, d, 2, 0.0);
        mp_real_get_dfloat(e, &er, 2, x);
        if (!(e[0] == 0.0 && er >= 0x1p-1074))
            TEST_FUNCTION_FAIL("subnormal\n");
        d[0] = DBL_MAX; d[1] = 0x1p969;
        mp_real_set_dfloat(x, d, 2, 0.0);
        mp_real_get_dfloat(e, &er, 2, x);
        if (!(e[0] == DBL_MAX && e[1] == 0x1p969 && er == 0.0))
            TEST_FUNCTION_FAIL("largest\n");
        /* DBL_MAX + 2^970 rounds to 2^1024 in the first component */
        d[1] = 0x1p970;
        mp_real_set_dfloat(x, d, 2, 0.0);
        mp_real_get_dfloat(e, &er, 2, x);
        if (!(e[0] == 0.0 && er == HUGE_VAL))
            TEST_FUNCTION_FAIL("overflow\n");
        mp_real_clear(x);
    }

    TEST_FUNCTION_END(state);
}
