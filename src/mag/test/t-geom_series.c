/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"
#include "mag.h"

/* Defined in t-binpow_uiui.c, t-exp_tail.c, t-geom_series.c, t-pow_ui.c and
 * t-pow_fmpz.c */
#ifndef arf_pow_binexp_fmpz
#define arf_pow_binexp_fmpz arf_pow_binexp_fmpz
void
arf_pow_binexp_fmpz(arf_t y, const arf_t b, const fmpz_t e,
    slong prec, arf_rnd_t rnd)
{
    slong i, wp, bits;

    if (fmpz_is_zero(e))
    {
        arf_set_ui(y, UWORD(1));
        return;
    }

    if (fmpz_sgn(e) < 0)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_neg(f, e);
        arf_pow_binexp_fmpz(y, b, f, prec + 2,
            (rnd == ARF_RND_FLOOR || rnd == ARF_RND_DOWN)
            ? ARF_RND_UP : ARF_RND_DOWN);
        arf_ui_div(y, UWORD(1), y, prec, rnd);
        fmpz_clear(f);
        return;
    }

    if (y == b)
    {
        arf_t t;
        arf_init(t);
        arf_set(t, b);
        arf_pow_binexp_fmpz(y, t, e, prec, rnd);
        arf_clear(t);
        return;
    }

    arf_set(y, b);

    bits = fmpz_bits(e);
    wp = ARF_PREC_ADD(prec, bits);

    for (i = bits - 2; i >= 0; i--)
    {
        arf_mul(y, y, y, wp, rnd);
        if (fmpz_tstbit(e, i))
            arf_mul(y, y, b, wp, rnd);
    }
}
#endif

/* Defined in t-binpow_uiui.c, t-exp_tail.c, t-geom_series.c, t-pow_ui.c and
 * t-pow_fmpz.c */
#ifndef arf_pow_binexp_ui
#define arf_pow_binexp_ui arf_pow_binexp_ui
void
arf_pow_binexp_ui(arf_t y, const arf_t b, ulong e, slong prec, arf_rnd_t rnd)
{
    fmpz_t f;
    fmpz_init_set_ui(f, e);
    arf_pow_binexp_fmpz(y, b, f, prec, rnd);
    fmpz_clear(f);
}
#endif

/* Defined in t-binpow_uiui.c, t-exp_tail.c, t-geom_series.c, t-pow_ui.c and
 * t-pow_fmpz.c */
#ifndef arf_pow_binexp_si
#define arf_pow_binexp_si arf_pow_binexp_si
void
arf_pow_binexp_si(arf_t y, const arf_t b, slong e, slong prec, arf_rnd_t rnd)
{
    fmpz_t f;
    fmpz_init(f);
    fmpz_set_si(f, e);
    arf_pow_binexp_fmpz(y, b, f, prec, rnd);
    fmpz_clear(f);
}
#endif

TEST_FUNCTION_START(mag_geom_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, t, y, z;
        mag_t xb, yb;
        ulong N, k;

        arf_init(x);
        arf_init(t);
        arf_init(y);
        arf_init(z);
        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 6);
        mag_randtest_special(yb, state, 6);
        N = n_randint(state, 100);

        mag_geom_series(yb, xb, N);

        arf_set_mag(x, xb);
        arf_set_mag(y, yb);

        arf_pow_binexp_ui(t, x, N, MAG_BITS, ARF_RND_DOWN);
        arf_set(z, t);

        for (k = 1; k < 50; k++)
        {
            arf_mul(t, t, x, MAG_BITS, ARF_RND_DOWN);
            arf_add(z, z, t, MAG_BITS, ARF_RND_DOWN);
        }

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(arf_cmpabs(z, y) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("N = %wu\n\n", N);
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            flint_abort();
        }

        mag_geom_series(xb, xb, N);

        if (!mag_equal(xb, yb))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(t);
        arf_clear(y);
        arf_clear(z);
        mag_clear(xb);
        mag_clear(yb);
    }

    TEST_FUNCTION_END(state);
}
