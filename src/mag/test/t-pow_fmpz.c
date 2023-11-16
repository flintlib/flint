/*
    Copyright (C) 2014 Fredrik Johansson

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

TEST_FUNCTION_START(mag_pow_fmpz, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, w;
        mag_t xb, yb;
        fmpz_t e;

        arf_init(x);
        arf_init(y);
        arf_init(w);
        fmpz_init(e);

        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 80);
        mag_randtest_special(yb, state, 80);
        fmpz_randtest(e, state, 200);

        arf_set_mag(x, xb);

        arf_pow_binexp_fmpz(y, x, e, 2 * MAG_BITS, ARF_RND_UP);
        if (arf_is_nan(y))
            arf_pos_inf(y);

        mag_pow_fmpz(yb, xb, e);
        arf_set_mag(w, yb);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(arf_cmpabs(y, w) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("e = "); fmpz_print(e); flint_printf("\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("w = "); arf_print(w); flint_printf("\n\n");
            flint_abort();
        }

        mag_pow_fmpz(xb, xb, e);

        if (!mag_equal(xb, yb))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_printf("e = "); fmpz_print(e); flint_printf("\n\n");
            mag_print(xb); flint_printf("\n\n");
            mag_print(yb); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(w);
        fmpz_clear(e);

        mag_clear(xb);
        mag_clear(yb);
    }

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, w;
        mag_t xb, yb;
        fmpz_t e;

        arf_init(x);
        arf_init(y);
        arf_init(w);
        fmpz_init(e);

        mag_init(xb);
        mag_init(yb);

        mag_randtest_special(xb, state, 80);
        mag_randtest_special(yb, state, 80);
        fmpz_randtest(e, state, 200);

        arf_set_mag(x, xb);

        arf_pow_binexp_fmpz(y, x, e, 2 * MAG_BITS, ARF_RND_DOWN);
        if (arf_is_nan(y))
            arf_pos_inf(y);

        mag_pow_fmpz_lower(yb, xb, e);
        arf_set_mag(w, yb);

        MAG_CHECK_BITS(xb)
        MAG_CHECK_BITS(yb)

        if (!(arf_cmpabs(w, y) <= 0))
        {
            flint_printf("FAIL\n\n");
            flint_printf("e = "); fmpz_print(e); flint_printf("\n\n");
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("w = "); arf_print(w); flint_printf("\n\n");
            flint_abort();
        }

        mag_pow_fmpz_lower(xb, xb, e);

        if (!mag_equal(xb, yb))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_printf("e = "); fmpz_print(e); flint_printf("\n\n");
            mag_print(xb); flint_printf("\n\n");
            mag_print(yb); flint_printf("\n\n");
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(w);
        fmpz_clear(e);

        mag_clear(xb);
        mag_clear(yb);
    }

    TEST_FUNCTION_END(state);
}
