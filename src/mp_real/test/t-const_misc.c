/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"

static const mp_real_const_func funcs[] = {
    mp_real_const_e,
    mp_real_const_log10,
    mp_real_const_catalan,
    mp_real_const_zeta3,
    mp_real_const_zeta5,
    mp_real_const_gamma_1_3,
    mp_real_const_gamma_1_4,
};

static const char * names[] = { "e", "log10", "catalan", "zeta3", "zeta5",
    "gamma(1/3)", "gamma(1/4)" };

/* 120 digits each */
static const char * refs[] = {
    "2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992",
    "2.30258509299404568401799145468436420760110148862877297603332790096757260967735248023599720508959829834196778404228624863",
    "0.915965594177219015054603514932384110774149374281672134266498119621763019776254769479356512926115106248574422619196199579",
    "1.20205690315959428539973816151144999076498629234049888179227155534183820578631309018645587360933525814619915779526071942",
    "1.03692775514336992633136548645703416805708091950191281197419267790380358978628148456004310655713333637962034146655660904",
    "2.67893853470774763365569294097467764412868937795730110095042832759041761016774381954098288904118878941915904920007226334",
    "3.62560990822190831193068515586767200299516768288006546743337799956991924353872912161836013672338430036147175139242071997",
};

#define NCONST 7

TEST_FUNCTION_START(mp_real_const_misc, state)
{
    slong iter, i;
    arb_t ref[NCONST];

    for (i = 0; i < NCONST; i++)
    {
        arb_init(ref[i]);
        arb_set_str(ref[i], refs[i], 420);
        /* the strings are exact to 120 digits */
        arb_add_error_2exp_si(ref[i], -390);
    }

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        slong n = 2 + n_randint(state, (iter % 10 == 0) ? 120 : 12);
        slong m = n + 1 + n_randint(state, 8);
        mp_real_t x, y;
        arb_t a, b;

        i = n_randint(state, NCONST);

        flint_set_num_threads(1 + n_randint(state, 4));

        mp_real_init(x);
        mp_real_init(y);
        arb_init(a);
        arb_init(b);

        funcs[i](x, n, 0);
        funcs[i](y, m, n_randint(state, 2));
        mp_real_get_arb(a, x);
        mp_real_get_arb(b, y);

        /* the reference digits */
        if (!arb_overlaps(a, ref[i]))
        {
            flint_printf("FAIL: %s does not match the reference\n", names[i]);
            flint_printf("n = %wd\n", n);
            flint_printf("x = "); arb_printn(a, 60, 0); flint_printf("\n");
            flint_abort();
        }

        /* two precisions agree */
        if (!arb_overlaps(a, b))
        {
            flint_printf("FAIL: %s at two precisions\n", names[i]);
            flint_printf("n = %wd, m = %wd\n", n, m);
            flint_printf("x = "); arb_printn(a, 60, 0); flint_printf("\n");
            flint_printf("y = "); arb_printn(b, 60, 0); flint_printf("\n");
            flint_abort();
        }

        /* and to the advertised accuracy: about FLINT_BITS (n - 1) bits,
           the top limb of a ball being possibly short, less the few
           bits accumulated by the splitting (measured: at most 11 for
           Catalan's constant and zeta(3) up to n = 2000) */
        if (arb_rel_accuracy_bits(a) < FLINT_BITS * (n - 1) - 24)
        {
            flint_printf("FAIL: %s is inaccurate\n", names[i]);
            flint_printf("n = %wd, accuracy = %wd\n", n,
                arb_rel_accuracy_bits(a));
            flint_abort();
        }

        mp_real_clear(x);
        mp_real_clear(y);
        arb_clear(a);
        arb_clear(b);
    }

    /* mp_real_root_ui and mp_real_rroot_ui against arb_root_ui */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        /* k up to 2^40 (every 13th: the series coefficients and their
           common denominator beyond a word) */
        ulong k = 1 + n_randint(state, (iter % 13 == 0) ? (UWORD(1) << 40) - 1
            : (iter % 7 == 0) ? 3000 : 20);
        /* many trailing zero bits: the generic denominator
           k^(r-1) (r-1)! then vanishes modulo the word */
        if (iter % 13 == 0)
            k = FLINT_MAX(2, (k >> 20) << (2 + n_randint(state, 19)));
        slong n = 2 + n_randint(state, (iter % 10 == 0 && k >= 4) ? 200 : 30);
        int recip = n_randint(state, 2);
        mp_real_t x, y;
        arb_t a, b, c;

        mp_real_init(x);
        mp_real_init(y);
        arb_init(a);
        arb_init(b);
        arb_init(c);

        mp_real_const_log2(x, n, 0);
        mp_real_mul_ui(x, x, 1 + n_randint(state, 1000), n);
        mp_real_mul_2exp_si(x, x, (slong) n_randint(state, 400) - 200);

        if (recip)
            mp_real_rroot_ui(y, x, k, n);
        else
            mp_real_root_ui(y, x, k, n);
        mp_real_get_arb(a, x);
        mp_real_get_arb(b, y);
        arb_root_ui(c, a, k, FLINT_BITS * n + 64);
        if (recip)
            arb_inv(c, c, FLINT_BITS * n + 64);

        if (!arb_overlaps(b, c))
        {
            flint_printf("FAIL: mp_real_root_ui k = %wu, n = %wd\n", k, n);
            flint_printf("y = "); arb_printn(b, 40, 0); flint_printf("\n");
            flint_printf("c = "); arb_printn(c, 40, 0); flint_printf("\n");
            flint_abort();
        }

        /* the square and cube roots bound the propagated radius and
           their Newton error by limb magnitudes and can lose a limb and
           a half on an inexact operand; the iteration for k >= 4 keeps
           everything within a few bits */
        if (arb_rel_accuracy_bits(b)
            < FLINT_MIN(arb_rel_accuracy_bits(a), FLINT_BITS * (n - 1))
                - ((k <= 3) ? 100 : 32))
        {
            flint_printf("FAIL: mp_real_root_ui accuracy, k = %wu, n = %wd\n",
                k, n);
            flint_printf("in  %wd, out %wd\n", arb_rel_accuracy_bits(a),
                arb_rel_accuracy_bits(b));
            flint_abort();
        }

        mp_real_clear(x);
        mp_real_clear(y);
        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
    }

    for (i = 0; i < NCONST; i++)
        arb_clear(ref[i]);

    flint_set_num_threads(1);
    TEST_FUNCTION_END(state);
}
