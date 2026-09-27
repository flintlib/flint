/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "ulong_extras.h"

#define ARRAY_LEN(x) ((slong) (sizeof(x) / sizeof((x)[0])))

/* compare flint_mpn_powmod_preinvn with fmpz_powm, e padded to en limbs */
static void
_check(const fmpz_t a, const fmpz_t e, const fmpz_t d, mp_size_t pad)
{
    fmpz_t r1, r2;
    nn_ptr dn, dinv, as, rs, es;
    mp_size_t n = fmpz_size(d), en = FLINT_MAX(fmpz_size(e), 1) + pad;
    flint_bitcnt_t norm;

    fmpz_init(r1);
    fmpz_init(r2);

    dn = flint_malloc((4 * n + en) * sizeof(ulong));
    dinv = dn + n;
    as = dinv + n;
    rs = as + n;
    es = rs + n;

    fmpz_get_ui_array(dn, n, d);
    fmpz_get_ui_array(as, n, a);
    fmpz_get_ui_array(es, en, e);

    norm = flint_clz(dn[n - 1]);

    if (norm)
    {
        mpn_lshift(dn, dn, n, norm);
        mpn_lshift(as, as, n, norm);
    }

    flint_mpn_preinvn(dinv, dn, n);
    flint_mpn_powmod_preinvn(rs, as, es, en, n, dn, dinv, norm);

    /* the result must carry the shift exactly, and be reduced */
    if (norm && (rs[0] & ((UWORD(1) << norm) - 1)) != 0)
    {
        flint_printf("FAIL: result not a shifted value\n");
        goto fail;
    }

    if (norm)
        mpn_rshift(rs, rs, n, norm);

    fmpz_set_ui_array(r2, rs, n);
    fmpz_powm(r1, a, e, d);

    if (!fmpz_equal(r1, r2) || fmpz_cmp(r2, d) >= 0)
    {
        flint_printf("FAIL: wrong result\n");
        goto fail;
    }

    flint_free(dn);
    fmpz_clear(r1);
    fmpz_clear(r2);
    return;

fail:
    flint_printf("n = %wd, norm = %wu, en = %wd, ebits = %wu\n",
            (slong) n, (ulong) norm, (slong) en, (ulong) fmpz_bits(e));
    flint_printf("a = "); fmpz_print(a); flint_printf("\n");
    flint_printf("e = "); fmpz_print(e); flint_printf("\n");
    flint_printf("d = "); fmpz_print(d); flint_printf("\n");
    flint_printf("expected = "); fmpz_print(r1); flint_printf("\n");
    flint_printf("got      = "); fmpz_print(r2); flint_printf("\n");
    fflush(stdout);
    flint_abort();
}

/* an n-limb modulus of bit length n B - norm, at least 2, of various shapes */
static void
_rand_modulus(fmpz_t d, flint_rand_t state, mp_size_t n, flint_bitcnt_t norm)
{
    flint_bitcnt_t bits = n * FLINT_BITS - norm;

    switch (n_randint(state, 8))
    {
        case 0:     /* all ones */
            fmpz_one(d);
            fmpz_mul_2exp(d, d, bits);
            fmpz_sub_ui(d, d, 1);
            break;
        case 1:     /* a power of two, just past it, or just below the next */
            fmpz_one(d);
            fmpz_mul_2exp(d, d, bits - 1);
            fmpz_add_ui(d, d, n_randint(state, 3));
            break;
        default:
            fmpz_randbits(d, state, bits);
            fmpz_abs(d, d);
            fmpz_setbit(d, bits - 1);
            if (n_randint(state, 2))
                fmpz_setbit(d, 0);
    }

    if (fmpz_cmp_ui(d, 2) < 0)
        fmpz_set_ui(d, 2);
}

/* a base reduced modulo d, favouring the small ones that have a path */
static void
_rand_base(fmpz_t a, flint_rand_t state, const fmpz_t d)
{
    switch (n_randint(state, 12))
    {
        case 0: fmpz_zero(a); break;
        case 1: fmpz_one(a); break;
        case 2: fmpz_set_ui(a, 2); break;
        case 3: fmpz_set_ui(a, 3); break;
        case 4: fmpz_set_ui(a, n_randint(state, 1000)); break;
        case 5: fmpz_set_ui(a, n_randint(state, UWORD(1) << (FLINT_BITS / 3))); break;
        case 6: fmpz_set_ui(a, n_randlimb(state)); break;
        case 7: fmpz_sub_ui(a, d, 1); break;
        case 8: fmpz_randbits(a, state, FLINT_BITS + 1 + n_randint(state, FLINT_BITS)); break;
        default: fmpz_randm(a, state, d);
    }

    fmpz_abs(a, a);
    fmpz_mod(a, a, d);
}

/* an exponent of at most maxbits bits, with the special shapes included */
static void
_rand_exponent(fmpz_t e, flint_rand_t state, flint_bitcnt_t maxbits)
{
    flint_bitcnt_t k = n_randint(state, maxbits + 1);

    switch (n_randint(state, 8))
    {
        case 0:
            fmpz_set_ui(e, n_randint(state, 4));
            break;
        case 1:     /* a single set bit */
            fmpz_one(e);
            fmpz_mul_2exp(e, e, k);
            break;
        case 2:     /* all set bits */
            fmpz_one(e);
            fmpz_mul_2exp(e, e, k);
            fmpz_sub_ui(e, e, 1);
            break;
        case 3:
            fmpz_randbits(e, state, n_randint(state, FLINT_BITS + 1));
            break;
        default:
            fmpz_randbits(e, state, k);
    }

    fmpz_abs(e, e);
}

TEST_FUNCTION_START(flint_mpn_powmod_preinvn, state)
{
    fmpz_t a, d, e;
    slong iter;

    fmpz_init(a);
    fmpz_init(d);
    fmpz_init(e);

    /*
        A sweep through every size and exponent length at which the choice
        of algorithm changes, for bases of each kind and moduli with and
        without a shift, so that every path and every boundary is taken.
    */
    {
        static const flint_bitcnt_t ebits_tab[] = { 1, 2, 7, 8, 15, 16, 17,
            31, 32, 33, 63, 64, 65, 127, 128, 129, 191, 192, 255, 256, 257,
            511, 512, 513, 1023, 1024, 1025, 2047, 2048, 2049 };
        const ulong base_tab[] = { 0, 1, 2, 3, 5,
            (UWORD(1) << (FLINT_BITS / 3)) - 1,     /* b^3 fits, just */
            (UWORD(1) << (FLINT_BITS / 3)) + 1,     /* b^3 does not */
            UWORD_MAX };
        mp_size_t n;
        slong ni, bi, ei;

        /* the tables stop at 22 limbs, and only go past 257 bits up to 13 */
        for (n = 1; n <= 23; n++)
        for (ni = 0; ni < 2; ni++)
        for (bi = 0; bi <= ARRAY_LEN(base_tab); bi++)
        {
            /* no shift, or a random nonzero one leaving at least 2 bits */
            flint_bitcnt_t norm = (ni == 0) ? 0 :
                1 + n_randint(state, n == 1 ? FLINT_BITS - 2 : FLINT_BITS - 1);

            fmpz_randbits(d, state, n * FLINT_BITS - norm);
            fmpz_abs(d, d);
            fmpz_setbit(d, n * FLINT_BITS - norm - 1);
            fmpz_setbit(d, 0);

            if (bi == ARRAY_LEN(base_tab))
                fmpz_randm(a, state, d);            /* a full-size base */
            else
            {
                fmpz_set_ui(a, base_tab[bi]);
                fmpz_mod(a, a, d);
            }

            for (ei = 0; ei < ARRAY_LEN(ebits_tab); ei++)
            {
                if (n > 13 && ebits_tab[ei] > 257)
                    break;

                fmpz_randbits(e, state, ebits_tab[ei]);
                fmpz_abs(e, e);
                fmpz_setbit(e, ebits_tab[ei] - 1);
                _check(a, e, d, 0);
            }
        }
    }

    /* random sizes up to well past the tables, with all shapes mixed */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        mp_size_t n;
        flint_bitcnt_t norm, maxbits;
        ulong r = n_randint(state, 100);

        if (r < 75)
            n = 1 + n_randint(state, 24);
        else if (r < 99)
            n = 25 + n_randint(state, 24);
        else
            n = 49 + n_randint(state, 112);

        norm = n_randint(state, 4) == 0 ? 0 :
                    n_randint(state, n == 1 ? FLINT_BITS - 1 : FLINT_BITS);

        _rand_modulus(d, state, n, norm);
        _rand_base(a, state, d);

        /* now and then long enough for the widest window */
        if (n <= 10 && n_randint(state, 50) == 0)
            maxbits = 28000 + n_randint(state, 4000);
        else if (n <= 48)
            maxbits = FLINT_MIN(4 * n * FLINT_BITS, 2048);
        else
            maxbits = 512;

        _rand_exponent(e, state, maxbits);

        _check(a, e, d, n_randint(state, 4) == 0 ? n_randint(state, 3) : 0);
    }

    fmpz_clear(a);
    fmpz_clear(d);
    fmpz_clear(e);

    TEST_FUNCTION_END(state);
}
