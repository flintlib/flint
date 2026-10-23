/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fft_small.h"

/* z[j] = sum_i a[i] v[i + j] for 0 <= j < zn, v of length an + zn - 1,
   as the range [an - 1, an - 1 + zn) of rev(a)*v */
static void _mul_tr_ref(ulong* z, ulong zn, const ulong* a, ulong an,
                        const ulong* v, nmod_t mod)
{
    ulong vn = an + zn - 1;
    ulong* ar = _nmod_vec_init(an);
    ulong* full = _nmod_vec_init(an + vn - 1);

    _nmod_poly_reverse(ar, a, an, an);
    _nmod_poly_mul_KS(full, v, vn, ar, an, mod);
    _nmod_vec_set(z, full + an - 1, zn);

    _nmod_vec_clear(ar);
    _nmod_vec_clear(full);
}

static void check(const ulong* z, const ulong* zref, ulong an, ulong zn,
                  nmod_t mod, const char* where)
{
    for (ulong i = 0; i < zn; i++)
    {
        if (z[i] != zref[i])
        {
            flint_printf("FAIL: %s\nmul_tr error at index %wu\n", where, i);
            flint_printf("an=%wu, zn=%wu\nmod: %wu\n", an, zn, mod.n);
            flint_printf("got %wu, want %wu\n", z[i], zref[i]);
            fflush(stdout);
            flint_abort();
        }
    }
}

static void test_one(ulong an, ulong zn, nmod_t mod, mpn_ctx_t R,
                     flint_rand_t state, const char* where)
{
    ulong vn = an + zn - 1;
    ulong* a = _nmod_vec_init(an);
    ulong* v = _nmod_vec_init(vn);
    ulong* z = _nmod_vec_init(zn);
    ulong* zref = _nmod_vec_init(zn);

    _nmod_vec_randtest(a, state, an, mod);
    _nmod_vec_randtest(v, state, vn, mod);

    _nmod_poly_mul_tr_mpn_ctx(z, zn, a, an, v, vn, mod, R);
    _mul_tr_ref(zref, zn, a, an, v, mod);
    check(z, zref, an, zn, mod, where);

    /* transposition: <u, mul_tr(a, v)> = <a*u, v> for u of length zn */
    {
        ulong* u = _nmod_vec_init(zn);
        ulong* au = _nmod_vec_init(vn);

        _nmod_vec_randtest(u, state, zn, mod);
        if (an >= zn)
            _nmod_poly_mul_KS(au, a, an, u, zn, mod);
        else
            _nmod_poly_mul_KS(au, u, zn, a, an, mod);

        if (_nmod_vec_dot(u, z, zn, mod, _nmod_vec_dot_params(zn, mod))
                != _nmod_vec_dot(au, v, vn, mod, _nmod_vec_dot_params(vn, mod)))
        {
            flint_printf("FAIL: %s\nmul_tr is not the transpose of mul\n", where);
            flint_printf("an=%wu, zn=%wu\nmod: %wu\n", an, zn, mod.n);
            fflush(stdout);
            flint_abort();
        }

        _nmod_vec_clear(u);
        _nmod_vec_clear(au);
    }

    /* the same output from a shorter v, read as padded with zeros */
    if (vn > 1)
    {
        ulong cut = 1 + n_randint(state, vn - 1);
        _nmod_vec_zero(v + vn - cut, cut);
        _nmod_poly_mul_tr_mpn_ctx(z, zn, a, an, v, vn - cut, mod, R);
        _mul_tr_ref(zref, zn, a, an, v, mod);
        check(z, zref, an, zn, mod, where);
    }

    _nmod_vec_clear(a);
    _nmod_vec_clear(v);
    _nmod_vec_clear(z);
    _nmod_vec_clear(zref);
}

TEST_FUNCTION_START(_nmod_poly_mul_tr_mpn_ctx, state)
{
    mpn_ctx_t R;
    nmod_t mod;

    mpn_ctx_init(R, UWORD(0x0003f00000000001));

    /* small everything, exhaustively in the shape */
    {
        ulong nn[4] = {UWORD(2), UWORD(3), UWORD(4611686018427387847),
                       UWORD(18446744073709551557)};

        for (int k = 0; k < 4; k++)
        {
            nmod_init(&mod, nn[k]);

            for (ulong an = 1; an <= 20; an++)
                for (ulong zn = 1; zn <= 21; zn++)
                    test_one(an, zn, mod, R, state, "small");
        }
    }

    /* random shapes and moduli, across the block and depth boundaries */
    {
        for (slong reps = 0; reps < 60 * flint_test_multiplier(); reps++)
        {
            flint_bitcnt_t bits = 2 + n_randint(state, FLINT_BITS - 1);
            ulong n = n_randtest_bits(state, bits);
            ulong an, zn;

            n = FLINT_MAX(n, UWORD(2));
            nmod_init(&mod, n);

            an = 1 + n_randint(state, 1500);
            zn = 1 + n_randint(state, 1500);

            test_one(an, zn, mod, R, state, "random");
        }
    }

    /* the shapes the transposition is degenerate at */
    {
        nmod_init(&mod, UWORD(0xffffffff00000001));

        test_one(1, 1, mod, R, state, "1 x 1");
        test_one(1, 1000, mod, R, state, "an = 1");
        test_one(1000, 1, mod, R, state, "zn = 1");
        test_one(BLK_SZ, BLK_SZ, mod, R, state, "block aligned");
        test_one(BLK_SZ + 1, BLK_SZ, mod, R, state, "block aligned + 1");
    }

    /* the single-prime branch: mod.n is one of the context's fft primes */
    {
        nmod_init(&mod, UWORD(0x0003f00000000001));

        test_one(300, 401, mod, R, state, "fft prime");
        test_one(2000, 2001, mod, R, state, "fft prime");
    }

    /* the direct fft branch: mod.n prime, mod.n - 1 with high 2-valuation,
       and operands long enough for _nmod_poly_should_directly_fft */
    {
        nmod_init(&mod, UWORD(7340033));  /* 7*2^20 + 1 */

        for (slong reps = 0; reps < 3; reps++)
        {
            flint_set_num_threads(1 + reps);
            test_one(1700 + 31*reps, 1701 + 22*reps, mod, R, state, "direct fft");
        }

        /* one side shorter than the direct-fft gate */
        test_one(1700, 50, mod, R, state, "direct fft, short zn");
        test_one(30, 2500, mod, R, state, "direct fft, short an");

        flint_set_num_threads(1);
    }

    /* threaded */
    {
        nmod_init(&mod, UWORD(18446744073709551557));

        for (slong reps = 1; reps <= 3; reps++)
        {
            flint_set_num_threads(reps);
            test_one(2000, 3001, mod, R, state, "threaded");
        }

        flint_set_num_threads(1);
    }

    mpn_ctx_clear(R);

    TEST_FUNCTION_END(state);
}
