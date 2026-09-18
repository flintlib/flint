/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "gmpcompat.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "gr.h"
#include "mpn_mod.h"
#include "profiler.h"

/* Inputs per measurement, so that a timing is not one lucky element. */
#define NX 16

/* mpz_jacobi is declared pure, so a discarded result gets optimized away. */
volatile slong sink;

/*
    A modular square root written in plain GMP: the same algorithm as
    mpn_mod_sqrt, with nothing but mpz arithmetic, so that the comparison
    measures the mpn_mod layer rather than the algorithm.
*/
static int
gmp_sqrtmod(mpz_t rop, const mpz_t a, const mpz_t p)
{
    mpz_t q, e, b, c, r, t, u;
    slong s, i, j, m;
    ulong k;

    if (flint_mpz_cmp_ui(a, 1) <= 0)
    {
        mpz_set(rop, a);
        return 1;
    }

    if (mpz_jacobi(a, p) == -1)
        return 0;

    mpz_init(q); mpz_init(e); mpz_init(b);
    mpz_init(c); mpz_init(r); mpz_init(t); mpz_init(u);

    flint_mpz_sub_ui(q, p, 1);
    s = mpz_scan1(q, 0);
    mpz_tdiv_q_2exp(q, q, s);

    if (s == 1)
    {
        flint_mpz_add_ui(e, p, 1);
        mpz_tdiv_q_2exp(e, e, 2);
        mpz_powm(rop, a, e, p);
    }
    else if (s == 2)
    {
        flint_mpz_sub_ui(e, p, 5);
        mpz_tdiv_q_2exp(e, e, 3);

        mpz_mul_2exp(t, a, 1); mpz_mod(t, t, p);        /* t = 2a */
        mpz_powm(b, t, e, p);
        mpz_mul(u, b, b); mpz_mod(u, u, p);
        mpz_mul(u, u, t); mpz_mod(u, u, p);             /* u = 2 a b^2 */
        flint_mpz_sub_ui(u, u, 1);
        mpz_mul(rop, b, a); mpz_mod(rop, rop, p);
        mpz_mul(rop, rop, u); mpz_mod(rop, rop, p);
    }
    else
    {
        flint_mpz_sub_ui(e, q, 1);
        mpz_tdiv_q_2exp(e, e, 1);

        mpz_powm(u, a, e, p);                           /* u = a^((q-1)/2) */
        mpz_mul(r, u, a); mpz_mod(r, r, p);
        mpz_mul(t, r, u); mpz_mod(t, t, p);             /* t = a^q */

        for (k = 3; ; k += 2)
        {
            flint_mpz_set_ui(e, k);

            if (mpz_jacobi(e, p) == -1)
                break;
        }

        mpz_powm(c, e, q, p);

        m = s;

        while (flint_mpz_cmp_ui(t, 1) != 0)
        {
            mpz_set(u, t);

            for (i = 1; i < m; i++)
            {
                mpz_mul(u, u, u); mpz_mod(u, u, p);

                if (flint_mpz_cmp_ui(u, 1) == 0)
                    break;
            }

            mpz_set(b, c);

            for (j = 0; j < m - i - 1; j++)
            {
                mpz_mul(b, b, b); mpz_mod(b, b, p);
            }

            m = i;
            mpz_mul(c, b, b); mpz_mod(c, c, p);
            mpz_mul(t, t, c); mpz_mod(t, t, p);
            mpz_mul(r, r, b); mpz_mod(r, r, p);
        }

        mpz_set(rop, r);
    }

    mpz_clear(q); mpz_clear(e); mpz_clear(b);
    mpz_clear(c); mpz_clear(r); mpz_clear(t); mpz_clear(u);

    return 1;
}

/* A prime of the given size with p - 1 divisible by exactly 2^s. */
static void
random_prime_val2(fmpz_t p, flint_rand_t state, flint_bitcnt_t bits, slong s)
{
    fmpz_t q;

    fmpz_init(q);

    while (1)
    {
        fmpz_randbits(q, state, bits - s);
        fmpz_abs(q, q);
        fmpz_setbit(q, 0);
        fmpz_setbit(q, bits - s - 1);
        fmpz_mul_2exp(p, q, s);
        fmpz_add_ui(p, p, 1);

        if (fmpz_bits(p) == bits && fmpz_is_probabprime(p))
            break;
    }

    fmpz_clear(q);
}

int main(void)
{
    flint_rand_t state;
    fmpz_t p, y;
    fmpz * x;
    gr_ctx_t ctx, ctx2;
    mpz_t pz, rz, xz[NX];
    nn_ptr xm, rm;
    slong bits_tab[] = { 128, 192, 256, 384, 512, 1024 };
    slong s_tab[] = { 1, 2, 3, 16, 64 };
    slong bi, si, i, nlimbs;
    flint_bitcnt_t bits;
    slong s;
    double t1, t2, t3, __;

    flint_rand_init(state);
    fmpz_init(p);
    fmpz_init(y);
    mpz_init(pz);
    mpz_init(rz);

    for (i = 0; i < NX; i++)
        mpz_init(xz[i]);

    x = _fmpz_vec_init(NX);

    flint_printf("Square roots modulo a prime p with p - 1 = q 2^s, q odd.\n");
    flint_printf("Times are microseconds per call, averaged over %wd inputs.\n\n", (slong) NX);

    flint_printf("                     is_square                      sqrt\n");
    flint_printf(" bits    s    mpn_mod fmpz_mod     gmp     mpn_mod fmpz_mod     gmp\n");

    for (bi = 0; bi < 6; bi++)
    {
        bits = bits_tab[bi];

        for (si = 0; si < 5; si++)
        {
            s = s_tab[si];

            if (s + 8 >= (slong) bits)
                continue;

            random_prime_val2(p, state, bits, s);

            GR_MUST_SUCCEED(gr_ctx_init_mpn_mod(ctx, p));
            gr_ctx_init_fmpz_mod(ctx2, p);
            GR_MUST_SUCCEED(gr_ctx_set_is_field(ctx, T_TRUE));
            GR_MUST_SUCCEED(gr_ctx_set_is_field(ctx2, T_TRUE));

            nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
            xm = flint_malloc((NX + 1) * nlimbs * sizeof(ulong));
            rm = xm + NX * nlimbs;

            fmpz_get_mpz(pz, p);

            /* squares, so that every implementation takes the long route */
            for (i = 0; i < NX; i++)
            {
                fmpz_randm(y, state, p);
                fmpz_mul(x + i, y, y);
                fmpz_mod(x + i, x + i, p);
                GR_MUST_SUCCEED(mpn_mod_set_fmpz(xm + i * nlimbs, x + i, ctx));
                fmpz_get_mpz(xz[i], x + i);
            }

            /* the three had better agree, up to the sign of the root */
            for (i = 0; i < NX; i++)
            {
                fmpz_t a, b;

                fmpz_init(a);
                fmpz_init(b);

                if (mpn_mod_sqrt(rm, xm + i * nlimbs, ctx) != GR_SUCCESS)
                    flint_throw(FLINT_ERROR, "mpn_mod_sqrt failed on a square\n");

                GR_MUST_SUCCEED(mpn_mod_get_fmpz(a, rm, ctx));
                fmpz_mul(b, a, a);
                fmpz_mod(b, b, p);

                if (!fmpz_equal(b, x + i))
                    flint_throw(FLINT_ERROR, "mpn_mod_sqrt returned a wrong root\n");

                if (!gmp_sqrtmod(rz, xz[i], pz))
                    flint_throw(FLINT_ERROR, "gmp_sqrtmod failed on a square\n");

                mpz_mul(rz, rz, rz);
                mpz_mod(rz, rz, pz);
                fmpz_set_mpz(b, rz);

                if (!fmpz_equal(b, x + i))
                    flint_throw(FLINT_ERROR, "gmp_sqrtmod returned a wrong root\n");

                fmpz_clear(a);
                fmpz_clear(b);
            }

            TIMEIT_START
            for (i = 0; i < NX; i++)
                sink += mpn_mod_is_square(xm + i * nlimbs, ctx);
            TIMEIT_STOP_VALUES(t1, __);

            TIMEIT_START
            for (i = 0; i < NX; i++)
                sink += gr_is_square(x + i, ctx2);
            TIMEIT_STOP_VALUES(t2, __);

            TIMEIT_START
            for (i = 0; i < NX; i++)
                sink += mpz_jacobi(xz[i], pz);
            TIMEIT_STOP_VALUES(t3, __);

            flint_printf("%5wd %4wd   %8.3f %8.3f %8.3f",
                (slong) bits, s, 1e6 * t1 / NX, 1e6 * t2 / NX, 1e6 * t3 / NX);

            TIMEIT_START
            for (i = 0; i < NX; i++)
                sink += mpn_mod_sqrt(rm, xm + i * nlimbs, ctx);
            TIMEIT_STOP_VALUES(t1, __);

            TIMEIT_START
            for (i = 0; i < NX; i++)
                sink += fmpz_sqrtmod(y, x + i, p);
            TIMEIT_STOP_VALUES(t2, __);

            TIMEIT_START
            for (i = 0; i < NX; i++)
                sink += gmp_sqrtmod(rz, xz[i], pz);
            TIMEIT_STOP_VALUES(t3, __);

            flint_printf("    %8.3f %8.3f %8.3f\n",
                1e6 * t1 / NX, 1e6 * t2 / NX, 1e6 * t3 / NX);

            (void) __;

            flint_free(xm);
            gr_ctx_clear(ctx);
            gr_ctx_clear(ctx2);
        }
    }

    _fmpz_vec_clear(x, NX);
    mpz_clear(pz);
    mpz_clear(rz);

    for (i = 0; i < NX; i++)
        mpz_clear(xz[i]);
    fmpz_clear(p);
    fmpz_clear(y);
    flint_rand_clear(state);

    return 0;
}
