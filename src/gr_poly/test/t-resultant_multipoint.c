/*
    Copyright (C) 2026 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"

/* a prime p = m*2^16 + 1 of about `bits` bits, for which the evaluation can
   use a DFT or zero if none exist.*/

static ulong _fft_prime(flint_rand_t state, int bits)
{
    ulong hi, lo, m;

    FLINT_ASSERT(bits > 17);

    hi = (1UL << (bits - 16)) - 1;
    lo = 1UL << (bits - 17);

    for (m = hi - n_randint(state, hi - lo); m > lo; m--)
    {
        ulong p = (m << 16) + 1;
        if (n_is_prime(p))
            return p;
    }

    return 0;
}

/* replaces the leading coefficient in y by a product of nroots linear
   factors, so that specialising x at a root of it drops the degree in y */
static int
_split_leading_coeff(gr_poly_t f, flint_rand_t state, slong nroots,
                     gr_ctx_t ctx, gr_ctx_t cctx)
{
    int status = GR_SUCCESS;
    gr_poly_t lc, lin;
    slong i;

    if (f->length == 0)
        return status;

    gr_poly_init(lc, cctx);
    gr_poly_init(lin, cctx);

    status |= gr_poly_one(lc, cctx);
    status |= gr_poly_set_coeff_ui(lin, 1, 1, cctx);

    for (i = 0; i < nroots; i++)
    {
        status |= gr_poly_set_coeff_ui(lin, 0, n_randint(state, 64), cctx);
        status |= gr_poly_mul(lc, lc, lin, cctx);
    }

    status |= gr_poly_set_coeff_scalar(f, f->length - 1, lc, ctx);

    gr_poly_clear(lc, cctx);
    gr_poly_clear(lin, cctx);

    return status;
}

/* multiplies the leading coefficient in y by x - root, so that the
   specialisation at x = root drops in degree */
static int
_root_in_leading_coeff(gr_poly_t f, ulong root, gr_ctx_t ctx, gr_ctx_t cctx)
{
    int status = GR_SUCCESS;
    gr_poly_t lc, lin;

    if (f->length == 0)
        return status;

    gr_poly_init(lc, cctx);
    gr_poly_init(lin, cctx);

    status |= gr_poly_set_coeff_si(lin, 0, -(slong) root, cctx);
    status |= gr_poly_set_coeff_ui(lin, 1, 1, cctx);

    status |= gr_poly_mul(lc, (gr_poly_struct *)
        GR_ENTRY(f->coeffs, f->length - 1, ctx->sizeof_elem), lin, cctx);
    status |= gr_poly_set_coeff_scalar(f, f->length - 1, lc, ctx);

    gr_poly_clear(lc, cctx);
    gr_poly_clear(lin, cctx);

    return status;
}

/* random bivariate polynomial of length at most leny in y whose
   coefficients have length at most lenx in x */
static int
_gr_poly_randtest_bivariate(gr_poly_t f, flint_rand_t state, slong leny,
                            slong lenx, gr_ctx_t ctx, gr_ctx_t cctx)
{
    int status = GR_SUCCESS;
    gr_poly_t c;
    slong i;

    gr_poly_init(c, cctx);
    status |= gr_poly_zero(f, ctx);

    for (i = 0; i < leny; i++)
    {
        status |= gr_poly_randtest(c, state, lenx, cctx);
        status |= gr_poly_set_coeff_scalar(f, i, c, ctx);
    }

    gr_poly_clear(c, cctx);

    return status;
}

TEST_FUNCTION_START(gr_poly_resultant_multipoint, state)
{
    slong iter;

    /* Compare with the division-free Sylvester determinant over Z/pZ[x][y] */
    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx, ctx;
        gr_poly_t f, g;
        gr_ptr x, y;
        ulong p;
        int status = GR_SUCCESS;
        int s1;

        /* the algorithm needs a prime field with enough elements; small
           primes also exercise the case where a specialisation of g drops
           in degree */
        if (n_randint(state, 2))
            p = n_randprime(state, 4 + n_randint(state, 4), 1);
        else
            p = n_randtest_prime(state, 0);

        gr_ctx_init_nmod(cctx, p);
        gr_ctx_init_gr_poly(ctx, cctx);

        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);

        status |= _gr_poly_randtest_bivariate(f, state, 1 + n_randint(state, 8),
                                              1 + n_randint(state, 8), ctx, cctx);
        status |= _gr_poly_randtest_bivariate(g, state, 1 + n_randint(state, 8),
                                              1 + n_randint(state, 8), ctx, cctx);

        s1 = gr_poly_resultant_multipoint(x, f, g, ctx);
        status |= gr_poly_resultant_sylvester(y, f, g, ctx);

        // if (s1 != GR_SUCCESS)
        //     flint_printf("p = %lu fail\n", p);
        // else
        //     flint_printf("p = %lu good\n", p);

        // gr_poly_resultant_multipoint is allowed to fail, if the field is too small to find a suitable
        // primitive root for the geometric progression
        // gr_poly_resultant_sylvester should never FAIL
        // So we condition the result only when gr_poly_resultant_multipoint succeeds
        if ((status != GR_SUCCESS) || 
            ((s1 == GR_SUCCESS) && (gr_equal(x, y, ctx) != T_TRUE)))
        {
            flint_printf("FAIL (vs sylvester):\n");
            gr_ctx_println(ctx);
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n\n");
            flint_printf("x = "); gr_println(x, ctx);
            flint_printf("y = "); gr_println(y, ctx);
            fflush(stdout);
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    /* Same, over prime fields supporting a DFT of the required length, which
       is the other evaluation path. Coefficients in y are often zero here, so
       that the degenerate cases of that path are covered too. */
    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx, ctx;
        gr_poly_t f, g;
        gr_ptr x, y;
        ulong p;
        int status = GR_SUCCESS;
        int s1;

        p = _fft_prime(state, 24 + n_randint(state, 27));

        if (p == 0)
            continue;

        gr_ctx_init_nmod(cctx, p);
        gr_ctx_init_gr_poly(ctx, cctx);

        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);

        status |= _gr_poly_randtest_bivariate(f, state, 1 + n_randint(state, 7),
                                              1 + n_randint(state, 7), ctx, cctx);
        status |= _gr_poly_randtest_bivariate(g, state, 1 + n_randint(state, 7),
                                              1 + n_randint(state, 7), ctx, cctx);

        if (n_randint(state, 3) == 0)
        {
            /* a common factor in y makes the resultant zero */
            gr_poly_t h;
            gr_poly_init(h, ctx);
            status |= _gr_poly_randtest_bivariate(h, state, 2 + n_randint(state, 2),
                                                  1 + n_randint(state, 3), ctx, cctx);
            status |= gr_poly_mul(f, f, h, ctx);
            status |= gr_poly_mul(g, g, h, ctx);
            gr_poly_clear(h, ctx);
        }

        s1 = gr_poly_resultant_multipoint(x, f, g, ctx);
        status |= gr_poly_resultant_sylvester(y, f, g, ctx);

        // gr_poly_resultant_multipoint is allowed to fail, if the field is too small to find a suitable
        // primitive root for the geometric progression
        // gr_poly_resultant_sylvester should never FAIL
        // So we condition the result only when gr_poly_resultant_multipoint succeeds
        if ((status != GR_SUCCESS) || 
            ((s1 == GR_SUCCESS) && (gr_equal(x, y, ctx) != T_TRUE)))
        {
            flint_printf("FAIL (vs sylvester, DFT evaluation):\n");
            gr_ctx_println(ctx);
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n\n");
            flint_printf("x = "); gr_println(x, ctx);
            flint_printf("y = "); gr_println(y, ctx);
            fflush(stdout);
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    /* Enough points that the evaluation is split into several blocks, which
       both paths handle separately: the first prime below admits no DFT and
       so goes through the geometric progression, the second one does admit a
       DFT. Being short in y and long in x makes the points numerous while
       each costs little, so the blocking is reached cheaply; the reference is
       the subresultant algorithm, the Sylvester determinant being far too
       slow at this size. */
    for (iter = 0; iter < 2 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx, ctx;
        gr_poly_t f, g;
        gr_ptr x, y;
        ulong p;
        int status = GR_SUCCESS;
        int s1;

        /* large enough for the blocked evaluation to be split over the
           thread pool */
        flint_set_num_threads(1 + n_randint(state, 8));

        if (iter % 2)
        {
            p = _fft_prime(state, 50);
        }
        else
        {
            /* p = 3 mod 4 leaves p - 1 with a 2-valuation of one, too small
               for any transform, so the geometric progression is used */
            p = n_randprime(state, 50, 1);

            while (p % 4 != 3)
                p = n_randprime(state, 50, 1);
        }

        if (p == 0)
            continue;

        gr_ctx_init_nmod(cctx, p);
        gr_ctx_init_gr_poly(ctx, cctx);

        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);

        status |= _gr_poly_randtest_bivariate(f, state, 6, 2400, ctx, cctx);
        status |= _gr_poly_randtest_bivariate(g, state, 6, 2400, ctx, cctx);

        s1 = gr_poly_resultant_multipoint(x, f, g, ctx);
        status |= gr_poly_resultant_subresultant(y, f, g, ctx);

        // gr_poly_resultant_multipoint is allowed to fail, if the field is too small to find a suitable
        // primitive root for the geometric progression
        // gr_poly_resultant_sylvester should never FAIL
        // So we condition the result only when gr_poly_resultant_multipoint succeeds
        if ((status != GR_SUCCESS) || 
            ((s1 == GR_SUCCESS) && (gr_equal(x, y, ctx) != T_TRUE)))
        {
            flint_printf("FAIL (vs subresultant, blocked evaluation):\n");
            gr_ctx_println(ctx);
            flint_printf("status = %d, multipoint returned %d\n", status, s1);
            fflush(stdout);
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    /* The leading coefficients in y of f and of g vanish at many points, and
       at shared ones, so that the specialisations keep dropping in degree,
       one of them or both at once. The algorithm has to account for every
       such point rather than avoid it, so it must succeed whenever the field
       is large enough for the progression: p >= 2 npoints + 1. */
    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx, ctx;
        gr_poly_t f, g, c;
        gr_ptr x, y;
        ulong p;
        slong lenf, leng, blenf, bleng, npoints, i;
        const gr_poly_struct * fc, * gc;
        int status = GR_SUCCESS;
        int s1;

        p = n_randprime(state, 8 + n_randint(state, 8), 1);

        gr_ctx_init_nmod(cctx, p);
        gr_ctx_init_gr_poly(ctx, cctx);

        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);
        gr_poly_init(c, cctx);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);

        status |= _gr_poly_randtest_bivariate(f, state, 2 + n_randint(state, 5),
                                              1 + n_randint(state, 5), ctx, cctx);
        status |= _gr_poly_randtest_bivariate(g, state, 2 + n_randint(state, 5),
                                              1 + n_randint(state, 5), ctx, cctx);

        status |= _split_leading_coeff(f, state, 1 + n_randint(state, 3), ctx, cctx);

        switch (n_randint(state, 3))
        {
            case 0:
                /* only f can drop in degree */
                if (g->length != 0)
                {
                    status |= gr_poly_set_coeff_ui(c, 0, 1 + n_randint(state, p - 1), cctx);
                    status |= gr_poly_set_coeff_scalar(g, g->length - 1, c, ctx);
                }
                break;
            case 1:
                /* both can, and x = 1 is a root of each: the progression
                   starts there, so both drop at once at that point */
                status |= _split_leading_coeff(g, state, 1 + n_randint(state, 3), ctx, cctx);
                status |= _root_in_leading_coeff(f, 1, ctx, cctx);
                status |= _root_in_leading_coeff(g, 1, ctx, cctx);
                break;
            default:
                break;
        }

        s1 = gr_poly_resultant_multipoint(x, f, g, ctx);
        status |= gr_poly_resultant_sylvester(y, f, g, ctx);

        /* the bound on deg_x of the resultant that the algorithm uses */
        if (f->length >= g->length)
        {
            fc = (const gr_poly_struct *) f->coeffs; lenf = f->length;
            gc = (const gr_poly_struct *) g->coeffs; leng = g->length;
        }
        else
        {
            fc = (const gr_poly_struct *) g->coeffs; lenf = g->length;
            gc = (const gr_poly_struct *) f->coeffs; leng = f->length;
        }

        blenf = 0;
        for (i = 0; i < lenf; i++)
            blenf = FLINT_MAX(blenf, fc[i].length);

        bleng = 0;
        for (i = 0; i < leng; i++)
            bleng = FLINT_MAX(bleng, gc[i].length);

        npoints = (leng - 1) * (blenf - 1) + (lenf - 1) * (bleng - 1) + 1;

        if ((status != GR_SUCCESS) ||
            ((s1 == GR_SUCCESS) && (gr_equal(x, y, ctx) != T_TRUE)) ||
            ((s1 != GR_SUCCESS) && leng >= 2 && (p - 1) / 2 >= (ulong) npoints))
        {
            flint_printf("FAIL (vs sylvester, degree drop in y):\n");
            gr_ctx_println(ctx);
            flint_printf("npoints = %wd, status = %d\n", npoints, s1);
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n\n");
            flint_printf("x = "); gr_println(x, ctx);
            flint_printf("y = "); gr_println(y, ctx);
            fflush(stdout);
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_poly_clear(c, cctx);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    /* Check res(f h, g) == res(f, g) res(h, g) over a large prime field */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx, ctx;
        gr_poly_t f, fh, g, h;
        gr_ptr x, y, z, yz;
        int status = GR_SUCCESS;

        gr_ctx_init_nmod(cctx, n_randprime(state, FLINT_BITS - 2, 1));
        gr_ctx_init_gr_poly(ctx, cctx);

        gr_poly_init(f, ctx);
        gr_poly_init(fh, ctx);
        gr_poly_init(g, ctx);
        gr_poly_init(h, ctx);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);
        z = gr_heap_init(ctx);
        yz = gr_heap_init(ctx);

        status |= _gr_poly_randtest_bivariate(f, state, 1 + n_randint(state, 6),
                                              1 + n_randint(state, 6), ctx, cctx);
        status |= _gr_poly_randtest_bivariate(g, state, 1 + n_randint(state, 6),
                                              1 + n_randint(state, 6), ctx, cctx);
        status |= _gr_poly_randtest_bivariate(h, state, 1 + n_randint(state, 6),
                                              1 + n_randint(state, 6), ctx, cctx);

        status |= gr_poly_mul(fh, f, h, ctx);

        status |= gr_poly_resultant_multipoint(x, fh, g, ctx);
        status |= gr_poly_resultant_multipoint(y, f, g, ctx);
        status |= gr_poly_resultant_multipoint(z, h, g, ctx);
        status |= gr_mul(yz, y, z, ctx);

        if (status == GR_SUCCESS && gr_equal(x, yz, ctx) == T_FALSE)
        {
            flint_printf("FAIL (res(f h, g) == res(f, g) res(h, g)):\n");
            gr_ctx_println(ctx);
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n\n");
            flint_printf("h = "); gr_poly_print(h, ctx); flint_printf("\n\n");
            flint_printf("x = "); gr_println(x, ctx);
            flint_printf("yz = "); gr_println(yz, ctx);
            fflush(stdout);
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(fh, ctx);
        gr_poly_clear(g, ctx);
        gr_poly_clear(h, ctx);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_heap_clear(z, ctx);
        gr_heap_clear(yz, ctx);
        gr_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    /* Check res(f, g) == (-1)^(deg f deg g) res(g, f), with larger degrees
       so that the generic resultant actually dispatches to this algorithm */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx, ctx;
        gr_poly_t f, g;
        gr_ptr x, y;
        int status = GR_SUCCESS;

        /* large enough for the evaluation and the univariate resultants to
           be split over the thread pool */
        flint_set_num_threads(1 + n_randint(state, 8));

        gr_ctx_init_nmod(cctx, n_randprime(state, FLINT_BITS - 2, 1));
        gr_ctx_init_gr_poly(ctx, cctx);

        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);
        x = gr_heap_init(ctx);
        y = gr_heap_init(ctx);

        status |= _gr_poly_randtest_bivariate(f, state, 1 + n_randint(state, 20),
                                              1 + n_randint(state, 20), ctx, cctx);
        status |= _gr_poly_randtest_bivariate(g, state, 1 + n_randint(state, 20),
                                              1 + n_randint(state, 20), ctx, cctx);

        status |= gr_poly_resultant(x, f, g, ctx);
        status |= gr_poly_resultant(y, g, f, ctx);

        if (((f->length - 1) * (g->length - 1)) % 2)
            status |= gr_neg(y, y, ctx);

        if (status != GR_SUCCESS)
        {
            flint_printf("FAIL: did not succeed over Z/pZ[x][y]\n\n");
            gr_ctx_println(ctx);
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        if (gr_equal(x, y, ctx) == T_FALSE)
        {
            flint_printf("FAIL (res(f, g) == (-1)^(deg f deg g) res(g, f)):\n");
            gr_ctx_println(ctx);
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n\n");
            flint_printf("x = "); gr_println(x, ctx);
            flint_printf("y = "); gr_println(y, ctx);
            fflush(stdout);
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_heap_clear(x, ctx);
        gr_heap_clear(y, ctx);
        gr_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
