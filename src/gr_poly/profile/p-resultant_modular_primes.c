/*
    Copyright (C) 2026 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* The multimodular resultant over Z reduces at word-size primes and runs the
   multipoint algorithm at each of them. Which primes are used is free, and the
   choice trades two effects against each other:

     - a prime p for which p - 1 has a large power of two as a factor, and
       which is small enough for the fft_small transforms, lets the multipoint
       algorithm evaluate with DFTs instead of the Bluestein products of the
       geometric method, which is several times faster;

     - such a prime has at most 50 bits, against the 62 bits of the primes
       chosen without that constraint, so 62/50 = 1.24 times as many of them
       are needed to reach the same bound.

   This program measures the product of the two, at a fixed input, by running
   the same multimodular scheme twice over the same coefficient bound and the
   same CRT, changing only the list of primes. */

#include <stdlib.h>
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "nmod.h"
#include "ulong_extras.h"
#include "gr.h"
#include "gr_poly.h"
#include "profiler.h"
#if FLINT_HAVE_FFT_SMALL
# include "fft_small.h"
#endif

/* fft_small works with doubles and accepts no modulus above 50 bits */
#define FFT_PRIME_BITS 50

/* as in resultant_modular.c */
#define MODULAR_PRIME_BITS (FLINT_BITS - 2)

/* ||res||_oo <= (sum_i ||A_i||_1^2)^(n/2) (sum_j ||B_j||_1^2)^(m/2), the bound
   used by the multimodular code */
static void
_bound(fmpz_t bound, const fmpz_poly_struct * A, slong lenA,
       const fmpz_poly_struct * B, slong lenB)
{
    fmpz_t sa, sb, t, u;
    slong i, k;

    fmpz_init(sa);
    fmpz_init(sb);
    fmpz_init(t);
    fmpz_init(u);

    for (i = 0; i < lenA; i++)
    {
        fmpz_zero(t);
        for (k = 0; k < A[i].length; k++)
        {
            fmpz_abs(u, A[i].coeffs + k);
            fmpz_add(t, t, u);
        }
        fmpz_addmul(sa, t, t);
    }

    for (i = 0; i < lenB; i++)
    {
        fmpz_zero(t);
        for (k = 0; k < B[i].length; k++)
        {
            fmpz_abs(u, B[i].coeffs + k);
            fmpz_add(t, t, u);
        }
        fmpz_addmul(sb, t, t);
    }

    fmpz_pow_ui(sa, sa, lenB - 1);
    fmpz_pow_ui(sb, sb, lenA - 1);
    fmpz_mul(bound, sa, sb);
    fmpz_sqrt(bound, bound);
    fmpz_add_ui(bound, bound, 1);

    fmpz_clear(sa);
    fmpz_clear(sb);
    fmpz_clear(t);
    fmpz_clear(u);
}

/* deg_x(res) <= deg_y(B) deg_x(A) + deg_y(A) deg_x(B), which is also the
   number of points the multipoint algorithm evaluates at */
static slong
_outlen(const fmpz_poly_struct * A, slong lenA,
        const fmpz_poly_struct * B, slong lenB)
{
    slong i, dxA = 0, dxB = 0;

    for (i = 0; i < lenA; i++)
        dxA = FLINT_MAX(dxA, A[i].length);
    for (i = 0; i < lenB; i++)
        dxB = FLINT_MAX(dxB, B[i].length);

    return (lenB - 1) * (dxA - 1) + (lenA - 1) * (dxB - 1) + 1;
}

/* the primes just above 2^62, as the multimodular code picks them */
static void
_generic_primes(nn_ptr primes, slong num, const fmpz_t l)
{
    ulong p = UWORD(1) << MODULAR_PRIME_BITS;
    slong n = 0;

    while (n < num)
    {
        p = n_nextprime(p, 0);

        if (fmpz_fdiv_ui(l, p) == 0)
            continue;

        primes[n++] = p;
    }
}

/* the largest primes p = m 2^k + 1 below 2^50, so that the multipoint
   algorithm can evaluate with a transform of depth up to k. Returns the number
   found, which is less than num when the exponent leaves too few candidates. */
static slong
_fft_primes(nn_ptr primes, slong num, flint_bitcnt_t k, const fmpz_t l)
{
    ulong m;
    slong n = 0;

    /* stopping at 2^(49 - k) keeps every prime above 2^49, so that num of
       them are enough for a bound of 49 num bits */
    for (m = (UWORD(1) << (FFT_PRIME_BITS - k)) - 1;
         m > (UWORD(1) << (FFT_PRIME_BITS - 1 - k)) && n < num; m--)
    {
        ulong p = (m << k) + 1;

        if (!n_is_prime(p))
            continue;
        if (!fft_small_mulmod_satisfies_bounds(p))
            continue;
        if (fmpz_fdiv_ui(l, p) == 0)
            continue;

        primes[n++] = p;
    }

    return n;
}

/* res_y(A, B) mod p, written to (out, outlen) with zero padding */
static void
_image(nn_ptr out, slong outlen, const fmpz_poly_struct * A, slong lenA,
       const fmpz_poly_struct * B, slong lenB, ulong p)
{
    gr_ctx_t cctx, ctx;
    gr_poly_t c, f, g;
    gr_ptr r;
    const gr_poly_struct * rr;
    slong i, k;

    gr_ctx_init_nmod(cctx, p);
    GR_MUST_SUCCEED(gr_ctx_set_is_field(cctx, T_TRUE));
    gr_ctx_init_gr_poly(ctx, cctx);

    gr_poly_init(c, cctx);
    gr_poly_init(f, ctx);
    gr_poly_init(g, ctx);
    r = gr_heap_init(ctx);

    for (i = 0; i < lenA; i++)
    {
        gr_poly_fit_length(c, A[i].length, cctx);
        for (k = 0; k < A[i].length; k++)
            ((nn_ptr) c->coeffs)[k] = fmpz_fdiv_ui(A[i].coeffs + k, p);
        _gr_poly_set_length(c, A[i].length, cctx);
        _gr_poly_normalise(c, cctx);
        GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(f, i, c, ctx));
    }

    for (i = 0; i < lenB; i++)
    {
        gr_poly_fit_length(c, B[i].length, cctx);
        for (k = 0; k < B[i].length; k++)
            ((nn_ptr) c->coeffs)[k] = fmpz_fdiv_ui(B[i].coeffs + k, p);
        _gr_poly_set_length(c, B[i].length, cctx);
        _gr_poly_normalise(c, cctx);
        GR_MUST_SUCCEED(gr_poly_set_coeff_scalar(g, i, c, ctx));
    }

    GR_MUST_SUCCEED(gr_poly_resultant(r, f, g, ctx));

    rr = r;
    for (k = 0; k < outlen; k++)
        out[k] = (k < rr->length) ? ((nn_srcptr) rr->coeffs)[k] : 0;

    gr_heap_clear(r, ctx);
    gr_poly_clear(f, ctx);
    gr_poly_clear(g, ctx);
    gr_poly_clear(c, cctx);
    gr_ctx_clear(ctx);
    gr_ctx_clear(cctx);
}

static void
_crt(fmpz * out, slong len, nn_srcptr residues, nn_srcptr primes, slong num)
{
    fmpz_comb_t comb;
    fmpz_comb_temp_t temp;
    nn_ptr r;
    slong i, k;

    r = flint_malloc(num * sizeof(ulong));

    fmpz_comb_init(comb, primes, num);
    fmpz_comb_temp_init(temp, comb);

    for (k = 0; k < len; k++)
    {
        for (i = 0; i < num; i++)
            r[i] = residues[i * len + k];

        fmpz_multi_CRT_ui(out + k, r, comb, temp, 1);
    }

    fmpz_comb_temp_clear(temp);
    fmpz_comb_clear(comb);
    flint_free(r);
}

/* the images and the reconstruction, the part of the multimodular scheme that
   the choice of primes changes */
static void
_multimodular(fmpz * res, slong outlen, const fmpz_poly_struct * A, slong lenA,
              const fmpz_poly_struct * B, slong lenB,
              nn_srcptr primes, slong num)
{
    nn_ptr residues = flint_malloc(num * outlen * sizeof(ulong));
    slong j;

    for (j = 0; j < num; j++)
        _image(residues + j * outlen, outlen, A, lenA, B, lenB, primes[j]);

    _crt(res, outlen, residues, primes, num);

    flint_free(residues);
}

/* dense bivariate polynomial of length leny in y, with coefficients of length
   lenx in x and `bits` bits */
static void
_randtest_bivariate(fmpz_poly_struct * A, slong leny, slong lenx,
                    flint_bitcnt_t bits, flint_rand_t state)
{
    slong i;

    for (i = 0; i < leny; i++)
    {
        fmpz_poly_randtest_not_zero(A + i, state, lenx, bits);
        fmpz_poly_set_coeff_ui(A + i, lenx - 1, 1 + n_randint(state, 1000));
    }
}

static void
_run(slong leny, slong lenx, flint_bitcnt_t bits, flint_rand_t state)
{
    fmpz_poly_struct * A, * B;
    fmpz * rg, * rf;
    fmpz_t bound, l, u;
    nn_ptr pg, pf;
    double t1, t2, i1, i2, FLINT_SET_BUT_UNUSED(tt);
    slong i, outlen, ng, nf, got;
    flint_bitcnt_t bound_bits, k;

    A = flint_malloc(leny * sizeof(fmpz_poly_struct));
    B = flint_malloc(leny * sizeof(fmpz_poly_struct));
    for (i = 0; i < leny; i++)
    {
        fmpz_poly_init(A + i);
        fmpz_poly_init(B + i);
    }

    _randtest_bivariate(A, leny, lenx, bits, state);
    _randtest_bivariate(B, leny, lenx, bits, state);

    outlen = _outlen(A, leny, B, leny);

    fmpz_init(bound);
    fmpz_init(l);
    fmpz_init(u);

    _bound(bound, A, leny, B, leny);
    bound_bits = fmpz_bits(bound) + 2;

    _fmpz_vec_content(l, A[leny - 1].coeffs, A[leny - 1].length);
    _fmpz_vec_content(u, B[leny - 1].coeffs, B[leny - 1].length);
    fmpz_mul(l, l, u);

    /* the transform has to reach outlen points, and one more level lets the
       algorithm choose between the two interpolations */
    k = n_clog2(outlen) + 1;

    ng = (bound_bits + MODULAR_PRIME_BITS - 1) / MODULAR_PRIME_BITS;
    nf = (bound_bits + FFT_PRIME_BITS - 2) / (FFT_PRIME_BITS - 1);

    /* the choice the multimodular code makes at these sizes */
    flint_printf("%4wd %6wd %8wd %5wu %6wu %4wd %4wd  %s", leny, lenx, outlen,
                 bits, bound_bits, ng, nf,
                 (_gr_poly_resultant_multipoint_cutoff(leny, leny, outlen)
                     && outlen >= 2048) ? "fft" : "gen");
    fflush(stdout);

    if (k >= FFT_PRIME_BITS - 8)
    {
        flint_printf("   too many points for a 50-bit FFT prime\n");
        goto cleanup;
    }

    pg = flint_malloc(ng * sizeof(ulong));
    pf = flint_malloc(nf * sizeof(ulong));

    _generic_primes(pg, ng, l);
    got = _fft_primes(pf, nf, k, l);

    if (got < nf)
    {
        flint_printf("   only %wd FFT primes of exponent %wu\n", got, k);
        flint_free(pg);
        flint_free(pf);
        goto cleanup;
    }

    rg = _fmpz_vec_init(outlen);
    rf = _fmpz_vec_init(outlen);

    /* one image at each kind of prime, which isolates the gain of the DFT
       evaluation from the cost of the extra primes */
    {
        nn_ptr scratch = flint_malloc(outlen * sizeof(ulong));

        TIMEIT_START;
        _image(scratch, outlen, A, leny, B, leny, pg[0]);
        TIMEIT_STOP_VALUES(tt, i1);

        TIMEIT_START;
        _image(scratch, outlen, A, leny, B, leny, pf[0]);
        TIMEIT_STOP_VALUES(tt, i2);

        flint_free(scratch);
    }

    flint_printf("  %8.2es %8.2es %6.2fx", i1, i2, i1 / i2);
    fflush(stdout);

    TIMEIT_START;
    _multimodular(rg, outlen, A, leny, B, leny, pg, ng);
    TIMEIT_STOP_VALUES(tt, t1);

    TIMEIT_START;
    _multimodular(rf, outlen, A, leny, B, leny, pf, nf);
    TIMEIT_STOP_VALUES(tt, t2);

    if (!_fmpz_vec_equal(rg, rf, outlen))
    {
        flint_printf("\nFAIL: the two prime sets disagree\n");
        flint_abort();
    }

    flint_printf("  %8.2es %8.2es %6.2fx\n", t1, t2, t1 / t2);

    _fmpz_vec_clear(rg, outlen);
    _fmpz_vec_clear(rf, outlen);
    flint_free(pg);
    flint_free(pf);

cleanup:

    fmpz_clear(bound);
    fmpz_clear(l);
    fmpz_clear(u);

    for (i = 0; i < leny; i++)
    {
        fmpz_poly_clear(A + i);
        fmpz_poly_clear(B + i);
    }
    flint_free(A);
    flint_free(B);
}

int main(int argc, char ** argv)
{
    flint_rand_t state;

#if !FLINT_HAVE_FFT_SMALL
    flint_printf("fft_small is not available, so no prime admits a DFT\n");
    return 0;
#else

    flint_rand_init(state);

    flint_printf("res_y(f, g) over Z[x][y], with f and g of length leny in y\n");
    flint_printf("and coefficients of length lenx in x and `bits` bits.\n\n");
    flint_printf("The same multimodular scheme is run over two sets of primes:\n");
    flint_printf("`generic` are the primes just above 2^62, as the library picks\n");
    flint_printf("them, and `fft` are the primes p = m 2^k + 1 just below 2^50,\n");
    flint_printf("at which the multipoint algorithm evaluates with a DFT.\n\n");

    flint_printf("`pick` is the choice the multimodular code makes at these\n");
    flint_printf("sizes, and should follow the sign of the last ratio.\n\n");

    flint_printf("leny   lenx  npoints  bits    bnd  #gen #fft pick"
                 "   one image at each prime"
                 "     the whole reconstruction\n");
    flint_printf("                                                 "
                 "  generic       fft  ratio"
                 "    generic       fft  ratio\n");

    if (argc >= 4)
    {
        _run(atol(argv[1]), atol(argv[2]), atol(argv[3]), state);
    }
    else
    {
        /* shapes with a large degree in x, where the evaluation is a large
           enough share of an image for the choice of prime to matter */
        static const slong shape[][2] = {
            {4, 8192}, {6, 8192}, {8, 256}, {8, 2048}, {16, 256},
            {16, 1024}, {16, 8192}, {24, 1024}, {32, 256}, {32, 1024}
        };
        slong i;

        for (i = 0; i < 10; i++)
            _run(shape[i][0], shape[i][1], 16, state);
    }

    flint_rand_clear(state);
    return 0;
#endif
}
