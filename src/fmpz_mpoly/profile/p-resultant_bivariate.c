/*
    Copyright (C) 2026 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* The dense bivariate resultant against the sparse algorithm that
   nmod_mpoly_resultant and fmpz_mpoly_resultant used to run unconditionally,
   a subresultant PRS over the univariate form. Both are timed directly, so
   the table shows the crossover whatever the cutoffs currently are; the
   `used` column says which of the two the dispatch picks, and `box` is the
   number of entries the dense representation has to materialise, which for a
   sparse input is much more than the input itself.

   The output is a markdown table. */

#include <stdlib.h>
#include "ulong_extras.h"
#include "n_poly.h"
#include "nmod_mpoly.h"
#include "nmod_mpoly_factor.h"
#include "fmpz_poly.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_factor.h"
#include "profiler.h"

/* mirrors the conditions in nmod_mpoly/resultant.c and fmpz_mpoly/resultant.c */
static int
_dispatch_takes_dense(slong lenA, slong blenA, slong termsA,
                      slong lenB, slong blenB, slong termsB,
                      double max_area, double max_points, double free_area)
{
    double areaA = (double) lenA * blenA;
    double areaB = (double) lenB * blenB;
    slong npoints;

    if (areaA > max_area || areaB > max_area)
        return 0;

    if ((double) (lenB - 1) * (blenA - 1) + (double) (lenA - 1) * (blenB - 1)
            + 1.0 > max_points)
        return 0;

    if (areaA > free_area && areaA > 64.0 * termsA)
        return 0;

    if (areaB > free_area && areaB > 64.0 * termsB)
        return 0;

    if (FLINT_MIN(termsA, termsB) < 5 && lenA * lenB < 256)
        return 0;

    npoints = (lenB - 1) * (blenA - 1) + (lenA - 1) * (blenB - 1) + 1;

    return n_bpoly_mod_resultant_cutoff(FLINT_MAX(lenA, lenB),
                                        FLINT_MIN(lenA, lenB), npoints);
}

/* random bivariate polynomial of degree dy in the variable the resultant
   eliminates and dx in the other, with about `terms` monomials. The two
   corners are always present, so that the degrees are the announced ones. */
static void
_rand_nmod(nmod_mpoly_t A, flint_rand_t state, slong dy, slong dx,
           slong terms, ulong p, const nmod_mpoly_ctx_t ctx)
{
    ulong bounds[2], exp[2];

    bounds[0] = dy + 1;
    bounds[1] = dx + 1;

    nmod_mpoly_randtest_bounds(A, state, terms, bounds, ctx);

    exp[0] = dy; exp[1] = dx;
    nmod_mpoly_set_coeff_ui_ui(A, 1 + n_randint(state, p - 1), exp, ctx);
    exp[0] = 0; exp[1] = 0;
    nmod_mpoly_set_coeff_ui_ui(A, 1 + n_randint(state, p - 1), exp, ctx);
}

static void
_rand_fmpz(fmpz_mpoly_t A, flint_rand_t state, slong dy, slong dx,
           slong terms, flint_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx)
{
    ulong bounds[2], exp[2];

    bounds[0] = dy + 1;
    bounds[1] = dx + 1;

    fmpz_mpoly_randtest_bounds(A, state, terms, bits, bounds, ctx);

    exp[0] = dy; exp[1] = dx;
    fmpz_mpoly_set_coeff_ui_ui(A, 1 + n_randint(state, 1000), exp, ctx);
    exp[0] = 0; exp[1] = 0;
    fmpz_mpoly_set_coeff_ui_ui(A, 1 + n_randint(state, 1000), exp, ctx);
}

static void
_run_nmod(flint_rand_t state, slong dy, slong dx, double density, ulong p)
{
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_t A, B, R, S;
    nmod_mpoly_univar_t Ax, Bx;
    n_bpoly_t Ab, Bb, Rb;
    n_poly_t rp;
    double t1, t2, FLINT_SET_BUT_UNUSED(tt);
    slong terms = (slong) (density * (dy + 1) * (dx + 1)) + 1;
    int dense;

    nmod_mpoly_ctx_init(ctx, 2, ORD_LEX, p);
    nmod_mpoly_init(A, ctx);
    nmod_mpoly_init(B, ctx);
    nmod_mpoly_init(R, ctx);
    nmod_mpoly_init(S, ctx);
    nmod_mpoly_univar_init(Ax, ctx);
    nmod_mpoly_univar_init(Bx, ctx);
    n_bpoly_init(Ab);
    n_bpoly_init(Bb);
    n_bpoly_init(Rb);
    n_poly_init(rp);

    _rand_nmod(A, state, dy, dx, terms, p, ctx);
    _rand_nmod(B, state, dy, dx, terms, p, ctx);

    TIMEIT_START;
    nmod_mpoly_to_univar(Ax, A, 0, ctx);
    nmod_mpoly_to_univar(Bx, B, 0, ctx);
    if (!nmod_mpoly_univar_resultant(S, Ax, Bx, ctx))
        flint_abort();
    TIMEIT_STOP_VALUES(tt, t1);

    TIMEIT_START;
    nmod_mpoly_get_bpoly(Ab, A, 0, 1, ctx);
    nmod_mpoly_get_bpoly(Bb, B, 0, 1, ctx);
    if (!n_bpoly_mod_resultant(rp, Ab, Bb, ctx->mod))
        flint_abort();
    n_bpoly_fit_length(Rb, 1);
    n_poly_set(Rb->coeffs + 0, rp);
    Rb->length = n_poly_is_zero(rp) ? 0 : 1;
    nmod_mpoly_set_bpoly(R, A->bits, Rb, 0, 1, ctx);
    TIMEIT_STOP_VALUES(tt, t2);

    if (!nmod_mpoly_equal(R, S, ctx))
    {
        flint_printf("\nFAIL: the two algorithms disagree\n");
        flint_abort();
    }

    dense = _dispatch_takes_dense(dy + 1, dx + 1, A->length,
                                  dy + 1, dx + 1, B->length,
                                  4294967296.0, 4294967296.0, 67108864.0);

    flint_printf("| %wd | %wd | %.4f | %wd | %wd | %.3es | %.3es | %.2fx | %s |\n",
        dy, dx, density, FLINT_MIN(A->length, B->length),
        (dy + 1) * (dx + 1), t1, t2, t1 / t2, dense ? "dense" : "sparse");
    fflush(stdout);

    n_poly_clear(rp);
    n_bpoly_clear(Ab);
    n_bpoly_clear(Bb);
    n_bpoly_clear(Rb);
    nmod_mpoly_univar_clear(Ax, ctx);
    nmod_mpoly_univar_clear(Bx, ctx);
    nmod_mpoly_clear(A, ctx);
    nmod_mpoly_clear(B, ctx);
    nmod_mpoly_clear(R, ctx);
    nmod_mpoly_clear(S, ctx);
    nmod_mpoly_ctx_clear(ctx);
}

static void
_run_fmpz(flint_rand_t state, slong dy, slong dx, double density,
          flint_bitcnt_t bits)
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t A, B, R, S;
    fmpz_mpoly_univar_t Ax, Bx;
    fmpz_bpoly_t Ab, Bb, Rb;
    fmpz_poly_t rp;
    double t1, t2, FLINT_SET_BUT_UNUSED(tt);
    slong terms = (slong) (density * (dy + 1) * (dx + 1)) + 1;
    int dense;

    fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
    fmpz_mpoly_init(A, ctx);
    fmpz_mpoly_init(B, ctx);
    fmpz_mpoly_init(R, ctx);
    fmpz_mpoly_init(S, ctx);
    fmpz_mpoly_univar_init(Ax, ctx);
    fmpz_mpoly_univar_init(Bx, ctx);
    fmpz_bpoly_init(Ab);
    fmpz_bpoly_init(Bb);
    fmpz_bpoly_init(Rb);
    fmpz_poly_init(rp);

    _rand_fmpz(A, state, dy, dx, terms, bits, ctx);
    _rand_fmpz(B, state, dy, dx, terms, bits, ctx);

    TIMEIT_START;
    fmpz_mpoly_to_univar(Ax, A, 0, ctx);
    fmpz_mpoly_to_univar(Bx, B, 0, ctx);
    if (!fmpz_mpoly_univar_resultant(S, Ax, Bx, ctx))
        flint_abort();
    TIMEIT_STOP_VALUES(tt, t1);

    TIMEIT_START;
    fmpz_mpoly_get_bpoly(Ab, A, 0, 1, ctx);
    fmpz_mpoly_get_bpoly(Bb, B, 0, 1, ctx);
    if (!fmpz_bpoly_resultant(rp, Ab, Bb, 1))
        flint_abort();
    fmpz_bpoly_fit_length(Rb, 1);
    fmpz_poly_set(Rb->coeffs + 0, rp);
    Rb->length = fmpz_poly_is_zero(rp) ? 0 : 1;
    fmpz_mpoly_set_fmpz_bpoly(R, A->bits, Rb, 0, 1, ctx);
    TIMEIT_STOP_VALUES(tt, t2);

    if (!fmpz_mpoly_equal(R, S, ctx))
    {
        flint_printf("\nFAIL: the two algorithms disagree\n");
        flint_abort();
    }

    dense = _dispatch_takes_dense(dy + 1, dx + 1, A->length,
                                  dy + 1, dx + 1, B->length,
                                  1073741824.0, 268435456.0, 8388608.0);

    flint_printf("| %wd | %wd | %.4f | %wu | %wd | %wd | %.3es | %.3es | %.2fx | %s |\n",
        dy, dx, density, bits, FLINT_MIN(A->length, B->length),
        (dy + 1) * (dx + 1), t1, t2, t1 / t2, dense ? "dense" : "sparse");
    fflush(stdout);

    fmpz_poly_clear(rp);
    fmpz_bpoly_clear(Ab);
    fmpz_bpoly_clear(Bb);
    fmpz_bpoly_clear(Rb);
    fmpz_mpoly_univar_clear(Ax, ctx);
    fmpz_mpoly_univar_clear(Bx, ctx);
    fmpz_mpoly_clear(A, ctx);
    fmpz_mpoly_clear(B, ctx);
    fmpz_mpoly_clear(R, ctx);
    fmpz_mpoly_clear(S, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

static void _nmod_header(void)
{
    flint_printf("| dy | dx | density | terms | box | sparse | dense | speedup | used |\n");
    flint_printf("|---:|---:|--------:|------:|----:|-------:|------:|--------:|:-----|\n");
}

static void _fmpz_header(void)
{
    flint_printf("| dy | dx | density | bits | terms | box | sparse | dense | speedup | used |\n");
    flint_printf("|---:|---:|--------:|-----:|------:|----:|-------:|------:|--------:|:-----|\n");
}

int main(int argc, char ** argv)
{
    flint_rand_t state;
    ulong p = n_nextprime(UWORD(1) << 49, 1);
    slong i, j, k;
    static const double densities[] = { 1.0, 0.25, 0.0625, 0.015625, 0.00390625 };
    static const slong shapes[][2] = {
        {4, 16}, {8, 8}, {8, 64}, {16, 16}, {16, 128}, {24, 24}, {32, 32}, {64, 64}
    };

    /* Over Z the sparse reference grows steeply in both the degrees and the
       coefficient size, to hours at the larger combinations, so the shapes
       are ordered by cost and only the cheaper ones are run at the larger
       coefficients. Any single row can still be had from the command line:

           p-resultant_bivariate <dy> <dx> <density> <bits | n>            */
    static const slong zshapes[][2] = {
        {4, 16}, {8, 8}, {16, 16}, {8, 64}, {24, 24}
    };
    static const flint_bitcnt_t bitsv[] = { 32, 100, 1000, 10000 };
    static const slong znshapes[] = { 5, 4, 2, 2 };
    static const slong zfirstdens[] = { 0, 0, 0, 1 };

    flint_rand_init(state);

    if (argc >= 5)
    {
        if (argv[4][0] == 'n')
            _run_nmod(state, atol(argv[1]), atol(argv[2]), atof(argv[3]), p);
        else
            _run_fmpz(state, atol(argv[1]), atol(argv[2]), atof(argv[3]),
                      atol(argv[4]));
        flint_rand_clear(state);
        return 0;
    }

    flint_printf("res_y(f, g) for bivariate f and g of degree dy in y and dx in x,\n");
    flint_printf("with `density` of the terms of the box present. `box` is the number\n");
    flint_printf("of entries the dense representation materialises, against `terms`\n");
    flint_printf("in the input (the smaller of the two). `used` is what the dispatch picks.\n\n");

    flint_printf("### nmod_mpoly, 50-bit prime\n\n");
    _nmod_header();

    for (i = 0; i < 8; i++)
        for (j = 0; j < 5; j++)
            _run_nmod(state, shapes[i][0], shapes[i][1], densities[j], p);

    for (k = 0; k < 4; k++)
    {
        flint_printf("\n### fmpz_mpoly, %wu-bit coefficients\n\n", bitsv[k]);
        _fmpz_header();

        for (i = 0; i < znshapes[k]; i++)
            for (j = zfirstdens[k]; j < 5; j++)
                _run_fmpz(state, zshapes[i][0], zshapes[i][1], densities[j],
                          bitsv[k]);
    }

    flint_rand_clear(state);
    return 0;
}
