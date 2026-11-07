/*
    Copyright (C) 2026 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Resultants of multivariate polynomials: the Kronecker substitution into the
   dense bivariate algorithm (nmod_mpoly_resultant_kronecker,
   fmpz_mpoly_resultant_kronecker) against the subresultant PRS on the
   univariate forms. Both are timed directly, so the table shows the
   crossover whatever the cutoffs currently are; the `used` column says which
   of the two nmod_mpoly_resultant or fmpz_mpoly_resultant picks, and
   `points` is the number of coefficients of the dense box of the resultant,
   which is the number of evaluation points.

   The output is a markdown table. Any single row can be had from the command
   line:

       p-resultant_multivariate <nvars> <dy> <dx> <density | terms> <bits | n>

   for inputs of degree dy in the variable eliminated and dx in each of the
   others, where a fourth argument above 1 is a number of terms. */

#include <stdlib.h>
#include "ulong_extras.h"
#include "fmpz.h"
#include "mpoly.h"
#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"
#include "profiler.h"

#define MAXVARS 8

/* exponent vector number k of the box (dy + 1) x (dx + 1)^(nvars - 1), the
   variable eliminated being the first */
static void
_box_exp(ulong * exp, slong k, slong nvars, slong dy, slong dx)
{
    slong j;

    exp[0] = k % (dy + 1);
    k /= dy + 1;
    for (j = 1; j < nvars; j++)
    {
        exp[j] = k % (dx + 1);
        k /= dx + 1;
    }
}

static slong
_box_size(slong nvars, slong dy, slong dx)
{
    slong j, box = dy + 1;

    for (j = 1; j < nvars; j++)
        box *= dx + 1;

    return box;
}

/* the monomials of a random polynomial with the announced degrees: the top
   corner and the constant term, then either each monomial of the box with
   probability `size`, or `size` monomials in all if it is above 1 */
static slong
_rand_monomials(ulong * exps, flint_rand_t state, slong nvars, slong dy,
                slong dx, double size)
{
    slong box = _box_size(nvars, dy, dx);
    slong j, k, len = 0;

    for (j = 0; j < nvars; j++)
        exps[j] = (j == 0) ? dy : dx;
    len = 1;
    for (j = 0; j < nvars; j++)
        exps[nvars + j] = 0;
    len = 2;

    if (size <= 1.0)
    {
        for (k = 1; k < box - 1; k++)
            if (n_randint(state, 1000000) < size*1000000)
                _box_exp(exps + nvars*len++, k, nvars, dy, dx);
    }
    else
    {
        for (k = 2; k < FLINT_MIN((slong) size, box); k++)
            _box_exp(exps + nvars*len++, 1 + n_randint(state, box - 2),
                     nvars, dy, dx);
    }

    return len;
}

static void
_rand_nmod(nmod_mpoly_t A, flint_rand_t state, slong nvars, slong dy,
           slong dx, double size, const nmod_mpoly_ctx_t ctx)
{
    slong box = _box_size(nvars, dy, dx);
    ulong * exps = flint_malloc(nvars*(box + 2)*sizeof(ulong));
    slong i, len = _rand_monomials(exps, state, nvars, dy, dx, size);

    nmod_mpoly_zero(A, ctx);
    for (i = 0; i < len; i++)
        nmod_mpoly_push_term_ui_ui(A, 1 + n_randint(state,
                                ctx->mod.n - 1), exps + nvars*i, ctx);
    nmod_mpoly_sort_terms(A, ctx);
    nmod_mpoly_combine_like_terms(A, ctx);

    flint_free(exps);
}

static void
_rand_fmpz(fmpz_mpoly_t A, flint_rand_t state, slong nvars, slong dy,
           slong dx, double size, flint_bitcnt_t bits,
           const fmpz_mpoly_ctx_t ctx)
{
    slong box = _box_size(nvars, dy, dx);
    ulong * exps = flint_malloc(nvars*(box + 2)*sizeof(ulong));
    slong i, len = _rand_monomials(exps, state, nvars, dy, dx, size);
    fmpz_t c;

    fmpz_init(c);
    fmpz_mpoly_zero(A, ctx);
    for (i = 0; i < len; i++)
    {
        do
            fmpz_randtest(c, state, bits);
        while (fmpz_is_zero(c));
        fmpz_mpoly_push_term_fmpz_ui(A, c, exps + nvars*i, ctx);
    }
    fmpz_mpoly_sort_terms(A, ctx);
    fmpz_mpoly_combine_like_terms(A, ctx);

    fmpz_clear(c);
    flint_free(exps);
}

static void
_run_nmod(flint_rand_t state, slong nvars, slong dy, slong dx, double size,
          ulong p)
{
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_t A, B, R, S;
    nmod_mpoly_univar_t Ax, Bx;
    mpoly_resultant_kronecker_t K;
    double t1, t2, FLINT_SET_BUT_UNUSED(tt);
    int used = 0;

    nmod_mpoly_ctx_init(ctx, nvars, ORD_LEX, p);
    nmod_mpoly_init(A, ctx);
    nmod_mpoly_init(B, ctx);
    nmod_mpoly_init(R, ctx);
    nmod_mpoly_init(S, ctx);
    nmod_mpoly_univar_init(Ax, ctx);
    nmod_mpoly_univar_init(Bx, ctx);
    mpoly_resultant_kronecker_init(K, nvars);

    _rand_nmod(A, state, nvars, dy, dx, size, ctx);
    _rand_nmod(B, state, nvars, dy, dx, size, ctx);

    TIMEIT_START;
    nmod_mpoly_to_univar(Ax, A, 0, ctx);
    nmod_mpoly_to_univar(Bx, B, 0, ctx);
    if (!nmod_mpoly_univar_resultant(S, Ax, Bx, ctx))
        flint_abort();
    TIMEIT_STOP_VALUES(tt, t1);

    TIMEIT_START;
    if (!nmod_mpoly_resultant_kronecker(R, A, B, 0, ctx))
        flint_abort();
    TIMEIT_STOP_VALUES(tt, t2);

    if (!nmod_mpoly_equal(R, S, ctx))
    {
        flint_printf("\nFAIL: the two algorithms disagree\n");
        flint_abort();
    }

    if (mpoly_resultant_kronecker_set(K, A->exps, A->length, A->bits,
                    B->exps, B->length, B->bits, 0, 1e300, ctx->minfo))
        used = _nmod_mpoly_resultant_kronecker_cutoff(K, A->length, B->length);

    flint_printf("| %wd | %wd | %wd | %wd | %.0f | %.3es | %.3es | %.2fx | %s |\n",
        nvars, dy, dx, FLINT_MIN(A->length, B->length), K->npoints,
        t1, t2, t1/t2, used ? "dense" : "sparse");
    fflush(stdout);

    mpoly_resultant_kronecker_clear(K);
    nmod_mpoly_univar_clear(Ax, ctx);
    nmod_mpoly_univar_clear(Bx, ctx);
    nmod_mpoly_clear(A, ctx);
    nmod_mpoly_clear(B, ctx);
    nmod_mpoly_clear(R, ctx);
    nmod_mpoly_clear(S, ctx);
    nmod_mpoly_ctx_clear(ctx);
}

static void
_run_fmpz(flint_rand_t state, slong nvars, slong dy, slong dx, double size,
          flint_bitcnt_t bits)
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t A, B, R, S;
    fmpz_mpoly_univar_t Ax, Bx;
    mpoly_resultant_kronecker_t K;
    double t1, t2, FLINT_SET_BUT_UNUSED(tt);
    int used = 0;

    fmpz_mpoly_ctx_init(ctx, nvars, ORD_LEX);
    fmpz_mpoly_init(A, ctx);
    fmpz_mpoly_init(B, ctx);
    fmpz_mpoly_init(R, ctx);
    fmpz_mpoly_init(S, ctx);
    fmpz_mpoly_univar_init(Ax, ctx);
    fmpz_mpoly_univar_init(Bx, ctx);
    mpoly_resultant_kronecker_init(K, nvars);

    _rand_fmpz(A, state, nvars, dy, dx, size, bits, ctx);
    _rand_fmpz(B, state, nvars, dy, dx, size, bits, ctx);

    TIMEIT_START;
    fmpz_mpoly_to_univar(Ax, A, 0, ctx);
    fmpz_mpoly_to_univar(Bx, B, 0, ctx);
    if (!fmpz_mpoly_univar_resultant(S, Ax, Bx, ctx))
        flint_abort();
    TIMEIT_STOP_VALUES(tt, t1);

    TIMEIT_START;
    if (!fmpz_mpoly_resultant_kronecker(R, A, B, 0, ctx))
        flint_abort();
    TIMEIT_STOP_VALUES(tt, t2);

    if (!fmpz_mpoly_equal(R, S, ctx))
    {
        flint_printf("\nFAIL: the two algorithms disagree\n");
        flint_abort();
    }

    if (mpoly_resultant_kronecker_set(K, A->exps, A->length, A->bits,
                    B->exps, B->length, B->bits, 0, 1e300, ctx->minfo))
        used = _fmpz_mpoly_resultant_kronecker_cutoff(K, A->length, B->length);

    flint_printf("| %wd | %wd | %wd | %wu | %wd | %.0f | %.3es | %.3es | %.2fx | %s |\n",
        nvars, dy, dx, bits, FLINT_MIN(A->length, B->length), K->npoints,
        t1, t2, t1/t2, used ? "dense" : "sparse");
    fflush(stdout);

    mpoly_resultant_kronecker_clear(K);
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
    flint_printf("| vars | dy | dx | terms | points | sparse | dense | speedup | used |\n");
    flint_printf("|-----:|---:|---:|------:|-------:|-------:|------:|--------:|:-----|\n");
}

static void _fmpz_header(void)
{
    flint_printf("| vars | dy | dx | bits | terms | points | sparse | dense | speedup | used |\n");
    flint_printf("|-----:|---:|---:|-----:|------:|-------:|-------:|------:|--------:|:-----|\n");
}

int main(int argc, char ** argv)
{
    flint_rand_t state;
    ulong p = n_nextprime(UWORD(1) << 49, 1);
    slong i, j, k;
    static const double sizes[] = { 1.0, 0.1, 8 };

    /* nvars, dy, dx, ordered by cost */
    static const slong shapes[][3] = {
        {3, 4, 4}, {3, 8, 4}, {3, 8, 8}, {3, 16, 8}, {4, 4, 3}, {4, 6, 4},
        {5, 3, 2}, {5, 4, 3}
    };

    /* Over Z the subresultant PRS takes minutes from three variables of
       degree 8 on, so fewer shapes are run there, fewer still at larger
       coefficients */
    static const slong zshapes[][3] = {
        {3, 4, 4}, {4, 3, 3}, {3, 6, 6}, {5, 3, 2}, {3, 8, 8}
    };
    static const flint_bitcnt_t bitsv[] = { 32, 100, 1000 };
    static const slong znshapes[] = { 5, 4, 3 };

    flint_rand_init(state);

    if (argc >= 6)
    {
        if (argv[5][0] == 'n')
            _run_nmod(state, atol(argv[1]), atol(argv[2]), atol(argv[3]),
                      atof(argv[4]), p);
        else
            _run_fmpz(state, atol(argv[1]), atol(argv[2]), atol(argv[3]),
                      atof(argv[4]), atol(argv[5]));
        flint_rand_clear(state);
        return 0;
    }

    flint_printf("res_y(f, g) for f and g of degree dy in y and dx in each of the other\n");
    flint_printf("variables, with each monomial of the box present with probability\n");
    flint_printf("`size`, or `size` monomials if it is above 1. `points` is the number of\n");
    flint_printf("coefficients of the dense box of the resultant. `used` is what the\n");
    flint_printf("dispatch picks.\n\n");

    flint_printf("### nmod_mpoly, 50-bit prime\n\n");
    _nmod_header();

    for (i = 0; i < 8; i++)
        for (j = 0; j < 3; j++)
            _run_nmod(state, shapes[i][0], shapes[i][1], shapes[i][2],
                      sizes[j], p);

    for (k = 0; k < 3; k++)
    {
        flint_printf("\n### fmpz_mpoly, %wu-bit coefficients\n\n", bitsv[k]);
        _fmpz_header();

        for (i = 0; i < znshapes[k]; i++)
            for (j = 0; j < 3; j++)
                _run_fmpz(state, zshapes[i][0], zshapes[i][1], zshapes[i][2],
                          sizes[j], bitsv[k]);
    }

    flint_rand_clear(state);
    return 0;
}
