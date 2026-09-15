/*
    Copyright (C) 2026 Mael Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* nmod_mpoly_resultant and fmpz_mpoly_resultant against the sparse algorithm
   they used to run unconditionally, which is a subresultant PRS over the
   univariate form. Bivariate inputs of a given degree in each variable are
   generated with a given fraction of the terms of their bounding box, since
   that fraction is what decides whether converting to the dense bivariate
   representation pays.

   The output is a markdown table. */

#include <stdlib.h>
#include "ulong_extras.h"
#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"
#include "profiler.h"

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
    double t1, t2, FLINT_SET_BUT_UNUSED(tt);
    slong terms = (slong) (density * (dy + 1) * (dx + 1)) + 1;

    nmod_mpoly_ctx_init(ctx, 2, ORD_LEX, p);
    nmod_mpoly_init(A, ctx);
    nmod_mpoly_init(B, ctx);
    nmod_mpoly_init(R, ctx);
    nmod_mpoly_init(S, ctx);
    nmod_mpoly_univar_init(Ax, ctx);
    nmod_mpoly_univar_init(Bx, ctx);

    _rand_nmod(A, state, dy, dx, terms, p, ctx);
    _rand_nmod(B, state, dy, dx, terms, p, ctx);

    /* the sparse algorithm, as nmod_mpoly_resultant ran it before */
    TIMEIT_START;
    nmod_mpoly_to_univar(Ax, A, 0, ctx);
    nmod_mpoly_to_univar(Bx, B, 0, ctx);
    if (!nmod_mpoly_univar_resultant(S, Ax, Bx, ctx))
        flint_abort();
    TIMEIT_STOP_VALUES(tt, t1);

    TIMEIT_START;
    if (!nmod_mpoly_resultant(R, A, B, 0, ctx))
        flint_abort();
    TIMEIT_STOP_VALUES(tt, t2);

    if (!nmod_mpoly_equal(R, S, ctx))
    {
        flint_printf("\nFAIL: the two algorithms disagree\n");
        flint_abort();
    }

    flint_printf("| %wd | %wd | %.3f | %wd | %.3es | %.3es | %.2fx |\n",
        dy, dx, density, A->length, t1, t2, t1 / t2);
    fflush(stdout);

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
    double t1, t2, FLINT_SET_BUT_UNUSED(tt);
    slong terms = (slong) (density * (dy + 1) * (dx + 1)) + 1;

    fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
    fmpz_mpoly_init(A, ctx);
    fmpz_mpoly_init(B, ctx);
    fmpz_mpoly_init(R, ctx);
    fmpz_mpoly_init(S, ctx);
    fmpz_mpoly_univar_init(Ax, ctx);
    fmpz_mpoly_univar_init(Bx, ctx);

    _rand_fmpz(A, state, dy, dx, terms, bits, ctx);
    _rand_fmpz(B, state, dy, dx, terms, bits, ctx);

    TIMEIT_START;
    fmpz_mpoly_to_univar(Ax, A, 0, ctx);
    fmpz_mpoly_to_univar(Bx, B, 0, ctx);
    if (!fmpz_mpoly_univar_resultant(S, Ax, Bx, ctx))
        flint_abort();
    TIMEIT_STOP_VALUES(tt, t1);

    TIMEIT_START;
    if (!fmpz_mpoly_resultant(R, A, B, 0, ctx))
        flint_abort();
    TIMEIT_STOP_VALUES(tt, t2);

    if (!fmpz_mpoly_equal(R, S, ctx))
    {
        flint_printf("\nFAIL: the two algorithms disagree\n");
        flint_abort();
    }

    flint_printf("| %wd | %wd | %.3f | %wu | %wd | %.3es | %.3es | %.2fx |\n",
        dy, dx, density, bits, A->length, t1, t2, t1 / t2);
    fflush(stdout);

    fmpz_mpoly_univar_clear(Ax, ctx);
    fmpz_mpoly_univar_clear(Bx, ctx);
    fmpz_mpoly_clear(A, ctx);
    fmpz_mpoly_clear(B, ctx);
    fmpz_mpoly_clear(R, ctx);
    fmpz_mpoly_clear(S, ctx);
    fmpz_mpoly_ctx_clear(ctx);
}

int main(int argc, char ** argv)
{
    flint_rand_t state;
    ulong p = n_nextprime(UWORD(1) << 49, 1);
    slong i, j;
    static const double densities[] = { 1.0, 0.25, 0.0625, 0.015625, 0.00390625 };
    static const slong shapes[][2] = {
        {4, 16}, {4, 256}, {8, 8}, {8, 64}, {16, 16},
        {16, 128}, {24, 24}, {32, 32}, {48, 48}, {64, 64}
    };

    flint_rand_init(state);

    if (argc >= 4)
    {
        _run_nmod(state, atol(argv[1]), atol(argv[2]), atof(argv[3]), p);
        _run_fmpz(state, atol(argv[1]), atol(argv[2]), atof(argv[3]),
                  argc >= 5 ? atol(argv[4]) : 32);
        flint_rand_clear(state);
        return 0;
    }

    flint_printf("res_y(f, g) for bivariate f and g of degree dy in y and dx\n");
    flint_printf("in x, with `density` of the terms of the box present.\n");
    flint_printf("`sparse` is the subresultant PRS over the univariate form,\n");
    flint_printf("`dense` is what the resultant functions now dispatch to.\n\n");

    flint_printf("### nmod_mpoly_resultant, 50-bit prime\n\n");
    flint_printf("| dy | dx | density | terms | sparse | dense | speedup |\n");
    flint_printf("|---:|---:|--------:|------:|-------:|------:|--------:|\n");

    for (i = 0; i < 10; i++)
        for (j = 0; j < 5; j++)
            _run_nmod(state, shapes[i][0], shapes[i][1], densities[j], p);

    flint_printf("\n### fmpz_mpoly_resultant\n\n");
    flint_printf("| dy | dx | density | bits | terms | sparse | dense | speedup |\n");
    flint_printf("|---:|---:|--------:|-----:|------:|-------:|------:|--------:|\n");

    for (i = 0; i < 7; i++)
        for (j = 0; j < 5; j++)
            _run_fmpz(state, shapes[i][0], shapes[i][1], densities[j], 32);

    flint_rand_clear(state);
    return 0;
}
