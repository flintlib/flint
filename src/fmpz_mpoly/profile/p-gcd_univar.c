/*
    Copyright (C) 2026 Max Horn

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* usage
make profile MOD=fmpz_mpoly && ./build/fmpz_mpoly/profile/p-gcd_univar [univar|control]

Times fmpz_mpoly_gcd on pairs where one input is a monomial times a univariate
polynomial ("univar"), and on general pairs as a control ("control"). Run the
two separately when comparing against another build, so that the long univar
section cannot skew the control through thermal throttling.
*/

#include <string.h>
#include "profiler.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

/*
    A polynomial with `terms` terms, of degree less than 3 in the first variable
    and total degree 3 in `used` further ones, in a ring that has more variables
    than that. This is the shape a numerator has when fractions over a multivariate
    ring are kept in lowest terms and the ring carries more variables than any
    one element uses.
*/
static void _make_numerator(
    fmpz_mpoly_t A,
    slong terms,
    slong used,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    ulong * e = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, ulong);

    FLINT_ASSERT(used + 1 <= ctx->minfo->nvars);

    for (i = 0; i < ctx->minfo->nvars; i++)
        e[i] = 0;

    fmpz_mpoly_zero(A, ctx);
    for (i = 0; i < terms; i++)
    {
        slong u = 1 + (i % used);
        slong v = 1 + ((i/used) % used);
        slong w = 1 + ((i/(used*used)) % used);

        e[0] = i % 3;
        e[u] += 1;
        e[v] += 1;
        e[w] += 1;
        fmpz_mpoly_push_term_si_ui(A, (i & 1) ? 1 + (i % 7) : -1 - (i % 5), e, ctx);
        e[0] = 0;
        e[u] = 0;
        e[v] = 0;
        e[w] = 0;
    }
    fmpz_mpoly_sort_terms(A, ctx);
    fmpz_mpoly_combine_like_terms(A, ctx);

    flint_free(e);
}

/* x1^deg - 1, the shape of the denominators such fractions have */
static void _make_denominator(
    fmpz_mpoly_t B,
    slong deg,
    const fmpz_mpoly_ctx_t ctx)
{
    ulong * e = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, ulong);
    slong i;

    for (i = 0; i < ctx->minfo->nvars; i++)
        e[i] = 0;

    fmpz_mpoly_zero(B, ctx);
    e[0] = deg;
    fmpz_mpoly_push_term_si_ui(B, 1, e, ctx);
    e[0] = 0;
    fmpz_mpoly_push_term_si_ui(B, -1, e, ctx);
    fmpz_mpoly_sort_terms(B, ctx);

    flint_free(e);
}

/* microseconds per gcd */
static double _time_gcd(
    const fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_t G;
    timeit_t timer;
    slong reps = 1;

    fmpz_mpoly_init(G, ctx);

    /* aim for at least 50ms of work */
    do {
        slong i;
        reps *= 2;
        timeit_start(timer);
        for (i = 0; i < reps; i++)
            fmpz_mpoly_gcd(G, A, B, ctx);
        timeit_stop(timer);
    } while (timer->wall < 50 && reps < WORD(1) << 24);

    fmpz_mpoly_clear(G, ctx);

    return 1000.0*timer->wall/reps;
}

int main(int argc, char ** argv)
{
    static const slong nvars_set[] = {13, 52};
    static const slong terms_set[] = {28, 500};
    /* the last degree is past UNIVAR_DIVISOR_MAX_DEG, so it falls back */
    static const slong deg_set[] = {2, 12, 1024};
    slong nvars, terms, deg, i, in, it, id;
    int do_univar = (argc < 2) || (strcmp(argv[1], "univar") == 0);
    int do_control = (argc < 2) || (strcmp(argv[1], "control") == 0);

    if (do_univar)
    {

    flint_printf("--- one input a monomial times a univariate polynomial ---\n");
    flint_printf("A: `terms` terms, degree < 3 in x1, cubic in 12 further variables\n");
    flint_printf("B: x1^deg - 1\n");
    flint_printf("cofactor: the same, with the factor x1 - 1 forced into A too\n\n");
    flint_printf("%5s %6s %5s %9s %10s  %s\n",
                            "nvars", "terms", "deg", "coprime", "cofactor", "gcd");

    for (in = 0; in < 2; in++)
    for (it = 0; it < 2; it++)
    for (id = 0; id < 3; id++)
    {
        nvars = nvars_set[in];
        terms = terms_set[it];
        deg = deg_set[id];
        {
            fmpz_mpoly_ctx_t ctx;
            fmpz_mpoly_t A, B, G, f;
            double t1, t2;
            slong alen;

            fmpz_mpoly_ctx_init(ctx, nvars, ORD_LEX);
            fmpz_mpoly_init(A, ctx);
            fmpz_mpoly_init(B, ctx);
            fmpz_mpoly_init(G, ctx);
            fmpz_mpoly_init(f, ctx);

            _make_numerator(A, terms, 12, ctx);
            _make_denominator(B, deg, ctx);
            alen = A->length;
            t1 = _time_gcd(A, B, ctx);

            /* now force x_0 - 1, a factor of B, into A as well */
            _make_denominator(f, 1, ctx);
            fmpz_mpoly_mul(A, A, f, ctx);
            t2 = _time_gcd(A, B, ctx);
            fmpz_mpoly_gcd(G, A, B, ctx);

            flint_printf("%5wd %6wd %5wd %9.2f %10.2f  ",
                                              nvars, alen, deg, t1, t2);
            fmpz_mpoly_print_pretty(G, NULL, ctx);
            flint_printf("\n");

            fmpz_mpoly_clear(f, ctx);
            fmpz_mpoly_clear(G, ctx);
            fmpz_mpoly_clear(B, ctx);
            fmpz_mpoly_clear(A, ctx);
            fmpz_mpoly_ctx_clear(ctx);
        }
    }

    }

    if (!do_control)
    {
        flint_cleanup_master();
        return 0;
    }

    flint_printf("--- control: general pairs sharing a factor ---\n");
    flint_printf("%5s %6s %10s  %s\n", "nvars", "terms", "us", "gcd terms");

    for (nvars = 2; nvars <= 6; nvars++)
    for (terms = 5; terms <= 20; terms += 15)
    {
        flint_rand_t state;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t A, B, G, a, b, t;
        double total = 0;

        flint_rand_init(state);
        flint_rand_set_seed(state, 12345, 6789);

        fmpz_mpoly_ctx_init(ctx, nvars, ORD_LEX);
        fmpz_mpoly_init(A, ctx); fmpz_mpoly_init(B, ctx); fmpz_mpoly_init(G, ctx);
        fmpz_mpoly_init(a, ctx); fmpz_mpoly_init(b, ctx); fmpz_mpoly_init(t, ctx);

        for (i = 0; i < 5; i++)
        {
            do {
                fmpz_mpoly_randtest_bound(t, state, 1 + terms/4, 40, 3, ctx);
            } while (fmpz_mpoly_is_zero(t, ctx));
            fmpz_mpoly_randtest_bound(a, state, terms, 40, 3, ctx);
            fmpz_mpoly_randtest_bound(b, state, terms, 40, 3, ctx);
            fmpz_mpoly_mul(A, a, t, ctx);
            fmpz_mpoly_mul(B, b, t, ctx);
            if (fmpz_mpoly_is_zero(A, ctx) || fmpz_mpoly_is_zero(B, ctx))
                continue;
            total += _time_gcd(A, B, ctx);
        }

        fmpz_mpoly_gcd(G, A, B, ctx);
        flint_printf("%5wd %6wd %10.2f  %wd\n", nvars, terms, total, G->length);

        fmpz_mpoly_clear(t, ctx); fmpz_mpoly_clear(b, ctx); fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(G, ctx); fmpz_mpoly_clear(B, ctx); fmpz_mpoly_clear(A, ctx);
        fmpz_mpoly_ctx_clear(ctx);
        flint_rand_clear(state);
    }

    flint_cleanup_master();
    return 0;
}
