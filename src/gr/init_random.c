/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include "mpoly.h"
#include "gr.h"

/* For random composite rings, some base rings that don't require
   memory allocation. */
static gr_ctx_struct _gr_some_base_rings[11];

static gr_ctx_struct *
_gr_random_base_ring(flint_rand_t state)
{
    gr_ctx_init_fmpz(_gr_some_base_rings + 0);
    gr_ctx_init_nmod(_gr_some_base_rings + 1, 1);
    gr_ctx_init_nmod(_gr_some_base_rings + 2, 2);
    gr_ctx_init_nmod(_gr_some_base_rings + 3, 11);
    gr_ctx_init_nmod(_gr_some_base_rings + 4, 12);
    gr_ctx_init_nmod(_gr_some_base_rings + 5, 257);
    gr_ctx_init_fmpq(_gr_some_base_rings + 6);
    gr_ctx_init_real_arb(_gr_some_base_rings + 7, 64);
    gr_ctx_init_real_arb(_gr_some_base_rings + 8, 256);
    gr_ctx_init_complex_acb(_gr_some_base_rings + 9, 64);
    gr_ctx_init_complex_acb(_gr_some_base_rings + 10, 256);

    return _gr_some_base_rings + n_randint(state, sizeof(_gr_some_base_rings) / sizeof(gr_ctx_struct));
}

void
gr_ctx_init_random_ring_composite(gr_ctx_t ctx, flint_rand_t state)
{
    gr_ctx_struct * base_ring = _gr_random_base_ring(state);

    switch (n_randint(state, 5))
    {
        case 0:
            gr_ctx_init_gr_poly(ctx, base_ring);
            break;
        case 1:
            gr_ctx_init_gr_mpoly(ctx, base_ring, n_randint(state, 3), mpoly_ordering_randtest(state));
            break;
        case 2:
            gr_ctx_init_gr_poly(ctx, base_ring);
/*
    this currently breaks some tests
            gr_ctx_init_gr_series(ctx, base_ring, n_randint(state, 6));
*/
            break;
        case 3:
            gr_ctx_init_gr_poly(ctx, base_ring);
/*
    this currently breaks some tests
            gr_ctx_init_gr_series_mod(ctx, base_ring, n_randint(state, 6));
            break;
*/
        case 4:
            gr_ctx_init_vector_space_gr_vec(ctx, base_ring, n_randint(state, 4));
            break;
/*
    this will break tests that currently assume commutativity

        case 5:
            gr_ctx_init_matrix_ring(ctx, base_ring, n_randint(state, 4));
            break;
*/
    }
}

void
gr_ctx_init_random_ring_integers_mod(gr_ctx_t ctx, flint_rand_t state)
{
    fmpz_t t;

    switch (n_randint(state, 4))
    {
        case 0:
            gr_ctx_init_nmod8(ctx, n_randtest(state) % 255 + 1);
            break;
        case 1:
            gr_ctx_init_nmod32(ctx, n_randtest(state) % UWORD(4294967295) + 1);
            break;
        case 2:
            gr_ctx_init_nmod(ctx, n_randtest_not_zero(state));
            break;
        case 3:
            fmpz_init(t);
            fmpz_randtest_not_zero(t, state, 100);
            fmpz_abs(t, t);
            gr_ctx_init_fmpz_mod(ctx, t);
            fmpz_clear(t);
            break;
    }
}

void
gr_ctx_init_random_ring_finite_field(gr_ctx_t ctx, flint_rand_t state)
{
    fmpz_t t;
    fmpz_init(t);

    switch (n_randint(state, 3))
    {
        case 0:
            fmpz_set_ui(t, n_randtest_prime(state, 0));
            gr_ctx_init_fq_nmod(ctx, t, 1 + n_randint(state, 4), NULL);
            break;

        case 1:
            fmpz_set_ui(t, n_randprime(state, 4, 0));
            gr_ctx_init_fq_zech(ctx, t, 1 + n_randint(state, 3), NULL);
            break;

        case 2:
            fmpz_randprime(t, state, 2 + n_randint(state, 100), 0);
            gr_ctx_init_fq(ctx, t, 1 + n_randint(state, 4), NULL);
            break;
    }

    fmpz_clear(t);
}

void
gr_ctx_init_random_ring_number_field(gr_ctx_t ctx, flint_rand_t state)
{
    fmpz_poly_t g;
    fmpq_poly_t f;

    fmpz_poly_init(g);
    fmpq_poly_init(f);

    do
    {
        fmpz_poly_randtest_irreducible(g, state, 2 + n_randint(state, 5), 1 + n_randint(state, 10));
    } while (g->length < 2);

    fmpq_poly_set_fmpz_poly(f, g);
    fmpq_poly_scalar_div_ui(f, f, 1 + n_randtest(state) % 256);

    gr_ctx_init_nf(ctx, f);

    fmpz_poly_clear(g);
    fmpq_poly_clear(f);
}

void
gr_ctx_init_random_ring_real_complex_ball(gr_ctx_t ctx, flint_rand_t state)
{
    if (n_randint(state, 2))
        gr_ctx_init_real_arb(ctx, 2 + n_randint(state, 200));
    else
        gr_ctx_init_complex_acb(ctx, 2 + n_randint(state, 200));
}

void
gr_ctx_init_random_ring_real_complex_exact(gr_ctx_t ctx, flint_rand_t state)
{
    switch (n_randint(state, 4))
    {
        case 0:
            gr_ctx_init_real_ca(ctx);
            break;
        case 1:
            gr_ctx_init_complex_ca(ctx);
            break;
        case 2:
            gr_ctx_init_real_algebraic_ca(ctx);
            break;
        case 3:
            gr_ctx_init_complex_algebraic_ca(ctx);
            break;
/*
    slow -- we'll want to enable these selectively
        case 4:
            gr_ctx_init_real_qqbar(ctx);
            break;
        case 5:
            gr_ctx_init_complex_qqbar(ctx);
            break;
*/
    }
}

void gr_ctx_init_random(gr_ctx_t ctx, flint_rand_t state)
{
    switch (n_randint(state, 11))
    {
        case 0:
        case 1:
        case 2:
            gr_ctx_init_fmpz(ctx);
            break;
        case 3:
            gr_ctx_init_fmpq(ctx);
            break;
        case 4:
            gr_ctx_init_fmpzi(ctx);
            break;
        case 5:
            gr_ctx_init_random_ring_integers_mod(ctx, state);
            break;
        case 6:
            gr_ctx_init_random_ring_finite_field(ctx, state);
            break;
        case 7:
            gr_ctx_init_random_ring_number_field(ctx, state);
            break;
        case 8:
            gr_ctx_init_random_ring_real_complex_ball(ctx, state);
            break;
        case 9:
            gr_ctx_init_random_ring_real_complex_exact(ctx, state);
            break;
        case 10:
            gr_ctx_init_random_ring_composite(ctx, state);
            break;
    }

/*
    flint_printf("  ");
    gr_ctx_println(ctx);
*/
}
