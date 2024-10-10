/*
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fq_nmod.h"
#include "fq.h"

/* Assumes that everything fits inside a small fmpz */
static void
fq_ctx_init_set_clear_small_fq_nmod_ctx(fq_ctx_t c1, fq_nmod_ctx_t c2)
{
    fmpz_mod_ctx_struct * ctxp = c1->ctxp;
    fmpz_mod_poly_struct * m1 = c1->modulus;
    fmpz_mod_poly_struct * i1 = c1->inv;
    nmod_poly_struct * m2 = c2->modulus;
    nmod_poly_struct * i2 = c2->inv;

    /* Init c1->ctxp */
    *ctxp->n = c2->mod.n;
    ctxp->add_fxn = _fmpz_mod_add1;
    ctxp->sub_fxn = _fmpz_mod_sub1;
    ctxp->mul_fxn = _fmpz_mod_mul1;
    ctxp->mod = c2->mod;
    ctxp->n_limbs[0] = 0;
    ctxp->n_limbs[1] = 0;
    ctxp->n_limbs[2] = 0;
    ctxp->ninv_huge = NULL;

    c1->sparse_modulus = c2->sparse_modulus;
    c1->is_conway = c2->is_conway;

    c1->a = (fmpz *) c2->a;
    c1->j = c2->j;
    c1->len = c2->len;

    /* Init c1->modulus and c1->inv */
    m1->coeffs = (fmpz *) m2->coeffs;
    m1->alloc = m2->alloc;
    m1->length = m2->length;
    i1->coeffs = (fmpz *) i2->coeffs;
    i1->alloc = i2->alloc;
    i1->length = i2->length;

    c1->var = c2->var;
}

int
_fq_ctx_init_conway_ui(fq_ctx_t ctx, ulong p, slong d, const char * var)
{
    fq_nmod_ctx_t ctx2;

    if (_fq_nmod_ctx_init_conway_ui(ctx2, p, d, var))
    {
        fq_ctx_init_set_clear_small_fq_nmod_ctx(ctx, ctx2);
        return 1;
    }
    else
        return 0;
}

void
fq_ctx_init_conway_ui(fq_ctx_t ctx, ulong p, slong d, const char *var)
{
    int result;

    result = _fq_ctx_init_conway_ui(ctx, p, d, var);

    if (!result)
        flint_throw(FLINT_ERROR, "Exception (fq_ctx_init_conway_ui).  The "
                "polynomial for (p, d) = (%wu, %wd) is not present "
                "in the database.\n", p, d);

    ctx->is_conway = 1;
}

int
_fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char * var)
{
    if (*p < 2 || *p > 109987) /* NOTE: This works. */
        return 0;
    else
        return _fq_ctx_init_conway_ui(ctx, *p, d, var);
}

void
fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char * var)
{
    if (*p < 2 || *p > 109987) /* NOTE: This works. */
        flint_throw(FLINT_ERROR, "Exception (fq_ctx_init_conway).  The "
                "polynomial for (p, d) = (%{fmpz}, %wd) is not present "
                "in the database.\n", p, d);
    else
        fq_ctx_init_conway_ui(ctx, *p, d, var);
}

void
fq_ctx_init(fq_ctx_t ctx, const fmpz_t p, slong d, const char * var)
{
    flint_rand_t state;
    fmpz_mod_poly_t poly;
    fmpz_mod_ctx_t ctxp;

    if (_fq_ctx_init_conway(ctx, p, d, var))
    {
        ctx->is_conway = 1;
	    return;
    }

    ctx->is_conway = 0;

    flint_rand_init(state);

    fmpz_mod_ctx_init(ctxp, p);
    fmpz_mod_poly_init2(poly, d + 1, ctxp);
    fmpz_mod_poly_randtest_sparse_irreducible(poly, state, d + 1, ctxp);

    fq_ctx_init_modulus(ctx, poly, ctxp, var);

    fmpz_mod_poly_clear(poly, ctxp);
    fmpz_mod_ctx_clear(ctxp);
    flint_rand_clear(state);
}

/* FIXME: Test super big primes? */
void
fq_ctx_init_randtest(fq_ctx_t ctx, flint_rand_t state, int type)
{
    fq_nmod_ctx_t nmod_ctx;

    fq_nmod_ctx_init_randtest(nmod_ctx, state, type);

    fq_ctx_init_set_clear_small_fq_nmod_ctx(ctx, nmod_ctx);
}

void
fq_ctx_init_randtest_reducible(fq_ctx_t ctx, flint_rand_t state, int type)
{
    fmpz_mod_ctx_t ctxp;
    fmpz_mod_poly_t mod;
    ulong prime;
    slong deg;

    /* Big prime < 2^30, big degree <= 30 */
    /* Small prime < 2^10, small degree <= 15 */
    switch (type)
    {
        case 0:
            prime = n_randprime(state, 2 + n_randint(state, 29), 1);
            deg = 1 + n_randint(state, 30);
            break;
        case 1:
            prime = n_randprime(state, 2 + n_randint(state, 29), 1);
            deg = 1 + n_randint(state, 15);
            break;
        case 2:
            prime = n_randprime(state, 2 + n_randint(state, 9), 1);
            deg = 1 + n_randint(state, 30);
            break;
        case 3:
            prime = n_randprime(state, 2 + n_randint(state, 9), 1);
            deg = 1 + n_randint(state, 15);
            break;
        default: FLINT_UNREACHABLE;
    }

    fmpz_mod_ctx_init_ui(ctxp, prime);
    fmpz_mod_poly_init(mod, ctxp);
    fmpz_mod_poly_randtest_monic(mod, state, deg + 1, ctxp);
    fq_ctx_init_modulus(ctx, mod, ctxp, "a");

    fmpz_mod_poly_clear(mod, ctxp);
    fmpz_mod_ctx_clear(ctxp);
}
