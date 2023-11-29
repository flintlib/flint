/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "nmod_poly.h"
#include "fq_nmod.h"
#include "fq_zech.h"

void
fq_zech_ctx_init(fq_zech_ctx_t ctx, const fmpz_t p, slong d, const char *var)
{
    if (!_fq_zech_ctx_init_conway(ctx, p, d, var))
        fq_zech_ctx_init_random(ctx, p, d, var);
}

void
fq_zech_ctx_init_conway(fq_zech_ctx_t ctx, const fmpz_t p, slong d,
                        const char *var)
{
    fq_nmod_ctx_struct * fq_nmod_ctx;

    fq_nmod_ctx = flint_malloc(sizeof(fq_nmod_ctx_struct));

    fq_nmod_ctx_init_conway(fq_nmod_ctx, p, d, var);
    fq_zech_ctx_init_fq_nmod_ctx(ctx, fq_nmod_ctx);
    ctx->owns_fq_nmod_ctx = 1;
    ctx->is_conway = 1;
}

int
_fq_zech_ctx_init_conway(fq_zech_ctx_t ctx, const fmpz_t p, slong d,
                         const char *var)
{
    int result;
    fq_nmod_ctx_struct * fq_nmod_ctx;

    fq_nmod_ctx = flint_malloc(sizeof(fq_nmod_ctx_struct));

    result = _fq_nmod_ctx_init_conway(fq_nmod_ctx, p, d, var);
    if (!result)
    {
        flint_free(fq_nmod_ctx);
	ctx->is_conway = 0;
        return result;
    } else
        ctx->is_conway = 1;

    fq_zech_ctx_init_fq_nmod_ctx(ctx, fq_nmod_ctx);
    ctx->owns_fq_nmod_ctx = 1;
    return result;
}

void
fq_zech_ctx_init_random(fq_zech_ctx_t ctx, const fmpz_t p, slong d,
                        const char *var)
{
    fq_nmod_ctx_struct * fq_nmod_ctx;
    flint_rand_t state;
    nmod_poly_t poly;

    fq_nmod_ctx = flint_malloc(sizeof(fq_nmod_ctx_struct));

    flint_randinit(state);

    nmod_poly_init2(poly, fmpz_get_ui(p), d + 1);
    nmod_poly_randtest_monic_primitive(poly, state, d + 1);

    fq_nmod_ctx_init_modulus(fq_nmod_ctx, poly, var);

    nmod_poly_clear(poly);
    flint_randclear(state);

    fq_zech_ctx_init_fq_nmod_ctx(ctx, fq_nmod_ctx);
    ctx->owns_fq_nmod_ctx = 1;
    ctx->is_conway = 0;
}

int
fq_zech_ctx_init_modulus_check(fq_zech_ctx_t ctx,
                         const nmod_poly_t modulus,
                         const char *var)
{
    int primitive;
    fq_nmod_ctx_struct * fq_nmod_ctx;

    fq_nmod_ctx = flint_malloc(sizeof(fq_nmod_ctx_struct));

    fq_nmod_ctx_init_modulus(fq_nmod_ctx, modulus, var);
    primitive = fq_zech_ctx_init_fq_nmod_ctx_check(ctx, fq_nmod_ctx);
    ctx->owns_fq_nmod_ctx = 1;
    return primitive;
}

void
fq_zech_ctx_init_modulus(fq_zech_ctx_t ctx,
		   const nmod_poly_t modulus,
		   const char *var)
{
   fq_zech_ctx_init_modulus_check(ctx, modulus, var);
}

int
fq_zech_ctx_init_fq_nmod_ctx_check(fq_zech_ctx_t ctx,
                             fq_nmod_ctx_t fq_nmod_ctx)
{
    ulong i, n;
    fq_nmod_t r, gen;
    slong up, q;
    fmpz_t result, order;
    mp_limb_t j, nz, result_ui;
    mp_limb_t *n_reverse_table;

    ctx->fq_nmod_ctx = fq_nmod_ctx;
    ctx->owns_fq_nmod_ctx = 0;

    fmpz_init(order);
    fq_nmod_ctx_order(order, fq_nmod_ctx);

    if (fmpz_bits(order) > FLINT_BITS)
    {
        flint_throw(FLINT_ERROR, "Exception (fq_zech_ctx_init_fq_nmod_ctx). Requires q < 2^FLINT_BITS\n");
    }

    q = fmpz_get_ui(order);
    up = fmpz_get_ui(fq_nmod_ctx_prime(fq_nmod_ctx));

    ctx->p = up;
    ctx->ppre = n_precompute_inverse(ctx->p);
    ctx->qm1 = q - 1;

    if (up == 2)
    {
        ctx->qm1o2 = 0;
    }
    else
    {
        ctx->qm1o2 = ctx->qm1 / 2;
    }

    ctx->qm1opm1 = ctx->qm1 / (up - 1);

/* 1. The field may not be defined with a Conway polynomial
 * 2. need to ensure prime_root is the norm of the generator
 * 3. so we take prime_root = (-1)^d * a_0, where d is the degree
 *    of the minimum polynomial P of the generator, and a_0 is the constant term of
 *    the generator.
 * 4. this is because if P(t) = (t-x_0)...(t-x_{d-1}), then the constant term of
 * P is the product of the x_i (ie the norm) and is equal to (-1)^d * a_0
 */
    ctx->prime_root = (fq_nmod_ctx_degree(fq_nmod_ctx) & 1) ?
        ctx->p - fq_nmod_ctx->a[0] : fq_nmod_ctx->a[0];

    ctx->zech_log_table = (mp_limb_t *) flint_malloc(q * sizeof(mp_limb_t));
    ctx->prime_field_table = (mp_limb_t *) flint_malloc(up * sizeof(mp_limb_t));
    n_reverse_table = (mp_limb_t *) flint_malloc(q * sizeof(mp_limb_t));
    ctx->eval_table = (mp_limb_t *) flint_malloc(q * sizeof(mp_limb_t));

    ctx->zech_log_table[ctx->qm1] = 0;
    ctx->prime_field_table[0] = ctx->qm1;
    for (i = 0; i < q; i++)
        n_reverse_table[i] = ctx->qm1;
    ctx->eval_table[ctx->qm1] = 0;

    fq_nmod_init(r, ctx->fq_nmod_ctx);
    fq_nmod_init(gen, ctx->fq_nmod_ctx);
    fq_nmod_one(r, ctx->fq_nmod_ctx);
    fq_nmod_gen(gen, ctx->fq_nmod_ctx);

    fmpz_init(result);

    for (i = 0; i < ctx->qm1; i++)
    {
        nmod_poly_evaluate_fmpz(result, r, fq_nmod_ctx_prime(fq_nmod_ctx));
        result_ui = fmpz_get_ui(result);
        if (n_reverse_table[result_ui] != ctx->qm1)
        {   /* clean up... */
            fq_nmod_clear(r, fq_nmod_ctx);
            fq_nmod_clear(gen, fq_nmod_ctx);
            flint_free(n_reverse_table);
            fmpz_clear(result);
            fmpz_clear(order);
            fq_zech_ctx_clear(ctx);
            return 0; /* failure: modulus not primitive */
        }
        n_reverse_table[result_ui] = i;
        ctx->eval_table[i] = result_ui;
        if (r->length == 1)
        {
            ctx->prime_field_table[result_ui] = i;
        }
        fq_nmod_mul(r, r, gen, fq_nmod_ctx);
    }

    i = 1;
    for (i = 0; i < q; i++)
    {
        j = n_reverse_table[i];
        n = i;
        if (n % up == up - 1)
        {
            nz = n - up + 1;
        }
        else
        {
            nz = n + 1;
        }
        ctx->zech_log_table[j] = n_reverse_table[nz];
    }

    fq_nmod_clear(r, fq_nmod_ctx);
    fq_nmod_clear(gen, fq_nmod_ctx);
    flint_free(n_reverse_table);
    fmpz_clear(result);
    fmpz_clear(order);

    return 1; /* success */
}

void
fq_zech_ctx_init_fq_nmod_ctx(fq_zech_ctx_t ctx, fq_nmod_ctx_t fq_nmod_ctx)
{
    if (!fq_zech_ctx_init_fq_nmod_ctx_check(ctx, fq_nmod_ctx))
    {
        flint_throw(FLINT_ERROR, "(fq_zech_ctx_init_fq_nmod_ctx): Polynomial is not primitive.\n");
    }
}
