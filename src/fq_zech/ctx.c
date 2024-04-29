/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "fmpz.h"
#include "fq_nmod.h"
#include "fq_zech.h"

static ulong
_nmod_poly_evaluate_ui(mp_srcptr poly, slong len, ulong xd)
{
    slong ix;
    ulong res;

    if (len == 0)
        return 0;

    res = poly[len - 1];

    for (ix = len - 2; ix >= 0; ix--)
    {
        res *= xd;
        res += poly[ix];
    }

    return res;
}

int
_fq_zech_ctx_init_conway_ui(fq_zech_ctx_t ctx, ulong p, slong d, const char * var)
{
    int ret;
    ulong conway_poly[410]; /* Largest degree in database is 409 */
    fq_nmod_ctx_struct * ctx_fq_nmod;
    nmod_poly_struct tmp;

    ret = _nmod_poly_conway(conway_poly, p, d);

    if (ret != 1 && ret != 0)
        FLINT_UNREACHABLE;

    if (!ret)
        return 0; /* We do not change ctx->is_conway */

    nmod_poly_init(&tmp, p);
    tmp.coeffs = conway_poly;
    tmp.length = d + 1;

    ctx_fq_nmod = flint_malloc(sizeof(fq_nmod_ctx_struct));
    fq_nmod_ctx_init_modulus(ctx_fq_nmod, &tmp, var);

    ctx->is_conway = 1;
    fq_zech_ctx_init_fq_nmod_ctx(ctx, ctx_fq_nmod);
    ctx->owns_fq_nmod_ctx = 1;

    /* No need to clear tmp */

    return 1;
}

void
fq_zech_ctx_init_conway_ui(fq_zech_ctx_t ctx, ulong p, slong d, const char * var)
{
    if (!_fq_zech_ctx_init_conway_ui(ctx, p, d, var))
        flint_throw(FLINT_ERROR,
                "Exception (fq_zech_ctx_init_conway_ui).  "
                "The polynomial for (p, d) = (%wu, %wd) is not present in the "
                "database.\n", p, d);
}

void
fq_zech_ctx_init_ui(fq_zech_ctx_t ctx, ulong p, slong d, const char *var)
{
    if (!_fq_zech_ctx_init_conway_ui(ctx, p, d, var))
        fq_zech_ctx_init_random_ui(ctx, p, d, var);
}

void
fq_zech_ctx_init_random_ui(fq_zech_ctx_t ctx, ulong p, slong d, const char * var)
{
    fq_nmod_ctx_struct * fq_nmod_ctx;
    flint_rand_t state;
    nmod_poly_t poly;
    ulong tmp_coeffs[FLINT_BITS]; /* p^d < 2^{FLINT_BITS} */

    fq_nmod_ctx = flint_malloc(sizeof(fq_nmod_ctx_struct));

    flint_randinit(state);
    nmod_poly_init(poly, p);

    poly->coeffs = tmp_coeffs;
    poly->alloc = FLINT_BITS; /* Ensure no reallocation in nmod_poly_randtest */
    poly->length = d + 1;

    nmod_poly_randtest_monic_primitive(poly, state, d + 1);

    fq_nmod_ctx_init_modulus(fq_nmod_ctx, poly, var);

    /* No need to clear poly */
    /* No need to clear state */

    fq_zech_ctx_init_fq_nmod_ctx(ctx, fq_nmod_ctx);

    ctx->owns_fq_nmod_ctx = 1;
    ctx->is_conway = 0;
}

int
fq_zech_ctx_init_modulus_check(fq_zech_ctx_t ctx, const nmod_poly_t modulus, const char * var)
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
fq_zech_ctx_init_modulus(fq_zech_ctx_t ctx, const nmod_poly_t modulus, const char * var)
{
   fq_zech_ctx_init_modulus_check(ctx, modulus, var);
}

int
fq_zech_ctx_init_fq_nmod_ctx_check(fq_zech_ctx_t ctx, fq_nmod_ctx_t ctx2)
{
    ulong i, n, q, up;
    fq_nmod_t r, gen;
    ulong result;
    ulong j, nz;
    mp_ptr n_reverse_table;

    ctx->fq_nmod_ctx = ctx2;
    ctx->owns_fq_nmod_ctx = 0;

    up = fq_nmod_ctx_prime(ctx2);

    q = _n_pow_check(up, fq_nmod_ctx_degree(ctx2));

    if (!q)
        flint_throw(FLINT_ERROR, "Exception (fq_zech_ctx_init_fq_nmod_ctx). "
                "Requires q < 2^FLINT_BITS\n");

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
    ctx->prime_root = (fq_nmod_ctx_degree(ctx2) & 1) ? up - ctx2->a[0] : ctx2->a[0];

    /* NOTE: Only malloc once for tables */
    ctx->zech_log_table = flint_malloc(sizeof(ulong) * (2 * q + up));
    ctx->prime_field_table = ctx->zech_log_table + q;
    ctx->eval_table = ctx->prime_field_table + up;

    n_reverse_table = flint_malloc(sizeof(ulong) * q);

    ctx->zech_log_table[ctx->qm1] = 0;
    ctx->prime_field_table[0] = ctx->qm1;
    for (i = 0; i < q; i++)
        n_reverse_table[i] = ctx->qm1;
    ctx->eval_table[ctx->qm1] = 0;

    fq_nmod_init(r, ctx->fq_nmod_ctx);
    fq_nmod_init(gen, ctx->fq_nmod_ctx);
    fq_nmod_one(r, ctx->fq_nmod_ctx);
    fq_nmod_gen(gen, ctx->fq_nmod_ctx);

    for (i = 0; i < ctx->qm1; i++)
    {
        result = _nmod_poly_evaluate_ui(r->coeffs, r->length, up);

        if (n_reverse_table[result] != ctx->qm1)
        {
            fq_nmod_clear(r, ctx2);
            fq_nmod_clear(gen, ctx2);
            flint_free(n_reverse_table);
            fq_zech_ctx_clear(ctx);
            return 0; /* failure: modulus not primitive */
        }

        n_reverse_table[result] = i;
        ctx->eval_table[i] = result;

        if (r->length == 1)
            ctx->prime_field_table[result] = i;

        fq_nmod_mul(r, r, gen, ctx2);
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

    fq_nmod_clear(r, ctx2);
    fq_nmod_clear(gen, ctx2);
    flint_free(n_reverse_table);

    return 1; /* success */
}

void
fq_zech_ctx_init_fq_nmod_ctx(fq_zech_ctx_t ctx, fq_nmod_ctx_t fq_nmod_ctx)
{
    if (!fq_zech_ctx_init_fq_nmod_ctx_check(ctx, fq_nmod_ctx))
        flint_throw(FLINT_ERROR, "(fq_zech_ctx_init_fq_nmod_ctx): "
                "Polynomial is not primitive.\n");
}

void
fq_zech_ctx_init_randtest(fq_zech_ctx_t ctx, flint_rand_t state, int type)
{
    slong deg;
    ulong prime;

    /* Big prime < 2^5, big degree <= 4 */
    /* Small prime < 2^4, small degree <= 3 */
    switch (type)
    {
        case 0:
            prime = n_randprime(state, 2 + n_randint(state, 4), 1);
            deg = 1 + n_randint(state, 4);
            break;
        case 1:
            prime = n_randprime(state, 2 + n_randint(state, 4), 1);
            deg = 1 + n_randint(state, 3);
            break;
        case 2:
            prime = n_randprime(state, 2 + n_randint(state, 3), 1);
            deg = 1 + n_randint(state, 4);
            break;
        case 3:
            prime = n_randprime(state, 2 + n_randint(state, 3), 1);
            deg = 1 + n_randint(state, 3);
            break;
        default: FLINT_UNREACHABLE;
    }

    fq_zech_ctx_init_random_ui(ctx, prime, deg, "a");

    ctx->owns_fq_nmod_ctx = 1;
}

void
fq_zech_ctx_init_randtest_reducible(fq_zech_ctx_t ctx, flint_rand_t state, int type)
{
    fq_zech_ctx_init_randtest(ctx, state, type);
}

void
fq_zech_ctx_clear(fq_zech_ctx_t ctx)
{
    /* NOTE: Only zech_log_table of the three tables was assigned by malloc */
    flint_free(ctx->zech_log_table);

    if (ctx->owns_fq_nmod_ctx)
    {
        fq_nmod_ctx_clear(ctx->fq_nmod_ctx);
        flint_free(ctx->fq_nmod_ctx);
    }
}

void
fq_zech_ctx_order(fmpz_t f, const fq_zech_ctx_t ctx)
{
    fq_nmod_ctx_order(f, ctx->fq_nmod_ctx);
}

/* Deprecated functions ******************************************************/

void fq_zech_ctx_init(fq_zech_ctx_t ctx, fmpz_t p, slong d, const char * var) { fq_zech_ctx_init_ui(ctx, fmpz_get_ui(p), d, var); }
int _fq_zech_ctx_init_conway(fq_zech_ctx_t ctx, fmpz_t p, slong d, const char * var) { return _fq_zech_ctx_init_conway_ui(ctx, fmpz_get_ui(p), d, var); }
void fq_zech_ctx_init_conway(fq_zech_ctx_t ctx, fmpz_t p, slong d, const char * var) { fq_zech_ctx_init_conway_ui(ctx, fmpz_get_ui(p), d, var); }
void fq_zech_ctx_init_random(fq_zech_ctx_t ctx, fmpz_t p, slong d, const char * var) { fq_zech_ctx_init_random_ui(ctx, fmpz_get_ui(p), d, var); }
