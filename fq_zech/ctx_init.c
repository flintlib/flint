/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Mike Hansen

******************************************************************************/
#include <stdio.h>
#include <string.h>

#include "flint.h"
#include "fq_zech.h"

void
fq_zech_ctx_init(fq_zech_ctx_t ctx, const fmpz_t p, slong d, const char *var)
{
    fq_nmod_ctx_struct * fq_nmod_ctx;

    fq_nmod_ctx = flint_malloc(sizeof(fq_nmod_ctx_struct));

    fq_nmod_ctx_init(fq_nmod_ctx, p, d, var);
    fq_zech_ctx_init_fq_nmod_ctx(ctx, fq_nmod_ctx);
    ctx->owns_fq_nmod_ctx = 1;
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
        return result;
    }
    
    fq_zech_ctx_init_fq_nmod_ctx(ctx, fq_nmod_ctx);
    ctx->owns_fq_nmod_ctx = 1;
    return result;
}


void
fq_zech_ctx_init_modulus(fq_zech_ctx_t ctx,
                         const nmod_poly_t modulus,
                         const char *var)
{
    fq_nmod_ctx_struct * fq_nmod_ctx;

    fq_nmod_ctx = flint_malloc(sizeof(fq_nmod_ctx_struct));

    fq_nmod_ctx_init_modulus(fq_nmod_ctx, modulus, var);
    fq_zech_ctx_init_fq_nmod_ctx(ctx, fq_nmod_ctx);
    ctx->owns_fq_nmod_ctx = 1;
}


void
fq_zech_ctx_init_fq_nmod_ctx(fq_zech_ctx_t ctx,
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
        flint_printf("Exception (fq_zech_ctx_init_nmod_ctx). Requires q < 2^FLINT_BITS\n");
        abort();
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

    ctx->prime_root = n_primitive_root_prime(ctx->p);

    ctx->zech_log_table = (mp_limb_t *) flint_malloc(q * sizeof(mp_limb_t));
    ctx->prime_field_table = (mp_limb_t *) flint_malloc(up * sizeof(mp_limb_t));
    n_reverse_table = (mp_limb_t *) flint_malloc(q * sizeof(mp_limb_t));
    ctx->eval_table = (mp_limb_t *) flint_malloc(q * sizeof(mp_limb_t));

    ctx->zech_log_table[ctx->qm1] = 0;
    ctx->prime_field_table[0] = ctx->qm1;
    n_reverse_table[0] = ctx->qm1;
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
}
