/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

/*
    Include routines for vectors over \code{fmpz_poly_struct}, 
    for use in the classical multiplication routine in the 
    $X$-direction.
 */

static fmpz_poly_struct *
__vec_init(slong len)
{
    slong i;
    fmpz_poly_struct *v;

    v = flint_malloc(len * sizeof(fmpz_poly_struct));
    for (i = 0; i < len; i++)
        fmpz_poly_init(v + i);
    return v;
}

static fmpz_poly_struct *
__vec_init2(slong len, slong n)
{
    slong i;
    fmpz_poly_struct *v;

    v = flint_malloc(len * sizeof(fmpz_poly_struct));
    for (i = 0; i < len; i++)
        fmpz_poly_init2(v + i, n);
    return v;
}

static void
__vec_clear(fmpz_poly_struct * v, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
        fmpz_poly_clear(v + i);
    flint_free(v);
}

static void
__scalar_addmul(fmpz_poly_struct * rop,
                const fmpz_poly_struct * op, slong len, const fmpz_poly_t x)
{
    slong i;

    if (fmpz_poly_is_zero(x))
    {
        return;
    }
    else if (fmpz_poly_is_one(x))
    {
        for (i = 0; i < len; i++)
            fmpz_poly_add(rop + i, rop + i, op + i);
    }
    else
    {
        fmpz_poly_t t;

        fmpz_poly_init(t);
        for (i = 0; i < len; i++)
        {
            fmpz_poly_mul(t, op + i, x);
            fmpz_poly_add(rop + i, rop + i, t);
        }
        fmpz_poly_clear(t);
    }
}

static void
__scalar_mul(fmpz_poly_struct * rop,
             const fmpz_poly_struct * op, slong len, const fmpz_poly_t x)
{
    slong i;

    if (fmpz_poly_is_zero(x))
    {
        for (i = 0; i < len; i++)
            fmpz_poly_zero(rop + i);
    }
    else if (fmpz_poly_is_one(x))
    {
        for (i = 0; i < len; i++)
            fmpz_poly_set(rop + i, op + i);
    }
    else
    {
        for (i = 0; i < len; i++)
            fmpz_poly_mul(rop + i, op + i, x);
    }
}

static void
__mul(fmpz_poly_struct * rop,
      fmpz_poly_struct * op1, slong len1, fmpz_poly_struct * op2, slong len2)
{
    if (len1 == 1 && len2 == 1)
    {
        fmpz_poly_mul(rop, op1, op2);
    }
    else
    {
        slong i;

        __scalar_mul(rop, op1, len1, op2);

        __scalar_mul(rop + len1, op2 + 1, len2 - 1, op1 + len1 - 1);

        for (i = 0; i < len1 - 1; i++)
            __scalar_addmul(rop + i + 1, op2 + 1, len2 - 1, op1 + i);
    }
}

void
_TEMPLATE(T, poly_mul_reorder) (TEMPLATE(T, struct) * rop,
                                const TEMPLATE(T, struct) * op1, slong len1,
                                const TEMPLATE(T, struct) * op2, slong len2,
                                const TEMPLATE(T, ctx_t) ctx)
{
    const slong d = TEMPLATE(T, ctx_degree) (ctx);

    fmpz_poly_struct *f, *g, *h;
    slong i, j, k, len;

    f = __vec_init(2 * d - 1);
    g = __vec_init2(d, len1);
    h = __vec_init2(d, len2);

    /* Convert (op1, len1) to (g, d) */
    for (i = 0; i < len1; i++)
        for (j = 0; j < fmpz_poly_length(op1 + i); j++)
            fmpz_set((g + j)->coeffs + i, (op1 + i)->coeffs + j);

    /* Convert (op2, len2) to (h, d) */
    for (i = 0; i < len2; i++)
        for (j = 0; j < fmpz_poly_length(op2 + i); j++)
            fmpz_set((h + j)->coeffs + i, (op2 + i)->coeffs + j);

    for (j = 0; j < d; j++)
    {
        _fmpz_poly_set_length(g + j, len1);
        _fmpz_poly_set_length(h + j, len2);
        _fmpz_poly_normalise(g + j);
        _fmpz_poly_normalise(h + j);
    }

    __mul(f, g, d, h, d);

    /* Normalise (f, len) */
    len = 2 * d - 1;
    while ((len) && fmpz_poly_is_zero(f + (len - 1)))
        len--;

    /* Reduce (f, j) using polynomial operations */
    if (len > d)
    {
        for (i = len - 1; i >= d; i--)
        {
            for (k = ctx->len - 2; k >= 0; k--)
            {
                fmpz_poly_scalar_submul_fmpz(f + ctx->j[k] + i - d, f + i,
                                             ctx->a + k);
            }
            fmpz_poly_zero(f + i);
        }
    }
    for (j = 0; j < FLINT_MIN(d, len); j++)
        fmpz_poly_scalar_mod_fmpz(f + j, f + j, TEMPLATE(T, ctx_prime) (ctx));

    /* Convert (f, d) to (rop, len1 + len2 - 1) */
    for (i = 0; i < len1 + len2 - 1; i++)
    {
        fmpz_poly_fit_length(rop + i, d);
        _fmpz_vec_zero((rop + i)->coeffs, d);
    }
    for (j = 0; j < d; j++)
        for (i = 0; i < fmpz_poly_length(f + j); i++)
            fmpz_set((rop + i)->coeffs + j, (f + j)->coeffs + i);
    for (i = 0; i < len1 + len2 - 1; i++)
    {
        _fmpz_poly_set_length(rop + i, d);
        _fmpz_poly_normalise(rop + i);
    }

    __vec_clear(f, 2 * d - 1);
    __vec_clear(g, d);
    __vec_clear(h, d);
}

void
TEMPLATE(T, poly_mul_reorder) (TEMPLATE(T, poly_t) rop,
                               const TEMPLATE(T, poly_t) op1,
                               const TEMPLATE(T, poly_t) op2,
                               const TEMPLATE(T, ctx_t) ctx)
{
    const slong len = op1->length + op2->length - 1;

    if (op1->length == 0 || op2->length == 0)
    {
        TEMPLATE(T, poly_zero) (rop, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length) (rop, len, ctx);
        _TEMPLATE(T, poly_mul_reorder) (rop->coeffs, op1->coeffs, op1->length,
                                        op2->coeffs, op2->length, ctx);
        _TEMPLATE(T, poly_set_length) (rop, len, ctx);
    }
}


#endif
