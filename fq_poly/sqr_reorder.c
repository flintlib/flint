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

    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "fq_poly.h"

/*
    Include routines for vectors over \code{fmpz_poly_struct}, 
    for use in the classical multiplication routine in the 
    $X$-direction.
 */

static fmpz_poly_struct * __vec_init(long len)
{
    long i;
    fmpz_poly_struct *v;

    v = flint_malloc(len * sizeof(fmpz_poly_struct));
    for (i = 0; i < len; i++)
        fmpz_poly_init(v + i);
    return v;
}

static fmpz_poly_struct * __vec_init2(long len, long n)
{
    long i;
    fmpz_poly_struct *v;

    v = flint_malloc(len * sizeof(fmpz_poly_struct));
    for (i = 0; i < len; i++)
        fmpz_poly_init2(v + i, n);
    return v;
}

static void __vec_clear(fmpz_poly_struct *v, long len)
{
    long i;

    for (i = 0; i < len; i++)
        fmpz_poly_clear(v + i);
    flint_free(v);
}

static void __scalar_addmul(fmpz_poly_struct *rop, 
    const fmpz_poly_struct *op, long len, const fmpz_poly_t x)
{
    long i;

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

static void __scalar_mul(fmpz_poly_struct *rop, 
    const fmpz_poly_struct *op, long len, const fmpz_poly_t x)
{
    long i;

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

static void __sqr(fmpz_poly_struct *rop, fmpz_poly_struct *op, long len)
{
    if (len == 1)
    {
        fmpz_poly_sqr(rop, op);
    }
    else
    {
        long i;
        fmpz_poly_t t;

        fmpz_poly_init(t);

        __scalar_mul(rop, op, len, op);

        __scalar_mul(rop + len, op + 1, len - 1, op + len - 1);

        for (i = 1; i < len - 1; i++)
            __scalar_addmul(rop + i + 1, op + 1, i - 1, op + i);

        for (i = 1; i < 2 * len - 2; i++)
            fmpz_poly_add(rop + i, rop + i, rop + i);

        for (i = 1; i < len - 1; i++)
        {
            fmpz_poly_sqr(t, op + i);
            fmpz_poly_add(rop + 2 * i, rop + 2 * i, t);
        }
        fmpz_poly_clear(t);
    }
}

void _fq_poly_sqr_reorder(fq_struct *rop, 
                          const fq_struct *op, long len, const fq_ctx_t ctx)
{
    const long d = fq_ctx_degree(ctx);

    fmpz_poly_struct *f, *g;
    long i, j, k, lenF;

    f = __vec_init(2 * d - 1);
    g = __vec_init2(d, len);

    /* Convert (op, len) to (g, d) */
    for (i = 0; i < len; i++)
        for (j = 0; j < fmpz_poly_length(op + i); j++)
            fmpz_set( (g + j)->coeffs + i, (op + i)->coeffs + j );

    for (j = 0; j < d; j++)
    {
        _fmpz_poly_set_length(g + j, len);
        _fmpz_poly_normalise(g + j);
    }

    __sqr(f, g, d);

    /* Normalise (f, len) */
    lenF = 2 * d - 1;
    while ((lenF) && fmpz_poly_is_zero(f + (lenF - 1)))
        lenF--;

    /* Reduce (f, j) using polynomial operations */
    if (lenF > d)
    {
        for (i = lenF - 1; i >= d; i--)
        {
            for (k = ctx->len - 2; k >= 0; k--)
            {
                fmpz_poly_scalar_submul_fmpz(f + ctx->j[k] + i - d, f + i, ctx->a + k);
            }
            fmpz_poly_zero(f + i);
        }
    }
    for (j = 0; j < FLINT_MIN(d, lenF); j++)
        fmpz_poly_scalar_mod_fmpz(f + j, f + j, fq_ctx_prime(ctx));

    /* Convert (f, d) to (rop, 2 * len - 1) */
    for (i = 0; i < 2 * len - 1; i++)
    {
        fmpz_poly_fit_length(rop + i, d);
        _fmpz_vec_zero((rop + i)->coeffs, d);
    }
    for (j = 0; j < d; j++)
        for (i = 0; i < fmpz_poly_length(f + j); i++)
            fmpz_set( (rop + i)->coeffs + j, (f + j)->coeffs + i );
    for (i = 0; i < 2 * len - 1; i++)
    {
        _fmpz_poly_set_length(rop + i, d);
        _fmpz_poly_normalise(rop + i);
    }

    __vec_clear(f, 2 * d - 1);
    __vec_clear(g, d);
}

void fq_poly_sqr_reorder(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)
{
    const long len = 2 * op->length - 1;

    if (op->length == 0)
    {
        fq_poly_zero(rop);
    }
    else
    {
        fq_poly_fit_length(rop, len);
        _fq_poly_sqr_reorder(rop->coeffs, op->coeffs, op->length, ctx);
        _fq_poly_set_length(rop, len);
    }
}
