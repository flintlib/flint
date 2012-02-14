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

    Copyright (C) 2011 Sebastian Pancratz
 
******************************************************************************/

#ifndef QADIC_H
#define QADIC_H

#undef ulong /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long

#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"
#include "padic.h"
#include "padic_poly.h"

/* Data types and context ****************************************************/

#define qadic_t padic_poly_t

#define qadic_struct padic_poly_struct

typedef struct
{
    padic_ctx_struct pctx;

    fmpz *a;
    long *j;
    long len;

    char *var;
}
qadic_ctx_struct;

typedef qadic_ctx_struct qadic_ctx_t[1];

void qadic_ctx_init_conway(qadic_ctx_t ctx, 
                           const fmpz_t p, long d, long N, const char *var, 
                           enum padic_print_mode mode);

void qadic_ctx_clear(qadic_ctx_t ctx);

static __inline__ long qadic_ctx_degree(const qadic_ctx_t ctx)
{
    return ctx->j[ctx->len - 1];
}

/* Memory management *********************************************************/

static __inline__ void qadic_init(qadic_t x)
{
    padic_poly_init(x);
}

static __inline__ void qadic_clear(qadic_t x)
{
    padic_poly_clear(x);
}

static __inline__ void
_fmpz_mod_poly_reduce(fmpz *R, long lenR, 
                      const fmpz *a, const long *j, long len, const fmpz_t p)
{
    const long d = j[len - 1];
    long i, k;

    FMPZ_VEC_NORM(R, lenR);

    _fmpz_vec_scalar_mod_fmpz(R, R, lenR, p);

    for (i = lenR - 1; i >= d; i--)
    {
        for (k = len - 2; k >= 0; k--)
        {
            const long t = j[k] + i - d;

            fmpz_submul(R + t, R + i, a + k);
            fmpz_mod(R + t, R + t, p);
        }
        fmpz_zero(R + i);
    }
}

static __inline__ void qadic_reduce(qadic_t x, const qadic_ctx_t ctx)
{
    if (x->val >= (&ctx->pctx)->N)
    {
        padic_poly_zero(x);
    }
    else
    {
        const long d = ctx->j[ctx->len - 1];

        if (x->length > d)
        {
            fmpz_t pow;
            int alloc;

            alloc = _padic_ctx_pow_ui(pow, (&ctx->pctx)->N - x->val, &ctx->pctx);

            _fmpz_mod_poly_reduce(x->coeffs, x->length, ctx->a, ctx->j, ctx->len, pow);
            _padic_poly_set_length(x, d);
            _padic_poly_normalise(x);
            padic_poly_canonicalise(x, (&ctx->pctx)->p);

            if (alloc)
                fmpz_clear(pow);
        }
    }
}

void qadic_scalar_mod_ppow(qadic_t rop, const qadic_t op, long e, 
                           const qadic_ctx_t ctx);

/* Randomisation *************************************************************/

static __inline__ void 
qadic_randtest(qadic_t x, flint_rand_t state, const qadic_ctx_t ctx)
{
    padic_poly_randtest(x, state, qadic_ctx_degree(ctx), &ctx->pctx);
}

static __inline__ void 
qadic_randtest_not_zero(qadic_t x, flint_rand_t state, const qadic_ctx_t ctx)
{
    padic_poly_randtest_not_zero(x, state, qadic_ctx_degree(ctx), 
                                 &ctx->pctx);
}

static __inline__ void 
qadic_randtest_val(qadic_t x, flint_rand_t state, long val, 
                   const qadic_ctx_t ctx)
{
    padic_poly_randtest_val(x, state, val, qadic_ctx_degree(ctx), &ctx->pctx);
}

static __inline__ void 
qadic_randtest_int(qadic_t x, flint_rand_t state, const qadic_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;

    if (N <= 0)
    {
        padic_poly_zero(x);
    }
    else
    {
        padic_poly_randtest_val(x, state, n_randint(state, N), 
                                qadic_ctx_degree(ctx), &ctx->pctx);
    }
}

/* Assignments and conversions ***********************************************/

static __inline__ void qadic_set(qadic_t x, const qadic_t y)
{
    padic_poly_set(x, y);
}

static __inline__ void qadic_zero(qadic_t x)
{
    padic_poly_zero(x);
}

static __inline__ void qadic_one(qadic_t x, const qadic_ctx_t ctx)
{
    padic_poly_one(x, &ctx->pctx);
}

/* Comparison ****************************************************************/

static __inline__ int 
qadic_is_zero(const qadic_t x)
{
    return padic_poly_is_zero(x);
}

static __inline__ int 
qadic_is_one(const qadic_t x, const qadic_ctx_t ctx)
{
    return padic_poly_is_one(x, &ctx->pctx);
}

static __inline__ int 
qadic_equal(const qadic_t x, const qadic_t y)
{
    return padic_poly_equal(x, y);
}

/* Basic arithmetic **********************************************************/

static __inline__ void 
qadic_add(qadic_t x, const qadic_t y, const qadic_t z, const qadic_ctx_t ctx)
{
    padic_poly_add(x, y, z, &ctx->pctx);
}

static __inline__ void 
qadic_sub(qadic_t x, const qadic_t y, const qadic_t z, const qadic_ctx_t ctx)
{
    padic_poly_sub(x, y, z, &ctx->pctx);
}

static __inline__ void 
qadic_neg(qadic_t x, const qadic_t y, const qadic_ctx_t ctx)
{
    padic_poly_neg(x, y, &ctx->pctx);
}

static __inline__ void 
qadic_mul(qadic_t x, const qadic_t y, const qadic_t z, const qadic_ctx_t ctx)
{
    padic_poly_mul(x, y, z, &ctx->pctx);
    qadic_reduce(x, ctx);
}

void _qadic_inv(fmpz *rop, const fmpz *op, long len, 
                const fmpz *a, const long *j, long lena, 
                const fmpz_t p, long N);

void qadic_inv(qadic_t x, const qadic_t y, const qadic_ctx_t ctx);

void _qadic_pow(fmpz *rop, const fmpz *op, long len, const fmpz_t e, 
                const fmpz *a, const long *j, long lena, 
                const fmpz_t p);

void qadic_pow(qadic_t x, const qadic_t y, const fmpz_t e, const qadic_ctx_t ctx);

/* Special functions *********************************************************/

void qadic_sigma(qadic_t rop, const qadic_t op, long e, const qadic_ctx_t ctx);

void _qadic_teichmuller(fmpz *rop, const fmpz *op, long len, 
                        const fmpz *a, const long *j, long lena, 
                        const fmpz_t p, long N);

void qadic_teichmuller(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

/* Output ********************************************************************/

int qadic_fprint_pretty(FILE *file, const qadic_t op, const qadic_ctx_t ctx);

static __inline__ int 
qadic_print_pretty(const qadic_t op, const qadic_ctx_t ctx)
{
    return qadic_fprint_pretty(stdout, op, ctx);
}

#endif

