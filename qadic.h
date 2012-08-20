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

    Copyright (C) 2011, 2012 Sebastian Pancratz
 
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

#define qadic_val(op) ((op)->val)

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

static __inline__ void qadic_ctx_print(const qadic_ctx_t ctx)
{
    long i, k;

    printf("p    = "), fmpz_print((&ctx->pctx)->p), printf("\n");
    printf("d    = %ld\n", ctx->j[ctx->len - 1]);
    printf("N    = %ld\n", (&ctx->pctx)->N);
    printf("f(X) = ");
    fmpz_print(ctx->a + 0);
    for (k = 1; k < ctx->len; k++)
    {
        i = ctx->j[k];
        printf(" + ");
        if (fmpz_is_one(ctx->a + k))
        {
            if (i == 1)
                printf("X");
            else
                printf("X^%ld", i);
        }
        else
        {
            fmpz_print(ctx->a + k);
            if (i == 1)
                printf("*X");
            else
                printf("*X^%ld", i);
        }
    }
    printf("\n");
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

/*
    TODO:  Consider renaming this function, prefix for the "qadic" module.
 */

static __inline__ void
_fmpz_poly_reduce(fmpz *R, long lenR, const fmpz *a, const long *j, long len)
{
    const long d = j[len - 1];
    long i, k;

    FMPZ_VEC_NORM(R, lenR);

    for (i = lenR - 1; i >= d; i--)
    {
        for (k = len - 2; k >= 0; k--)
        {
            fmpz_submul(R + j[k] + i - d, R + i, a + k);
        }
        fmpz_zero(R + i);
    }
}

/*
    TODO:  Consider renaming this function, prefix for the "qadic" module.
 */

static __inline__ void 
_fmpz_mod_poly_reduce(fmpz *R, long lenR, 
                      const fmpz *a, const long *j, long len, const fmpz_t p)
{
    const long d = j[len - 1];

    if (lenR > d)
    {
        _fmpz_poly_reduce(R, lenR, a, j, len);
        _fmpz_vec_scalar_mod_fmpz(R, R, d, p);
    }
    else
    {
        _fmpz_vec_scalar_mod_fmpz(R, R, lenR, p);
    }
}

static __inline__ void qadic_reduce(qadic_t x, const qadic_ctx_t ctx)
{
    const long N = (&ctx->pctx)->N;
    const long d = ctx->j[ctx->len - 1];

    if (x->length == 0 || x->val >= N)
    {
        padic_poly_zero(x);
    }
    else
    {
        fmpz_t pow;
        int alloc;

        alloc = _padic_ctx_pow_ui(pow, (&ctx->pctx)->N - x->val, &ctx->pctx);

        _fmpz_mod_poly_reduce(x->coeffs, x->length, ctx->a, ctx->j, ctx->len, pow);
        _padic_poly_set_length(x, FLINT_MIN(x->length, d));
        _padic_poly_normalise(x);
        padic_poly_canonicalise(x, (&ctx->pctx)->p);

        if (alloc)
            fmpz_clear(pow);
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

static __inline__ void qadic_zero(qadic_t x)
{
    padic_poly_zero(x);
}

static __inline__ void qadic_one(qadic_t x, const qadic_ctx_t ctx)
{
    padic_poly_one(x, &ctx->pctx);
}

static __inline__ void qadic_gen(qadic_t x, const qadic_ctx_t ctx)
{
    const long d = qadic_ctx_degree(ctx);

    if (d > 1)
    {
        if ((&ctx->pctx)->N > 0)
        {
            padic_poly_fit_length(x, 2);
            fmpz_zero(x->coeffs + 0);
            fmpz_one(x->coeffs + 1);
            _padic_poly_set_length(x, 2);
            x->val = 0;
        }
        else
        {
            padic_poly_zero(x);
        }
    }
    else
    {
        printf("Exception (qadic_gen).  Extension degree d = 1.\n");
        abort();
    }
}

static __inline__ 
void qadic_set_ui(qadic_t rop, ulong op, const qadic_ctx_t ctx)
{
    padic_poly_set_ui(rop, op, &ctx->pctx);
}

static __inline__ int 
qadic_get_padic(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    if (op->length > 0)
    {
        if (_fmpz_vec_is_zero(op->coeffs + 1, op->length - 1))
        {
            fmpz_set(padic_unit(rop), op->coeffs + 0);
            padic_val(rop) = op->val;
            _padic_canonicalise(rop, &ctx->pctx);
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        padic_zero(rop);
        return 1;
    }
}

static __inline__ void qadic_set(qadic_t x, const qadic_t y)
{
    padic_poly_set(x, y);
}

void qadic_set_fmpz_poly(qadic_t rop, const fmpz_poly_t op, 
                         const qadic_ctx_t ctx);

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

void qadic_mul(qadic_t x, const qadic_t y, const qadic_t z,
                          const qadic_ctx_t ctx);

void _qadic_inv(fmpz *rop, const fmpz *op, long len, 
                const fmpz *a, const long *j, long lena, 
                const fmpz_t p, long N);

void qadic_inv(qadic_t x, const qadic_t y, const qadic_ctx_t ctx);

void _qadic_pow(fmpz *rop, const fmpz *op, long len, const fmpz_t e, 
                const fmpz *a, const long *j, long lena, 
                const fmpz_t p);

void qadic_pow(qadic_t x, const qadic_t y, const fmpz_t e, const qadic_ctx_t ctx);

/* Special functions *********************************************************/

void _qadic_exp_rectangular(fmpz *rop, const fmpz *op, long v, long len, 
                            const fmpz *a, const long *j, long lena, 
                            const fmpz_t p, long N, const fmpz_t pN);

int qadic_exp_rectangular(qadic_t rop, const qadic_t op, 
                          const qadic_ctx_t ctx);

void _qadic_exp_balanced(fmpz *rop, const fmpz *op, long v, long len, 
                         const fmpz *a, const long *j, long lena, 
                         const fmpz_t p, long N, const fmpz_t pN);

int qadic_exp_balanced(qadic_t rop, const qadic_t op, 
                       const qadic_ctx_t ctx);

void _qadic_exp(fmpz *rop, const fmpz *op, long v, long len, 
                           const fmpz *a, const long *j, long lena, 
                           const fmpz_t p, long N, const fmpz_t pN);

int qadic_exp(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_log_rectangular(fmpz *z, const fmpz *y, long v, long len, 
                            const fmpz *a, const long *j, long lena, 
                            const fmpz_t p, long N, const fmpz_t pN);

int qadic_log_rectangular(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_log_balanced(fmpz *z, const fmpz *y, long len, 
                         const fmpz *a, const long *j, long lena, 
                         const fmpz_t p, long N, const fmpz_t pN);

int qadic_log_balanced(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_log(fmpz *z, const fmpz *y, long v, long len, 
                const fmpz *a, const long *j, long lena, 
                const fmpz_t p, long N, const fmpz_t pN);

int qadic_log(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void qadic_frobenius(qadic_t rop, const qadic_t op, long e, const qadic_ctx_t ctx);

void _qadic_teichmuller(fmpz *rop, const fmpz *op, long len, 
                        const fmpz *a, const long *j, long lena, 
                        const fmpz_t p, long N);

void qadic_teichmuller(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_trace(fmpz_t rop, const fmpz *op, long len, 
                  const fmpz *a, const long *j, long lena, const fmpz_t pN);

void qadic_trace(padic_t rop, const qadic_t op, const qadic_ctx_t ctx);

void _qadic_norm(fmpz_t rop, const fmpz *op, long len, 
                 const fmpz *a, const long *j, long lena, 
                 const fmpz_t p, long N);

void qadic_norm(padic_t rop, const qadic_t op, const qadic_ctx_t ctx);
void qadic_norm_analytic(padic_t rop, const qadic_t op, const qadic_ctx_t ctx);
void qadic_norm_resultant(padic_t rop, const qadic_t op, const qadic_ctx_t ctx);

int qadic_sqrt(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx);

/* Output ********************************************************************/

int qadic_fprint_pretty(FILE *file, const qadic_t op, const qadic_ctx_t ctx);

static __inline__ int 
qadic_print_pretty(const qadic_t op, const qadic_ctx_t ctx)
{
    return qadic_fprint_pretty(stdout, op, ctx);
}

#endif

