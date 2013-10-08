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

    Copyright (C) 2011 Sebastian Pancratz, 2012 Andres Goens
 
******************************************************************************/

#ifndef FQ_H
#define FQ_H

#undef ulong                /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long

#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

/* Data types and context ****************************************************/

typedef fmpz_poly_t fq_t;
typedef fmpz_poly_struct fq_struct;

typedef struct
{
    fmpz p;

    fmpz *a;
    long *j;
    long len;

    char *var;
}
fq_ctx_struct;

typedef fq_ctx_struct fq_ctx_t[1];

void fq_ctx_init_conway(fq_ctx_t ctx,
                        const fmpz_t p, long d, const char *var);

void fq_ctx_clear(fq_ctx_t ctx);

static __inline__ long fq_ctx_degree(const fq_ctx_t ctx)
{
    return ctx->j[ctx->len - 1];
}

#define fq_ctx_prime(ctx)  (&((ctx)->p))

static __inline__ void fq_ctx_print(const fq_ctx_t ctx)
{
    long i, k;

    printf("p = "), fmpz_print(fq_ctx_prime(ctx)), printf("\n");
    printf("d = %ld\n", ctx->j[ctx->len - 1]);
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

/* Memory managment  *********************************************************/

static __inline__ void fq_init(fq_t rop)
{
    fmpz_poly_init(rop);
}

static __inline__ void fq_init2(fq_t rop, const fq_ctx_t ctx)
{
    fmpz_poly_init2(rop, fq_ctx_degree(ctx));
}

static __inline__ void fq_clear(fq_t rop)
{
    fmpz_poly_clear(rop);
}

static __inline__ 
void _fq_reduce(fmpz *R, long lenR, 
                const fmpz *a, const long *j, long len, const fmpz_t p)
{
    const long d = j[len - 1];

    FMPZ_VEC_NORM(R, lenR);

    if (lenR > d)
    {
        long i, k;

        for (i = lenR - 1; i >= d; i--)
        {
            for (k = len - 2; k >= 0; k--)
            {
                fmpz_submul(R + j[k] + i - d, R + i, a + k);
            }
            fmpz_zero(R + i);
        }
    }

    _fmpz_vec_scalar_mod_fmpz(R, R, FLINT_MIN(d, lenR), p);
}

static __inline__ void fq_reduce(fq_t rop, const fq_ctx_t ctx)
{
    _fq_reduce(rop->coeffs, rop->length, 
               ctx->a, ctx->j, ctx->len, fq_ctx_prime(ctx));
    _fmpz_poly_normalise(rop);
}

/* Basic arithmetic **********************************************************/

void fq_add(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx);

void fq_sub(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx);

void fq_neg(fq_t rop, const fq_t op1, const fq_ctx_t ctx);

void fq_mul(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx);

void fq_mul_fmpz(fq_t rop, const fq_t op, const fmpz_t x, const fq_ctx_t ctx);

void fq_mul_si(fq_t rop, const fq_t op, long x, const fq_ctx_t ctx);

void fq_mul_ui(fq_t rop, const fq_t op, ulong x, const fq_ctx_t ctx);

void fq_sqr(fq_t rop, const fq_t op, const fq_ctx_t ctx);

void fq_inv(fq_t rop, const fq_t op1, const fq_ctx_t ctx);

void _fq_pow(fmpz *rop, const fmpz *op, long len, const fmpz_t e, 
             const fmpz *a, const long *j, long lena, const fmpz_t p);

void fq_pow(fq_t rop, const fq_t op1, const fmpz_t e, const fq_ctx_t ctx);

/* Randomisation *************************************************************/

void fq_randtest(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);

void fq_randtest_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);

/* Comparison ****************************************************************/

static __inline__ int fq_equal(const fq_t op1, const fq_t op2)
{
    return fmpz_poly_equal(op1, op2);
}

static __inline__ int fq_is_zero(const fq_t op)
{
    return fmpz_poly_is_zero(op);
}

static __inline__ int fq_is_one(const fq_t op)
{
    return fmpz_poly_is_one(op);
}

/* Assignments and conversions ***********************************************/

static __inline__ void fq_set(fq_t rop, const fq_t op)
{
    fmpz_poly_set(rop, op);
}

static __inline__ void fq_swap(fq_t op1, fq_t op2)
{
    fmpz_poly_swap(op1, op2);
}

static __inline__ void fq_zero(fq_t rop)
{
    fmpz_poly_zero(rop);
}

static __inline__ void fq_one(fq_t rop)
{
    fmpz_poly_one(rop);
}

/* Output ********************************************************************/

static __inline__ 
int fq_fprint_pretty(FILE * file, const fq_t op, const fq_ctx_t ctx)
{
    return fmpz_poly_fprint_pretty(file, op, ctx->var);
}

static __inline__ 
int fq_print_pretty(const fq_t op, const fq_ctx_t ctx)
{
    return fmpz_poly_print_pretty(op, ctx->var);
}

/* Special functions *********************************************************/

void _fq_trace(fmpz_t rop, const fmpz *op, long len, 
               const fmpz *a, const long *j, long lena, const fmpz_t p);

void fq_trace(fmpz_t rop, const fq_t op, const fq_ctx_t ctx);

void _fq_frobenius(fmpz *rop, const fmpz *op, long len, long e, 
                   const fmpz *a, const long *j, long lena, 
                   const fmpz_t p);

void fq_frobenius(fq_t rop, const fq_t op, long e, const fq_ctx_t ctx);

void _fq_norm(fmpz_t rop, const fmpz *op, long len, 
                          const fmpz *a, const long *j, long lena, 
                          const fmpz_t p, long N);

void fq_norm(fmpz_t rop, const fq_t op, const fq_ctx_t ctx);

#endif

