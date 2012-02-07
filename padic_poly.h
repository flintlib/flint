/*============================================================================

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

#ifndef PADIC_POLY_H
#define PADIC_POLY_H

#include <limits.h>

#include <mpir.h>
#include "fmpz.h"
#include "fmpq.h"
#include "padic.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

/*  Type definitions  ********************************************************/

typedef struct
{
    fmpz *coeffs;
    long alloc;
    long length;
    long val;
} padic_poly_struct;

typedef padic_poly_struct padic_poly_t[1];

/*  Helper functions  ********************************************************/

/*
    Returns the minimum $p$-adic valuation of \code{(vec, len)}, 
    assuming this fits into a \code{signed long}.

    If \code{len} is zero, returns $0$.
 */
static __inline__ long _fmpz_vec_ord_p(const fmpz *vec, long len, const fmpz_t p)
{
    if (len == 0)
    {
        return 0;
    }
    else
    {
        fmpz_t t;
        long i, min = LONG_MAX, v;

        fmpz_init(t);
        for (i = 0; (min > 0) && (i < len); i++)
        {
            if (!fmpz_is_zero(vec + i))
            {
                v   = fmpz_remove(t, vec + i, p);
                min = FLINT_MIN(min, v);
            }
        }
        fmpz_clear(t);
        return (min < LONG_MAX) ? min : 0;
    }
}

/*  Memory management  *******************************************************/

void padic_poly_init(padic_poly_t poly);

void padic_poly_init2(padic_poly_t poly, long alloc);

void padic_poly_clear(padic_poly_t poly);

void padic_poly_realloc(padic_poly_t f, long alloc, const fmpz_t p);

void padic_poly_fit_length(padic_poly_t f, long len);

static __inline__
void _padic_poly_set_length(padic_poly_t poly, long len)
{
    if (poly->length > len)
    {
        long i;

        for (i = len; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i); 
    }
    poly->length = len;
}

void _padic_poly_normalise(padic_poly_t f);

void _padic_poly_canonicalise(fmpz *poly, long *v, long len, const fmpz_t p);

void padic_poly_canonicalise(padic_poly_t poly, const fmpz_t p);

void padic_poly_reduce(padic_poly_t f, const padic_ctx_t ctx);

static __inline__ 
void padic_poly_truncate(padic_poly_t poly, long n, const fmpz_t p)
{
    if (poly->length > n)
    {
        long i;

        for (i = n; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i);
        poly->length = n;
        _padic_poly_normalise(poly);
        padic_poly_canonicalise(poly, p);
    }
}

/*  Polynomial parameters  ***************************************************/

static __inline__ long padic_poly_degree(const padic_poly_t poly)
{
    return poly->length - 1;
}

static __inline__ long padic_poly_length(const padic_poly_t poly)
{
    return poly->length;
}

static __inline__ long padic_poly_val(const padic_poly_t poly)
{
    return poly->val;
}

/*  Randomisation  ***********************************************************/

void padic_poly_randtest(padic_poly_t f, flint_rand_t state, 
                         long len, const padic_ctx_t ctx);

void padic_poly_randtest_not_zero(padic_poly_t f, flint_rand_t state,
                                  long len, const padic_ctx_t ctx);

void padic_poly_randtest_val(padic_poly_t f, flint_rand_t state, 
                             long val, long len, const padic_ctx_t ctx);

/*  Assignment and basic manipulation  ***************************************/

void padic_poly_set(padic_poly_t f, const padic_poly_t g);

void padic_poly_set_padic(padic_poly_t poly, const padic_t x);

void padic_poly_set_si(padic_poly_t poly, long x, const padic_ctx_t ctx);

void padic_poly_set_ui(padic_poly_t poly, ulong x, const padic_ctx_t ctx);

void padic_poly_set_fmpz(padic_poly_t poly, const fmpz_t x, 
                         const padic_ctx_t ctx);

void padic_poly_set_fmpq(padic_poly_t poly, const fmpq_t x, 
                         const padic_ctx_t ctx);

void padic_poly_set_fmpz_poly(padic_poly_t rop, const fmpz_poly_t op, 
                              const padic_ctx_t ctx);

void padic_poly_set_fmpq_poly(padic_poly_t rop, 
                              const fmpq_poly_t op, const padic_ctx_t ctx);

int padic_poly_get_fmpz_poly(fmpz_poly_t rop, const padic_poly_t op, 
                             const padic_ctx_t ctx);

void padic_poly_get_fmpq_poly(fmpq_poly_t rop, 
                              const padic_poly_t op, const padic_ctx_t ctx);

static __inline__ void padic_poly_zero(padic_poly_t poly)
{
    _padic_poly_set_length(poly, 0);
    poly->val = 0;
}

static __inline__ void 
padic_poly_one(padic_poly_t poly, const padic_ctx_t ctx)
{
    if (ctx->N > 0)
    {
        padic_poly_fit_length(poly, 1);
        fmpz_set_ui(poly->coeffs, 1);
        _padic_poly_set_length(poly, 1);
        poly->val = 0;
    }
    else
    {
        padic_poly_zero(poly);
    }
}

void padic_poly_swap(padic_poly_t poly1, padic_poly_t poly2);

/*  Getting and setting coefficients  ****************************************/

void padic_poly_get_coeff_padic(padic_t c, const padic_poly_t poly, long n, 
                                           const padic_ctx_t ctx);

void padic_poly_set_coeff_padic(padic_poly_t f, long n, const padic_t c, 
                                                const padic_ctx_t ctx);

/*  Comparison  **************************************************************/

int padic_poly_equal(const padic_poly_t f, const padic_poly_t g);

static __inline__ 
int padic_poly_is_zero(const padic_poly_t poly)
{
    return poly->length == 0;
}

static __inline__ 
int padic_poly_is_one(const padic_poly_t poly, const padic_ctx_t ctx)
{
    if (ctx->N > 0)
    {
        return (poly->length == 1 && fmpz_is_one(poly->coeffs));
    }
    else
    {
        return (poly->length == 0);
    }
}

/*  Addition and subtraction  ************************************************/

void _padic_poly_add(fmpz *rop, long *rval, const fmpz *op1, long val1, long len1, 
                                            const fmpz *op2, long val2, long len2, 
                                            const padic_ctx_t ctx);

void padic_poly_add(padic_poly_t f, 
                    const padic_poly_t g, const padic_poly_t h, 
                    const padic_ctx_t ctx);

void _padic_poly_sub(fmpz *rop, long *rval, const fmpz *op1, long val1, long len1, 
                                            const fmpz *op2, long val2, long len2, 
                                            const padic_ctx_t ctx);

void padic_poly_sub(padic_poly_t f, 
                    const padic_poly_t g, const padic_poly_t h, 
                    const padic_ctx_t ctx);

void padic_poly_neg(padic_poly_t f, const padic_poly_t g, 
                    const padic_ctx_t ctx);

/*  Scalar multiplication and division  **************************************/

void _padic_poly_scalar_mul_padic(fmpz *rop, long *rval, 
                                  const fmpz *op, long val, long len, 
                                  const padic_t c, const padic_ctx_t ctx);

void padic_poly_scalar_mul_padic(padic_poly_t rop, const padic_poly_t op, 
                                 const padic_t c, const padic_ctx_t ctx);

/*  Multiplication  **********************************************************/

void _padic_poly_mul(fmpz *rop, long *rval, 
                     const fmpz *op1, long val1, long len1, 
                     const fmpz *op2, long val2, long len2, 
                     const padic_ctx_t ctx);

void padic_poly_mul(padic_poly_t f, 
                    const padic_poly_t g, const padic_poly_t h, 
                    const padic_ctx_t ctx);

/*  Powering  ****************************************************************/

void _padic_poly_pow(fmpz *rop, long *rval, 
                     const fmpz *op, long val, long len, ulong e,
                     const padic_ctx_t ctx);

void padic_poly_pow(padic_poly_t rop, const padic_poly_t op, ulong e, 
                    const padic_ctx_t ctx);

/*  Series inversion  ********************************************************/

void padic_poly_inv_series(padic_poly_t Qinv, const padic_poly_t Q, long n, 
                           const padic_ctx_t ctx);

/*  Derivative  **************************************************************/

void _padic_poly_derivative(fmpz *rop, long *rval, 
                            const fmpz *op, long val, long len, 
                            const padic_ctx_t ctx);

void padic_poly_derivative(padic_poly_t rop, 
                           const padic_poly_t op, const padic_ctx_t ctx);

/*  Shifting  ****************************************************************/

void padic_poly_shift_left(padic_poly_t rop, const padic_poly_t op, long n);

void padic_poly_shift_right(padic_poly_t rop, const padic_poly_t op, long n, 
                            const padic_ctx_t ctx);

/*  Evaluation  **************************************************************/

void _padic_poly_evaluate_padic(fmpz_t u, long *v, 
                                const fmpz *poly, long val, long len, 
                                const fmpz_t a, long b, const padic_ctx_t ctx);

void padic_poly_evaluate_padic(padic_t y, const padic_poly_t poly, 
                                          const padic_t x, const padic_ctx_t ctx);

/*  Composition  *************************************************************/

void _padic_poly_compose(fmpz *rop, long *rval, 
                         const fmpz *op1, long val1, long len1, 
                         const fmpz *op2, long val2, long len2, 
                         const padic_ctx_t ctx);

void padic_poly_compose(padic_poly_t rop, 
                        const padic_poly_t op1, const padic_poly_t op2, 
                        const padic_ctx_t ctx);

void _padic_poly_compose_pow(fmpz *rop, long *rval, 
                             const fmpz *op, long val, long len, long k, 
                             const padic_ctx_t ctx);

void padic_poly_compose_pow(padic_poly_t rop, const padic_poly_t op, long k, 
                            const padic_ctx_t ctx);

/*  Input and output  ********************************************************/

static __inline__ int padic_poly_debug(const padic_poly_t poly)
{
    printf("(alloc = %ld, length = %ld, val = %ld, vec = ", 
        poly->alloc, poly->length, poly->val);
    if (poly->coeffs)
    {
        printf("{");
        _fmpz_vec_print(poly->coeffs, poly->alloc);
        printf("}");
    }
    else
    {
        printf("NULL");
    }
    printf(")");
    fflush(stdout);

    return 1;
}

int _padic_poly_fprint(FILE *file, const fmpz *poly, long val, long len, 
                       const padic_ctx_t ctx);

int padic_poly_fprint(FILE *file, const padic_poly_t poly, 
                      const padic_ctx_t ctx);

static __inline__ 
int _padic_poly_print(const fmpz *poly, long val, long len, 
                      const padic_ctx_t ctx)
{
    return _padic_poly_fprint(stdout, poly, val, len, ctx);
}

static __inline__ 
int padic_poly_print(const padic_poly_t poly, const padic_ctx_t ctx)
{
    return padic_poly_fprint(stdout, poly, ctx);
}

int _padic_poly_fprint_pretty(FILE *file, 
                              const fmpz *poly, long val, long len, 
                              const char *var, 
                              const padic_ctx_t ctx);

int padic_poly_fprint_pretty(FILE *file, 
                             const padic_poly_t poly, const char *var, 
                             const padic_ctx_t ctx);

static __inline__ 
int _padic_poly_print_pretty(FILE *file, 
                             const fmpz *poly, long val, long len, 
                             const char *var, 
                             const padic_ctx_t ctx)
{
    return _padic_poly_fprint_pretty(stdout, poly, val, len, var, ctx);
}

static __inline__ 
int padic_poly_print_pretty(const padic_poly_t poly, const char *var, 
                            const padic_ctx_t ctx)
{
    return padic_poly_fprint_pretty(stdout, poly, var, ctx);
}

#endif

