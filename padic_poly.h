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

    Copyright (C) 2011, 2012 Sebastian Pancratz
 
******************************************************************************/

#ifndef PADIC_POLY_H
#define PADIC_POLY_H

#ifdef PADIC_POLY_INLINES_C
#define PADIC_POLY_INLINE FLINT_DLL
#else
#define PADIC_POLY_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx/* interferes with system includes */
#include <limits.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "fmpz.h"
#include "fmpq.h"
#include "padic.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif


/*  Type definitions  ********************************************************/

typedef struct
{
    fmpz *coeffs;
    slong alloc;
    slong length;
    slong val;
    slong N;
} padic_poly_struct;

typedef padic_poly_struct padic_poly_t[1];

/*  Helper functions  ********************************************************/

/*
    Returns the minimum $p$-adic valuation of \code{(vec, len)}, 
    assuming this fits into a \code{signed long}.

    If \code{len} is zero, returns $0$.
 */
PADIC_POLY_INLINE slong _fmpz_vec_ord_p(const fmpz *vec, slong len, const fmpz_t p)
{
    if (len == 0)
    {
        return 0;
    }
    else
    {
        fmpz_t t;
        slong i, min = WORD_MAX, v;

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
        return (min < WORD_MAX) ? min : 0;
    }
}

/*  Memory management  *******************************************************/

FLINT_DLL void padic_poly_init(padic_poly_t poly);

FLINT_DLL void padic_poly_init2(padic_poly_t poly, slong alloc, slong prec);

FLINT_DLL void padic_poly_clear(padic_poly_t poly);

FLINT_DLL void padic_poly_realloc(padic_poly_t f, slong alloc, const fmpz_t p);

FLINT_DLL void padic_poly_fit_length(padic_poly_t f, slong len);

PADIC_POLY_INLINE
void _padic_poly_set_length(padic_poly_t poly, slong len)
{
    if (poly->length > len)
    {
        slong i;

        for (i = len; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i); 
    }
    poly->length = len;
}

FLINT_DLL void _padic_poly_normalise(padic_poly_t f);

FLINT_DLL void _padic_poly_canonicalise(fmpz *poly, slong *v, slong len, const fmpz_t p);

FLINT_DLL void padic_poly_canonicalise(padic_poly_t poly, const fmpz_t p);

FLINT_DLL void padic_poly_reduce(padic_poly_t f, const padic_ctx_t ctx);

PADIC_POLY_INLINE 
void padic_poly_truncate(padic_poly_t poly, slong n, const fmpz_t p)
{
    if (poly->length > n)
    {
        slong i;

        for (i = n; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i);
        poly->length = n;
        _padic_poly_normalise(poly);
        padic_poly_canonicalise(poly, p);
    }
}

/*  Polynomial parameters  ***************************************************/

PADIC_POLY_INLINE slong padic_poly_degree(const padic_poly_t poly)
{
    return poly->length - 1;
}

PADIC_POLY_INLINE slong padic_poly_length(const padic_poly_t poly)
{
    return poly->length;
}

PADIC_POLY_INLINE slong padic_poly_val(const padic_poly_t poly)
{
    return poly->val;
}

#define padic_poly_val(poly)   ((poly)->val)

#define padic_poly_prec(poly)  ((poly)->N)

/*  Randomisation  ***********************************************************/

FLINT_DLL void padic_poly_randtest(padic_poly_t f, flint_rand_t state, 
                         slong len, const padic_ctx_t ctx);

FLINT_DLL void padic_poly_randtest_not_zero(padic_poly_t f, flint_rand_t state,
                                  slong len, const padic_ctx_t ctx);

FLINT_DLL void padic_poly_randtest_val(padic_poly_t f, flint_rand_t state, 
                             slong val, slong len, const padic_ctx_t ctx);

/*  Assignment and basic manipulation  ***************************************/

FLINT_DLL void padic_poly_set(padic_poly_t f, const padic_poly_t g, 
                    const padic_ctx_t ctx);

FLINT_DLL void padic_poly_set_padic(padic_poly_t poly, const padic_t x, 
                          const padic_ctx_t ctx);

FLINT_DLL void padic_poly_set_si(padic_poly_t poly, slong x, const padic_ctx_t ctx);

FLINT_DLL void padic_poly_set_ui(padic_poly_t poly, ulong x, const padic_ctx_t ctx);

FLINT_DLL void padic_poly_set_fmpz(padic_poly_t poly, const fmpz_t x, 
                         const padic_ctx_t ctx);

FLINT_DLL void padic_poly_set_fmpq(padic_poly_t poly, const fmpq_t x, 
                         const padic_ctx_t ctx);

FLINT_DLL void padic_poly_set_fmpz_poly(padic_poly_t rop, const fmpz_poly_t op, 
                              const padic_ctx_t ctx);

FLINT_DLL void padic_poly_set_fmpq_poly(padic_poly_t rop, 
                              const fmpq_poly_t op, const padic_ctx_t ctx);

FLINT_DLL int padic_poly_get_fmpz_poly(fmpz_poly_t rop, const padic_poly_t op, 
                             const padic_ctx_t ctx);

FLINT_DLL void padic_poly_get_fmpq_poly(fmpq_poly_t rop, 
                              const padic_poly_t op, const padic_ctx_t ctx);

PADIC_POLY_INLINE void padic_poly_zero(padic_poly_t poly)
{
    _padic_poly_set_length(poly, 0);
    poly->val = 0;
}

PADIC_POLY_INLINE void padic_poly_one(padic_poly_t poly)
{
    if (padic_poly_prec(poly) > 0)
    {
        padic_poly_fit_length(poly, 1);
        fmpz_one(poly->coeffs);
        _padic_poly_set_length(poly, 1);
        poly->val = 0;
    }
    else
    {
        padic_poly_zero(poly);
    }
}

FLINT_DLL void padic_poly_swap(padic_poly_t poly1, padic_poly_t poly2);

/*  Getting and setting coefficients  ****************************************/

FLINT_DLL void padic_poly_get_coeff_padic(padic_t c, const padic_poly_t poly, slong n, 
                                           const padic_ctx_t ctx);

FLINT_DLL void padic_poly_set_coeff_padic(padic_poly_t f, slong n, const padic_t c, 
                                                const padic_ctx_t ctx);

/*  Comparison  **************************************************************/

FLINT_DLL int padic_poly_equal(const padic_poly_t f, const padic_poly_t g);

PADIC_POLY_INLINE int padic_poly_is_zero(const padic_poly_t poly)
{
    return poly->length == 0;
}

PADIC_POLY_INLINE int padic_poly_is_one(const padic_poly_t poly)
{
    return (poly->length == 1) && fmpz_is_one(poly->coeffs) && 
           (poly->val == 0);
}

/*  Addition and subtraction  ************************************************/

FLINT_DLL void _padic_poly_add(fmpz *rop, slong *rval, slong N, 
                     const fmpz *op1, slong val1, slong len1, slong N1, 
                     const fmpz *op2, slong val2, slong len2, slong N2, 
                     const padic_ctx_t ctx);

FLINT_DLL void padic_poly_add(padic_poly_t f, 
                    const padic_poly_t g, const padic_poly_t h, 
                    const padic_ctx_t ctx);

FLINT_DLL void _padic_poly_sub(fmpz *rop, slong *rval, slong N, 
                     const fmpz *op1, slong val1, slong len1, slong N1, 
                     const fmpz *op2, slong val2, slong len2, slong N2, 
                     const padic_ctx_t ctx);

FLINT_DLL void padic_poly_sub(padic_poly_t f, 
                    const padic_poly_t g, const padic_poly_t h, 
                    const padic_ctx_t ctx);

FLINT_DLL void padic_poly_neg(padic_poly_t f, const padic_poly_t g, 
                    const padic_ctx_t ctx);

/*  Scalar multiplication and division  **************************************/

FLINT_DLL void _padic_poly_scalar_mul_padic(fmpz *rop, slong *rval, slong N, 
                                  const fmpz *op, slong val, slong len, 
                                  const padic_t c, const padic_ctx_t ctx);

FLINT_DLL void padic_poly_scalar_mul_padic(padic_poly_t rop, const padic_poly_t op, 
                                 const padic_t c, const padic_ctx_t ctx);

/*  Multiplication  **********************************************************/

FLINT_DLL void _padic_poly_mul(fmpz *rop, slong *rval, slong N, 
                     const fmpz *op1, slong val1, slong len1, 
                     const fmpz *op2, slong val2, slong len2, 
                     const padic_ctx_t ctx);

FLINT_DLL void padic_poly_mul(padic_poly_t f, 
                    const padic_poly_t g, const padic_poly_t h, 
                    const padic_ctx_t ctx);

/*  Powering  ****************************************************************/

FLINT_DLL void _padic_poly_pow(fmpz *rop, slong *rval, slong N, 
                     const fmpz *op, slong val, slong len, ulong e,
                     const padic_ctx_t ctx);

FLINT_DLL void padic_poly_pow(padic_poly_t rop, const padic_poly_t op, ulong e, 
                    const padic_ctx_t ctx);

/*  Series inversion  ********************************************************/

FLINT_DLL void padic_poly_inv_series(padic_poly_t Qinv, const padic_poly_t Q, slong n, 
                           const padic_ctx_t ctx);

/*  Derivative  **************************************************************/

FLINT_DLL void _padic_poly_derivative(fmpz *rop, slong *rval, slong N, 
                            const fmpz *op, slong val, slong len, 
                            const padic_ctx_t ctx);

FLINT_DLL void padic_poly_derivative(padic_poly_t rop, 
                           const padic_poly_t op, const padic_ctx_t ctx);

/*  Shifting  ****************************************************************/

FLINT_DLL void padic_poly_shift_left(padic_poly_t rop, const padic_poly_t op, slong n, 
                           const padic_ctx_t ctx);

FLINT_DLL void padic_poly_shift_right(padic_poly_t rop, const padic_poly_t op, slong n, 
                            const padic_ctx_t ctx);

/*  Evaluation  **************************************************************/

FLINT_DLL void _padic_poly_evaluate_padic(fmpz_t u, slong *v, slong N, 
                                const fmpz *poly, slong val, slong len, 
                                const fmpz_t a, slong b, const padic_ctx_t ctx);

FLINT_DLL void padic_poly_evaluate_padic(padic_t y, const padic_poly_t poly, 
                                          const padic_t x, const padic_ctx_t ctx);

/*  Composition  *************************************************************/

FLINT_DLL void _padic_poly_compose(fmpz *rop, slong *rval, slong N, 
                         const fmpz *op1, slong val1, slong len1, 
                         const fmpz *op2, slong val2, slong len2, 
                         const padic_ctx_t ctx);

FLINT_DLL void padic_poly_compose(padic_poly_t rop, 
                        const padic_poly_t op1, const padic_poly_t op2, 
                        const padic_ctx_t ctx);

FLINT_DLL void _padic_poly_compose_pow(fmpz *rop, slong *rval, slong N, 
                             const fmpz *op, slong val, slong len, slong k, 
                             const padic_ctx_t ctx);

FLINT_DLL void padic_poly_compose_pow(padic_poly_t rop, const padic_poly_t op, slong k, 
                            const padic_ctx_t ctx);

/*  Input and output  ********************************************************/

PADIC_POLY_INLINE int padic_poly_debug(const padic_poly_t poly)
{
    flint_printf("(alloc = %wd, length = %wd, val = %wd, N = %wd, vec = ", 
        poly->alloc, poly->length, poly->val, poly->N);
    if (poly->coeffs)
    {
        flint_printf("{");
        _fmpz_vec_print(poly->coeffs, poly->alloc);
        flint_printf("}");
    }
    else
    {
        flint_printf("NULL");
    }
    flint_printf(")");

    return 1;
}

FLINT_DLL int _padic_poly_fprint(FILE *file, const fmpz *poly, slong val, slong len, 
                       const padic_ctx_t ctx);

FLINT_DLL int padic_poly_fprint(FILE *file, const padic_poly_t poly, 
                      const padic_ctx_t ctx);

PADIC_POLY_INLINE 
int _padic_poly_print(const fmpz *poly, slong val, slong len, 
                      const padic_ctx_t ctx)
{
    return _padic_poly_fprint(stdout, poly, val, len, ctx);
}

PADIC_POLY_INLINE 
int padic_poly_print(const padic_poly_t poly, const padic_ctx_t ctx)
{
    return padic_poly_fprint(stdout, poly, ctx);
}

FLINT_DLL int _padic_poly_fprint_pretty(FILE *file, 
                              const fmpz *poly, slong val, slong len, 
                              const char *var, 
                              const padic_ctx_t ctx);

FLINT_DLL int padic_poly_fprint_pretty(FILE *file, 
                             const padic_poly_t poly, const char *var, 
                             const padic_ctx_t ctx);

PADIC_POLY_INLINE 
int _padic_poly_print_pretty(FILE *file, 
                             const fmpz *poly, slong val, slong len, 
                             const char *var, 
                             const padic_ctx_t ctx)
{
    return _padic_poly_fprint_pretty(stdout, poly, val, len, var, ctx);
}

PADIC_POLY_INLINE 
int padic_poly_print_pretty(const padic_poly_t poly, const char *var, 
                            const padic_ctx_t ctx)
{
    return padic_poly_fprint_pretty(stdout, poly, var, ctx);
}

/*  Testing  *****************************************************************/

FLINT_DLL int _padic_poly_is_canonical(const fmpz *op, slong val, slong len, 
                             const padic_ctx_t ctx);

FLINT_DLL int padic_poly_is_canonical(const padic_poly_t op, const padic_ctx_t ctx);

FLINT_DLL int _padic_poly_is_reduced(const fmpz *op, slong val, slong len, slong N, 
                           const padic_ctx_t ctx);

FLINT_DLL int padic_poly_is_reduced(const padic_poly_t op, const padic_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

