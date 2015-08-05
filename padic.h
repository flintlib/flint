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

#ifndef PADIC_H
#define PADIC_H

#ifdef PADIC_INLINES_C
#define PADIC_INLINE FLINT_DLL
#else
#define PADIC_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong

#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

#ifdef __cplusplus
 extern "C" {
#endif

#define PADIC_DEFAULT_PREC WORD(20)

#define PADIC_TEST_PREC_MIN WORD(-100)
#define PADIC_TEST_PREC_MAX  WORD(100)

typedef struct {
    fmpz u;
    slong v;
    slong N;
} padic_struct;

typedef padic_struct padic_t[1];

#define padic_val(x)   ((x)->v)
#define padic_prec(x)  ((x)->N)

PADIC_INLINE
fmpz * padic_unit(const padic_t x)
{
   return (fmpz *)(&(x->u));
}

PADIC_INLINE
slong padic_get_val(const padic_t x)
{
   return x->v;
}

PADIC_INLINE
slong padic_get_prec(const padic_t x)
{
   return x->N;
}

enum padic_print_mode
{
    PADIC_TERSE, 
    PADIC_SERIES, 
    PADIC_VAL_UNIT
};

typedef struct {

    fmpz_t p;

    double pinv;

    fmpz *pow;
    slong min;
    slong max;

    enum padic_print_mode mode;

} padic_ctx_struct;

typedef padic_ctx_struct padic_ctx_t[1];

typedef struct {
    slong n;
    fmpz *pow;
} padic_inv_struct;

typedef padic_inv_struct padic_inv_t[1];

/* Context *******************************************************************/

FLINT_DLL void padic_ctx_init(padic_ctx_t ctx, const fmpz_t p, slong min, slong max, 
                    enum padic_print_mode mode);

FLINT_DLL void padic_ctx_clear(padic_ctx_t ctx);

PADIC_INLINE 
int _padic_ctx_pow_ui(fmpz_t rop, ulong e, const padic_ctx_t ctx)
{
    if (ctx->min <= (slong) e && (slong) e < ctx->max)
    {
        *rop   = *(ctx->pow + (e - ctx->min));
        return 0;
    }
    else
    {
        slong l = (slong) e;
        if (l < 0)
        {
            flint_printf("Exception (_padic_ctx_pow_ui). Power too large.\n");
            flint_printf("e = %wu\n", e);
            flint_printf("l = %wd\n", l);
            abort();
        }

        fmpz_init(rop);
        fmpz_pow_ui(rop, ctx->p, e);
        return 1;
    }
}

PADIC_INLINE 
void padic_ctx_pow_ui(fmpz_t rop, ulong e, const padic_ctx_t ctx)
{
    if (ctx->min <= (slong) e && (slong) e < ctx->max)
        fmpz_set(rop, ctx->pow + (e - ctx->min));
    else
    {
        slong l = (slong) e;
        if (l < 0)
        {
            flint_printf("Exception (padic_ctx_pow_ui). Power too large.\n");
            flint_printf("e = %wu\n", e);
            flint_printf("l = %wd\n", l);
            abort();
        }

        fmpz_pow_ui(rop, ctx->p, e);
    }
}

/* Memory management *********************************************************/

FLINT_DLL void padic_init(padic_t rop);

FLINT_DLL void padic_init2(padic_t rop, slong N);

FLINT_DLL void padic_clear(padic_t rop);

PADIC_INLINE void _padic_canonicalise(padic_t rop, const padic_ctx_t ctx)
{
    if (!fmpz_is_zero(padic_unit(rop)))
    {
        padic_val(rop) += _fmpz_remove(padic_unit(rop), ctx->p, ctx->pinv);
    }
    else
    {
        padic_val(rop) = 0;
    }
}

FLINT_DLL void _padic_reduce(padic_t rop, const padic_ctx_t ctx);

FLINT_DLL void padic_reduce(padic_t rop, const padic_ctx_t ctx);

/* Randomisation *************************************************************/

FLINT_DLL void padic_randtest(padic_t rop, flint_rand_t state, const padic_ctx_t ctx);

FLINT_DLL void padic_randtest_not_zero(padic_t rop, flint_rand_t state, 
                             const padic_ctx_t ctx);

FLINT_DLL void padic_randtest_int(padic_t rop, flint_rand_t state, 
                        const padic_ctx_t ctx);

/* Assignments and conversions ***********************************************/

FLINT_DLL void padic_set(padic_t rop, const padic_t op, const padic_ctx_t ctx);

FLINT_DLL void padic_set_si(padic_t rop, slong op, const padic_ctx_t ctx);

FLINT_DLL void padic_set_ui(padic_t rop, ulong op, const padic_ctx_t ctx);

FLINT_DLL void padic_set_fmpz(padic_t rop, const fmpz_t op, const padic_ctx_t ctx);

FLINT_DLL void padic_set_fmpq(padic_t rop, const fmpq_t op, const padic_ctx_t ctx);

FLINT_DLL void padic_set_mpz(padic_t rop, const mpz_t op, const padic_ctx_t ctx);

FLINT_DLL void padic_set_mpq(padic_t rop, const mpq_t op, const padic_ctx_t ctx);

FLINT_DLL void padic_get_fmpz(fmpz_t rop, const padic_t op, const padic_ctx_t ctx);

FLINT_DLL void padic_get_fmpq(fmpq_t rop, const padic_t op, const padic_ctx_t ctx);

FLINT_DLL void padic_get_mpz(mpz_t rop, const padic_t op, const padic_ctx_t ctx);

FLINT_DLL void padic_get_mpq(mpq_t rop, const padic_t op, const padic_ctx_t ctx);

PADIC_INLINE void padic_swap(padic_t op1, padic_t op2)
{
    slong t;

    fmpz_swap(padic_unit(op1), padic_unit(op2));

    t              = padic_val(op1);
    padic_val(op1) = padic_val(op2);
    padic_val(op2) = t;

    t               = padic_prec(op1);
    padic_prec(op1) = padic_prec(op2);
    padic_prec(op2) = t;
}

PADIC_INLINE void padic_zero(padic_t rop)
{
    fmpz_zero(padic_unit(rop));
    padic_val(rop) = 0;
}

PADIC_INLINE void padic_one(padic_t rop)
{
    if (padic_prec(rop) > 0)
    {
        fmpz_one(padic_unit(rop));
        padic_val(rop) = 0;
    }
    else
    {
        padic_zero(rop);
    }
}

/* Comparison ****************************************************************/

PADIC_INLINE int padic_is_zero(const padic_t op)
{
    return fmpz_is_zero(padic_unit(op));
}

PADIC_INLINE int padic_is_one(const padic_t op)
{
    return fmpz_is_one(padic_unit(op)) && (padic_val(op) == 0);
}

PADIC_INLINE int padic_equal(const padic_t op1, const padic_t op2)
{
    return (padic_val(op1) == padic_val(op2)) && 
           (fmpz_equal(padic_unit(op1), padic_unit(op2)));
}

/* Arithmetic operations *****************************************************/

slong * _padic_lifts_exps(slong *n, slong N);

FLINT_DLL void _padic_lifts_pows(fmpz *pow, const slong *a, slong n, const fmpz_t p);

FLINT_DLL void padic_add(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx);

FLINT_DLL void padic_sub(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx);

FLINT_DLL void padic_neg(padic_t rop, const padic_t op, const padic_ctx_t ctx);

FLINT_DLL void padic_mul(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx);

FLINT_DLL void padic_shift(padic_t rop, const padic_t op, slong v, const padic_ctx_t ctx);

FLINT_DLL void padic_div(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx);

FLINT_DLL void _padic_inv_precompute(padic_inv_t S, const fmpz_t p, slong N);

FLINT_DLL void _padic_inv_clear(padic_inv_t S);

FLINT_DLL void _padic_inv_precomp(fmpz_t rop, const fmpz_t op, const padic_inv_t S);

FLINT_DLL void _padic_inv(fmpz_t rop, const fmpz_t op, const fmpz_t p, slong N);

FLINT_DLL void padic_inv(padic_t rop, const padic_t op, const padic_ctx_t ctx);

FLINT_DLL int padic_sqrt(padic_t rop, const padic_t op, const padic_ctx_t ctx);

FLINT_DLL void padic_pow_si(padic_t rop, const padic_t op, slong e, 
                  const padic_ctx_t ctx);

/* Exponential ***************************************************************/

FLINT_DLL slong _padic_exp_bound(slong v, slong N, const fmpz_t p);

FLINT_DLL void _padic_exp(fmpz_t rop, const fmpz_t u, slong v, const fmpz_t p, slong N);
FLINT_DLL void _padic_exp_rectangular(fmpz_t rop, const fmpz_t u, slong v, const fmpz_t p, slong N);
FLINT_DLL void _padic_exp_balanced(fmpz_t rop, const fmpz_t u, slong v, const fmpz_t p, slong N);

FLINT_DLL int padic_exp(padic_t rop, const padic_t op, const padic_ctx_t ctx);
FLINT_DLL int padic_exp_rectangular(padic_t rop, const padic_t op, const padic_ctx_t ctx);
FLINT_DLL int padic_exp_balanced(padic_t rop, const padic_t op, const padic_ctx_t ctx);

/* Logarithm *****************************************************************/

FLINT_DLL slong _padic_log_bound(slong v, slong N, const fmpz_t p);

FLINT_DLL void _padic_log(fmpz_t z, const fmpz_t y, slong v, const fmpz_t p, slong N);
FLINT_DLL void _padic_log_rectangular(fmpz_t z, const fmpz_t y, slong v, const fmpz_t p, slong N);
FLINT_DLL void _padic_log_satoh(fmpz_t z, const fmpz_t y, slong v, const fmpz_t p, slong N);
FLINT_DLL void _padic_log_balanced(fmpz_t z, const fmpz_t y, slong v, const fmpz_t p, slong N);

FLINT_DLL int padic_log(padic_t rop, const padic_t op, const padic_ctx_t ctx);
FLINT_DLL int padic_log_rectangular(padic_t rop, const padic_t op, const padic_ctx_t ctx);
FLINT_DLL int padic_log_satoh(padic_t rop, const padic_t op, const padic_ctx_t ctx);
FLINT_DLL int padic_log_balanced(padic_t rop, const padic_t op, const padic_ctx_t ctx);

/* Special functions *********************************************************/

FLINT_DLL void _padic_teichmuller(fmpz_t rop, const fmpz_t op, const fmpz_t p, slong N);

FLINT_DLL void padic_teichmuller(padic_t rop, const padic_t op, const padic_ctx_t ctx);

FLINT_DLL ulong padic_val_fac_ui_2(ulong N);

FLINT_DLL ulong padic_val_fac_ui(ulong N, const fmpz_t p);

FLINT_DLL void padic_val_fac(fmpz_t rop, const fmpz_t op, const fmpz_t p);

/* Input and output **********************************************************/

char * padic_get_str(char * str, const padic_t op, const padic_ctx_t ctx);

FLINT_DLL int _padic_fprint(FILE * file, const fmpz_t u, slong v, const padic_ctx_t ctx);

FLINT_DLL int padic_fprint(FILE * file, const padic_t op, const padic_ctx_t ctx);

PADIC_INLINE 
int _padic_print(const fmpz_t u, slong v, const padic_ctx_t ctx)
{
    return _padic_fprint(stdout, u, v, ctx);
}

PADIC_INLINE int padic_print(const padic_t op, const padic_ctx_t ctx)
{
    return padic_fprint(stdout, op, ctx);
}

PADIC_INLINE void padic_debug(const padic_t op)
{
    flint_printf("(");
    fmpz_print(padic_unit(op)); 
    flint_printf(" %wd %wd)", padic_val(op), padic_prec(op));
}

#ifdef __cplusplus
}
#endif

#endif

