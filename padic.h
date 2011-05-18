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

#ifndef PADIC_H
#define PADIC_H

#include <stdlib.h>
#include <stdio.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

typedef long padic_t[2];

#define padic_unit(x)  (x)
#define padic_val(x)   ((x)[1])

enum padic_print_mode
{
    PADIC_TERSE, 
    PADIC_SERIES, 
    PADIC_VAL_UNIT
};

typedef struct {

    fmpz_t p;
    long N;

    double pinv;

    fmpz *pow;
    long min;
    long max;

    enum padic_print_mode mode;

} padic_ctx_struct;

typedef padic_ctx_struct padic_ctx_t[1];

/* Context *******************************************************************/

void padic_ctx_init(padic_ctx_t ctx, const fmpz_t p, long N, 
                    enum padic_print_mode mode);

void padic_ctx_clear(padic_ctx_t ctx);

static __inline__ 
void _padic_ctx_pow_ui(fmpz *rop, int *alloc, ulong e, const padic_ctx_t ctx)
{
    if (ctx->min <= e && e < ctx->max)
    {
        *alloc = 0;
        *rop   = *(ctx->pow + (e - ctx->min));
    }
    else
    {
        *alloc = 1;
        *rop   = 0L;
        fmpz_pow_ui(rop, ctx->p, e);
    }
}

/* Memory management *********************************************************/

void padic_init(padic_t rop, const padic_ctx_t ctx);

void padic_clear(padic_t rop, const padic_ctx_t ctx);

void padic_normalise(padic_t rop, const padic_ctx_t ctx);

void _padic_reduce_unit(padic_t rop, const padic_ctx_t ctx);

/* Randomisation *************************************************************/

void padic_randtest(padic_t rop, flint_rand_t state, const padic_ctx_t ctx);

void padic_randtest_not_zero(padic_t rop, flint_rand_t state, 
                             const padic_ctx_t ctx);

/* Assignment ****************************************************************/

void padic_set(padic_t rop, const padic_t op, const padic_ctx_t ctx);

void padic_set_si(padic_t rop, long op, const padic_ctx_t ctx);

void padic_set_ui(padic_t rop, ulong op, const padic_ctx_t ctx);

void padic_set_fmpz(padic_t rop, const fmpz_t op, const padic_ctx_t ctx);

void padic_set_mpz(padic_t rop, const mpz_t op, const padic_ctx_t ctx);

void padic_set_mpq(padic_t rop, const mpq_t op, const padic_ctx_t ctx);

void padic_get_fmpz(fmpz_t rop, const padic_t op, const padic_ctx_t ctx);

void padic_get_mpz(mpz_t rop, const padic_t op, const padic_ctx_t ctx);

void padic_get_mpq(mpq_t rop, const padic_t op, const padic_ctx_t ctx);

static __inline__
void padic_swap(padic_t op1, padic_t op2, const padic_ctx_t ctx)
{
    long t;

    fmpz_swap(padic_unit(op1), padic_unit(op2));
    t              = padic_val(op1);
    padic_val(op1) = padic_val(op2);
    padic_val(op2) = t;
}

static __inline__ 
void padic_zero(padic_t rop, const padic_ctx_t ctx)
{
    fmpz_zero(padic_unit(rop));
    padic_val(rop) = 0;
}

static __inline__ 
void padic_one(padic_t rop, const padic_ctx_t ctx)
{
    if (ctx->N > 0)
    {
        fmpz_set_ui(padic_unit(rop), 1);
        padic_val(rop) = 0;
    }
    else
        padic_zero(rop, ctx);
}

/* Arithmetic operations *****************************************************/

void padic_add(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx);

void padic_sub(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx);

void padic_neg(padic_t rop, const padic_t op, const padic_ctx_t ctx);

void padic_mul(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx);

void padic_div(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx);

void padic_shift(padic_t rop, const padic_t op, long v, const padic_ctx_t ctx);

void _padic_inv_naive(fmpz_t rop, const fmpz_t op, const fmpz_t p, long N);

void padic_inv_naive(padic_t rop, const padic_t op, const padic_ctx_t ctx);

void _padic_inv_hensel(fmpz_t rop, const fmpz_t op, const fmpz_t p, long N);

void padic_inv_hensel(padic_t rop, const padic_t op, const padic_ctx_t ctx);

static __inline__ 
void _padic_inv(fmpz_t rop, const fmpz_t op, const fmpz_t p, long N)
{
    _padic_inv_hensel(rop, op, p, N);
}

static __inline__ 
void padic_inv(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    padic_inv_hensel(rop, op, ctx);
}

int padic_sqrt(padic_t rop, const padic_t op, const padic_ctx_t ctx);

void padic_pow_si(padic_t rop, const padic_t op, long e, 
                  const padic_ctx_t ctx);

int padic_exp(padic_t rop, const padic_t op, const padic_ctx_t ctx);

/* Comparison ****************************************************************/

static __inline__
int padic_is_zero(const padic_t op, const padic_ctx_t ctx)
{
    return fmpz_is_zero(padic_unit(op)) || (padic_val(op) >= ctx->N);
}

static __inline__
int padic_is_one(const padic_t op, const padic_ctx_t ctx)
{
    if (ctx->N > 0)
        return (padic_val(op) == 0) && fmpz_is_one(padic_unit(op));
    else
        return padic_val(op) >= ctx->N;
}

static __inline__
int padic_equal(const padic_t op1, const padic_t op2, const padic_ctx_t ctx)
{
    /* Quick guess */
    if (   (padic_val(op1) == padic_val(op2))
        && (fmpz_equal(padic_unit(op1), padic_unit(op2))))
        return 1;

    if (padic_is_zero(op1, ctx))
        return padic_is_zero(op2, ctx);

    if (padic_is_zero(op2, ctx))
        return 0;

    if (padic_val(op1) == padic_val(op2))
    {
        fmpz_t pow;
        int alloc, ans;

        _padic_ctx_pow_ui(pow, &alloc, ctx->N - padic_val(op1), ctx);

        if (   (fmpz_sgn(padic_unit(op1)) > 0) 
            && (fmpz_sgn(padic_unit(op2)) > 0)
            && (fmpz_cmp(padic_unit(op1), pow) < 0) 
            && (fmpz_cmp(padic_unit(op2), pow) < 0))
        {
            ans = fmpz_equal(padic_unit(op1), padic_unit(op2));
        }
        else
        {
            fmpz_t u1, u2;

            fmpz_init(u1);
            fmpz_init(u2);
            fmpz_mod(u1, padic_unit(op1), pow);
            fmpz_mod(u2, padic_unit(op2), pow);
            ans = fmpz_equal(u1, u2);
            fmpz_clear(u1);
            fmpz_clear(u2);
        }

        if (alloc)
            fmpz_clear(pow);

        return ans;
    }
    else
        return 0;
}

/* Special functions *********************************************************/

void padic_teichmuller(padic_t rop, const padic_t op, const padic_ctx_t ctx);

static 
void padic_val_factorial(fmpz_t rop, const fmpz_t op, const fmpz_t p)
{
    fmpz_t t, q, pow;

    if (fmpz_sgn(op) <= 0)
    {
        printf("Exception (padic_val_factorial).  op is non-positive.\n");
        abort();
    }

    fmpz_init(t);
    fmpz_init(q);
    fmpz_init(pow);
    fmpz_set_ui(pow, 1);

    do 
    {
        fmpz_mul(pow, pow, p);
        fmpz_fdiv_q(q, op, pow);
        fmpz_add(t, t, q);
    }
    while (!fmpz_is_zero(q));

    fmpz_swap(rop, t);
    fmpz_clear(t);
    fmpz_clear(q);
    fmpz_clear(pow);
}

/* Input and output **********************************************************/

char * padic_get_str(const padic_t op, const padic_ctx_t ctx);

int padic_fprint(FILE * file, const padic_t op, const padic_ctx_t ctx);

static __inline__
int padic_print(const padic_t op, const padic_ctx_t ctx)
{
    return padic_fprint(stdout, op, ctx);
}

static __inline__
void padic_debug(const padic_t op, const padic_ctx_t ctx)
{
    printf("[u = "), fmpz_print(op), printf(", v = %ld]", op[1]);
}

#endif

