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
#include "ulong_extras.h"

typedef long padic_t[2];

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
    enum padic_print_mode mode;
} padic_ctx_struct;

typedef padic_ctx_struct padic_ctx_t[1];

/* Context *******************************************************************/

void padic_ctx_init(padic_ctx_t ctx, const fmpz_t p, long N, 
                    enum padic_print_mode mode);

void padic_ctx_clear(padic_ctx_t ctx);

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

static __inline__
void padic_swap(padic_t op1, padic_t op2, const padic_ctx_t ctx)
{
    long t;

    fmpz_swap(op1, op2);
    t      = op1[1];
    op1[1] = op2[1];
    op2[1] = t;
}

static __inline__ 
void padic_zero(padic_t rop, const padic_ctx_t ctx)
{
    fmpz_zero(rop);
    rop[1] = 0;
}

static __inline__ 
void padic_one(padic_t rop, const padic_ctx_t ctx)
{
    fmpz_set_ui(rop, 1);
    rop[1] = 0;
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

void padic_inv(padic_t rop, const padic_t op, const padic_ctx_t ctx);

void padic_pow_si(padic_t rop, const padic_t op, long e, 
                  const padic_ctx_t ctx);

/* Comparison ****************************************************************/

static __inline__
int padic_is_zero(const padic_t op, const padic_ctx_t ctx)
{
    return fmpz_is_zero(op);
}

static __inline__
int padic_is_one(const padic_t op, const padic_ctx_t ctx)
{
    return (op[1] == 0) && fmpz_is_one(op);
}

static __inline__
int padic_equal(const padic_t op1, const padic_t op2, const padic_ctx_t ctx)
{
    return (op1[1] == op2[1]) && fmpz_equal(op1, op2);
}

/* Input and output **********************************************************/

char * padic_get_str(const padic_t op, const padic_ctx_t ctx);

int padic_set_str(padic_t rop, const char *str, const padic_ctx_t ctx);

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

