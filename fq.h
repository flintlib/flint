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

#undef ulong                    /* interferes with system includes */
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
#include "qadic.h"

/* Data types and context ****************************************************/

/*
#define fq_t padic_poly_t
#define fq_struct padic_poly_struct
*/

typedef qadic_t fq_t;
typedef qadic_struct fq_struct;
typedef qadic_ctx_struct fq_ctx_t[1];

void fq_ctx_init_conway(fq_ctx_t ctx,
                        const fmpz_t p, long d, const char *var,
                        enum padic_print_mode mode);

void fq_ctx_clear(fq_ctx_t ctx);

static __inline__ long
fq_ctx_dim(const fq_ctx_t ctx)
{
    return ctx->j[ctx->len - 1];
}

static __inline__ void
fq_ctx_print(const fq_ctx_t ctx)
{
    printf("p = "), fmpz_print((&ctx->pctx)->p), printf("\n");
    printf("d = %ld\n", ctx->j[ctx->len - 1]);
}

/* Basic arithmetic **********************************************************/

void fq_add(fq_t x, const fq_t y, const fq_t z, const fq_ctx_t ctx);

void fq_sub(fq_t x, const fq_t y, const fq_t z, const fq_ctx_t ctx);

void fq_neg(fq_t x, const fq_t y, const fq_ctx_t ctx);

void fq_mul(fq_t x, const fq_t y, const fq_t z, const fq_ctx_t ctx);

void fq_inv(qadic_t x, const qadic_t y, const qadic_ctx_t ctx);

void fq_pow(qadic_t x, const qadic_t y, const fmpz_t e, const qadic_ctx_t ctx);

/* Memory managment  *********************************************************/

static __inline__ void
fq_init(fq_t x)
{
    padic_poly_init(x);
}

static __inline__ void
fq_clear(fq_t x)
{
    padic_poly_clear(x);
}

/* Randomisation *************************************************************/

void fq_randtest(fq_t x, flint_rand_t state, const fq_ctx_t ctx);


void fq_randtest_not_zero(fq_t x, flint_rand_t state, const fq_ctx_t ctx);


void fq_randtest_val(fq_t x, flint_rand_t state, long val, const fq_ctx_t ctx);

void fq_randtest_int(fq_t x, flint_rand_t state, const fq_ctx_t ctx);


/* Comparison ****************************************************************/

int fq_equal(const fq_t x, const fq_t y);

static __inline__ int
fq_is_zero(const fq_t poly)
{
    return poly->length == 0;
}

static __inline__ int
fq_is_one(const fq_t poly)
{
    return (poly->length == 1 && fmpz_is_one(poly->coeffs));
}


/* Assignments and conversions ***********************************************/

void fq_set(fq_t x, const fq_t y);

static __inline__ void
fq_zero(fq_t x)
{
    padic_poly_zero(x);
}

static __inline__ void
fq_one(fq_t x, const fq_ctx_t ctx)
{
    padic_poly_one(x, &ctx->pctx);
}


/* Output ********************************************************************/

int fq_fprint_pretty(FILE * file, const fq_t op, const fq_ctx_t ctx);

int fq_print_pretty(const fq_t op, const fq_ctx_t ctx);


#endif
