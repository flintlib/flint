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

    Copyright (C) 2013 Mike Hansen
 
******************************************************************************/

#ifndef FQ_ZECH_H
#define FQ_ZECH_H

#include "fq_nmod.h"

/* Data types and context ****************************************************/
#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    mp_limb_t value;
} fq_zech_struct;

typedef fq_zech_struct fq_zech_t[1];

typedef struct
{
    mp_limb_t qm1;              /* q - 1 */
    mp_limb_t qm1o2;            /* (q - 1) / 2 or 1 when p == 2 */
    mp_limb_t qm1opm1;          /* (q - 1) / (p - 1) */
    mp_limb_t p;
    double ppre;
    mp_limb_t prime_root;       /* primitive root for prime subfield */
    mp_limb_t *zech_log_table;
    mp_limb_t *prime_field_table;
    mp_limb_t *eval_table;

    fq_nmod_ctx_struct *fq_nmod_ctx;
    int owns_fq_nmod_ctx;

} fq_zech_ctx_struct;

typedef fq_zech_ctx_struct fq_zech_ctx_t[1];

void fq_zech_ctx_init(fq_zech_ctx_t ctx, const fmpz_t p, slong d, const char *var);

void fq_zech_ctx_init_fq_nmod_ctx(fq_zech_ctx_t ctx, fq_nmod_ctx_t ctxn);

int _fq_zech_ctx_init_conway(fq_zech_ctx_t ctx, const fmpz_t p, slong d, const char *var);

void fq_zech_ctx_init_conway(fq_zech_ctx_t ctx, const fmpz_t p, slong d, const char *var);

void fq_zech_ctx_init_modulus(fq_zech_ctx_t ctx,
                              const nmod_poly_t modulus,
                              const char *var);

void fq_zech_ctx_randtest(fq_zech_ctx_t ctx, flint_rand_t state);

void fq_zech_ctx_clear(fq_zech_ctx_t ctx);

static __inline__ slong
fq_zech_ctx_degree(const fq_zech_ctx_t ctx)
{
    return fq_nmod_ctx_degree(ctx->fq_nmod_ctx);
}

static __inline__ void
fq_zech_ctx_order(fmpz_t f, const fq_zech_ctx_t ctx)
{
    fq_nmod_ctx_order(f, ctx->fq_nmod_ctx);
}

static __inline__ mp_limb_t
fq_zech_ctx_order_ui(const fq_zech_ctx_t ctx)
{
    return ctx->qm1 + 1;
}

#define fq_zech_ctx_prime(ctx)  fq_nmod_ctx_prime(ctx->fq_nmod_ctx)

static __inline__ int
fq_zech_ctx_fprint(FILE * file, const fq_zech_ctx_t ctx)
{
    int r;
    r = flint_fprintf(file, "Zech Representation:\n");
    if (r <= 0)
        return r;
    return fq_nmod_ctx_fprint(file, ctx->fq_nmod_ctx);
}

static __inline__ void
fq_zech_ctx_print(const fq_zech_ctx_t ctx)
{
    fq_zech_ctx_fprint(stdout, ctx);
}

/* Memory managment  *********************************************************/

static __inline__ void
fq_zech_init(fq_zech_t rop, const fq_zech_ctx_t ctx)
{
    rop->value = ctx->qm1;
}

static __inline__ void
fq_zech_init2(fq_zech_t rop, const fq_zech_ctx_t ctx)
{
    rop->value = ctx->qm1;
}

static __inline__ void
fq_zech_clear(fq_zech_t rop, const fq_zech_ctx_t ctx)
{
}

static __inline__ void
fq_zech_reduce(fq_zech_t rop, const fq_zech_ctx_t ctx)
{
    mp_limb_t order = fq_zech_ctx_order_ui(ctx);
    if (rop->value >= order)
    {
        rop->value -= order;
    }
}

/* Basic arithmetic **********************************************************/

void fq_zech_add(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2,
                 const fq_zech_ctx_t ctx);

void fq_zech_sub(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2,
                 const fq_zech_ctx_t ctx);

void fq_zech_sub_one(fq_zech_t rop, const fq_zech_t op1,
                     const fq_zech_ctx_t ctx);

void fq_zech_neg(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx);

void fq_zech_mul(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2,
                 const fq_zech_ctx_t ctx);

void fq_zech_mul_fmpz(fq_zech_t rop, const fq_zech_t op, const fmpz_t x,
                      const fq_zech_ctx_t ctx);

void fq_zech_mul_si(fq_zech_t rop, const fq_zech_t op, slong x,
                    const fq_zech_ctx_t ctx);

void fq_zech_mul_ui(fq_zech_t rop, const fq_zech_t op, ulong x,
                    const fq_zech_ctx_t ctx);

void fq_zech_sqr(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx);

void fq_zech_inv(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx);

void _fq_zech_pow(fmpz * rop, const fmpz * op, slong len, const fmpz_t e,
                  const fmpz * a, const slong *j, slong lena, const fmpz_t p);

void fq_zech_pow(fq_zech_t rop, const fq_zech_t op1, const fmpz_t e,
                 const fq_zech_ctx_t ctx);

void fq_zech_pow_ui(fq_zech_t rop, const fq_zech_t op1, const ulong e,
                    const fq_zech_ctx_t ctx);

void
fq_zech_pth_root(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx);

/* Randomisation *************************************************************/

void fq_zech_randtest(fq_zech_t rop, flint_rand_t state,
                      const fq_zech_ctx_t ctx);

void fq_zech_randtest_not_zero(fq_zech_t rop, flint_rand_t state,
                               const fq_zech_ctx_t ctx);

/* Comparison ****************************************************************/

static __inline__ int
fq_zech_equal(const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
{
    return op1->value == op2->value;
}

static __inline__ int
fq_zech_is_zero(const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    return op->value == ctx->qm1;
}

static __inline__ int
fq_zech_is_one(const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    return op->value == 0;
}

/* Assignments and conversions ***********************************************/

static __inline__ void
fq_zech_set(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    rop->value = op->value;
}

void
fq_zech_set_fmpz(fq_zech_t rop, const fmpz_t x, const fq_zech_ctx_t ctx);

static __inline__ void
fq_zech_set_ui(fq_zech_t rop, const ulong x, const fq_zech_ctx_t ctx)
{
    fmpz_t xx;
    fmpz_init_set_ui(xx, x);
    fq_zech_set_fmpz(rop, xx, ctx);
    fmpz_clear(xx);
}

static __inline__ void
fq_zech_swap(fq_zech_t op1, fq_zech_t op2, const fq_zech_ctx_t ctx)
{
    slong temp;
    temp = op2->value;
    op2->value = op1->value;
    op1->value = temp;
}

static __inline__ void
fq_zech_zero(fq_zech_t rop, const fq_zech_ctx_t ctx)
{
    rop->value = ctx->qm1;
}

static __inline__ void
fq_zech_one(fq_zech_t rop, const fq_zech_ctx_t ctx)
{
    rop->value = 0;
}

static __inline__ void
fq_zech_gen(fq_zech_t rop, const fq_zech_ctx_t ctx)
{
    rop->value = 1;
}

void
fq_zech_set_fq_nmod(fq_zech_t rop, const fq_nmod_t op, const fq_zech_ctx_t ctx);

void
fq_zech_get_fq_nmod(fq_nmod_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx);


/* Output ********************************************************************/
static __inline__ int
fq_zech_fprint_pretty(FILE * file, const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    return flint_fprintf(file, "%s^%wd", ctx->fq_nmod_ctx->var, op->value);
}

static __inline__ void
fq_zech_print_pretty(const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    fq_zech_fprint_pretty(stdout, op, ctx);
}

static __inline__ int
fq_zech_fprint(FILE * file, const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    return flint_fprintf(file, "%wd", op->value);
}

static __inline__ void
fq_zech_print(const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    fq_zech_fprint(stdout, op, ctx);
}

char *
fq_zech_get_str(const fq_zech_t op, const fq_zech_ctx_t ctx);

char *
fq_zech_get_str_pretty(const fq_zech_t op, const fq_zech_ctx_t ctx);


/* Special functions *********************************************************/

void fq_zech_trace(fmpz_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx);

void fq_zech_frobenius(fq_zech_t rop, const fq_zech_t op, slong e,
                       const fq_zech_ctx_t ctx);

void fq_zech_norm(fmpz_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx);

/* Bit packing ******************************************************/

void
fq_zech_bit_pack(fmpz_t f, const fq_zech_t op, mp_bitcnt_t bit_size,
                 const fq_zech_ctx_t ctx);

void
fq_zech_bit_unpack(fq_zech_t rop, const fmpz_t f, mp_bitcnt_t bit_size,
                   const fq_zech_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
