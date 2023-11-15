/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_ZECH_H
#define FQ_ZECH_H

#ifdef FQ_ZECH_INLINES_C
#define FQ_ZECH_INLINE
#ifndef FQ_TEMPLATES_INLINE
# define FQ_TEMPLATES_INLINE
#endif
#else
#define FQ_ZECH_INLINE static inline
#ifndef FQ_TEMPLATES_INLINE
# define FQ_TEMPLATES_INLINE static inline
#endif
#endif

#include "fq_zech_types.h"

/* Data types and context ****************************************************/
#ifdef __cplusplus
extern "C" {
#endif

void fq_zech_ctx_init(fq_zech_ctx_t ctx, const fmpz_t p, slong d, const char *var);

int fq_zech_ctx_init_fq_nmod_ctx_check(fq_zech_ctx_t ctx, fq_nmod_ctx_t ctxn);

void fq_zech_ctx_init_fq_nmod_ctx(fq_zech_ctx_t ctx, fq_nmod_ctx_t ctxn);

int _fq_zech_ctx_init_conway(fq_zech_ctx_t ctx, const fmpz_t p, slong d, const char *var);

void fq_zech_ctx_init_conway(fq_zech_ctx_t ctx, const fmpz_t p, slong d, const char *var);

void fq_zech_ctx_init_random(fq_zech_ctx_t ctx, const fmpz_t p, slong d, const char *var);

void fq_zech_ctx_init_modulus(fq_zech_ctx_t ctx,
                              const nmod_poly_t modulus,
                              const char *var);

int fq_zech_ctx_init_modulus_check(fq_zech_ctx_t ctx,
                              const nmod_poly_t modulus,
                              const char *var);

void fq_zech_ctx_randtest(fq_zech_ctx_t ctx, flint_rand_t state);

void fq_zech_ctx_randtest_reducible(fq_zech_ctx_t ctx, flint_rand_t state);

void fq_zech_ctx_clear(fq_zech_ctx_t ctx);

const nmod_poly_struct * fq_zech_ctx_modulus(const fq_zech_ctx_t ctx);

slong fq_zech_ctx_degree(const fq_zech_ctx_t ctx);

void fq_zech_ctx_order(fmpz_t f, const fq_zech_ctx_t ctx);

FQ_ZECH_INLINE mp_limb_t
fq_zech_ctx_order_ui(const fq_zech_ctx_t ctx)
{
    return ctx->qm1 + 1;
}

#ifdef FLINT_HAVE_FILE
int fq_zech_ctx_fprint(FILE * file, const fq_zech_ctx_t ctx);
#endif

void fq_zech_ctx_print(const fq_zech_ctx_t ctx);

/* Memory management  *********************************************************/

FQ_ZECH_INLINE void
fq_zech_init(fq_zech_t rop, const fq_zech_ctx_t ctx)
{
    rop->value = ctx->qm1;
}

FQ_ZECH_INLINE void
fq_zech_init2(fq_zech_t rop, const fq_zech_ctx_t ctx)
{
    rop->value = ctx->qm1;
}

FQ_ZECH_INLINE void
fq_zech_clear(fq_zech_t rop, const fq_zech_ctx_t ctx)
{
}

FQ_ZECH_INLINE void
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

/* Roots *********************************************************************/

int fq_zech_sqrt(fq_zech_t rop, const fq_zech_t op1,
		                                      const fq_zech_ctx_t ctx);

void fq_zech_pth_root(fq_zech_t rop,
		                 const fq_zech_t op1, const fq_zech_ctx_t ctx);

int fq_zech_is_square(const fq_zech_t op1, const fq_zech_ctx_t ctx);

/* Randomisation *************************************************************/

void fq_zech_randtest(fq_zech_t rop, flint_rand_t state,
                      const fq_zech_ctx_t ctx);

void fq_zech_randtest_not_zero(fq_zech_t rop, flint_rand_t state,
                               const fq_zech_ctx_t ctx);

void fq_zech_rand(fq_zech_t rop, flint_rand_t state,
                                                      const fq_zech_ctx_t ctx);

void fq_zech_rand_not_zero(fq_zech_t rop, flint_rand_t state,
                                                      const fq_zech_ctx_t ctx);

/* Comparison ****************************************************************/

FQ_ZECH_INLINE int
fq_zech_equal(const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
{
    return op1->value == op2->value;
}

FQ_ZECH_INLINE int
fq_zech_is_zero(const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    return op->value == ctx->qm1;
}

FQ_ZECH_INLINE int
fq_zech_is_one(const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    return op->value == 0;
}

/* Assignments and conversions ***********************************************/

FQ_ZECH_INLINE void
fq_zech_set(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    rop->value = op->value;
}

void fq_zech_set_fmpz(fq_zech_t rop, const fmpz_t x, const fq_zech_ctx_t ctx);

void fq_zech_set_si(fq_zech_t rop, const slong x, const fq_zech_ctx_t ctx);
void fq_zech_set_ui(fq_zech_t rop, const ulong x, const fq_zech_ctx_t ctx);

FQ_ZECH_INLINE void
fq_zech_swap(fq_zech_t op1, fq_zech_t op2, const fq_zech_ctx_t ctx)
{
    slong temp;
    temp = op2->value;
    op2->value = op1->value;
    op1->value = temp;
}

FQ_ZECH_INLINE void
fq_zech_zero(fq_zech_t rop, const fq_zech_ctx_t ctx)
{
    rop->value = ctx->qm1;
}

FQ_ZECH_INLINE void
fq_zech_one(fq_zech_t rop, const fq_zech_ctx_t ctx)
{
    rop->value = 0;
}

FQ_ZECH_INLINE void
fq_zech_gen(fq_zech_t rop, const fq_zech_ctx_t ctx)
{
    rop->value = 1;
}


int fq_zech_get_fmpz(fmpz_t a, const fq_zech_t op,
                                                      const fq_zech_ctx_t ctx);

void fq_zech_set_fq_nmod(fq_zech_t rop, const fq_nmod_t op,
                                                      const fq_zech_ctx_t ctx);

void fq_zech_get_fq_nmod(fq_nmod_t rop, const fq_zech_t op,
                                                      const fq_zech_ctx_t ctx);

void fq_zech_get_nmod_poly(nmod_poly_t a, const fq_zech_t b,
                                                      const fq_zech_ctx_t ctx);

void fq_zech_set_nmod_poly(fq_zech_t a, const nmod_poly_t b,
                                                      const fq_zech_ctx_t ctx);


/* Output ********************************************************************/

#ifdef FLINT_HAVE_FILE
int fq_zech_fprint(FILE * file, const fq_zech_t op, const fq_zech_ctx_t ctx);
int fq_zech_fprint_pretty(FILE * file, const fq_zech_t op, const fq_zech_ctx_t ctx);
#endif

void fq_zech_print(const fq_zech_t op, const fq_zech_ctx_t ctx);
void fq_zech_print_pretty(const fq_zech_t op, const fq_zech_ctx_t ctx);

char * fq_zech_get_str(const fq_zech_t op, const fq_zech_ctx_t ctx);
char * fq_zech_get_str_pretty(const fq_zech_t op, const fq_zech_ctx_t ctx);


/* Special functions *********************************************************/

void fq_zech_trace(fmpz_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx);

void fq_zech_frobenius(fq_zech_t rop, const fq_zech_t op, slong e,
                       const fq_zech_ctx_t ctx);

void fq_zech_norm(fmpz_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx);

/* Bit packing ******************************************************/

void fq_zech_bit_pack(fmpz_t f, const fq_zech_t op, flint_bitcnt_t bit_size,
                 const fq_zech_ctx_t ctx);

void fq_zech_bit_unpack(fq_zech_t rop, const fmpz_t f, flint_bitcnt_t bit_size,
                   const fq_zech_ctx_t ctx);

/* Inlines *******************************************************************/

void __fq_zech_ctx_prime(fmpz_t p, fq_zech_ctx_t ctx);

#ifdef T
#undef T
#endif

#define T fq_zech
#define CAP_T FQ_ZECH
#define B nmod
#include "fq_templates.h"
#undef B
#undef CAP_T
#undef T

#ifdef __cplusplus
}
#endif

#endif
