/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_H
#define FQ_H

#ifdef FQ_INLINES_C
#define FQ_INLINE
#define FQ_TEMPLATES_INLINE
#else
#define FQ_INLINE static inline
#define FQ_TEMPLATES_INLINE static inline
#endif

#include "fq_types.h"

/* Data types and context ****************************************************/

#ifdef __cplusplus
extern "C" {
#endif

void fq_ctx_init(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var);
int _fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var);
void fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var);
void fq_ctx_init_modulus(fq_ctx_t ctx, const fmpz_mod_poly_t modulus, const fmpz_mod_ctx_t ctxp, const char *var);

void fq_ctx_clear(fq_ctx_t ctx);

void fq_ctx_randtest(fq_ctx_t ctx, flint_rand_t state);
void fq_ctx_randtest_reducible(fq_ctx_t ctx, flint_rand_t state);

FQ_INLINE const fmpz_mod_poly_struct* fq_ctx_modulus(const fq_ctx_t ctx)
{
    return ctx->modulus;
}

FQ_INLINE slong fq_ctx_degree(const fq_ctx_t ctx)
{
    return ctx->modulus->length - 1;
}

const fmpz * fq_ctx_prime(const fq_ctx_t ctx);

void fq_ctx_order(fmpz_t f, const fq_ctx_t ctx);

#ifdef FLINT_HAVE_FILE
int fq_ctx_fprint(FILE * file, const fq_ctx_t ctx);
#endif

void fq_ctx_print(const fq_ctx_t ctx);

/* Memory management  *********************************************************/

void fq_init(fq_t rop, const fq_ctx_t ctx);
void fq_init2(fq_t rop, const fq_ctx_t ctx);

void fq_clear(fq_t rop, const fq_ctx_t ctx);

void _fq_sparse_reduce(fmpz * R, slong lenR, const fq_ctx_t ctx);
void _fq_dense_reduce(fmpz * R, slong lenR, const fq_ctx_t ctx);
void _fq_reduce(fmpz * R, slong lenR, const fq_ctx_t ctx);
void fq_reduce(fq_t rop, const fq_ctx_t ctx);

/* Basic arithmetic **********************************************************/

void fq_add(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx);
void fq_sub(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx);

void fq_sub_one(fq_t rop, const fq_t op1, const fq_ctx_t ctx);

void fq_neg(fq_t rop, const fq_t op1, const fq_ctx_t ctx);

void fq_mul(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx);
void fq_mul_fmpz(fq_t rop, const fq_t op, const fmpz_t x, const fq_ctx_t ctx);
void fq_mul_si(fq_t rop, const fq_t op, slong x, const fq_ctx_t ctx);
void fq_mul_ui(fq_t rop, const fq_t op, ulong x, const fq_ctx_t ctx);

void fq_sqr(fq_t rop, const fq_t op, const fq_ctx_t ctx);

void fq_inv(fq_t rop, const fq_t op1, const fq_ctx_t ctx);

void _fq_pow(fmpz *rop, const fmpz *op, slong len, const fmpz_t e, const fq_ctx_t ctx);
void fq_pow(fq_t rop, const fq_t op1, const fmpz_t e, const fq_ctx_t ctx);
void fq_pow_ui(fq_t rop, const fq_t op, const ulong e, const fq_ctx_t ctx);

/* Roots *********************************************************************/

int fq_sqrt(fq_t rop, const fq_t op, const fq_ctx_t ctx);

void fq_pth_root(fq_t rop, const fq_t op1, const fq_ctx_t ctx);

int fq_is_square(const fq_t op, const fq_ctx_t ctx);

/* Randomisation *************************************************************/

void fq_randtest(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);
void fq_randtest_dense(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);
void fq_randtest_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);

void fq_rand(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);
void fq_rand_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);

/* Comparison ****************************************************************/

int fq_equal(const fq_t op1, const fq_t op2, const fq_ctx_t ctx);

int fq_is_zero(const fq_t op, const fq_ctx_t ctx);
int fq_is_one(const fq_t op, const fq_ctx_t ctx);

/* Assignments and conversions ***********************************************/

void fq_set(fq_t rop, const fq_t op, const fq_ctx_t ctx);
void fq_set_ui(fq_t rop, const ulong x, const fq_ctx_t ctx);
void fq_set_si(fq_t rop, const slong x, const fq_ctx_t ctx);
void fq_set_fmpz(fq_t rop, const fmpz_t x, const fq_ctx_t ctx);
void fq_set_fmpz_poly(fq_t a, const fmpz_poly_t b, const fq_ctx_t ctx);
void fq_set_fmpz_mod_poly(fq_t a, const fmpz_mod_poly_t b, const fq_ctx_t ctx);

int fq_get_fmpz(fmpz_t a, const fq_t b, const fq_ctx_t ctx);
void fq_get_fmpz_poly(fmpz_poly_t a, const fq_t b, const fq_ctx_t ctx);
void fq_get_fmpz_mod_poly(fmpz_mod_poly_t a, const fq_t b, const fq_ctx_t ctx);

void fq_zero(fq_t rop, const fq_ctx_t ctx);
void fq_one(fq_t rop, const fq_ctx_t ctx);

void fq_swap(fq_t op1, fq_t op2, const fq_ctx_t ctx);

void fq_gen(fq_t rop, const fq_ctx_t ctx);

/* Output ********************************************************************/

#ifdef FLINT_HAVE_FILE
int fq_fprint(FILE * file, const fq_t op, const fq_ctx_t ctx);
void fq_print(const fq_t op, const fq_ctx_t ctx);
int fq_fprint_pretty(FILE * file, const fq_t op, const fq_ctx_t ctx);
int fq_print_pretty(const fq_t op, const fq_ctx_t ctx);
#endif

char * fq_get_str(const fq_t op, const fq_ctx_t ctx);
char * fq_get_str_pretty(const fq_t op, const fq_ctx_t ctx);

/* Special functions *********************************************************/

void _fq_trace(fmpz_t rop, const fmpz *op, slong len, const fq_ctx_t ctx);
void fq_trace(fmpz_t rop, const fq_t op, const fq_ctx_t ctx);

void _fq_frobenius(fmpz *rop, const fmpz *op, slong len, slong e, const fq_ctx_t ctx);
void fq_frobenius(fq_t rop, const fq_t op, slong e, const fq_ctx_t ctx);

void _fq_norm(fmpz_t rop, const fmpz *op, slong len, const fq_ctx_t ctx);
void fq_norm(fmpz_t rop, const fq_t op, const fq_ctx_t ctx);


/* Bit packing ******************************************************/

void fq_bit_pack(fmpz_t f, const fq_t op, flint_bitcnt_t bit_size, const fq_ctx_t ctx);

void fq_bit_unpack(fq_t rop, const fmpz_t f, flint_bitcnt_t bit_size, const fq_ctx_t ctx);

/* Inlines *******************************************************************/

void __fq_ctx_prime(fmpz_t p, fq_ctx_t ctx);

#ifdef T
#undef T
#endif

#define T fq
#define CAP_T FQ
#define B fmpz_mod
#include "fq_templates.h"
#undef B
#undef CAP_T
#undef T

#ifdef __cplusplus
}
#endif

#endif

