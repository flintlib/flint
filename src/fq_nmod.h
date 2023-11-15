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

#ifndef FQ_NMOD_H
#define FQ_NMOD_H

#ifdef FQ_NMOD_INLINES_C
#define FQ_NMOD_INLINE
#define FQ_TEMPLATES_INLINE
#else
#define FQ_NMOD_INLINE static inline
#define FQ_TEMPLATES_INLINE static inline
#endif

#include "fq_nmod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Context ********************************************************************/

void fq_nmod_ctx_init(fq_nmod_ctx_t ctx,
                      const fmpz_t p, slong d, const char *var);

int _fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx,
                             const fmpz_t p, slong d, const char *var);

void fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx,
                             const fmpz_t p, slong d, const char *var);

void fq_nmod_ctx_init_modulus(fq_nmod_ctx_t ctx,
                              const nmod_poly_t modulus,
                              const char *var);

void fq_nmod_ctx_randtest(fq_nmod_ctx_t ctx, flint_rand_t state);

void fq_nmod_ctx_randtest_reducible(fq_nmod_ctx_t ctx, flint_rand_t state);

void fq_nmod_ctx_clear(fq_nmod_ctx_t ctx);

FQ_NMOD_INLINE const nmod_poly_struct* fq_nmod_ctx_modulus(const fq_nmod_ctx_t ctx)
{
    return ctx->modulus;
}

FQ_NMOD_INLINE slong fq_nmod_ctx_degree(const fq_nmod_ctx_t ctx)
{
    return ctx->modulus->length - 1;
}

void fq_nmod_ctx_order(fmpz_t f, const fq_nmod_ctx_t ctx);

#ifdef FLINT_HAVE_FILE
int fq_nmod_ctx_fprint(FILE * file, const fq_nmod_ctx_t ctx);
#endif

void fq_nmod_ctx_print(const fq_nmod_ctx_t ctx);

/* Memory management  *********************************************************/

void fq_nmod_init(fq_nmod_t rop, const fq_nmod_ctx_t ctx);
void fq_nmod_init2(fq_nmod_t rop, const fq_nmod_ctx_t ctx);

void fq_nmod_clear(fq_nmod_t rop, const fq_nmod_ctx_t ctx);

void _fq_nmod_sparse_reduce(mp_limb_t *R, slong lenR, const fq_nmod_ctx_t ctx);
void _fq_nmod_dense_reduce(mp_limb_t* R, slong lenR, const fq_nmod_ctx_t ctx);
void _fq_nmod_reduce(mp_limb_t* R, slong lenR, const fq_nmod_ctx_t ctx);
void fq_nmod_reduce(fq_nmod_t rop, const fq_nmod_ctx_t ctx);

/* Basic arithmetic **********************************************************/

void fq_nmod_add(fq_nmod_t rop, const fq_nmod_t op1,
                                 const fq_nmod_t op2, const fq_nmod_ctx_t ctx);

void fq_nmod_sub(fq_nmod_t rop, const fq_nmod_t op1,
                                 const fq_nmod_t op2, const fq_nmod_ctx_t ctx);

void fq_nmod_sub_one(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx);

void fq_nmod_neg(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx);

void fq_nmod_mul(fq_nmod_t rop,
            const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx);

void fq_nmod_mul_fmpz(fq_nmod_t rop,
                  const fq_nmod_t op, const fmpz_t x, const fq_nmod_ctx_t ctx);

void fq_nmod_mul_si(fq_nmod_t rop,
                         const fq_nmod_t op, slong x, const fq_nmod_ctx_t ctx);

void fq_nmod_mul_ui(fq_nmod_t rop,
                         const fq_nmod_t op, ulong x, const fq_nmod_ctx_t ctx);

void fq_nmod_sqr(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx);

void fq_nmod_inv(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx);

void _fq_nmod_pow(mp_limb_t *rop, const mp_limb_t *op,
                           slong len, const fmpz_t e, const fq_nmod_ctx_t ctx);

void fq_nmod_pow(fq_nmod_t rop, const fq_nmod_t op1,
                                      const fmpz_t e, const fq_nmod_ctx_t ctx);

void fq_nmod_pow_ui(fq_nmod_t rop,
                  const fq_nmod_t op1, const ulong e, const fq_nmod_ctx_t ctx);

/* Roots ********************************************************************/

int fq_nmod_sqrt(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx);

void fq_nmod_pth_root(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx);

int fq_nmod_is_square(const fq_nmod_t op, const fq_nmod_ctx_t ctx);

/* Randomisation *************************************************************/

void fq_nmod_randtest(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx);

void fq_nmod_randtest_dense(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx);

void fq_nmod_randtest_not_zero(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx);

void fq_nmod_rand(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx);

void fq_nmod_rand_not_zero(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx);


/* Comparison ****************************************************************/

int fq_nmod_equal(const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx);
int fq_nmod_is_zero(const fq_nmod_t op, const fq_nmod_ctx_t ctx);
int fq_nmod_is_one(const fq_nmod_t op, const fq_nmod_ctx_t ctx);
int fq_nmod_cmp(const fq_nmod_t a, const fq_nmod_t b, const fq_nmod_ctx_t ctx);

/* Assignments and conversions ***********************************************/

void fq_nmod_set(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx);
void fq_nmod_set_si(fq_nmod_t rop, const slong x, const fq_nmod_ctx_t ctx);
void fq_nmod_set_ui(fq_nmod_t rop, const ulong x, const fq_nmod_ctx_t ctx);
void fq_nmod_set_fmpz(fq_nmod_t rop, const fmpz_t x, const fq_nmod_ctx_t ctx);
void fq_nmod_set_nmod_poly(fq_nmod_t a, const nmod_poly_t b, const fq_nmod_ctx_t ctx);

int fq_nmod_get_fmpz(fmpz_t a, const fq_nmod_t b, const fq_nmod_ctx_t ctx);
void fq_nmod_get_nmod_poly(nmod_poly_t a, const fq_nmod_t b, const fq_nmod_ctx_t ctx);

void fq_nmod_swap(fq_nmod_t op1, fq_nmod_t op2, const fq_nmod_ctx_t ctx);

void fq_nmod_zero(fq_nmod_t rop,  const fq_nmod_ctx_t ctx);
void fq_nmod_one(fq_nmod_t rop,  const fq_nmod_ctx_t ctx);

void fq_nmod_gen(fq_nmod_t rop, const fq_nmod_ctx_t ctx);

/* Output ********************************************************************/

#ifdef FLINT_HAVE_FILE
int fq_nmod_fprint(FILE * file, const fq_nmod_t op, const fq_nmod_ctx_t ctx);
int fq_nmod_fprint_pretty(FILE * file, const fq_nmod_t op, const fq_nmod_ctx_t ctx);
#endif

void fq_nmod_print(const fq_nmod_t op, const fq_nmod_ctx_t ctx);
void fq_nmod_print_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx);

char * fq_nmod_get_str(const fq_nmod_t op, const fq_nmod_ctx_t ctx);
char * fq_nmod_get_str_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx);

/* Special functions *********************************************************/

void _fq_nmod_trace(fmpz_t rop, const mp_limb_t *op, slong len,
                    const fq_nmod_ctx_t ctx);

void fq_nmod_trace(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx);

void _fq_nmod_frobenius(mp_limb_t *rop, const mp_limb_t *op, slong len, slong e,
                        const fq_nmod_ctx_t ctx);

void fq_nmod_frobenius(fq_nmod_t rop, const fq_nmod_t op, slong e, const fq_nmod_ctx_t ctx);

void _fq_nmod_norm(fmpz_t rop, const mp_limb_t *op, slong len,
                   const fq_nmod_ctx_t ctx);

void fq_nmod_norm(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx);

/* Bit packing ******************************************************/

void fq_nmod_bit_pack(fmpz_t f, const fq_nmod_t op, flint_bitcnt_t bit_size,
                 const fq_nmod_ctx_t ctx);

void fq_nmod_bit_unpack(fq_nmod_t rop, const fmpz_t f, flint_bitcnt_t bit_size,
                   const fq_nmod_ctx_t ctx);

/* Inlines *******************************************************************/

void __fq_nmod_ctx_prime(fmpz_t p, fq_nmod_ctx_t ctx);

#ifdef T
#undef T
#endif

#define T fq_nmod
#define CAP_T FQ_NMOD
#define B nmod
#include "fq_templates.h"
#undef B
#undef CAP_T
#undef T

#ifdef __cplusplus
}
#endif

#endif

