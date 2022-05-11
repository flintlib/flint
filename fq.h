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

#include "fq_mini.h"

/* Data types and context ****************************************************/

#ifdef __cplusplus
extern "C" {
#endif

FLINT_DLL void fq_ctx_init(fq_ctx_t ctx, const fmpz_t p, slong d,
                                                              const char *var);

FLINT_DLL int _fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d,
                                                              const char *var);

FLINT_DLL void fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d,
                                                              const char *var);

FLINT_DLL void fq_ctx_init_modulus(fq_ctx_t ctx, const fmpz_mod_poly_t modulus,
                                   const fmpz_mod_ctx_t ctxp, const char *var);

FLINT_DLL void fq_ctx_randtest(fq_ctx_t ctx, flint_rand_t state);

FLINT_DLL void fq_ctx_randtest_reducible(fq_ctx_t ctx, flint_rand_t state);

FLINT_DLL void fq_ctx_clear(fq_ctx_t ctx);

#if defined (H_STDIO)               \
  || defined (_H_STDIO)             \
  || defined (_STDIO_H)             \
  || defined (_STDIO_H_)            \
  || defined (__STDIO_H)            \
  || defined (__STDIO_H__)          \
  || defined (_STDIO_INCLUDED)      \
  || defined (__dj_include_stdio_h_)\
  || defined (_FILE_DEFINED)        \
  || defined (__STDIO__)            \
  || defined (_MSL_STDIO_H)         \
  || defined (_STDIO_H_INCLUDED)    \
  || defined (_ISO_STDIO_ISO_H)     \
  || defined (__STDIO_LOADED)       \
  || defined (_STDIO)

FLINT_DLL int fq_ctx_fprint(FILE * file, const fq_ctx_t ctx);

FQ_INLINE void fq_ctx_print(const fq_ctx_t ctx)
{
    fq_ctx_fprint(stdout, ctx);
}

#endif

/* Memory managment  *********************************************************/

FQ_INLINE void fq_init(fq_t rop, const fq_ctx_t ctx)
{
    fmpz_poly_init(rop);
}

FQ_INLINE void fq_init2(fq_t rop, const fq_ctx_t ctx)
{
    fmpz_poly_init2(rop, fq_ctx_degree(ctx));
}

FQ_INLINE void fq_clear(fq_t rop, const fq_ctx_t ctx)
{
    fmpz_poly_clear(rop);
}

/* Basic arithmetic **********************************************************/

FLINT_DLL void fq_add(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx);
FLINT_DLL void fq_sub(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx);

FLINT_DLL void fq_sub_one(fq_t rop, const fq_t op1, const fq_ctx_t ctx);

FLINT_DLL void fq_neg(fq_t rop, const fq_t op1, const fq_ctx_t ctx);

FLINT_DLL void fq_mul(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx);
FLINT_DLL void fq_mul_fmpz(fq_t rop, const fq_t op, const fmpz_t x, const fq_ctx_t ctx);
FLINT_DLL void fq_mul_si(fq_t rop, const fq_t op, slong x, const fq_ctx_t ctx);
FLINT_DLL void fq_mul_ui(fq_t rop, const fq_t op, ulong x, const fq_ctx_t ctx);

FLINT_DLL void fq_sqr(fq_t rop, const fq_t op, const fq_ctx_t ctx);

FLINT_DLL void _fq_pow(fmpz *rop, const fmpz *op, slong len, const fmpz_t e, const fq_ctx_t ctx);
FLINT_DLL void fq_pow(fq_t rop, const fq_t op1, const fmpz_t e, const fq_ctx_t ctx);
FLINT_DLL void fq_pow_ui(fq_t rop, const fq_t op, const ulong e, const fq_ctx_t ctx);

/* Roots *********************************************************************/

FLINT_DLL int fq_sqrt(fq_t rop, const fq_t op, const fq_ctx_t ctx);

FLINT_DLL void fq_pth_root(fq_t rop, const fq_t op1, const fq_ctx_t ctx);

FLINT_DLL int fq_is_square(const fq_t op, const fq_ctx_t ctx);

/* Randomisation *************************************************************/

FLINT_DLL void fq_randtest(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);

FLINT_DLL void fq_randtest_dense(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);

FLINT_DLL void fq_randtest_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);

FLINT_DLL void fq_rand(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);

FLINT_DLL void fq_rand_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx);

/* Assignments and conversions ***********************************************/

FQ_INLINE void fq_set_ui(fq_t rop, const ulong x, const fq_ctx_t ctx)
{
    fmpz_poly_set_ui(rop, x);
    fq_reduce(rop, ctx);
}

FQ_INLINE void fq_set_si(fq_t rop, const slong x, const fq_ctx_t ctx)
{
    fmpz_poly_set_si(rop, x);
    fq_reduce(rop, ctx);
}

FQ_INLINE void fq_one(fq_t rop, const fq_ctx_t ctx)
{
    fmpz_poly_one(rop);
}

FQ_INLINE void fq_gen(fq_t rop, const fq_ctx_t ctx)
{
    if (ctx->modulus->length == 2)
    {
        fmpz_poly_fit_length(rop, 1);
        fmpz_invmod(rop->coeffs, ctx->modulus->coeffs + 1, fq_ctx_prime(ctx));
        fmpz_neg(rop->coeffs, rop->coeffs);
        fmpz_mul(rop->coeffs, rop->coeffs, ctx->modulus->coeffs);
        fmpz_mod(rop->coeffs, rop->coeffs, fq_ctx_prime(ctx));
        _fmpz_poly_set_length(rop, !fmpz_is_zero(rop->coeffs));
    }
    else
    {
        fmpz_poly_zero(rop);
        fmpz_poly_set_coeff_ui(rop, 0, 0);
        fmpz_poly_set_coeff_ui(rop, 1, 1);
    }
}

FLINT_DLL int fq_get_fmpz(fmpz_t a, const fq_t b, const fq_ctx_t ctx);
FLINT_DLL void fq_get_fmpz_poly(fmpz_poly_t a, const fq_t b, const fq_ctx_t ctx);
FLINT_DLL void fq_get_fmpz_mod_poly(fmpz_mod_poly_t a, const fq_t b, const fq_ctx_t ctx);

FLINT_DLL void fq_set_fmpz_poly(fq_t a, const fmpz_poly_t b, const fq_ctx_t ctx);
FLINT_DLL void fq_set_fmpz_mod_poly(fq_t a, const fmpz_mod_poly_t b, const fq_ctx_t ctx);

/* Output ********************************************************************/

#if defined (FILE)                  \
  || defined (H_STDIO)              \
  || defined (_H_STDIO)             \
  || defined (_STDIO_H)             \
  || defined (_STDIO_H_)            \
  || defined (__STDIO_H)            \
  || defined (__STDIO_H__)          \
  || defined (_STDIO_INCLUDED)      \
  || defined (__dj_include_stdio_h_)\
  || defined (_FILE_DEFINED)        \
  || defined (__STDIO__)            \
  || defined (_MSL_STDIO_H)         \
  || defined (_STDIO_H_INCLUDED)    \
  || defined (_ISO_STDIO_ISO_H)     \
  || defined (__STDIO_LOADED)       \
  || defined (_STDIO)               \
  || defined (__DEFINED_FILE)

#include "fmpz_poly.h"

FQ_INLINE
int fq_fprint(FILE * file, const fq_t op, const fq_ctx_t ctx)
{
    return fmpz_poly_fprint(file, op);
}

FQ_INLINE
void fq_print(const fq_t op, const fq_ctx_t ctx)
{
    fmpz_poly_print(op);
}

FQ_INLINE 
int fq_fprint_pretty(FILE * file, const fq_t op, const fq_ctx_t ctx)
{
    return fmpz_poly_fprint_pretty(file, op, ctx->var);
}

FQ_INLINE 
int fq_print_pretty(const fq_t op, const fq_ctx_t ctx)
{
    return fmpz_poly_print_pretty(op, ctx->var);
}

FLINT_DLL char * fq_get_str(const fq_t op, const fq_ctx_t ctx);

FLINT_DLL char * fq_get_str_pretty(const fq_t op, const fq_ctx_t ctx);

#endif

/* Special functions *********************************************************/

FLINT_DLL void _fq_trace(fmpz_t rop, const fmpz *op, slong len, const fq_ctx_t ctx);
FLINT_DLL void fq_trace(fmpz_t rop, const fq_t op, const fq_ctx_t ctx);

FLINT_DLL void _fq_frobenius(fmpz *rop, const fmpz *op, slong len, slong e, const fq_ctx_t ctx);
FLINT_DLL void fq_frobenius(fq_t rop, const fq_t op, slong e, const fq_ctx_t ctx);

FLINT_DLL void _fq_norm(fmpz_t rop, const fmpz *op, slong len, const fq_ctx_t ctx);
FLINT_DLL void fq_norm(fmpz_t rop, const fq_t op, const fq_ctx_t ctx);

/* Bit packing ******************************************************/

FLINT_DLL void fq_bit_pack(fmpz_t f, const fq_t op, flint_bitcnt_t bit_size, const fq_ctx_t ctx);

FLINT_DLL void fq_bit_unpack(fq_t rop, const fmpz_t f, flint_bitcnt_t bit_size, const fq_ctx_t ctx);

/* Inlines *******************************************************************/

FLINT_DLL void __fq_ctx_prime(fmpz_t p, fq_ctx_t ctx);

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
