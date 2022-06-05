/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define ARITH_INLINES_C
#define DOUBLE_EXTRAS_INLINES_C
#define D_MAT_INLINES_C
#define FFT_INLINES_C
#define FMPQ_INLINES_C
#define FMPQ_MAT_INLINES_C
#define FMPQ_MPOLY_FACTOR_INLINES_C
#define FMPQ_MPOLY_INLINES_C
#define FMPQ_POLY_INLINES_C
#define FMPZ_FACTOR_INLINES_C
#define FMPZ_INLINES_C
#define FMPZ_MAT_INLINES_C
#define FMPZ_MOD_INLINES_C
#define FMPZ_MOD_MAT_INLINES_C
#define FMPZ_MOD_POLY_FACTOR_INLINES_C
#define FMPZ_MOD_POLY_INLINES_C
#define FMPZ_MPOLY_FACTOR_INLINES_C
#define FMPZ_MPOLY_INLINES_C
#define FMPZ_POLY_FACTOR_INLINES_C
#define FMPZ_POLY_INLINES_C
#define FMPZ_POLY_MAT_INLINES_C
#define FMPZ_POLY_Q_INLINES_C
#define FMPZ_VEC_INLINES_C
#define FQ_DEFAULT_INLINES_C
#define FQ_DEFAULT_MAT_INLINES_C
#define FQ_DEFAULT_POLY_FACTOR_INLINES_C
#define FQ_DEFAULT_POLY_INLINES_C
#define FQ_INLINES_C
#define FQ_MAT_INLINES_C
#define FQ_NMOD_INLINES_C
#define FQ_NMOD_MAT_INLINES_C
#define FQ_NMOD_MPOLY_FACTOR_INLINES_C
#define FQ_NMOD_MPOLY_INLINES_C
#define FQ_NMOD_POLY_FACTOR_INLINES_C
#define FQ_NMOD_POLY_INLINES_C
#define FQ_NMOD_VEC_INLINES_C
#define FQ_POLY_FACTOR_INLINES_C
#define FQ_POLY_INLINES_C
#define FQ_VEC_INLINES_C
#define FQ_ZECH_INLINES_C
#define FQ_ZECH_MAT_INLINES_C
#define FQ_ZECH_POLY_FACTOR_INLINES_C
#define FQ_ZECH_POLY_INLINES_C
#define FQ_ZECH_VEC_INLINES_C
#define LONG_EXTRAS_INLINES_C
#define MPF_MAT_INLINES_C
#define MPOLY_INLINES_C
#define NMOD_INLINES_C
#define NMOD_MAT_INLINES_C
#define NMOD_MPOLY_FACTOR_INLINES_C
#define NMOD_MPOLY_INLINES_C
#define NMOD_POLY_FACTOR_INLINES_C
#define NMOD_POLY_INLINES_C
#define NMOD_POLY_MAT_INLINES_C
#define NMOD_VEC_INLINES_C
#define PADIC_INLINES_C
#define PADIC_MAT_INLINES_C
#define PADIC_POLY_INLINES_C
#define PERM_INLINES_C
#define QADIC_INLINES_C
#define ULONG_EXTRAS_INLINES_C

#include <stdio.h> /* Ensure that FLINT_HAVE_FILE */
#include <stdarg.h> /* Ensure that FLINT_HAVE_VA_LIST */
#include "arith.h"
#include "double_extras.h"
#include "d_mat.h"
#include "fft.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "fmpq_mpoly_factor.h"
#include "fmpq_mpoly.h"
#include "fmpq_poly.h"
#include "fmpz_factor.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"
#include "fmpz_mod_poly_factor.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mpoly_factor.h"
#include "fmpz_mpoly.h"
#include "fmpz_poly_factor.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"
#include "fmpz_poly_q.h"
#include "fmpz_vec.h"
#include "fq_default.h"
#include "fq_default_mat.h"
#include "fq_default_poly_factor.h"
#include "fq_default_poly.h"
#include "fq.h"
#include "fq_mat.h"
#include "fq_nmod.h"
#include "fq_nmod_mat.h"
#include "fq_nmod_mpoly_factor.h"
#include "fq_nmod_mpoly.h"
#include "fq_nmod_poly_factor.h"
#include "fq_nmod_poly.h"
#include "fq_nmod_vec.h"
#include "fq_poly_factor.h"
#include "fq_poly.h"
#include "fq_vec.h"
#include "fq_zech.h"
#include "fq_zech_mat.h"
#include "fq_zech_poly_factor.h"
#include "fq_zech_poly.h"
#include "fq_zech_vec.h"
#include "long_extras.h"
#include "mpf_mat.h"
#include "mpoly.h"
#include "nmod.h"
#include "nmod_mat.h"
#include "nmod_mpoly_factor.h"
#include "nmod_mpoly.h"
#include "nmod_poly_factor.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "nmod_vec.h"
#include "padic.h"
#include "padic_mat.h"
#include "padic_poly.h"
#include "perm.h"
#include "qadic.h"
#include "ulong_extras.h"

/* fmpq */
void fmpq_numerator(fmpz_t n, const fmpq_t q)
{
   fmpz_set(n, fmpq_numref(q));
}

void fmpq_denominator(fmpz_t n, const fmpq_t q)
{
   fmpz_set(n, fmpq_denref(q));
}

fmpz * fmpq_numerator_ptr(fmpq_t q)
{
   return fmpq_numref(q);
}

fmpz * fmpq_denominator_ptr(fmpq_t q)
{
   return fmpq_denref(q);
}

int fmpq_equal_fmpz(fmpq_t q, fmpz_t n)
{
   return fmpz_equal(fmpq_numref(q), n) && q->den == WORD(1);
}

/* fmpz_factor */
void fmpz_factor_get_fmpz(fmpz_t z, const fmpz_factor_t factor, slong i)
{
    fmpz_set(z, factor->p + i);
}
 
/* fmpz */
fmpz * __new_fmpz()
{
    return flint_calloc(sizeof(fmpz), 1);
}

void __free_fmpz(fmpz * f)
{
   _fmpz_demote(f);
   flint_free(f);
}   

void __fmpz_set_si(fmpz_t f, slong val)
{
    if (val < COEFF_MIN || val > COEFF_MAX) /* val is large */
    {
        __mpz_struct *mpz_coeff = _fmpz_promote(f);
        flint_mpz_set_si(mpz_coeff, val);
    }
    else
    {
        _fmpz_demote(f);
        *f = val;               /* val is small */
    }
}

void __fmpz_set_ui(fmpz_t f, ulong val)
{
    if (val > COEFF_MAX)        /* val is large */
    {
        __mpz_struct *mpz_coeff = _fmpz_promote(f);
        flint_mpz_set_ui(mpz_coeff, val);
    }
    else
    {
        _fmpz_demote(f);
        *f = val;               /* val is small */
    }
}

void __fmpz_init(fmpz_t f)
{
	(*f) = WORD(0);
}

void __fmpz_init_set_ui(fmpz_t f, ulong g)
{
    if (g <= COEFF_MAX)
    {
        *f = g;
    }
    else
    {
        __mpz_struct *ptr;

        ptr = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(ptr);
        flint_mpz_set_ui(ptr, g);
    }
}

void __fmpz_clear(fmpz_t f)
{
	_fmpz_demote(f);
}

int __fmpz_lt(fmpz_t f, fmpz_t g)
{
   return fmpz_cmp(f, g) < 0;
}

int __fmpz_gt(fmpz_t f, fmpz_t g)
{
   return fmpz_cmp(f, g) > 0;
}

int __fmpz_lte(fmpz_t f, fmpz_t g)
{
   return fmpz_cmp(f, g) <= 0;
}

int __fmpz_gte(fmpz_t f, fmpz_t g)
{
   return fmpz_cmp(f, g) >= 0;
}

int __fmpz_eq(fmpz_t f, fmpz_t g)
{
   return fmpz_cmp(f, g) == 0;
}

int __fmpz_neq(fmpz_t f, fmpz_t g)
{
   return fmpz_cmp(f, g) != 0;
}

void __fmpz_init_set(fmpz_t f, const fmpz_t g)
{
    if (!COEFF_IS_MPZ(*g))
    {
        *f = *g;
    }
    else
    {
        __mpz_struct *ptr;

        ptr = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(ptr);
        mpz_set(ptr, COEFF_TO_PTR(*g));
    }
}

void __fmpz_neg(fmpz_t f1, const fmpz_t f2)
{
    if (!COEFF_IS_MPZ(*f2))     /* coeff is small */
    {
        fmpz t = -*f2;          /* Need to save value in case of aliasing */
        _fmpz_demote(f1);
        *f1 = t;
    }
    else                        /* coeff is large */
    {
        /* No need to retain value in promotion, as if aliased, both already large */
        __mpz_struct * mf1 = _fmpz_promote(f1);
        mpz_neg(mf1, COEFF_TO_PTR(*f2));
    }
}

/* fmpz_mod_mat */
void fmpz_mod_mat_get_entry(fmpz_t x, const fmpz_mod_mat_t mat, slong i, slong j)
{
  fmpz_set(x, fmpz_mod_mat_entry(mat, i, j));
}

/* fmpz_mod_poly_factor */
void fmpz_mod_poly_factor_get_fmpz_mod_poly(fmpz_mod_poly_t z,
                 fmpz_mod_poly_factor_t fac, slong i, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_set(z, fac->poly + i, ctx);
}

/* fmpz_mod_poly */
void fmpz_mod_poly_add_si(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                             slong c, const fmpz_mod_ctx_t ctx)
{
   fmpz_t d;

   fmpz_init(d);
   fmpz_set_si(d, c);

   if (fmpz_size(fmpz_mod_ctx_modulus(ctx)) > 1)
   {
      if (c < 0)
         fmpz_add(d, d, fmpz_mod_ctx_modulus(ctx));
   } else
      fmpz_mod(d, d, fmpz_mod_ctx_modulus(ctx));
      
   if (poly->length == 0)
      fmpz_mod_poly_set_fmpz(res, d, ctx);
   else
   {
      fmpz_mod_poly_set(res, poly, ctx);

      fmpz_add(res->coeffs + 0, res->coeffs + 0, d);
      if (fmpz_cmp(res->coeffs + 0, fmpz_mod_ctx_modulus(ctx)) >= 0)
         fmpz_sub(res->coeffs + 0, res->coeffs + 0, fmpz_mod_ctx_modulus(ctx));

      _fmpz_mod_poly_normalise(res);
   }

   fmpz_clear(d);
}

void fmpz_mod_poly_sub_si(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                             slong c, const fmpz_mod_ctx_t ctx)
{
   fmpz_t d;

   fmpz_init(d);
   fmpz_set_si(d, c);

   if (fmpz_size(fmpz_mod_ctx_modulus(ctx)) > 1)
   {
      if (c < 0)
         fmpz_add(d, d, fmpz_mod_ctx_modulus(ctx));
   } else
      fmpz_mod(d, d, fmpz_mod_ctx_modulus(ctx));
      
   if (poly->length == 0)
   {
      fmpz_sub(d, fmpz_mod_ctx_modulus(ctx), d);
      if (fmpz_cmp(d, fmpz_mod_ctx_modulus(ctx)) == 0)
         fmpz_zero(d);
      fmpz_mod_poly_set_fmpz(res, d, ctx);
   }
   else
   {
      fmpz_mod_poly_set(res, poly, ctx);

      fmpz_sub(res->coeffs + 0, res->coeffs + 0, d);
      if (fmpz_sgn(res->coeffs + 0) < 0)
         fmpz_add(res->coeffs + 0, res->coeffs + 0, fmpz_mod_ctx_modulus(ctx));

      _fmpz_mod_poly_normalise(res);
   }

   fmpz_clear(d);
}

void fmpz_mod_poly_si_sub(fmpz_mod_poly_t res, slong c,
                          const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
   fmpz_t d;

   fmpz_init(d);
   fmpz_set_si(d, c);

   if (fmpz_size(fmpz_mod_ctx_modulus(ctx)) > 1)
   {
      if (c < 0)
         fmpz_add(d, d, fmpz_mod_ctx_modulus(ctx));
   } else
      fmpz_mod(d, d, fmpz_mod_ctx_modulus(ctx));
      

   if (poly->length == 0)
      fmpz_mod_poly_set_fmpz(res, d, ctx);
   else
   {
      fmpz_mod_poly_neg(res, poly, ctx);

      fmpz_add(res->coeffs + 0, res->coeffs + 0, d);
      if (fmpz_cmp(res->coeffs + 0, fmpz_mod_ctx_modulus(ctx)) >= 0)
         fmpz_sub(res->coeffs + 0, res->coeffs + 0, fmpz_mod_ctx_modulus(ctx));

      _fmpz_mod_poly_normalise(res);
   }

   fmpz_clear(d);
}

void fmpz_mod_poly_add_fmpz(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                     const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
   fmpz_t d;

   fmpz_init(d);
   if (fmpz_sgn(c) < 0 || fmpz_cmp(c, fmpz_mod_ctx_modulus(ctx)) >= 0)
      fmpz_mod(d, c, fmpz_mod_ctx_modulus(ctx));
   else
      fmpz_set(d, c);

   if (poly->length == 0)
      fmpz_mod_poly_set_fmpz(res, d, ctx);
   else
   {
      fmpz_mod_poly_set(res, poly, ctx);

      fmpz_add(res->coeffs + 0, res->coeffs + 0, d);
      if (fmpz_cmp(res->coeffs + 0, fmpz_mod_ctx_modulus(ctx)) >= 0)
         fmpz_sub(res->coeffs + 0, res->coeffs + 0, fmpz_mod_ctx_modulus(ctx));

      _fmpz_mod_poly_normalise(res);
   }

   fmpz_clear(d);
}

void fmpz_mod_poly_sub_fmpz(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                      const fmpz_t c, const fmpz_mod_ctx_t ctx)
{
   fmpz_t d;

   fmpz_init(d);
   if (fmpz_sgn(c) < 0 || fmpz_cmp(c, fmpz_mod_ctx_modulus(ctx)) >= 0)
      fmpz_mod(d, c, fmpz_mod_ctx_modulus(ctx));
   else
      fmpz_set(d, c);

   if (poly->length == 0)
   {
      fmpz_sub(d, fmpz_mod_ctx_modulus(ctx), d);
      if (fmpz_cmp(d, fmpz_mod_ctx_modulus(ctx)) == 0)
         fmpz_zero(d);
      fmpz_mod_poly_set_fmpz(res, d, ctx);
   } else
   {
      fmpz_mod_poly_set(res, poly, ctx);

      fmpz_sub(res->coeffs + 0, res->coeffs + 0, d);
      if (fmpz_sgn(res->coeffs + 0) < 0)
         fmpz_add(res->coeffs + 0, res->coeffs + 0, fmpz_mod_ctx_modulus(ctx));

      _fmpz_mod_poly_normalise(res);
   }

   fmpz_clear(d);
}

void fmpz_mod_poly_fmpz_sub(fmpz_mod_poly_t res, const fmpz_t c,
                          const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
   fmpz_t d;

   fmpz_init(d);
   if (fmpz_sgn(c) < 0 || fmpz_cmp(c, fmpz_mod_ctx_modulus(ctx)) >= 0)
      fmpz_mod(d, c, fmpz_mod_ctx_modulus(ctx));
   else
      fmpz_set(d, c);


   if (poly->length == 0)
      fmpz_mod_poly_set_fmpz(res, d, ctx);
   else
   {
      fmpz_mod_poly_neg(res, poly, ctx);

      fmpz_add(res->coeffs + 0, res->coeffs + 0, d);
      if (fmpz_cmp(res->coeffs + 0, fmpz_mod_ctx_modulus(ctx)) >= 0)
         fmpz_sub(res->coeffs + 0, res->coeffs + 0, fmpz_mod_ctx_modulus(ctx));

      _fmpz_mod_poly_normalise(res);
   }

   fmpz_clear(d);
}

/* fmpz_poly_factor */
void fmpz_poly_factor_get_fmpz_poly(fmpz_poly_t z, const fmpz_poly_factor_t F, slong i)
{
    fmpz_poly_set(z, F->p + i);
}

void fmpz_poly_factor_get_fmpz(fmpz_t z, const fmpz_poly_factor_t F)
{
    fmpz_set(z, &(F->c));
}

/* fq */
void __fq_ctx_prime(fmpz_t p, fq_ctx_t ctx)
{
   fmpz_set(p, fq_ctx_prime(ctx));
}

/* fq_nmod */
void __fq_nmod_ctx_prime(fmpz_t p, fq_nmod_ctx_t ctx)
{
   fmpz_set(p, fq_nmod_ctx_prime(ctx));
}

/* fq_nmod_poly_factor */
void fq_nmod_poly_factor_get_poly(fq_nmod_poly_t z,
             const fq_nmod_poly_factor_t fac, slong i, const fq_nmod_ctx_t ctx)
{
    fq_nmod_poly_set(z, fac->poly + i, ctx);
}

/* fq_poly_factor */
void fq_poly_factor_get_poly(fq_poly_t z,
                       const fq_poly_factor_t fac, slong i, const fq_ctx_t ctx)
{
    fq_poly_set(z, fac->poly + i, ctx);
}

/* fq_zech */
void __fq_zech_ctx_prime(fmpz_t p, fq_zech_ctx_t ctx)
{
   fmpz_set(p, fq_zech_ctx_prime(ctx));
}

/* fq_zech_poly_factor */
void fq_zech_poly_factor_get_poly(fq_zech_poly_t z,
             const fq_zech_poly_factor_t fac, slong i, const fq_zech_ctx_t ctx)
{
   fq_zech_poly_set(z, fac->poly + i, ctx);
}

/* nmod_mat */
void nmod_mat_set_entry(nmod_mat_t mat, slong i, slong j, mp_limb_t x)
{
  nmod_mat_entry(mat, i, j) = x;
}

/* nmod_poly_factor */
void nmod_poly_factor_get_nmod_poly(nmod_poly_t z, nmod_poly_factor_t fac, slong i)
{
    nmod_poly_set(z, fac->p + i);
}

/* padic */
void padic_get_unit(fmpz_t f, padic_t p)
{
   fmpz_set(f, padic_unit(p));
}
