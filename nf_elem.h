/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 William Hart

******************************************************************************/

#ifndef NF_ELEM_H
#define NF_ELEM_H

#ifdef NF_ELEM_INLINES_C
#define NF_ELEM_INLINE ANTIC_DLL
#else
#define NF_ELEM_INLINE static __inline__
#endif

#include "gmp.h"
#include "flint.h"
#include "fmpq_poly.h"
#include "fmpq_mat.h"
#include "fmpz_mat.h"
#include "fmpz_mod_poly.h"
#include "nf.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct /* element of a linear number field */
{
   fmpz_t num;
   fmpz_t den;
} lnf_elem_struct;

typedef lnf_elem_struct lnf_elem_t[1];

typedef struct /* element of a quadratic number field */
{
   fmpz num[3]; /* extra coeff for delayed reduction */
   fmpz_t den;
} qnf_elem_struct;

typedef qnf_elem_struct qnf_elem_t[1];

typedef union /* element in a number field (specified by an nf_t) */
{
   fmpq_poly_t elem; /* general case */
   lnf_elem_t lelem; /* linear number field */
   qnf_elem_t qelem; /* quadratic number field */
} nf_elem_struct;

typedef nf_elem_struct nf_elem_t[1];

#define NF_ELEM_NUMREF(xxx) fmpq_poly_numref((xxx)->elem)
#define NF_ELEM_DENREF(xxx) fmpq_poly_denref((xxx)->elem)
#define LNF_ELEM_NUMREF(xxx) ((xxx)->lelem->num)
#define LNF_ELEM_DENREF(xxx) ((xxx)->lelem->den)
#define QNF_ELEM_NUMREF(xxx) ((xxx)->qelem->num)
#define QNF_ELEM_DENREF(xxx) ((xxx)->qelem->den)
#define NF_ELEM(xxx) (xxx)->elem
#define LNF_ELEM(xxx) (xxx)->lelem
#define QNF_ELEM(xxx) (xxx)->qelem

/******************************************************************************

    Initialisation

******************************************************************************/

ANTIC_DLL void nf_elem_init(nf_elem_t a, const nf_t nf);

ANTIC_DLL void nf_elem_clear(nf_elem_t a, const nf_t nf);

ANTIC_DLL void nf_elem_randtest(nf_elem_t a, flint_rand_t state,
                                              mp_bitcnt_t bits, const nf_t nf);

ANTIC_DLL void nf_elem_randtest_not_zero(nf_elem_t a, flint_rand_t state,
                                              mp_bitcnt_t bits, const nf_t nf);

NF_ELEM_INLINE
void nf_elem_canonicalise(nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR) {
      _fmpq_canonicalise(LNF_ELEM_NUMREF(a), LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
      _fmpq_poly_canonicalise(QNF_ELEM_NUMREF(a), QNF_ELEM_DENREF(a), 3);
   else
      fmpq_poly_canonicalise(NF_ELEM(a));
}

ANTIC_DLL void _nf_elem_reduce(nf_elem_t a, const nf_t nf);

ANTIC_DLL void nf_elem_reduce(nf_elem_t a, const nf_t nf);

ANTIC_DLL int _nf_elem_invertible_check(nf_elem_t a, const nf_t nf);

/******************************************************************************

    Comparison

******************************************************************************/

ANTIC_DLL int _nf_elem_equal(const nf_elem_t a, const nf_elem_t b, const nf_t nf);

ANTIC_DLL int nf_elem_equal(const nf_elem_t a, const nf_elem_t b, const nf_t nf);

NF_ELEM_INLINE
int nf_elem_is_zero(const nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      return fmpz_is_zero(LNF_ELEM_NUMREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const anum = QNF_ELEM_NUMREF(a);

      return fmpz_is_zero(anum) && fmpz_is_zero(anum + 1);
   } else
      return fmpq_poly_is_zero(a->elem);
}

NF_ELEM_INLINE
int nf_elem_is_one(const nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      return fmpz_is_one(LNF_ELEM_NUMREF(a)) && fmpz_is_one(LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      const fmpz * const anum = QNF_ELEM_NUMREF(a);

      return fmpz_is_one(anum) && fmpz_is_zero(anum + 1)
          && fmpz_is_one(QNF_ELEM_DENREF(a));
   } else
      return fmpq_poly_is_one(a->elem);
}

ANTIC_DLL int nf_elem_is_gen(const nf_elem_t a, const nf_t nf);

NF_ELEM_INLINE
int nf_elem_is_integer(const nf_elem_t a, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
        return fmpz_is_one(LNF_ELEM_DENREF(a));
    else if (nf->flag & NF_QUADRATIC)
        return fmpz_is_zero(QNF_ELEM_NUMREF(a) + 1) &&
               fmpz_is_one(QNF_ELEM_DENREF(a));
    else
        return NF_ELEM(a)->length <= 1 && fmpz_is_one(NF_ELEM_DENREF(a));
}

NF_ELEM_INLINE
int nf_elem_is_rational(const nf_elem_t a, const nf_t nf)
{
    if (nf->flag & NF_LINEAR) return 1;
    else if (nf->flag & NF_QUADRATIC)
        return fmpz_is_zero(QNF_ELEM_NUMREF(a) + 1);
    else
        return NF_ELEM(a)->length <= 1;
}

NF_ELEM_INLINE
int nf_elem_equal_si(const nf_elem_t a, const slong b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
        return fmpz_is_one(LNF_ELEM_DENREF(a)) &&
               fmpz_equal_si(LNF_ELEM_NUMREF(a), b);
    else if (nf->flag & NF_QUADRATIC)
        return fmpz_is_zero(QNF_ELEM_NUMREF(a) + 1) &&
               fmpz_is_one(QNF_ELEM_DENREF(a)) &&
               fmpz_equal_si(QNF_ELEM_NUMREF(a), b);
    else
    {
        if (b == 0) return fmpq_poly_is_zero(NF_ELEM(a));
        else return NF_ELEM(a)->length == 1 &&
                    fmpz_is_one(NF_ELEM_DENREF(a)) &&
                    fmpz_equal_si(NF_ELEM_NUMREF(a), b);
    }
}

NF_ELEM_INLINE
int nf_elem_equal_ui(const nf_elem_t a, const ulong b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
        return fmpz_is_one(LNF_ELEM_DENREF(a)) &&
               fmpz_equal_ui(LNF_ELEM_NUMREF(a), b);
    else if (nf->flag & NF_QUADRATIC)
        return fmpz_is_zero(QNF_ELEM_NUMREF(a) + 1) &&
               fmpz_is_one(QNF_ELEM_DENREF(a)) &&
               fmpz_equal_ui(QNF_ELEM_NUMREF(a), b);
    else
    {
        if (b == 0) return fmpq_poly_is_zero(NF_ELEM(a));
        else return NF_ELEM(a)->length == 1 &&
                    fmpz_is_one(NF_ELEM_DENREF(a)) &&
                    fmpz_equal_ui(NF_ELEM_NUMREF(a), b);
    }
}

NF_ELEM_INLINE
int nf_elem_equal_fmpz(const nf_elem_t a, const fmpz_t b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
        return fmpz_is_one(LNF_ELEM_DENREF(a)) && fmpz_equal(LNF_ELEM_NUMREF(a), b);
    else if (nf->flag & NF_QUADRATIC)
        return fmpz_is_zero(QNF_ELEM_NUMREF(a) + 1) &&
               fmpz_is_one(QNF_ELEM_DENREF(a)) &&
               fmpz_equal(QNF_ELEM_NUMREF(a), b);
    else
    {
        if (NF_ELEM(a)->length == 0)
            return fmpz_is_zero(b);
        else if (NF_ELEM(a)->length == 1)
            return fmpz_is_one(NF_ELEM_DENREF(a)) &&
                   fmpz_equal(NF_ELEM_NUMREF(a), b);
        else
            return 0;
    }
}

NF_ELEM_INLINE
int nf_elem_equal_fmpq(const nf_elem_t a, const fmpq_t b, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
        return fmpz_equal(LNF_ELEM_NUMREF(a), fmpq_numref(b)) &&
               fmpz_equal(LNF_ELEM_DENREF(a), fmpq_denref(b));
    else if (nf->flag & NF_QUADRATIC)
        return fmpz_is_zero(QNF_ELEM_NUMREF(a) + 1) &&
               fmpz_equal(QNF_ELEM_NUMREF(a), fmpq_numref(b)) &&
               fmpz_equal(QNF_ELEM_DENREF(a), fmpq_denref(b)) ;
    else
    {
        if (NF_ELEM(a)->length == 0)
            return fmpq_is_zero(b);
        else if (NF_ELEM(a)->length == 1)
            return fmpz_equal(NF_ELEM_NUMREF(a), fmpq_numref(b)) &&
                   fmpz_equal(NF_ELEM_DENREF(a), fmpq_denref(b));
        else
            return 0;
    }
}

/******************************************************************************

    I/O

******************************************************************************/

ANTIC_DLL void nf_elem_print_pretty(const nf_elem_t a,
                             const nf_t nf, const char * var);

ANTIC_DLL char * nf_elem_get_str_pretty(const nf_elem_t a,
                              const char * var, const nf_t nf);

/******************************************************************************

    Element creation

******************************************************************************/

NF_ELEM_INLINE
void nf_elem_zero(nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_zero(LNF_ELEM_NUMREF(a));
      fmpz_one(LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);

      fmpz_zero(anum);
      fmpz_zero(anum + 1);
      fmpz_one(QNF_ELEM_DENREF(a));
   } else
      fmpq_poly_zero(NF_ELEM(a));
}

NF_ELEM_INLINE
void nf_elem_one(nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_one(LNF_ELEM_NUMREF(a));
      fmpz_one(LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);

      fmpz_one(anum);
      fmpz_zero(anum + 1);
      fmpz_one(QNF_ELEM_DENREF(a));
   } else
      fmpq_poly_one(NF_ELEM(a));
}

NF_ELEM_INLINE
void nf_elem_gen(nf_elem_t a, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_neg(LNF_ELEM_NUMREF(a), nf->pol->coeffs);
      fmpz_set(LNF_ELEM_DENREF(a), nf->pol->coeffs + 1);
      _fmpq_canonicalise(LNF_ELEM_NUMREF(a), LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);

      fmpz_one(anum + 1);
      fmpz_zero(anum);
      fmpz_one(QNF_ELEM_DENREF(a));
   } else
   {
      fmpq_poly_zero(NF_ELEM(a));
      fmpq_poly_set_coeff_ui(NF_ELEM(a), 1, 1);
   }
}

NF_ELEM_INLINE
void nf_elem_set_si(nf_elem_t a, slong c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_set_si(LNF_ELEM_NUMREF(a), c);
      fmpz_one(LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);

      fmpz_set_si(anum, c);
      fmpz_zero(anum + 1);
      fmpz_one(QNF_ELEM_DENREF(a));
   } else
      fmpq_poly_set_si(NF_ELEM(a), c);
}

NF_ELEM_INLINE
void nf_elem_set_ui(nf_elem_t a, ulong c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_set_ui(LNF_ELEM_NUMREF(a), c);
      fmpz_one(LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);

      fmpz_set_ui(anum, c);
      fmpz_zero(anum + 1);
      fmpz_one(QNF_ELEM_DENREF(a));
   } else
      fmpq_poly_set_ui(NF_ELEM(a), c);
}

NF_ELEM_INLINE
void nf_elem_set_fmpz(nf_elem_t a, const fmpz_t c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_set(LNF_ELEM_NUMREF(a), c);
      fmpz_one(LNF_ELEM_DENREF(a));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);

      fmpz_set(anum, c);
      fmpz_zero(anum + 1);
      fmpz_one(QNF_ELEM_DENREF(a));
   } else
      fmpq_poly_set_fmpz(NF_ELEM(a), c);
}

NF_ELEM_INLINE
void nf_elem_set_fmpq(nf_elem_t a, const fmpq_t c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_set(LNF_ELEM_NUMREF(a), fmpq_numref(c));
      fmpz_set(LNF_ELEM_DENREF(a), fmpq_denref(c));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);

      fmpz_set(anum, fmpq_numref(c));
      fmpz_zero(anum + 1);
      fmpz_set(QNF_ELEM_DENREF(a), fmpq_denref(c));
   } else
      fmpq_poly_set_fmpq(NF_ELEM(a), c);
}

ANTIC_DLL void nf_elem_set_fmpq_poly(nf_elem_t a, const fmpq_poly_t pol, const nf_t nf);

/******************************************************************************

    Conversion

******************************************************************************/

ANTIC_DLL
void nf_elem_set_fmpz_mat_row(nf_elem_t b, const fmpz_mat_t M,
                                     const slong i, fmpz_t den, const nf_t nf);

ANTIC_DLL
void nf_elem_get_fmpz_mat_row(fmpz_mat_t M, const slong i, fmpz_t den,
                                             const nf_elem_t b, const nf_t nf);

FLINT_DLL
void nf_elem_get_fmpq_poly(fmpq_poly_t pol, const nf_elem_t a, const nf_t nf);

ANTIC_DLL
void _nf_elem_get_nmod_poly(nmod_poly_t pol, const nf_elem_t a, const nf_t nf);

ANTIC_DLL
void nf_elem_get_nmod_poly_den(nmod_poly_t pol,
                                    const nf_elem_t a, const nf_t nf, int den);

ANTIC_DLL
void nf_elem_get_nmod_poly(nmod_poly_t pol, const nf_elem_t a, const nf_t nf);

#if __FLINT_RELEASE >= 20700
ANTIC_DLL
void _nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol,
                   const nf_elem_t a, const nf_t nf, const fmpz_mod_ctx_t ctx);

ANTIC_DLL
void nf_elem_get_fmpz_mod_poly_den(fmpz_mod_poly_t pol,
          const nf_elem_t a, const nf_t nf, int den, const fmpz_mod_ctx_t ctx);

ANTIC_DLL
void nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol,
                   const nf_elem_t a, const nf_t nf, const fmpz_mod_ctx_t ctx);
#else
ANTIC_DLL
void _nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol, const nf_elem_t a, const nf_t nf);

ANTIC_DLL
void nf_elem_get_fmpz_mod_poly_den(fmpz_mod_poly_t pol, const nf_elem_t a, const nf_t nf, int den);

ANTIC_DLL
void nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol, const nf_elem_t a, const nf_t nf);
#endif

/******************************************************************************

    Basic manipulation

******************************************************************************/

NF_ELEM_INLINE
void nf_elem_get_den(fmpz_t d, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
     fmpz_set(d, LNF_ELEM_DENREF(b));
   } else if (nf->flag & NF_QUADRATIC)
   {
     fmpz_set(d, QNF_ELEM_DENREF(b));
   } else
   {
     fmpz_set(d, NF_ELEM_DENREF(b));
   }
}

NF_ELEM_INLINE
void nf_elem_set_den(nf_elem_t b, fmpz_t d, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
     fmpz_set(LNF_ELEM_DENREF(b), d);
   } else if (nf->flag & NF_QUADRATIC)
   {
     fmpz_set(QNF_ELEM_DENREF(b), d);
   } else
   {
     fmpz_set(NF_ELEM_DENREF(b), d);
   }
}

ANTIC_DLL
void nf_elem_get_coeff_fmpq(fmpq_t a, const nf_elem_t b, slong i, const nf_t nf);

ANTIC_DLL
void nf_elem_get_coeff_fmpz(fmpz_t a, const nf_elem_t b, slong i, const nf_t nf);

NF_ELEM_INLINE
int nf_elem_den_is_one(const nf_elem_t a, const nf_t nf)
{
    if (nf->flag & NF_LINEAR)
    {
        return fmpz_is_one(LNF_ELEM_DENREF(a));
    } else if (nf->flag & NF_QUADRATIC)
    {
        return fmpz_is_one(QNF_ELEM_DENREF(a));
    } else
    {
        return fmpz_is_one(NF_ELEM_DENREF(a));
    }
}

ANTIC_DLL
void _nf_elem_set_coeff_num_fmpz(nf_elem_t a, slong i, const fmpz_t b, const nf_t nf);

/******************************************************************************

    Arithmetic

******************************************************************************/

NF_ELEM_INLINE
void nf_elem_set(nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_set(LNF_ELEM_NUMREF(a), LNF_ELEM_NUMREF(b));
      fmpz_set(LNF_ELEM_DENREF(a), LNF_ELEM_DENREF(b));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      const fmpz * const bnum = QNF_ELEM_NUMREF(b);

      fmpz_set(anum, bnum);
      fmpz_set(anum + 1, bnum + 1);
      fmpz_set(QNF_ELEM_DENREF(a), QNF_ELEM_DENREF(b));
   } else
      fmpq_poly_set(NF_ELEM(a), NF_ELEM(b));
}

NF_ELEM_INLINE
void nf_elem_neg(nf_elem_t a, const nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_neg(LNF_ELEM_NUMREF(a), LNF_ELEM_NUMREF(b));
      fmpz_set(LNF_ELEM_DENREF(a), LNF_ELEM_DENREF(b));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      const fmpz * const bnum = QNF_ELEM_NUMREF(b);

      fmpz_neg(anum, bnum);
      fmpz_neg(anum + 1, bnum + 1);
      fmpz_set(QNF_ELEM_DENREF(a), QNF_ELEM_DENREF(b));
   } else
      fmpq_poly_neg(NF_ELEM(a), NF_ELEM(b));
}

NF_ELEM_INLINE
void nf_elem_swap(nf_elem_t a, nf_elem_t b, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
   {
      fmpz_swap(LNF_ELEM_NUMREF(a), LNF_ELEM_NUMREF(b));
      fmpz_swap(LNF_ELEM_DENREF(a), LNF_ELEM_DENREF(b));
   } else if (nf->flag & NF_QUADRATIC)
   {
      fmpz * const anum = QNF_ELEM_NUMREF(a);
      fmpz * const bnum = QNF_ELEM_NUMREF(b);

      fmpz_swap(anum, bnum);
      fmpz_swap(anum + 1, bnum + 1);
      fmpz_swap(anum + 2, bnum + 2);
      fmpz_swap(QNF_ELEM_DENREF(a), QNF_ELEM_DENREF(b));
   } else
      fmpq_poly_swap(NF_ELEM(a), NF_ELEM(b));
}

ANTIC_DLL void nf_elem_add_si(nf_elem_t a,
                                   const nf_elem_t b, slong c, const nf_t nf);

ANTIC_DLL void nf_elem_add_fmpz(nf_elem_t a,
                                  const nf_elem_t b, const fmpz_t c, const nf_t nf);

ANTIC_DLL void nf_elem_add_fmpq(nf_elem_t a,
                                  const nf_elem_t b, const fmpq_t c, const nf_t nf);

ANTIC_DLL void nf_elem_sub_si(nf_elem_t a,
                                   const nf_elem_t b, slong c, const nf_t nf);

ANTIC_DLL void nf_elem_sub_fmpz(nf_elem_t a,
                                  const nf_elem_t b, const fmpz_t c, const nf_t nf);

ANTIC_DLL void nf_elem_sub_fmpq(nf_elem_t a,
                                  const nf_elem_t b, const fmpq_t c, const nf_t nf);

ANTIC_DLL void nf_elem_si_sub(nf_elem_t a,
                                   slong c, const nf_elem_t b, const nf_t nf);

ANTIC_DLL void nf_elem_fmpz_sub(nf_elem_t a,
                                  const fmpz_t c, const nf_elem_t b, const nf_t nf);

ANTIC_DLL void nf_elem_fmpq_sub(nf_elem_t a,
                                  const fmpq_t c, const nf_elem_t b, const nf_t nf);

ANTIC_DLL void nf_elem_scalar_mul_si(nf_elem_t a, const nf_elem_t b,
                                                      slong c, const nf_t nf);

ANTIC_DLL void nf_elem_scalar_mul_fmpz(nf_elem_t a, const nf_elem_t b,
                                                     const fmpz_t c, const nf_t nf);

ANTIC_DLL void nf_elem_scalar_mul_fmpq(nf_elem_t a, const nf_elem_t b,
                                                     const fmpq_t c, const nf_t nf);
									
ANTIC_DLL void nf_elem_scalar_div_si(nf_elem_t a, const nf_elem_t b,
                                                      slong c, const nf_t nf);

ANTIC_DLL void nf_elem_scalar_div_fmpz(nf_elem_t a, const nf_elem_t b,
                                                     const fmpz_t c, const nf_t nf);

ANTIC_DLL void nf_elem_scalar_div_fmpq(nf_elem_t a, const nf_elem_t b,
                                                     const fmpq_t c, const nf_t nf);
									
ANTIC_DLL void _nf_elem_add_lf(nf_elem_t a, const nf_elem_t b,
                                   const nf_elem_t c, const nf_t nf, int can);

ANTIC_DLL void _nf_elem_sub_lf(nf_elem_t a, const nf_elem_t b,
                                   const nf_elem_t c, const nf_t nf, int can);

ANTIC_DLL void _nf_elem_add_qf(nf_elem_t a, const nf_elem_t b,
                                   const nf_elem_t c, const nf_t nf, int can);

ANTIC_DLL void _nf_elem_sub_qf(nf_elem_t a, const nf_elem_t b,
                                   const nf_elem_t c, const nf_t nf, int can);

ANTIC_DLL void nf_elem_add_qf(nf_elem_t a, const nf_elem_t b,
                                            const nf_elem_t c, const nf_t nf);

ANTIC_DLL void nf_elem_sub_qf(nf_elem_t a, const nf_elem_t b,
                                            const nf_elem_t c, const nf_t nf);

NF_ELEM_INLINE
void _nf_elem_add(nf_elem_t a, const nf_elem_t b,
                                              const nf_elem_t c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
      _nf_elem_add_lf(a, b, c, nf, 0);
   else if (nf->flag & NF_QUADRATIC)
      _nf_elem_add_qf(a, b, c, nf, 0);
   else
      fmpq_poly_add_can(NF_ELEM(a), NF_ELEM(b), NF_ELEM(c), 0);
}

NF_ELEM_INLINE
void _nf_elem_sub(nf_elem_t a, const nf_elem_t b,
                                              const nf_elem_t c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
      _nf_elem_sub_lf(a, b, c, nf, 0);
   else if (nf->flag & NF_QUADRATIC)
      _nf_elem_sub_qf(a, b, c, nf, 0);
   else
      fmpq_poly_sub_can(NF_ELEM(a), NF_ELEM(b), NF_ELEM(c), 0);
}

NF_ELEM_INLINE
void nf_elem_add(nf_elem_t a, const nf_elem_t b,
                                              const nf_elem_t c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
      _nf_elem_add_lf(a, b, c, nf, 1);
   else if (nf->flag & NF_QUADRATIC)
      nf_elem_add_qf(a, b, c, nf);
   else
      fmpq_poly_add_can(NF_ELEM(a), NF_ELEM(b), NF_ELEM(c), 1);
}

NF_ELEM_INLINE
void nf_elem_sub(nf_elem_t a, const nf_elem_t b,
                                              const nf_elem_t c, const nf_t nf)
{
   if (nf->flag & NF_LINEAR)
      _nf_elem_sub_lf(a, b, c, nf, 1);
   else if (nf->flag & NF_QUADRATIC)
      nf_elem_sub_qf(a, b, c, nf);
   else
      fmpq_poly_sub_can(NF_ELEM(a), NF_ELEM(b), NF_ELEM(c), 1);
}

FLINT_DLL void nf_elem_mul_gen(nf_elem_t a, const nf_elem_t b, const nf_t nf);

ANTIC_DLL void _nf_elem_mul(nf_elem_t a, const nf_elem_t b,
                                             const nf_elem_t c, const nf_t nf);

ANTIC_DLL void nf_elem_mul(nf_elem_t a, const nf_elem_t b,
                                             const nf_elem_t c, const nf_t nf);

ANTIC_DLL void _nf_elem_mul_red(nf_elem_t a, const nf_elem_t b,
                                    const nf_elem_t c, const nf_t nf, int red);

ANTIC_DLL void nf_elem_mul_red(nf_elem_t a, const nf_elem_t b,
                                    const nf_elem_t c, const nf_t nf, int red);

ANTIC_DLL void _nf_elem_inv(nf_elem_t a, const nf_elem_t b, const nf_t nf);

ANTIC_DLL void nf_elem_inv(nf_elem_t a, const nf_elem_t b, const nf_t nf);

ANTIC_DLL void _nf_elem_div(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf);

ANTIC_DLL void nf_elem_div(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf);

ANTIC_DLL void _nf_elem_pow(nf_elem_t res, const nf_elem_t b, ulong e, const nf_t nf);

ANTIC_DLL void nf_elem_pow(nf_elem_t res, const nf_elem_t a, ulong e, const nf_t nf);

ANTIC_DLL void _nf_elem_norm(fmpz_t rnum, fmpz_t rden, const nf_elem_t a, const nf_t nf);

ANTIC_DLL void nf_elem_norm(fmpq_t res, const nf_elem_t a, const nf_t nf);

ANTIC_DLL void _nf_elem_norm_div(fmpz_t rnum, fmpz_t rden, const nf_elem_t a,
                             const nf_t nf, const fmpz_t divisor, slong nbits);

ANTIC_DLL void nf_elem_norm_div(fmpq_t res, const nf_elem_t a, const nf_t nf,
                                            const fmpz_t divisor, slong nbits);

ANTIC_DLL void _nf_elem_trace(fmpz_t rnum, fmpz_t rden, const nf_elem_t a,
                                                                const nf_t nf);

ANTIC_DLL void nf_elem_trace(fmpq_t res, const nf_elem_t a, const nf_t nf);

ANTIC_DLL void nf_elem_rep_mat(fmpq_mat_t res, const nf_elem_t a, const nf_t nf);

ANTIC_DLL void nf_elem_rep_mat_fmpz_mat_den(fmpz_mat_t res, fmpz_t den, const nf_elem_t a, const nf_t nf);

/******************************************************************************

    Modular reduction

******************************************************************************/

ANTIC_DLL
void _nf_elem_mod_fmpz(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf, int sign);

ANTIC_DLL
void nf_elem_mod_fmpz_den(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf, int den);

ANTIC_DLL
void nf_elem_smod_fmpz_den(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf, int den);

ANTIC_DLL
void nf_elem_mod_fmpz(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf);

ANTIC_DLL
void nf_elem_smod_fmpz(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf);

ANTIC_DLL
void nf_elem_coprime_den(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf);

ANTIC_DLL
void nf_elem_coprime_den_signed(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf);

#ifdef __cplusplus
}
#endif

/******************************************************************************

    Helpers for compatibility with FLINT 2.6

******************************************************************************/
#if __FLINT_RELEASE >= 20700
#define FMPZ_MOD_POLY_FIT_LENGTH(POL, N, CTX) fmpz_mod_poly_fit_length(POL, N, CTX)
#define FMPZ_MOD(F, G, CTX, P) fmpz_mod(F, G, (CTX)->n)
#define FMPZ_MOD_POLY_ZERO(POL, CTX) fmpz_mod_poly_zero(POL, CTX)
#define FMPZ_MOD_POLY_SCALAR_DIV_FMPZ(RES, POL, X, CTX) fmpz_mod_poly_scalar_div_fmpz(RES, POL, X, CTX)
#else
#define FMPZ_MOD_POLY_FIT_LENGTH(POL, N, CTX) fmpz_mod_poly_fit_length(POL, N)
#define FMPZ_MOD(F, G, CTX, P) fmpz_mod(F, G, P)
#define FMPZ_MOD_POLY_ZERO(POL, CTX) fmpz_mod_poly_zero(POL)
#define FMPZ_MOD_POLY_SCALAR_DIV_FMPZ(RES, POL, X, CTX) fmpz_mod_poly_scalar_div_fmpz(RES, POL, X)
#endif

#endif
