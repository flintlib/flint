/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_DEFAULT_POLY_H
#define FQ_DEFAULT_POLY_H

#ifdef FQ_DEFAULT_POLY_INLINES_C
#define FQ_DEFAULT_POLY_INLINE FLINT_DLL
#else
#define FQ_DEFAULT_POLY_INLINE static __inline__
#endif

#include "ulong_extras.h"
#include "fmpz.h"
#include "fq.h"
#include "fq_nmod.h"
#include "fq_zech.h"
#include "fq_default.h"
#include "fq_poly.h"
#include "fq_nmod_poly.h"
#include "fq_zech_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/*  Type definitions *********************************************************/

typedef union fq_default_poly_struct
{
    fq_poly_t fq;
    fq_nmod_poly_t fq_nmod;
    fq_zech_poly_t fq_zech;
}
fq_default_poly_struct;

typedef fq_default_poly_struct fq_default_poly_t[1];

/*  Memory management ********************************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_init(fq_default_poly_t poly,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_init(poly->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_init(poly->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_init(poly->fq, ctx->ctx.fq);
   } 
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_init2(fq_default_poly_t poly,
                                       slong alloc, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_init2(poly->fq_zech, alloc, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_init2(poly->fq_nmod, alloc, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_init2(poly->fq, alloc, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_realloc(fq_default_poly_t poly,
                                       slong alloc, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_realloc(poly->fq_zech, alloc, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_realloc(poly->fq_nmod, alloc, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_realloc(poly->fq, alloc, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_truncate(fq_default_poly_t poly,
                                         slong len, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_truncate(poly->fq_zech, len, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_truncate(poly->fq_nmod, len, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_truncate(poly->fq, len, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_set_trunc(fq_default_poly_t poly1,
                fq_default_poly_t poly2, slong len, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_set_trunc(poly1->fq_zech,
		                        poly2->fq_zech, len, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_set_trunc(poly1->fq_nmod,
		                        poly2->fq_nmod, len, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_set_trunc(poly1->fq, poly2->fq, len, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_fit_length(fq_default_poly_t poly,
                                         slong len, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_fit_length(poly->fq_zech, len, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_fit_length(poly->fq_nmod, len, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_fit_length(poly->fq, len, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void _fq_default_poly_set_length(fq_default_poly_t poly,
                                         slong len, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      _fq_zech_poly_set_length(poly->fq_zech, len, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      _fq_nmod_poly_set_length(poly->fq_nmod, len, ctx->ctx.fq_nmod);
   } else
   {
      _fq_poly_set_length(poly->fq, len, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_clear(fq_default_poly_t poly,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_clear(poly->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_clear(poly->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_clear(poly->fq, ctx->ctx.fq);
   }
}

/*  Polynomial parameters  ***************************************************/

FQ_DEFAULT_POLY_INLINE slong
fq_default_poly_length(const fq_default_poly_t poly,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_length(poly->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_length(poly->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_length(poly->fq, ctx->ctx.fq);
}

FQ_DEFAULT_POLY_INLINE slong
fq_default_poly_degree(const fq_default_poly_t poly,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_degree(poly->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_degree(poly->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_degree(poly->fq, ctx->ctx.fq);
}

/*  Randomisation  ***********************************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_randtest(fq_default_poly_t f,
                     flint_rand_t state, slong len, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_randtest(f->fq_zech, state, len, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_randtest(f->fq_nmod, state, len, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_randtest(f->fq, state, len, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_randtest_not_zero(fq_default_poly_t f,
                     flint_rand_t state, slong len, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_randtest_not_zero(f->fq_zech, state, len, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_randtest_not_zero(f->fq_nmod, state, len, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_randtest_not_zero(f->fq, state, len, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_randtest_monic(fq_default_poly_t f,
                     flint_rand_t state, slong len, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_randtest_monic(f->fq_zech, state, len, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_randtest_monic(f->fq_nmod, state, len, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_randtest_monic(f->fq, state, len, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_randtest_irreducible(fq_default_poly_t f,
                     flint_rand_t state, slong len, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_randtest_irreducible(f->fq_zech,
                                                 state, len, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_randtest_irreducible(f->fq_nmod,
                                                 state, len, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_randtest_irreducible(f->fq, state, len, ctx->ctx.fq);
   }
}

/*  Assignment and basic manipulation  ***************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_set(fq_default_poly_t rop,
                        const fq_default_poly_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_set(rop->fq_zech, op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_set(rop->fq_nmod, op->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_set(rop->fq, op->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_set_fq_default(fq_default_poly_t poly,
                              const fq_default_t c, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_set_fq_zech(poly->fq_zech, c->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_set_fq_nmod(poly->fq_nmod, c->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_set_fq(poly->fq, c->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_swap(fq_default_poly_t op1,
                             fq_default_poly_t op2, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_swap(op1->fq_zech, op2->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_swap(op1->fq_nmod, op2->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_swap(op1->fq, op2->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void
fq_default_poly_zero(fq_default_poly_t poly, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_zero(poly->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_zero(poly->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_zero(poly->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_one(fq_default_poly_t poly,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_one(poly->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_one(poly->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_one(poly->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_gen(fq_default_poly_t f,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_gen(f->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_gen(f->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_gen(f->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_make_monic(fq_default_poly_t rop,
                        const fq_default_poly_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_make_monic(rop->fq_zech, op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_make_monic(rop->fq_nmod, op->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_make_monic(rop->fq, op->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_reverse(fq_default_poly_t res,
             const fq_default_poly_t poly, slong n, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_reverse(res->fq_zech, poly->fq_zech, n, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_reverse(res->fq_nmod, poly->fq_nmod, n, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_reverse(res->fq, poly->fq, n, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
ulong fq_default_poly_deflation(const fq_default_poly_t input,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_deflation(input->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_deflation(input->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_deflation(input->fq, ctx->ctx.fq);
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_deflate(fq_default_poly_t result,
    const fq_default_poly_t input, ulong deflation, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_deflate(result->fq_zech,
                                  input->fq_zech, deflation, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_deflate(result->fq_nmod,
                                  input->fq_nmod, deflation, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_deflate(result->fq, input->fq, deflation, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_inflate(fq_default_poly_t result,
    const fq_default_poly_t input, ulong inflation, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_inflate(result->fq_zech,
                                  input->fq_zech, inflation, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_inflate(result->fq_nmod,
                                  input->fq_nmod, inflation, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_inflate(result->fq, input->fq, inflation, ctx->ctx.fq);
   }
}

/*  Getting and setting coefficients  ****************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_get_coeff(fq_default_t x,
             const fq_default_poly_t poly, slong n, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_get_coeff(x->fq_zech, poly->fq_zech, n, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_get_coeff(x->fq_nmod, poly->fq_nmod, n, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_get_coeff(x->fq, poly->fq, n, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_set_coeff(fq_default_poly_t poly,
                     slong n, const fq_default_t x, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_set_coeff(poly->fq_zech, n, x->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_set_coeff(poly->fq_nmod, n, x->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_set_coeff(poly->fq, n, x->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void 
fq_default_poly_set_coeff_fmpz(fq_default_poly_t poly,
                           slong n, const fmpz_t x, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_set_coeff_fmpz(poly->fq_zech, n, x, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_set_coeff_fmpz(poly->fq_nmod, n, x, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_set_coeff_fmpz(poly->fq, n, x, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_set_nmod_poly(fq_default_poly_t rop,
                          const     nmod_poly_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_set_nmod_poly(rop->fq_zech, op, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_set_nmod_poly(rop->fq_nmod, op, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_set_nmod_poly(rop->fq, op, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_set_fmpz_mod_poly(fq_default_poly_t rop,
                          const fmpz_mod_poly_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_set_fmpz_mod_poly(rop->fq_zech, op, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_set_fmpz_mod_poly(rop->fq_nmod, op, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_set_fmpz_mod_poly(rop->fq, op, ctx->ctx.fq);
   }
}

FLINT_DLL
void fq_default_poly_set_fmpz_poly(fq_default_poly_t rop,
                             const fmpz_poly_t op, const fq_default_ctx_t ctx);

/*  Comparison  **************************************************************/

FQ_DEFAULT_POLY_INLINE
int fq_default_poly_equal(const fq_default_poly_t poly1,
                     const fq_default_poly_t poly2, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_equal(poly1->fq_zech, poly2->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_equal(poly1->fq_nmod, poly2->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_equal(poly1->fq, poly2->fq, ctx->ctx.fq);
}

FQ_DEFAULT_POLY_INLINE
int fq_default_poly_equal_trunc(const fq_default_poly_t poly1,
            const fq_default_poly_t poly2, slong n, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_equal_trunc(poly1->fq_zech,
                                          poly2->fq_zech, n, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_equal_trunc(poly1->fq_nmod,
                                          poly2->fq_nmod, n, ctx->ctx.fq_nmod);
   }
   return fq_poly_equal_trunc(poly1->fq, poly2->fq, n, ctx->ctx.fq);
}

FQ_DEFAULT_POLY_INLINE int
fq_default_poly_is_zero(const fq_default_poly_t poly,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_is_zero(poly->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_is_zero(poly->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_is_zero(poly->fq, ctx->ctx.fq);
}

FQ_DEFAULT_POLY_INLINE int
fq_default_poly_is_one(const fq_default_poly_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_is_one(op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_is_one(op->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_is_one(op->fq, ctx->ctx.fq);
}

FQ_DEFAULT_POLY_INLINE int
fq_default_poly_is_unit(const fq_default_poly_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_is_unit(op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_is_unit(op->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_is_unit(op->fq, ctx->ctx.fq);
}

FQ_DEFAULT_POLY_INLINE int
fq_default_poly_is_gen(const fq_default_poly_t poly,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_is_gen(poly->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_is_gen(poly->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_is_gen(poly->fq, ctx->ctx.fq);
}

FQ_DEFAULT_POLY_INLINE int
fq_default_poly_equal_fq_default(const fq_default_poly_t poly,
                              const fq_default_t c, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_equal_fq_zech(poly->fq_zech,
		                                 c->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_equal_fq_nmod(poly->fq_nmod,
		                                 c->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_equal_fq(poly->fq, c->fq, ctx->ctx.fq);
}

/*  Addition and subtraction  ************************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_add(fq_default_poly_t rop, const fq_default_poly_t op1,
                       const fq_default_poly_t op2, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_add(rop->fq_zech,
                                 op1->fq_zech, op2->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_add(rop->fq_nmod,
                                 op1->fq_nmod, op2->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_add(rop->fq, op1->fq, op2->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_add_si(fq_default_poly_t rop,
              const fq_default_poly_t op1, slong c, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_add_si(rop->fq_zech, op1->fq_zech, c, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_add_si(rop->fq_nmod, op1->fq_nmod, c, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_add_si(rop->fq, op1->fq, c, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_add_series(fq_default_poly_t rop,
		const fq_default_poly_t op1, const fq_default_poly_t op2,
                                           slong n, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_add_series(rop->fq_zech,
                              op1->fq_zech, op2->fq_zech, n, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_add_series(rop->fq_nmod,
                              op1->fq_nmod, op2->fq_nmod, n, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_add_series(rop->fq, op1->fq, op2->fq, n, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_sub(fq_default_poly_t rop,
                const fq_default_poly_t op1, const fq_default_poly_t op2,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_sub(rop->fq_zech,
                                 op1->fq_zech, op2->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_sub(rop->fq_nmod,
                                 op1->fq_nmod, op2->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_sub(rop->fq, op1->fq, op2->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_sub_series(fq_default_poly_t rop,
		const fq_default_poly_t op1, const fq_default_poly_t op2,
                                           slong n, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_sub_series(rop->fq_zech,
                              op1->fq_zech, op2->fq_zech, n, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_sub_series(rop->fq_nmod,
                              op1->fq_nmod, op2->fq_nmod, n, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_sub_series(rop->fq, op1->fq, op2->fq, n, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_neg(fq_default_poly_t rop,
                        const fq_default_poly_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_neg(rop->fq_zech, op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_neg(rop->fq_nmod, op->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_neg(rop->fq, op->fq, ctx->ctx.fq);
   }
}

/*  Scalar multiplication and division  **************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_scalar_mul_fq_default(fq_default_poly_t rop,
                 const fq_default_poly_t op, const fq_default_t x,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_scalar_mul_fq_zech(rop->fq_zech,
                                    op->fq_zech, x->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_scalar_mul_fq_nmod(rop->fq_nmod,
                                    op->fq_nmod, x->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_scalar_mul_fq(rop->fq, op->fq, x->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_scalar_div_fq_default(fq_default_poly_t rop,
                 const fq_default_poly_t op, const fq_default_t x,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_scalar_div_fq_zech(rop->fq_zech,
                                    op->fq_zech, x->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_scalar_div_fq_nmod(rop->fq_nmod,
                                    op->fq_nmod, x->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_scalar_div_fq(rop->fq, op->fq, x->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_scalar_addmul_fq_default(fq_default_poly_t rop,
                 const fq_default_poly_t op, const fq_default_t x,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_scalar_addmul_fq_zech(rop->fq_zech,
                                    op->fq_zech, x->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_scalar_addmul_fq_nmod(rop->fq_nmod,
                                    op->fq_nmod, x->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_scalar_addmul_fq(rop->fq, op->fq, x->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_scalar_submul_fq_default(fq_default_poly_t rop,
                 const fq_default_poly_t op, const fq_default_t x,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_scalar_submul_fq_zech(rop->fq_zech,
                                    op->fq_zech, x->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_scalar_submul_fq_nmod(rop->fq_nmod,
                                    op->fq_nmod, x->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_scalar_submul_fq(rop->fq, op->fq, x->fq, ctx->ctx.fq);
   }
}

/*  Multiplication  **********************************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_mul(fq_default_poly_t rop,
                 const fq_default_poly_t op1, const fq_default_poly_t op2,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_mul(rop->fq_zech,
                                 op1->fq_zech, op2->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_mul(rop->fq_nmod,
                                 op1->fq_nmod, op2->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_mul(rop->fq, op1->fq, op2->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_mullow(fq_default_poly_t rop,
                 const fq_default_poly_t op1, const fq_default_poly_t op2,
                                           slong n, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_mullow(rop->fq_zech,
                              op1->fq_zech, op2->fq_zech, n, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_mullow(rop->fq_nmod,
                              op1->fq_nmod, op2->fq_nmod, n, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_mullow(rop->fq, op1->fq, op2->fq, n, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_mulhigh(fq_default_poly_t rop,
                 const fq_default_poly_t op1, const fq_default_poly_t op2,
                                       slong start, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_mulhigh(rop->fq_zech,
                          op1->fq_zech, op2->fq_zech, start, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_mulhigh(rop->fq_nmod,
                          op1->fq_nmod, op2->fq_nmod, start, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_mulhigh(rop->fq, op1->fq, op2->fq, start, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_mulmod(fq_default_poly_t res,
                 const fq_default_poly_t poly1, const fq_default_poly_t poly2,
                         const fq_default_poly_t f, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_mulmod(res->fq_zech, poly1->fq_zech, poly2->fq_zech, f->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_mulmod(res->fq_nmod, poly1->fq_nmod, poly2->fq_nmod, f->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_mulmod(res->fq, poly1->fq, poly2->fq, f->fq, ctx->ctx.fq);
   }
}

/* Squaring ******************************************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_sqr(fq_default_poly_t rop,
                        const fq_default_poly_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_sqr(rop->fq_zech, op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_sqr(rop->fq_nmod, op->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_sqr(rop->fq, op->fq, ctx->ctx.fq);
   }
}

/*  Powering  ****************************************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_pow(fq_default_poly_t rop,
                                      const fq_default_poly_t op, ulong e,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_pow(rop->fq_zech, op->fq_zech, e, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_pow(rop->fq_nmod, op->fq_nmod, e, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_pow(rop->fq, op->fq, e, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_pow_trunc(fq_default_poly_t res,
                                    const fq_default_poly_t poly, ulong e,
				       slong trunc, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_pow_trunc(res->fq_zech,
                                    poly->fq_zech, e, trunc, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_pow_trunc(res->fq_nmod,
                                    poly->fq_nmod, e, trunc, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_pow_trunc(res->fq, poly->fq, e, trunc, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_powmod_fmpz_binexp(fq_default_poly_t res,
                             const fq_default_poly_t poly, const fmpz_t e,
                         const fq_default_poly_t f, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_powmod_fmpz_binexp(res->fq_zech,
                               poly->fq_zech, e, f->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_powmod_fmpz_binexp(res->fq_nmod,
                               poly->fq_nmod, e, f->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_powmod_fmpz_binexp(res->fq, poly->fq, e, f->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_powmod_ui_binexp(fq_default_poly_t res,
                         const fq_default_poly_t poly, ulong e,
                         const fq_default_poly_t f, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_powmod_ui_binexp(res->fq_zech,
                               poly->fq_zech, e, f->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_powmod_ui_binexp(res->fq_nmod,
                               poly->fq_nmod, e, f->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_powmod_ui_binexp(res->fq, poly->fq, e, f->fq, ctx->ctx.fq);
   }
}

/*  Shifting  ****************************************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_shift_left(fq_default_poly_t rop,
               const fq_default_poly_t op, slong n, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_shift_left(rop->fq_zech, op->fq_zech, n, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_shift_left(rop->fq_nmod, op->fq_nmod, n, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_shift_left(rop->fq, op->fq, n, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_shift_right(fq_default_poly_t rop,
               const fq_default_poly_t op, slong n, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_shift_right(rop->fq_zech, op->fq_zech, n, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_shift_right(rop->fq_nmod, op->fq_nmod, n, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_shift_right(rop->fq, op->fq, n, ctx->ctx.fq);
   }
}

/*  Norms  *******************************************************************/

FQ_DEFAULT_POLY_INLINE
slong fq_default_poly_hamming_weight(const fq_default_poly_t op,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_hamming_weight(op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_hamming_weight(op->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_hamming_weight(op->fq, ctx->ctx.fq);
}

/*  Greatest common divisor  *************************************************/

FQ_DEFAULT_POLY_INLINE void fq_default_poly_gcd(fq_default_poly_t rop,
                  const fq_default_poly_t op1, const fq_default_poly_t op2,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_gcd(rop->fq_zech,
		                 op1->fq_zech, op2->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_gcd(rop->fq_nmod,
		                 op1->fq_nmod, op2->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_gcd(rop->fq, op1->fq, op2->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void
fq_default_poly_xgcd(fq_default_poly_t G,
                       fq_default_poly_t S, fq_default_poly_t T,
                  const fq_default_poly_t A, const fq_default_poly_t B,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_xgcd(G->fq_zech,
             S->fq_zech, T->fq_zech, A->fq_zech, B->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_xgcd(G->fq_nmod,
             S->fq_nmod, T->fq_nmod, A->fq_nmod, B->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_xgcd(G->fq, S->fq, T->fq, A->fq, B->fq, ctx->ctx.fq);
   }
}

/*  Euclidean division  ******************************************************/

FQ_DEFAULT_POLY_INLINE ulong fq_default_poly_remove(fq_default_poly_t f,
                         const fq_default_poly_t g, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_remove(f->fq_zech, g->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_remove(f->fq_nmod, g->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_remove(f->fq, g->fq, ctx->ctx.fq);
}

FQ_DEFAULT_POLY_INLINE void
fq_default_poly_divrem(fq_default_poly_t Q, fq_default_poly_t R,
                   const fq_default_poly_t A, const fq_default_poly_t B,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_divrem(Q->fq_zech,
                         R->fq_zech, A->fq_zech, B->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_divrem(Q->fq_nmod,
                         R->fq_nmod, A->fq_nmod, B->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_divrem(Q->fq, R->fq, A->fq, B->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void
fq_default_poly_rem(fq_default_poly_t R,
             const fq_default_poly_t A, const fq_default_poly_t B,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_rem(R->fq_zech, A->fq_zech, B->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_rem(R->fq_nmod, A->fq_nmod, B->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_rem(R->fq, A->fq, B->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void
fq_default_poly_inv_series(fq_default_poly_t Qinv,
                               const fq_default_poly_t Q, slong n,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_inv_series(Qinv->fq_zech, Q->fq_zech, n, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_inv_series(Qinv->fq_nmod, Q->fq_nmod, n, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_inv_series(Qinv->fq, Q->fq, n, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE void fq_default_poly_div_series(fq_default_poly_t Q, 
               const fq_default_poly_t A, const fq_default_poly_t B, 
                                           slong n, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_div_series(Q->fq_zech,
                                  A->fq_zech, B->fq_zech, n, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_div_series(Q->fq_nmod,
                                  A->fq_nmod, B->fq_nmod, n, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_div_series(Q->fq, A->fq, B->fq, n, ctx->ctx.fq);
   }
}

/*  Divisibility testing  ****************************************************/

FQ_DEFAULT_POLY_INLINE int fq_default_poly_divides(fq_default_poly_t Q,
                      const fq_default_poly_t A, const fq_default_poly_t B,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_divides(Q->fq_zech,
                                     A->fq_zech, B->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_divides(Q->fq_nmod,
                                     A->fq_nmod, B->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_divides(Q->fq, A->fq, B->fq, ctx->ctx.fq);
}

/*  Derivative  **************************************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_derivative(fq_default_poly_t rop,
                        const fq_default_poly_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_derivative(rop->fq_zech, op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_derivative(rop->fq_nmod, op->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_derivative(rop->fq, op->fq, ctx->ctx.fq);
   }
}

/*  Evaluation  **************************************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_evaluate_fq_default(fq_default_t res,
                           const fq_default_poly_t f, const fq_default_t a,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_evaluate_fq_zech(res->fq_zech,
		                     f->fq_zech, a->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_evaluate_fq_nmod(res->fq_nmod,
		                     f->fq_nmod, a->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_evaluate_fq(res->fq, f->fq, a->fq, ctx->ctx.fq);
   }
}

/*  Composition  *************************************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_compose(fq_default_poly_t rop,
            const fq_default_poly_t op1, const fq_default_poly_t op2,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_compose(rop->fq_zech,
                                 op1->fq_zech, op2->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_compose(rop->fq_nmod,
                                 op1->fq_nmod, op2->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_compose(rop->fq, op1->fq, op2->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_POLY_INLINE
void fq_default_poly_compose_mod(fq_default_poly_t res,
            const fq_default_poly_t poly1, const fq_default_poly_t poly2,
                     const fq_default_poly_t poly3, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_poly_compose_mod(res->fq_zech,
             poly1->fq_zech, poly2->fq_zech, poly3->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_poly_compose_mod(res->fq_nmod,
             poly1->fq_nmod, poly2->fq_nmod, poly3->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_poly_compose_mod(res->fq,
                                 poly1->fq, poly2->fq, poly3->fq, ctx->ctx.fq);
   }
}

/*  Input and output  ********************************************************/

FQ_DEFAULT_POLY_INLINE int fq_default_poly_fprint_pretty(FILE * file,
                             const fq_default_poly_t poly, const char *x,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_fprint_pretty(file,
                                           poly->fq_zech, x, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_fprint_pretty(file,
                                           poly->fq_nmod, x, ctx->ctx.fq_nmod);
   }
   return fq_poly_fprint_pretty(file, poly->fq, x, ctx->ctx.fq);
}

FQ_DEFAULT_POLY_INLINE int fq_default_poly_fprint(FILE * file,
                      const fq_default_poly_t poly, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_fprint(file, poly->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_fprint(file, poly->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_fprint(file, poly->fq, ctx->ctx.fq);
}

FQ_DEFAULT_POLY_INLINE int
fq_default_poly_print(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_print(poly->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_print(poly->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_print(poly->fq, ctx->ctx.fq);
}


FQ_DEFAULT_POLY_INLINE int
fq_default_poly_print_pretty(const fq_default_poly_t poly,
                                     const char *x, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_print_pretty(poly->fq_zech, x, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_print_pretty(poly->fq_nmod, x, ctx->ctx.fq_nmod);
   }
   return fq_poly_print_pretty(poly->fq, x, ctx->ctx.fq);
}

FQ_DEFAULT_POLY_INLINE
char * fq_default_poly_get_str_pretty(const fq_default_poly_t poly,
                                     const char *x, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_get_str_pretty(poly->fq_zech, x, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_get_str_pretty(poly->fq_nmod, x, ctx->ctx.fq_nmod);
   }
   return fq_poly_get_str_pretty(poly->fq, x, ctx->ctx.fq);
}

FQ_DEFAULT_POLY_INLINE
char * fq_default_poly_get_str(const fq_default_poly_t poly,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_poly_get_str(poly->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_poly_get_str(poly->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_poly_get_str(poly->fq, ctx->ctx.fq);
}

#include "fq_default_mat.h"

/* Characteristic polynomial *************************************************/

FQ_DEFAULT_POLY_INLINE
void fq_default_mat_charpoly(fq_default_poly_t p,
                          const fq_default_mat_t M, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_mat_charpoly(p->fq_zech, M->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_mat_charpoly(p->fq_nmod, M->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_mat_charpoly(p->fq, M->fq, ctx->ctx.fq);
   }
}

/* Minimal polynomial ********************************************************/

FQ_DEFAULT_POLY_INLINE 
void fq_default_mat_minpoly(fq_default_poly_t p, 
                          const fq_default_mat_t X, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_mat_minpoly(p->fq_zech, X->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_mat_minpoly(p->fq_nmod, X->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_mat_minpoly(p->fq, X->fq, ctx->ctx.fq);
   }
}

#ifdef __cplusplus
}
#endif

#endif
