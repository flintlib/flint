/*
    Copyright (C) 2021 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_DEFAULT_H
#define FQ_DEFAULT_H

#ifdef FQ_DEFAULT_INLINES_C
#define FQ_DEFAULT_INLINE FLINT_DLL
#else
#define FQ_DEFAULT_INLINE static __inline__
#endif

#include "ulong_extras.h"
#include "fmpz.h"
#include "fq.h"
#include "fq_nmod.h"
#include "fq_zech.h"

/* Data types and context ****************************************************/

#ifdef __cplusplus
extern "C" {
#endif

typedef union fq_default_struct
{
   fq_t fq;
   fq_nmod_t fq_nmod;
   fq_zech_t fq_zech;
} fq_default_struct;

typedef fq_default_struct fq_default_t[1];

typedef struct
{
   int type;
   union ctx
   {
      fq_ctx_t fq;
      fq_nmod_ctx_t fq_nmod;
      fq_zech_ctx_t fq_zech;
   } ctx;
} fq_default_ctx_struct;

typedef fq_default_ctx_struct fq_default_ctx_t[1];

FQ_DEFAULT_INLINE void fq_default_ctx_init_type(fq_default_ctx_t ctx,
                            const fmpz_t p, slong d, const char *var, int type)
{
   int bits = fmpz_bits(p);

   if (type == 1 ||
       (type == 0 && bits*d <= 16 &&
	n_pow(fmpz_get_ui(p), d) < (UWORD(1) << 16)))
   {
      ctx->type = 1;
      fq_zech_ctx_init(ctx->ctx.fq_zech, p, d, var);
   } else if (type == 2 || (type == 0 && fmpz_abs_fits_ui(p)))
   {
      ctx->type = 2;
      fq_nmod_ctx_init(ctx->ctx.fq_nmod, p, d, var);
   } else
   {
      ctx->type = 3;
      fq_ctx_init(ctx->ctx.fq, p, d, var);
   }
}

FQ_DEFAULT_INLINE void fq_default_ctx_init(fq_default_ctx_t ctx,
                                      const fmpz_t p, slong d, const char *var)
{
   fq_default_ctx_init_type(ctx, p, d, var, 0);
}

FLINT_DLL
void fq_default_ctx_init_modulus_type(fq_default_ctx_t ctx,
                const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx,
                                                   const char * var, int type);
FLINT_DLL
void fq_default_ctx_init_modulus(fq_default_ctx_t ctx,
      const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var);

FLINT_DLL
void fq_default_ctx_init_modulus_nmod_type(fq_default_ctx_t ctx,
                        const nmod_poly_t modulus, const char * var, int type);

FLINT_DLL
void fq_default_ctx_init_modulus_nmod(fq_default_ctx_t ctx,
                                  const nmod_poly_t modulus, const char * var);

FQ_DEFAULT_INLINE void fq_default_ctx_clear(fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_ctx_clear(ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_ctx_clear(ctx->ctx.fq_nmod);
   } else
   {
      fq_ctx_clear(ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE int fq_default_ctx_type(const fq_default_ctx_t ctx)
{
   return ctx->type;
}

FQ_DEFAULT_INLINE slong fq_default_ctx_degree(const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_ctx_degree(ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_ctx_degree(ctx->ctx.fq_nmod);
   } else
   {
      return fq_ctx_degree(ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_ctx_prime(fmpz_t prime,
                                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fmpz_set(prime, fq_zech_ctx_prime(ctx->ctx.fq_zech));
   } else if (ctx->type == 2)
   {
      fmpz_set(prime, fq_nmod_ctx_prime(ctx->ctx.fq_nmod));
   } else
   {
      fmpz_set(prime, fq_ctx_prime(ctx->ctx.fq));
   }
}

FLINT_DLL void fq_default_ctx_modulus(fmpz_mod_poly_t p,
		                                   const fq_default_ctx_t ctx);

FQ_DEFAULT_INLINE void fq_default_ctx_order(fmpz_t f,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_ctx_order(f, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_ctx_order(f, ctx->ctx.fq_nmod);
   } else
   {
      fq_ctx_order(f, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE int fq_default_ctx_fprint(FILE * file,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_ctx_fprint(file, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_ctx_fprint(file, ctx->ctx.fq_nmod);
   }
   return fq_ctx_fprint(file, ctx->ctx.fq);
}

FQ_DEFAULT_INLINE void fq_default_ctx_print(const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_ctx_print(ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_ctx_print(ctx->ctx.fq_nmod);
   } else
   {
      fq_ctx_print(ctx->ctx.fq);
   }
}


/* Memory managment  *********************************************************/

FQ_DEFAULT_INLINE void fq_default_init(fq_default_t rop,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_init(rop->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_init(rop->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_init(rop->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_init2(fq_default_t rop,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_init2(rop->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_init2(rop->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_init2(rop->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_clear(fq_default_t rop,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_clear(rop->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_clear(rop->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_clear(rop->fq, ctx->ctx.fq);
   }
}

/* Predicates ****************************************************************/

FQ_DEFAULT_INLINE int fq_default_is_invertible(const fq_default_t op,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_is_invertible(op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_is_invertible(op->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_is_invertible(op->fq, ctx->ctx.fq);
}

/* Basic arithmetic **********************************************************/

FQ_DEFAULT_INLINE void fq_default_add(fq_default_t rop, const fq_default_t op1,
		            const fq_default_t op2, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_add(rop->fq_zech, op1->fq_zech, op2->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_add(rop->fq_nmod, op1->fq_nmod, op2->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_add(rop->fq, op1->fq, op2->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_sub(fq_default_t rop, const fq_default_t op1,
		            const fq_default_t op2, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_sub(rop->fq_zech, op1->fq_zech, op2->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_sub(rop->fq_nmod, op1->fq_nmod, op2->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_sub(rop->fq, op1->fq, op2->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_sub_one(fq_default_t rop,
		            const fq_default_t op1, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_sub_one(rop->fq_zech, op1->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_sub_one(rop->fq_nmod, op1->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_sub_one(rop->fq, op1->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_neg(fq_default_t rop,
		            const fq_default_t op1, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_neg(rop->fq_zech, op1->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_neg(rop->fq_nmod, op1->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_neg(rop->fq, op1->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_mul(fq_default_t rop, const fq_default_t op1,
                            const fq_default_t op2, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_mul(rop->fq_zech, op1->fq_zech, op2->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_mul(rop->fq_nmod, op1->fq_nmod, op2->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_mul(rop->fq, op1->fq, op2->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_mul_fmpz(fq_default_t rop,
             const fq_default_t op, const fmpz_t x, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_mul_fmpz(rop->fq_zech, op->fq_zech, x, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_mul_fmpz(rop->fq_nmod, op->fq_nmod, x, ctx->ctx.fq_nmod);
   } else
   {
      fq_mul_fmpz(rop->fq, op->fq, x, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_mul_si(fq_default_t rop,
		    const fq_default_t op, slong x, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_mul_si(rop->fq_zech, op->fq_zech, x, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_mul_si(rop->fq_nmod, op->fq_nmod, x, ctx->ctx.fq_nmod);
   } else
   {
      fq_mul_si(rop->fq, op->fq, x, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_mul_ui(fq_default_t rop,
		    const fq_default_t op, ulong x, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_mul_ui(rop->fq_zech, op->fq_zech, x, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_mul_ui(rop->fq_nmod, op->fq_nmod, x, ctx->ctx.fq_nmod);
   } else
   {
      fq_mul_ui(rop->fq, op->fq, x, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_sqr(fq_default_t rop,
                             const fq_default_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_sqr(rop->fq_zech, op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_sqr(rop->fq_nmod, op->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_sqr(rop->fq, op->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_inv(fq_default_t rop,
		            const fq_default_t op1, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_inv(rop->fq_zech, op1->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_inv(rop->fq_nmod, op1->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_inv(rop->fq, op1->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_div(fq_default_t rop, fq_default_t op1,
                                  fq_default_t op2, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_div(rop->fq_zech, op1->fq_zech, op2->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_div(rop->fq_nmod, op1->fq_nmod, op2->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_div(rop->fq, op1->fq, op2->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_pow(fq_default_t rop,
	    const fq_default_t op1, const fmpz_t e, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_pow(rop->fq_zech, op1->fq_zech, e, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_pow(rop->fq_nmod, op1->fq_nmod, e, ctx->ctx.fq_nmod);
   } else
   {
      fq_pow(rop->fq, op1->fq, e, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_pow_ui(fq_default_t rop,
              const fq_default_t op, const ulong e, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_pow_ui(rop->fq_zech, op->fq_zech, e, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_pow_ui(rop->fq_nmod, op->fq_nmod, e, ctx->ctx.fq_nmod);
   } else
   {
      fq_pow_ui(rop->fq, op->fq, e, ctx->ctx.fq);
   }
}

/* Roots *********************************************************************/

FQ_DEFAULT_INLINE int fq_default_sqrt(fq_default_t rop,
		             const fq_default_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_sqrt(rop->fq_zech, op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_sqrt(rop->fq_nmod, op->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_sqrt(rop->fq, op->fq, ctx->ctx.fq);
}

FQ_DEFAULT_INLINE void fq_default_pth_root(fq_default_t rop,
		            const fq_default_t op1, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_pth_root(rop->fq_zech, op1->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_pth_root(rop->fq_nmod, op1->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_pth_root(rop->fq, op1->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE int fq_default_is_square(const fq_default_t op,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_is_square(op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_is_square(op->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_is_square(op->fq, ctx->ctx.fq);
}

/* Randomisation *************************************************************/

FQ_DEFAULT_INLINE void fq_default_randtest(fq_default_t rop,
		                flint_rand_t state, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_randtest(rop->fq_zech, state, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_randtest(rop->fq_nmod, state, ctx->ctx.fq_nmod);
   } else
   {
      fq_randtest(rop->fq, state, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_randtest_not_zero(fq_default_t rop,
		                flint_rand_t state, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_randtest_not_zero(rop->fq_zech, state, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_randtest_not_zero(rop->fq_nmod, state, ctx->ctx.fq_nmod);
   } else
   {
      fq_randtest_not_zero(rop->fq, state, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_rand(fq_default_t rop,
		                flint_rand_t state, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_rand(rop->fq_zech, state, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_rand(rop->fq_nmod, state, ctx->ctx.fq_nmod);
   } else
   {
      fq_rand(rop->fq, state, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_rand_not_zero(fq_default_t rop,
		                flint_rand_t state, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_rand_not_zero(rop->fq_zech, state, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_rand_not_zero(rop->fq_nmod, state, ctx->ctx.fq_nmod);
   } else
   {
      fq_rand_not_zero(rop->fq, state, ctx->ctx.fq);
   }
}

/* Comparison ****************************************************************/

FQ_DEFAULT_INLINE int fq_default_equal(const fq_default_t op1,
		            const fq_default_t op2, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_equal(op1->fq_zech, op2->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_equal(op1->fq_nmod, op2->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_equal(op1->fq, op2->fq, ctx->ctx.fq);
}

FQ_DEFAULT_INLINE int fq_default_is_zero(const fq_default_t op,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_is_zero(op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_is_zero(op->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_is_zero(op->fq, ctx->ctx.fq);
}

FQ_DEFAULT_INLINE int fq_default_is_one(const fq_default_t op,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_is_one(op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_is_one(op->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_is_one(op->fq, ctx->ctx.fq);
}

/* Assignments and conversions ***********************************************/

FQ_DEFAULT_INLINE void fq_default_set(fq_default_t rop,
		             const fq_default_t op, const fq_default_ctx_t ctx)
{

   if (ctx->type == 1)
   {
      fq_zech_set(rop->fq_zech, op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_set(rop->fq_nmod, op->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_set(rop->fq, op->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_set_fmpz(fq_default_t rop,
		                    const fmpz_t x, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_set_fmpz(rop->fq_zech, x, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_set_fmpz(rop->fq_nmod, x, ctx->ctx.fq_nmod);
   } else
   {
      fq_set_fmpz(rop->fq, x, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_set_ui(fq_default_t rop,
		                     const ulong x, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_set_ui(rop->fq_zech, x, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_set_ui(rop->fq_nmod, x, ctx->ctx.fq_nmod);
   } else
   {
      fq_set_ui(rop->fq, x, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_set_si(fq_default_t rop,
		                     const slong x, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_set_si(rop->fq_zech, x, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_set_si(rop->fq_nmod, x, ctx->ctx.fq_nmod);
   } else
   {
      fq_set_si(rop->fq, x, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_zero(fq_default_t rop, 
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_zero(rop->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_zero(rop->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_zero(rop->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_one(fq_default_t rop,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_one(rop->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_one(rop->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_one(rop->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_swap(fq_default_t op1,
		                  fq_default_t op2, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_swap(op1->fq_zech, op2->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_swap(op1->fq_nmod, op2->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_swap(op1->fq, op2->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_gen(fq_default_t rop,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_gen(rop->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_gen(rop->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_gen(rop->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_get_nmod_poly(nmod_poly_t poly,
                             const fq_default_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_get_nmod_poly(poly, op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_get_nmod_poly(poly, op->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      flint_printf("Impossible conversion\n");
      flint_abort();
   }
}

FQ_DEFAULT_INLINE void fq_default_set_nmod_poly(fq_default_t op,
                            const nmod_poly_t poly, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_set_nmod_poly(op->fq_zech, poly, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_set_nmod_poly(op->fq_nmod, poly, ctx->ctx.fq_nmod);
   } else
   {
      flint_printf("Impossible conversion\n");
      flint_abort();
   }
}

FLINT_DLL void fq_default_get_fmpz_mod_poly(fmpz_mod_poly_t poly,
                            const fq_default_t op, const fq_default_ctx_t ctx);

FLINT_DLL void fq_default_set_fmpz_mod_poly(fq_default_t op,
                       const fmpz_mod_poly_t poly, const fq_default_ctx_t ctx);

FLINT_DLL void fq_default_get_fmpz_poly(fmpz_poly_t poly,
                            const fq_default_t op, const fq_default_ctx_t ctx);

FLINT_DLL void fq_default_set_fmpz_poly(fq_default_t op,
                           const fmpz_poly_t poly, const fq_default_ctx_t ctx);

FQ_DEFAULT_INLINE void fq_default_get_coeff_fmpz(fmpz_t c,
                          fq_default_t op, slong n, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      nmod_poly_t p;
      ulong c0;
      nmod_poly_init(p, fmpz_get_ui(fq_zech_ctx_prime(ctx->ctx.fq_zech)));
      fq_zech_get_nmod_poly(p, op->fq_zech, ctx->ctx.fq_zech);
      c0 = nmod_poly_get_coeff_ui(p, n);
      fmpz_set_ui(c, c0);
      nmod_poly_clear(p);
   } else if (ctx->type == 2)
   {
      ulong c0 = nmod_poly_get_coeff_ui((nmod_poly_struct *) op->fq_nmod, n);
      fmpz_set_ui(c, c0);
   } else
   {
      fmpz_mod_ctx_t mod_ctx;
      fmpz_mod_ctx_init(mod_ctx, fq_ctx_prime(ctx->ctx.fq));
      fmpz_mod_poly_get_coeff_fmpz(c,
		                  (fmpz_mod_poly_struct *) op->fq, n, mod_ctx);
      fmpz_mod_ctx_clear(mod_ctx);
   }
}

/* Output ********************************************************************/

FQ_DEFAULT_INLINE
int fq_default_fprint(FILE * file, const fq_default_t op,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_fprint(file, op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_fprint(file, op->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_fprint(file, op->fq, ctx->ctx.fq);
}

FQ_DEFAULT_INLINE
void fq_default_print(const fq_default_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_print(op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_print(op->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_print(op->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE 
int fq_default_fprint_pretty(FILE * file, const fq_default_t op,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_fprint_pretty(file, op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_fprint_pretty(file, op->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_fprint_pretty(file, op->fq, ctx->ctx.fq);
}

FQ_DEFAULT_INLINE 
void fq_default_print_pretty(const fq_default_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_print_pretty(op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_print_pretty(op->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_print_pretty(op->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE char * fq_default_get_str(const fq_default_t op,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_get_str(op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_get_str(op->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_get_str(op->fq, ctx->ctx.fq);
}

FQ_DEFAULT_INLINE char * fq_default_get_str_pretty(const fq_default_t op,
		                                    const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      return fq_zech_get_str_pretty(op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      return fq_nmod_get_str_pretty(op->fq_nmod, ctx->ctx.fq_nmod);
   }
   return fq_get_str_pretty(op->fq, ctx->ctx.fq);
}

/* Special functions *********************************************************/

FQ_DEFAULT_INLINE void fq_default_trace(fmpz_t rop,
		             const fq_default_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_trace(rop, op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_trace(rop, op->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_trace(rop, op->fq, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_frobenius(fq_default_t rop,
		    const fq_default_t op, slong e, const fq_default_ctx_t ctx)
{

   if (ctx->type == 1)
   {
      fq_zech_frobenius(rop->fq_zech, op->fq_zech, e, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_frobenius(rop->fq_nmod, op->fq_nmod, e, ctx->ctx.fq_nmod);
   } else
   {
      fq_frobenius(rop->fq, op->fq, e, ctx->ctx.fq);
   }
}

FQ_DEFAULT_INLINE void fq_default_norm(fmpz_t rop,
		             const fq_default_t op, const fq_default_ctx_t ctx)
{
   if (ctx->type == 1)
   {
      fq_zech_norm(rop, op->fq_zech, ctx->ctx.fq_zech);
   } else if (ctx->type == 2)
   {
      fq_nmod_norm(rop, op->fq_nmod, ctx->ctx.fq_nmod);
   } else
   {
      fq_norm(rop, op->fq, ctx->ctx.fq);
   }
}

#ifdef __cplusplus
}
#endif

#endif

