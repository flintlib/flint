/*
    Copyright (C) 2012 Andres Goens
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"
#include "nmod_types.h"
#include "fmpz_mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*  Memory management ********************************************************/

void TEMPLATE(T, poly_init)(TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_init2)(TEMPLATE(T, poly_t) poly, slong alloc,
                        const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_realloc)(TEMPLATE(T, poly_t) poly, slong alloc,
                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_truncate)(TEMPLATE(T, poly_t) poly, slong len,
                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_set_trunc)(TEMPLATE(T, poly_t) poly1, TEMPLATE(T, poly_t) poly2, slong len,
                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_fit_length)(TEMPLATE(T, poly_t) poly, slong len,
                             const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_clear)(TEMPLATE(T, poly_t) poly,
                        const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_normalise)(TEMPLATE(T, poly_t) poly,
                             const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_normalise2)(const TEMPLATE(T, struct) *poly, slong *length,
                              const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_set_length)(TEMPLATE(T, poly_t) poly, slong len,
                              const TEMPLATE(T, ctx_t) ctx);

/*  Polynomial parameters  ***************************************************/

FQ_POLY_TEMPLATES_INLINE slong
TEMPLATE(T, poly_length)(const TEMPLATE(T, poly_t) poly,
                         const TEMPLATE(T, ctx_t) ctx)
{
    return poly->length;
}

FQ_POLY_TEMPLATES_INLINE slong
TEMPLATE(T, poly_degree)(const TEMPLATE(T, poly_t) poly,
                         const TEMPLATE(T, ctx_t) ctx)
{
    return poly->length - 1;
}

FQ_POLY_TEMPLATES_INLINE TEMPLATE(T, struct) *
TEMPLATE(T, poly_lead)(const TEMPLATE(T, poly_t) poly,
                       const TEMPLATE(T, ctx_t) ctx)
{
    return poly->length > 0 ? poly->coeffs + (poly->length - 1) : NULL;
}

/*  Randomisation  ***********************************************************/

void TEMPLATE(T, poly_randtest)(TEMPLATE(T, poly_t) f, flint_rand_t state,
                           slong len, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_randtest_not_zero)(TEMPLATE(T, poly_t) f, flint_rand_t state,
                                    slong len, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_randtest_monic) (TEMPLATE(T, poly_t) f, flint_rand_t state,
                                  slong len, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_randtest_irreducible) (TEMPLATE(T, poly_t) f,
                                        flint_rand_t state, slong len,
                                        const TEMPLATE(T, ctx_t) ctx);


/*  Assignment and basic manipulation  ***************************************/

void _TEMPLATE(T, poly_set)(TEMPLATE(T, struct) *rop, const TEMPLATE(T, struct) *op,
                       slong len, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_set)(TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op,
                      const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE3(T, poly_set, T)(TEMPLATE(T, poly_t) poly, const TEMPLATE(T, t) c,
                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_set_fmpz_mod_poly)(TEMPLATE(T, poly_t) rop,
                                                   const fmpz_mod_poly_t op,
                                                   const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_set_nmod_poly)(TEMPLATE(T, poly_t) rop,
                                               const nmod_poly_t op,
                                               const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_swap)(TEMPLATE(T, poly_t) op1, TEMPLATE(T, poly_t) op2,
                       const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_zero)(TEMPLATE(T, struct) *rop, slong len, const TEMPLATE(T, ctx_t) ctx);
void TEMPLATE(T, poly_zero)(TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx);
void TEMPLATE(T, poly_one)(TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_gen)(TEMPLATE(T, poly_t) f, const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_make_monic)(TEMPLATE(T, struct) *rop,
                              const TEMPLATE(T, struct) *op, slong length,
                              const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_make_monic)(TEMPLATE(T, poly_t) rop,
                                  const TEMPLATE(T, poly_t) op,
                                  const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_reverse)(TEMPLATE(T, struct) * res,
                           const TEMPLATE(T, struct) * poly, slong len, slong n,
                           const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_reverse)(TEMPLATE(T, poly_t) res,
                          const TEMPLATE(T, poly_t) poly, slong n,
                          const TEMPLATE(T, ctx_t) ctx);

ulong TEMPLATE(T, poly_deflation)(const TEMPLATE(T, poly_t) input,
                            const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_deflate)(TEMPLATE(T, poly_t) result,
                          const TEMPLATE(T, poly_t) input, ulong deflation,
                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_inflate)(TEMPLATE(T, poly_t) result,
                          const TEMPLATE(T, poly_t) input, ulong inflation,
                          const TEMPLATE(T, ctx_t) ctx);

/*  Getting and setting coefficients  ****************************************/

void TEMPLATE(T, poly_get_coeff)(TEMPLATE(T, t) x, const TEMPLATE(T, poly_t) poly,
                            slong n, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_set_coeff)(TEMPLATE(T, poly_t) poly, slong n,
                                 const TEMPLATE(T, t) x,
                                 const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_set_coeff_fmpz)(TEMPLATE(T, poly_t) poly, slong n,
                                 const fmpz_t x, const TEMPLATE(T, ctx_t) ctx);

int TEMPLATE(T, poly_is_gen)(const TEMPLATE(T, poly_t) poly,
                                const TEMPLATE(T, ctx_t) ctx);

/*  Comparison  **************************************************************/

int TEMPLATE(T, poly_equal)(const TEMPLATE(T, poly_t) poly1,
                            const TEMPLATE(T, poly_t) poly2,
                            const TEMPLATE(T, ctx_t) ctx);

int TEMPLATE(T, poly_equal_trunc)(const TEMPLATE(T, poly_t) poly1,
                            const TEMPLATE(T, poly_t) poly2,
                            slong n, const TEMPLATE(T, ctx_t) ctx);

FQ_POLY_TEMPLATES_INLINE int
TEMPLATE(T, poly_is_zero)(const TEMPLATE(T, poly_t) poly,
                          const TEMPLATE(T, ctx_t) ctx)
{
    return (poly->length == 0);
}

int TEMPLATE(T, poly_is_one)(const TEMPLATE(T, poly_t) op, const TEMPLATE(T, ctx_t) ctx);
int TEMPLATE(T, poly_is_unit)(const TEMPLATE(T, poly_t) op, const TEMPLATE(T, ctx_t) ctx);
int TEMPLATE3(T, poly_equal, T)(const TEMPLATE(T, poly_t) poly, const TEMPLATE(T, t) c, const TEMPLATE(T, ctx_t) ctx);

/*  Addition and subtraction  ************************************************/

void _TEMPLATE(T, poly_add)(TEMPLATE(T, struct) *res,
                       const TEMPLATE(T, struct) *poly1, slong len1,
                       const TEMPLATE(T, struct) *poly2, slong len2,
                       const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_add)(TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op1,
                      const TEMPLATE(T, poly_t) op2,
                      const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_add_si)(TEMPLATE(T, poly_t) rop,
                                    const TEMPLATE(T, poly_t) op1, slong c,
                                                 const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_add_series)(TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op1,
                      const TEMPLATE(T, poly_t) op2,
                      slong n, const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_sub)(TEMPLATE(T, struct) *res,
                       const TEMPLATE(T, struct) *poly1, slong len1,
                       const TEMPLATE(T, struct) *poly2, slong len2,
                       const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_sub)(TEMPLATE(T, poly_t) rop,
                      const TEMPLATE(T, poly_t) op1,
                      const TEMPLATE(T, poly_t) op2,
                      const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_sub_series)(TEMPLATE(T, poly_t) rop,
                      const TEMPLATE(T, poly_t) op1,
                      const TEMPLATE(T, poly_t) op2,
                      slong n, const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_neg)(TEMPLATE(T, struct) *rop,
                       const TEMPLATE(T, struct) *op, slong len,
                       const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_neg)(TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op,
                      const TEMPLATE(T, ctx_t) ctx);

/*  Scalar multiplication and division  **************************************/

void _TEMPLATE3(T, poly_scalar_mul, T)(TEMPLATE(T, struct) *rop,
                                  const TEMPLATE(T, struct) *op, slong len,
                                  const TEMPLATE(T, t) x,
                                  const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE3(T, poly_scalar_mul, T)(TEMPLATE(T, poly_t) rop,
                                 const TEMPLATE(T, poly_t) op,
                                 const TEMPLATE(T, t) x,
                                 const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE3(T, poly_scalar_div, T)(TEMPLATE(T, struct) *rop,
                                  const TEMPLATE(T, struct) *op, slong len,
                                  const TEMPLATE(T, t) x,
                                  const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE3(T, poly_scalar_div, T)(TEMPLATE(T, poly_t) rop,
                                 const TEMPLATE(T, poly_t) op,
                                 const TEMPLATE(T, t) x,
                                 const TEMPLATE(T, ctx_t) ctx);

void  _TEMPLATE3(T, poly_scalar_addmul, T)(TEMPLATE(T, struct) *rop,
                                     const TEMPLATE(T, struct) *op, slong len,
                                     const TEMPLATE(T, t) x,
                                     const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE3(T, poly_scalar_addmul, T)(TEMPLATE(T, poly_t) rop,
                                    const TEMPLATE(T, poly_t) op,
                                    const TEMPLATE(T, t) x,
                                    const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE3(T, poly_scalar_submul, T)(TEMPLATE(T, struct) *rop,
                                     const TEMPLATE(T, struct) *op, slong len,
                                     const TEMPLATE(T, t) x,
                                     const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE3(T, poly_scalar_submul, T)(TEMPLATE(T, poly_t) rop,
                                    const TEMPLATE(T, poly_t) op,
                                    const TEMPLATE(T, t) x,
                                    const TEMPLATE(T, ctx_t) ctx);

/*  Multiplication  **********************************************************/

void _TEMPLATE(T, poly_mul_classical)(TEMPLATE(T, struct) *rop,
                                 const TEMPLATE(T, struct) *op1, slong len1,
                                 const TEMPLATE(T, struct) *op2, slong len2,
                                 const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_mul_classical)(TEMPLATE(T, poly_t) rop,
                                const TEMPLATE(T, poly_t) op1,
                                const TEMPLATE(T, poly_t) op2,
                                const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_mul_reorder)(TEMPLATE(T, struct) *rop,
                               const TEMPLATE(T, struct) *op1, slong len1,
                               const TEMPLATE(T, struct) *op2, slong len2,
                               const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_mul_reorder)(TEMPLATE(T, poly_t) rop,
                              const TEMPLATE(T, poly_t) op1,
                              const TEMPLATE(T, poly_t) op2,
                              const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_mul_univariate)(TEMPLATE(T, struct) *rop,
                          const TEMPLATE(T, struct) *op1, slong len1,
                          const TEMPLATE(T, struct) *op2, slong len2,
                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_mul_univariate)(TEMPLATE(T, poly_t) rop,
                         const TEMPLATE(T, poly_t) op1,
                         const TEMPLATE(T, poly_t) op2,
                         const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_mul_KS)(TEMPLATE(T, struct) *rop,
                          const TEMPLATE(T, struct) *op1, slong len1,
                          const TEMPLATE(T, struct) *op2, slong len2,
                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_mul_KS)(TEMPLATE(T, poly_t) rop,
                         const TEMPLATE(T, poly_t) op1,
                         const TEMPLATE(T, poly_t) op2,
                         const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_mul)(TEMPLATE(T, struct) *rop,
                       const TEMPLATE(T, struct) *op1, slong len1,
                       const TEMPLATE(T, struct) *op2, slong len2,
                       const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_mul)(TEMPLATE(T, poly_t) rop,
                      const TEMPLATE(T, poly_t) op1,
                      const TEMPLATE(T, poly_t) op2,
                      const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_mullow_classical)(TEMPLATE(T, struct) *rop,
                                    const TEMPLATE(T, struct) *op1, slong len1,
                                    const TEMPLATE(T, struct) *op2, slong len2,
                                    slong n,
                                    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_mullow_classical)(TEMPLATE(T, poly_t) rop,
                                   const TEMPLATE(T, poly_t) op1,
                                   const TEMPLATE(T, poly_t) op2,
                                   slong n, const
                                   TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_mullow_KS)(TEMPLATE(T, struct) *rop,
                             const TEMPLATE(T, struct) *op1, slong len1,
                             const TEMPLATE(T, struct) *op2, slong len2,
                             slong n,
                             const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_mullow_KS)(TEMPLATE(T, poly_t) rop,
                            const TEMPLATE(T, poly_t) op1,
                            const TEMPLATE(T, poly_t) op2,
                            slong n,
                            const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_mullow_univariate)(TEMPLATE(T, struct) *rop,
                             const TEMPLATE(T, struct) *op1, slong len1,
                             const TEMPLATE(T, struct) *op2, slong len2,
                             slong n,
                             const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_mullow_univariate)(TEMPLATE(T, poly_t) rop,
                            const TEMPLATE(T, poly_t) op1,
                            const TEMPLATE(T, poly_t) op2,
                            slong n,
                            const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_mullow)(TEMPLATE(T, struct) *rop,
                          const TEMPLATE(T, struct) *op1, slong len1,
                          const TEMPLATE(T, struct) *op2, slong len2,
                          slong n,
                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_mullow)(TEMPLATE(T, poly_t) rop,
                         const TEMPLATE(T, poly_t) op1,
                         const TEMPLATE(T, poly_t) op2,
                         slong n,
                         const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_mulhigh_classical)(
    TEMPLATE(T, struct)* rop,
    const TEMPLATE(T, struct)* op1, slong len1,
    const TEMPLATE(T, struct)* op2, slong len2,
    slong start,
    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_mulhigh_classical)(TEMPLATE(T, poly_t) rop,
                                    const TEMPLATE(T, poly_t) op1,
                                    const TEMPLATE(T, poly_t) op2,
                                    slong start,
                                    const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_mulhigh)(TEMPLATE(T, struct)* res,
                           const TEMPLATE(T, struct)* poly1, slong len1,
                           const TEMPLATE(T, struct)* poly2, slong len2,
                           slong n, TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_mulhigh)(TEMPLATE(T, poly_t) rop,
                          const TEMPLATE(T, poly_t) op1,
                          const TEMPLATE(T, poly_t) op2, slong start,
                          const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_mulmod)(TEMPLATE(T, struct) * res,
                          const TEMPLATE(T, struct) * poly1, slong len1,
                          const TEMPLATE(T, struct) * poly2, slong len2,
                          const TEMPLATE(T, struct) * f, slong lenf,
                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_mulmod)(TEMPLATE(T, poly_t) res,
                         const TEMPLATE(T, poly_t) poly1,
                         const TEMPLATE(T, poly_t) poly2,
                         const TEMPLATE(T, poly_t) f,
                         const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_mulmod_preinv)(TEMPLATE(T, struct) *res,
                                 const TEMPLATE(T, struct) *poly1, slong len1,
                                 const TEMPLATE(T, struct)  *poly2, slong len2,
                                 const TEMPLATE(T, struct) *f, slong lenf,
                                 const TEMPLATE(T, struct) *finv, slong lenfinv,
                                 const TEMPLATE(T, ctx_t) ctx);
void TEMPLATE(T, poly_mulmod_preinv)(TEMPLATE(T, poly_t) res,
                                const TEMPLATE(T, poly_t) poly1,
                                const TEMPLATE(T, poly_t) poly2,
                                const TEMPLATE(T, poly_t) f,
                                const TEMPLATE(T, poly_t) finv,
                                const TEMPLATE(T, ctx_t) ctx);

/* Squaring ******************************************************************/

void _TEMPLATE(T, poly_sqr_classical)(TEMPLATE(T, struct) *rop,
                                 const TEMPLATE(T, struct) *op, slong len,
                                 const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_sqr_classical)(TEMPLATE(T, poly_t) rop,
                                const TEMPLATE(T, poly_t) op,
                                const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_sqr_reorder)(TEMPLATE(T, struct) *rop,
                               const TEMPLATE(T, struct) *op, slong len,
                               const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_sqr_reorder)(TEMPLATE(T, poly_t) rop,
                              const TEMPLATE(T, poly_t) op,
                              const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_sqr_KS)(TEMPLATE(T, struct) *rop,
                          const TEMPLATE(T, struct) *op, slong len,
                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_sqr_KS)(TEMPLATE(T, poly_t) rop,
                         const TEMPLATE(T, poly_t) op,
                         const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_sqr)(TEMPLATE(T, struct) *rop,
                       const TEMPLATE(T, struct) *op, slong len,
                       const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_sqr)(TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op,
                      const TEMPLATE(T, ctx_t) ctx);

/*  Powering  ****************************************************************/

void _TEMPLATE(T, poly_pow)(TEMPLATE(T, struct) *rop,
                       const TEMPLATE(T, struct) *op, slong len,
                       ulong e,
                       const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_pow)(TEMPLATE(T, poly_t) rop,
                      const TEMPLATE(T, poly_t) op,
                      ulong e,
                      const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_pow_trunc_binexp) (TEMPLATE(T, struct) * res,
                                   const TEMPLATE(T, struct) * poly,  ulong e,
			            slong trunc, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_pow_trunc_binexp) (TEMPLATE(T, poly_t) res,
                                     const TEMPLATE(T, poly_t) poly, ulong e,
                                    slong trunc, const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_pow_trunc) (TEMPLATE(T, struct) * res,
                                   const TEMPLATE(T, struct) * poly,  ulong e,
                                    slong trunc, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_pow_trunc) (TEMPLATE(T, poly_t) res,
                                     const TEMPLATE(T, poly_t) poly, ulong e,                                                 slong trunc, const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_powmod_fmpz_binexp)(TEMPLATE(T, struct) * res,
                                      const TEMPLATE(T, struct) * poly,
                                      const fmpz_t e,
                                      const TEMPLATE(T, struct) * f, slong lenf,
                                      const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_powmod_fmpz_binexp)(TEMPLATE(T, poly_t) res,
                                     const TEMPLATE(T, poly_t) poly,
                                     const fmpz_t e,
                                     const TEMPLATE(T, poly_t) f,
                                     const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_powmod_fmpz_binexp_preinv)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly,
    const fmpz_t e,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) *finv, slong lenfinv,
    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_powmod_fmpz_binexp_preinv)(
    TEMPLATE(T, poly_t) res,
    const TEMPLATE(T, poly_t) poly, const fmpz_t e,
    const TEMPLATE(T, poly_t) f, const TEMPLATE(T, poly_t) finv,
    const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_powmod_ui_binexp)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly,
    ulong e,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_powmod_ui_binexp)(
    TEMPLATE(T, poly_t) res,
    const TEMPLATE(T, poly_t) poly, ulong e,
    const TEMPLATE(T, poly_t) f,
    const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_powmod_ui_binexp_preinv)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly,
    ulong e,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * finv, slong lenfinv,
    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_powmod_ui_binexp_preinv)(
    TEMPLATE(T, poly_t) res,
    const TEMPLATE(T, poly_t) poly,
    ulong e,
    const TEMPLATE(T, poly_t) f,
    const TEMPLATE(T, poly_t) finv,
    const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly,
    const fmpz_t e, ulong k,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * finv, slong lenfinv,
    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(
    TEMPLATE(T, poly_t) res,
    const TEMPLATE(T, poly_t) poly,
    const fmpz_t e, ulong k,
    const TEMPLATE(T, poly_t) f,
    const TEMPLATE(T, poly_t) finv,
    const TEMPLATE(T, ctx_t) ctx);


void _TEMPLATE(T, poly_powmod_x_fmpz_preinv)(
    TEMPLATE(T, struct) * res, const fmpz_t e,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct)* finv, slong lenfinv,
    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_powmod_x_fmpz_preinv)(
    TEMPLATE(T, poly_t) res, const fmpz_t e,
    const TEMPLATE(T, poly_t) f,
    const TEMPLATE(T, poly_t) finv,
    const TEMPLATE(T, ctx_t) ctx);

/*  Shifting  ****************************************************************/

void _TEMPLATE(T, poly_shift_left)(TEMPLATE(T, struct) *rop,
                              const TEMPLATE(T, struct) *op, slong len,
                              slong n,
                              const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_shift_left)(TEMPLATE(T, poly_t) rop,
                             const TEMPLATE(T, poly_t) op,
                             slong n,
                             const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_shift_right)(TEMPLATE(T, struct) *rop,
                               const TEMPLATE(T, struct) *op, slong len,
                               slong n,
                               const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_shift_right)(TEMPLATE(T, poly_t) rop,
                              const TEMPLATE(T, poly_t) op,
                              slong n,
                              const TEMPLATE(T, ctx_t) ctx);

/*  Norms  *******************************************************************/

slong _TEMPLATE(T, poly_hamming_weight)(const TEMPLATE(T, struct) *op, slong len,
                                  const TEMPLATE(T, ctx_t) ctx);

slong TEMPLATE(T, poly_hamming_weight)(const TEMPLATE(T, poly_t) op,
                                 const TEMPLATE(T, ctx_t) ctx);

/*  Greatest common divisor  *************************************************/

void TEMPLATE(T, poly_gcd_euclidean)(TEMPLATE(T, poly_t) rop,
                                const TEMPLATE(T, poly_t) op1,
                                const TEMPLATE(T, poly_t) op2,
                                const TEMPLATE(T, ctx_t) ctx);

slong _TEMPLATE(T, poly_gcd_euclidean)(TEMPLATE(T, struct)* G,
                                 const TEMPLATE(T, struct)* A, slong lenA,
                                 const TEMPLATE(T, struct)* B, slong lenB,
                                 const TEMPLATE(T, t) invB,
                                 const TEMPLATE(T, ctx_t) ctx);

slong _TEMPLATE(T, poly_gcd)(TEMPLATE(T, struct)* G,
                       const TEMPLATE(T, struct)* A, slong lenA,
                       const TEMPLATE(T, struct)* B, slong lenB,
                       const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_gcd)(TEMPLATE(T, poly_t) rop,
                      const TEMPLATE(T, poly_t) op1,
                      const TEMPLATE(T, poly_t) op2,
                      const TEMPLATE(T, ctx_t) ctx);

slong _TEMPLATE(T, poly_gcd_euclidean_f)(TEMPLATE(T, t) f, TEMPLATE(T, struct)* G,
                                   const TEMPLATE(T, struct)* A, slong lenA,
                                   const TEMPLATE(T, struct)* B, slong lenB,
                                   const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_gcd_euclidean_f)(TEMPLATE(T, t) f, TEMPLATE(T, poly_t) G,
                                  const TEMPLATE(T, poly_t) A,
                                  const TEMPLATE(T, poly_t) B,
                                  const TEMPLATE(T, ctx_t) ctx);

slong _TEMPLATE(T, poly_xgcd_euclidean_f)(TEMPLATE(T, t) f, TEMPLATE(T, struct) *G,
                                    TEMPLATE(T, struct) *S,
                                    TEMPLATE(T, struct) *T,
                                    const TEMPLATE(T, struct) *A, slong lenA,
                                    const TEMPLATE(T, struct) *B, slong lenB,
                                    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_xgcd_euclidean_f)(TEMPLATE(T, t) f, TEMPLATE(T, poly_t) G,
                                   TEMPLATE(T, poly_t) S,
                                   TEMPLATE(T, poly_t) T,
                                   const TEMPLATE(T, poly_t) A,
                                   const TEMPLATE(T, poly_t) B,
                                   const TEMPLATE(T, ctx_t) ctx);

slong _TEMPLATE(T, poly_xgcd)(TEMPLATE(T, struct) *G,
                                  TEMPLATE(T, struct) *S,
                                  TEMPLATE(T, struct) *T,
                                  const TEMPLATE(T, struct) *A, slong lenA,
                                  const TEMPLATE(T, struct) *B, slong lenB,
                                  const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_xgcd)(TEMPLATE(T, poly_t) G,
                                 TEMPLATE(T, poly_t) S, TEMPLATE(T, poly_t) T,
                                 const TEMPLATE(T, poly_t) A,
                                 const TEMPLATE(T, poly_t) B,
                                 const TEMPLATE(T, ctx_t) ctx);


/*  Euclidean division  ******************************************************/

ulong TEMPLATE(T, poly_remove)(TEMPLATE(T, poly_t) f,
                         const TEMPLATE(T, poly_t) g,
                         const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_div)(TEMPLATE(T, struct) *Q,
                               const TEMPLATE(T, struct) *A, slong lenA,
                               const TEMPLATE(T, struct) *B, slong lenB,
                               const TEMPLATE(T, t) invB,
                               const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_div)(TEMPLATE(T, poly_t) Q,
                               const TEMPLATE(T, poly_t) A,
                               const TEMPLATE(T, poly_t) B,
                               const TEMPLATE(T, ctx_t) ctx);

/* flint 2.x compatibility needed by Nemo */
void TEMPLATE(T, poly_div_basecase)(TEMPLATE(T, poly_t) Q,
                               const TEMPLATE(T, poly_t) A,
                               const TEMPLATE(T, poly_t) B,
                               const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_divrem)(TEMPLATE(T, struct) *Q,
                                   TEMPLATE(T, struct) *R,
                                   const TEMPLATE(T, struct) *A, slong lenA,
                                   const TEMPLATE(T, struct) *B, slong lenB,
                                   const TEMPLATE(T, t) invB,
                                   const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_divrem)(TEMPLATE(T, poly_t) Q, TEMPLATE(T, poly_t) R,
                                  const TEMPLATE(T, poly_t) A,
                                  const TEMPLATE(T, poly_t) B,
                                  const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_rem)(TEMPLATE(T, struct) *R,
                       const TEMPLATE(T, struct) *A, slong lenA,
                       const TEMPLATE(T, struct) *B, slong lenB,
                       const TEMPLATE(T, t) invB,
                       const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_rem)(TEMPLATE(T, poly_t) R,
                      const TEMPLATE(T, poly_t) A, const TEMPLATE(T, poly_t) B,
                      const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_inv_series_newton)(TEMPLATE(T, struct) * Qinv,
                                     const TEMPLATE(T, struct) * Q, slong n,
                                     const TEMPLATE(T, t) cinv,
                                     const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_inv_series_newton)(TEMPLATE(T, poly_t) Qinv,
                                    const TEMPLATE(T, poly_t) Q, slong n,
                                    const TEMPLATE(T, ctx_t) ctx);

FQ_POLY_TEMPLATES_INLINE void
_TEMPLATE(T, poly_inv_series)(TEMPLATE(T, struct) * Qinv,
                                     const TEMPLATE(T, struct) * Q, slong n,
                                     const TEMPLATE(T, t) cinv,
                                     const TEMPLATE(T, ctx_t) ctx)
{
   _TEMPLATE(T, poly_inv_series_newton) (Qinv, Q, n, cinv, ctx);
}

FQ_POLY_TEMPLATES_INLINE void
TEMPLATE(T, poly_inv_series)(TEMPLATE(T, poly_t) Qinv,
                                    const TEMPLATE(T, poly_t) Q, slong n,
                                    const TEMPLATE(T, ctx_t) ctx)
{
   TEMPLATE(T, poly_inv_series_newton) (Qinv, Q, n, ctx);
}

void _TEMPLATE(T, poly_div_series) (TEMPLATE(T, struct) * Q,
                          const TEMPLATE(T, struct) * A, slong Alen,
                          const TEMPLATE(T, struct) * B, slong Blen,
                          slong n, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_div_series)(TEMPLATE(T, poly_t) Q,
               const TEMPLATE(T, poly_t) A, const TEMPLATE(T, poly_t) B,
                                        slong n, const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_div_newton_n_preinv) (
    TEMPLATE(T, struct) *Q,
    const TEMPLATE(T, struct) *A, slong lenA,
    const TEMPLATE(T, struct)* B, slong lenB,
    const TEMPLATE(T, struct)* Binv, slong lenBinv,
    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_div_newton_n_preinv) (TEMPLATE(T, poly_t) Q,
                                       const TEMPLATE(T, poly_t) A,
                                       const TEMPLATE(T, poly_t) B,
                                       const TEMPLATE(T, poly_t) Binv,
                                       const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_divrem_newton_n_preinv) (
    TEMPLATE(T, struct)* Q, TEMPLATE(T, struct)* R,
    const TEMPLATE(T, struct)* A, slong lenA,
    const TEMPLATE(T, struct)* B, slong lenB,
    const TEMPLATE(T, struct)* Binv, slong lenBinv,
    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_divrem_newton_n_preinv)(
    TEMPLATE(T, poly_t) Q, TEMPLATE(T, poly_t) R,
    const TEMPLATE(T, poly_t) A, const TEMPLATE(T, poly_t) B,
    const TEMPLATE(T, poly_t) Binv, const TEMPLATE(T, ctx_t) ctx);


void _TEMPLATE(T, poly_divrem_f)(TEMPLATE(T, t) f,
                            TEMPLATE(T, struct)* Q, TEMPLATE(T, struct)* R,
                            const TEMPLATE(T, struct)* A, slong lenA,
                            const TEMPLATE(T, struct)* B, slong lenB,
                            const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_divrem_f)(TEMPLATE(T, t) f,
                           TEMPLATE(T, poly_t) Q, TEMPLATE(T, poly_t) R,
                           const TEMPLATE(T, poly_t) A,
                           const TEMPLATE(T, poly_t) B,
                           const TEMPLATE(T, ctx_t) ctx);

/*  Divisibility testing  ***************************************************/

int _TEMPLATE(T, poly_divides)(TEMPLATE(T, struct) *Q,
                           const TEMPLATE(T, struct) *A, slong lenA,
                           const TEMPLATE(T, struct) *B, slong lenB,
                           const TEMPLATE(T, t) invB,
                           const TEMPLATE(T, ctx_t) ctx);

int TEMPLATE(T, poly_divides)(TEMPLATE(T, poly_t) Q,
                          const TEMPLATE(T, poly_t) A, const TEMPLATE(T, poly_t) B,
                          const TEMPLATE(T, ctx_t) ctx);

/*  Derivative  **************************************************************/

void _TEMPLATE(T, poly_derivative)(TEMPLATE(T, struct) *rop,
                              const TEMPLATE(T, struct) *op, slong len,
                              const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_derivative)(TEMPLATE(T, poly_t) rop,
                             const TEMPLATE(T, poly_t) op,
                             const TEMPLATE(T, ctx_t) ctx);

/* Square root ***************************************************************/

void _TEMPLATE(T, poly_invsqrt_series)(TEMPLATE(T, struct) * g,
               const TEMPLATE(T, struct) * h, slong n, TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_invsqrt_series)(TEMPLATE(T, poly_t) g,
                 const TEMPLATE(T, poly_t) h, slong n, TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_sqrt_series)(TEMPLATE(T, struct) * g,
               const TEMPLATE(T, struct) * h, slong n, TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_sqrt_series)(TEMPLATE(T, poly_t) g,
                 const TEMPLATE(T, poly_t) h, slong n, TEMPLATE(T, ctx_t) ctx);

int _TEMPLATE(T, poly_sqrt)(TEMPLATE(T, struct) * s,
             const TEMPLATE(T, struct) * p, slong len, TEMPLATE(T, ctx_t) ctx);

int TEMPLATE(T, poly_sqrt)(TEMPLATE(T, poly_t) b,
                          const TEMPLATE(T, poly_t) a, TEMPLATE(T, ctx_t) ctx);

/*  Evaluation  **************************************************************/

void _TEMPLATE3(T, poly_evaluate, T)(TEMPLATE(T, t) rop,
                                     const TEMPLATE(T, struct) *op, slong len,
                                     const TEMPLATE(T, t) a,
                                     const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE3(T, poly_evaluate, T)(TEMPLATE(T, t) res,
                                    const TEMPLATE(T, poly_t) f,
                                    const TEMPLATE(T, t) a,
                                    const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE4(T, poly_evaluate, T, vec)(TEMPLATE(T, struct) * ys,
                                     const TEMPLATE(T, struct) *coeffs, slong len,
                                     const TEMPLATE(T, struct) *xs, slong n,
                                     const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE4(T, poly_evaluate, T, vec)(TEMPLATE(T, struct) * ys,
                                    const TEMPLATE(T, poly_t) poly,
                                    const TEMPLATE(T, struct) * xs, slong n,
                                    const TEMPLATE(T, ctx_t) ctx);

TEMPLATE(T, poly_struct) **
_TEMPLATE(T, poly_tree_alloc)(slong len, const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_tree_free)(TEMPLATE(T, poly_struct) ** tree, slong len,
                             const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_tree_build)(TEMPLATE(T, poly_struct) ** tree,
                              const TEMPLATE(T, struct) * roots,
                              slong len,
                              const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE4(T, poly_evaluate, T, vec_fast_precomp)
    (TEMPLATE(T, struct) * vs,
     const TEMPLATE(T, struct) * poly, slong plen,
     TEMPLATE(T, poly_struct) * const * tree, slong len,
     const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE4(T, poly_evaluate, T, vec_fast)(TEMPLATE(T, struct) * ys,
                                          const TEMPLATE(T, struct) * poly, slong plen,
                                          const TEMPLATE(T, struct) * xs, slong n,
                                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE4(T, poly_evaluate, T, vec_fast)(TEMPLATE(T, struct) * ys,
                                         const TEMPLATE(T, poly_t) poly,
                                         const TEMPLATE(T, struct) *xs, slong n,
                                         const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE4(T, poly_evaluate, T, vec_iter)(TEMPLATE(T, struct) * ys,
                                          const TEMPLATE(T, struct) * coeffs, slong len,
                                          const TEMPLATE(T, struct) * xs, slong n,
                                          const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE4(T, poly_evaluate, T, vec_iter)(TEMPLATE(T, struct) * ys,
                                         const TEMPLATE(T, poly_t) poly,
                                         const TEMPLATE(T, struct) * xs, slong n,
                                         const TEMPLATE(T, ctx_t) ctx);

/*  Composition  *************************************************************/

void _TEMPLATE(T, poly_compose)(TEMPLATE(T, struct) *rop,
                           const TEMPLATE(T, struct) *op1, slong len1,
                           const TEMPLATE(T, struct) *op2, slong len2,
                           const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_compose)(TEMPLATE(T, poly_t) rop,
                          const TEMPLATE(T, poly_t) op1,
                          const TEMPLATE(T, poly_t) op2,
                          const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_compose_mod)(TEMPLATE(T, struct) * res,
                     const TEMPLATE(T, struct) * f, slong lenf,
                     const TEMPLATE(T, struct) * g,
                     const TEMPLATE(T, struct) * h, slong lenh,
                     const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_compose_mod)(TEMPLATE(T, poly_t) res,
                              const TEMPLATE(T, poly_t) poly1,
                              const TEMPLATE(T, poly_t) poly2,
                              const TEMPLATE(T, poly_t) poly3,
                              const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_compose_mod_preinv)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * g,
    const TEMPLATE(T, struct) * h, slong lenh,
    const TEMPLATE(T, struct) * hinv, slong lenhinv,
    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_compose_mod_preinv)(TEMPLATE(T, poly_t) res,
                                     const TEMPLATE(T, poly_t) poly1,
                                     const TEMPLATE(T, poly_t) poly2,
                                     const TEMPLATE(T, poly_t) poly3,
                                     const TEMPLATE(T, poly_t) poly3inv,
                                     const TEMPLATE(T, ctx_t) ctx);


void _TEMPLATE(T, poly_compose_mod_horner)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * g,
    const TEMPLATE(T, struct) * h, slong lenh,
    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_compose_mod_horner)(TEMPLATE(T, poly_t) res,
                                     const TEMPLATE(T, poly_t) poly1,
                                     const TEMPLATE(T, poly_t) poly2,
                                     const TEMPLATE(T, poly_t) poly3,
                                     const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_compose_mod_horner_preinv)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * g,
    const TEMPLATE(T, struct) * h, slong lenh,
    const TEMPLATE(T, struct) * hinv, slong lenhinv,
    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_compose_mod_horner_preinv)(TEMPLATE(T, poly_t) res,
                                            const TEMPLATE(T, poly_t) poly1,
                                            const TEMPLATE(T, poly_t) poly2,
                                            const TEMPLATE(T, poly_t) poly3,
                                            const TEMPLATE(T, poly_t) poly3inv,
                                            const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_compose_mod_brent_kung)(TEMPLATE(T, poly_t) res,
                                         const TEMPLATE(T, poly_t) poly1,
                                         const TEMPLATE(T, poly_t) poly2,
                                         const TEMPLATE(T, poly_t) poly3,
                                         const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_compose_mod_brent_kung)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly1, slong len1,
    const TEMPLATE(T, struct) * poly2,
    const TEMPLATE(T, struct) * poly3, slong len3,
    const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_compose_mod_brent_kung_preinv)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly1, slong len1,
    const TEMPLATE(T, struct) * poly2,
    const TEMPLATE(T, struct) * poly3, slong len3,
    const TEMPLATE(T, struct) * poly3inv, slong len3inv,
    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_compose_mod_brent_kung_preinv)(
    TEMPLATE(T, poly_t) res,
    const TEMPLATE(T, poly_t) poly1,
    const TEMPLATE(T, poly_t) poly2,
    const TEMPLATE(T, poly_t) poly3,
    const TEMPLATE(T, poly_t) poly3inv,
    const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_reduce_matrix_mod_poly) (TEMPLATE(T, mat_t) A,
                                           const TEMPLATE(T, mat_t) B,
                                           const TEMPLATE(T, poly_t) f,
                                           const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_precompute_matrix) (
    TEMPLATE(T, mat_t A),
    const TEMPLATE(T, struct)* poly1,
    const TEMPLATE(T, struct)* poly2, slong len2,
    const TEMPLATE(T, struct)* poly2inv, slong len2inv,
    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_precompute_matrix) (TEMPLATE(T, mat_t A),
                                     const TEMPLATE(T, poly_t) poly1,
                                     const TEMPLATE(T, poly_t) poly2,
                                     const TEMPLATE(T, poly_t) poly2inv,
                                     const TEMPLATE(T, ctx_t) ctx);

void _TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv)(
    TEMPLATE(T, struct)* res,
    const TEMPLATE(T, struct)* poly1, slong len1,
    const TEMPLATE(T, mat_t) A,
    const TEMPLATE(T, struct)* poly3, slong len3,
    const TEMPLATE(T, struct)* poly3inv, slong len3inv,
    const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv)(
    TEMPLATE(T, poly_t) res,
    const TEMPLATE(T, poly_t) poly1,
    const TEMPLATE(T, mat_t) A,
    const TEMPLATE(T, poly_t) poly3,
    const TEMPLATE(T, poly_t) poly3inv,
    const TEMPLATE(T, ctx_t) ctx);

/*  Input and output  ********************************************************/

#ifdef FLINT_HAVE_FILE
int _TEMPLATE(T, poly_fprint_pretty)(FILE *file, const TEMPLATE(T, struct) *poly, slong len, const char *x, const TEMPLATE(T, ctx_t) ctx);
int TEMPLATE(T, poly_fprint_pretty)(FILE * file, const TEMPLATE(T, poly_t) poly, const char *x, const TEMPLATE(T, ctx_t) ctx);
int _TEMPLATE(T, poly_fprint)(FILE * file, const TEMPLATE(T, struct) *poly, slong len, const TEMPLATE(T, ctx_t) ctx);
int TEMPLATE(T, poly_fprint)(FILE * file, const TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx);
#endif

int _TEMPLATE(T, poly_print)(const TEMPLATE(T, struct) *poly, slong len, const TEMPLATE(T, ctx_t) ctx);
int TEMPLATE(T, poly_print)(const TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx);
int _TEMPLATE(T, poly_print_pretty)(const TEMPLATE(T, struct) *poly, slong len, const char *x, const TEMPLATE(T, ctx_t) ctx);
int TEMPLATE(T, poly_print_pretty)(const TEMPLATE(T, poly_t) poly, const char *x, const TEMPLATE(T, ctx_t) ctx);

char * _TEMPLATE(T, poly_get_str_pretty)(const TEMPLATE(T, struct) * poly, slong len,
                                  const char *x,
                                  const TEMPLATE(T, ctx_t) ctx);

char * TEMPLATE(T, poly_get_str_pretty)(const TEMPLATE(T, poly_t) poly,
                                 const char *x,
                                 const TEMPLATE(T, ctx_t) ctx);

char * _TEMPLATE(T, poly_get_str)(const TEMPLATE(T, struct) * poly, slong len,
                           const TEMPLATE(T, ctx_t) ctx);

char * TEMPLATE(T, poly_get_str)(const TEMPLATE(T, poly_t) poly,
                          const TEMPLATE(T, ctx_t) ctx);

/* Characteristic polynomial ************************************************/


void TEMPLATE(T, mat_charpoly_danilevsky) (TEMPLATE(T, poly_t) p,
                      const TEMPLATE(T, mat_t) A, const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, mat_charpoly)(TEMPLATE(T, poly_t) p,
                       const TEMPLATE(T, mat_t) M, const TEMPLATE(T, ctx_t) ctx);

/* Minimal polynomial ************************************************/


void TEMPLATE(T, mat_minpoly) (TEMPLATE(T, poly_t) p,
                      const TEMPLATE(T, mat_t) X, const TEMPLATE(T, ctx_t) ctx);

#ifdef __cplusplus
}
#endif

#endif
