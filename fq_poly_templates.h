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

    Copyright (C) 2012 Andres Goens
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#ifdef T

#ifdef __cplusplus
extern "C" {
#endif

/*  Type definitions *********************************************************/

typedef struct
{
    TEMPLATE(T, struct) *coeffs;
    slong alloc;
    slong length;
}
TEMPLATE(T, poly_struct);

typedef TEMPLATE(T, poly_struct) TEMPLATE(T, poly_t)[1];

/*  Memory management ********************************************************/

FLINT_DLL void TEMPLATE(T, poly_init)(TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_init2)(TEMPLATE(T, poly_t) poly, slong alloc,
                        const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_realloc)(TEMPLATE(T, poly_t) poly, slong alloc,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_truncate)(TEMPLATE(T, poly_t) poly, slong len,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_set_trunc)(TEMPLATE(T, poly_t) poly1, TEMPLATE(T, poly_t) poly2, slong len,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_fit_length)(TEMPLATE(T, poly_t) poly, slong len,
                             const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_clear)(TEMPLATE(T, poly_t) poly,
                        const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_normalise)(TEMPLATE(T, poly_t) poly,
                             const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_normalise2)(TEMPLATE(T, struct) *poly, slong *length,
                              const TEMPLATE(T, ctx_t) ctx);

FQ_POLY_TEMPLATES_INLINE void
_TEMPLATE(T, poly_set_length)(TEMPLATE(T, poly_t) poly, slong len,
                              const TEMPLATE(T, ctx_t) ctx)
{
    if (poly->length > len)
    {
        slong i;

        for (i = len; i < poly->length; i++)
            TEMPLATE(T, zero)(poly->coeffs + i, ctx);
    }
    poly->length = len;
}

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

FLINT_DLL void TEMPLATE(T, poly_randtest)(TEMPLATE(T, poly_t) f, flint_rand_t state,
                           slong len, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_randtest_not_zero)(TEMPLATE(T, poly_t) f, flint_rand_t state,
                                    slong len, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_randtest_monic) (TEMPLATE(T, poly_t) f, flint_rand_t state,
                                  slong len, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_randtest_irreducible) (TEMPLATE(T, poly_t) f,
                                        flint_rand_t state, slong len,
                                        const TEMPLATE(T, ctx_t) ctx);


/*  Assignment and basic manipulation  ***************************************/

FLINT_DLL void _TEMPLATE(T, poly_set)(TEMPLATE(T, struct) *rop, const TEMPLATE(T, struct) *op,
                       slong len, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_set)(TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op,
                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE3(T, poly_set, T)(TEMPLATE(T, poly_t) poly, const TEMPLATE(T, t) c,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_swap)(TEMPLATE(T, poly_t) op1, TEMPLATE(T, poly_t) op2,
                       const TEMPLATE(T, ctx_t) ctx);

FQ_POLY_TEMPLATES_INLINE void
_TEMPLATE(T, poly_zero)(TEMPLATE(T, struct) *rop, slong len,
                        const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        TEMPLATE(T, zero)(rop + i, ctx);
}

FQ_POLY_TEMPLATES_INLINE void
TEMPLATE(T, poly_zero)(TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx)
{
    _TEMPLATE(T, poly_set_length)(poly, 0, ctx);
}

FLINT_DLL void TEMPLATE(T, poly_one)(TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_gen)(TEMPLATE(T, poly_t) f, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_make_monic)(TEMPLATE(T, struct) *rop,
                              const TEMPLATE(T, struct) *op, slong length,
                              const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_make_monic)(TEMPLATE(T, poly_t) rop,
                                  const TEMPLATE(T, poly_t) op,
                                  const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_reverse)(TEMPLATE(T, struct) * res,
                           const TEMPLATE(T, struct) * poly, slong len, slong n,
                           const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_reverse)(TEMPLATE(T, poly_t) res,
                          const TEMPLATE(T, poly_t) poly, slong n,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL ulong TEMPLATE(T, poly_deflation)(const TEMPLATE(T, poly_t) input,
                            const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_deflate)(TEMPLATE(T, poly_t) result,
                          const TEMPLATE(T, poly_t) input, ulong deflation,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_inflate)(TEMPLATE(T, poly_t) result,
                          const TEMPLATE(T, poly_t) input, ulong inflation,
                          const TEMPLATE(T, ctx_t) ctx);

/*  Getting and setting coefficients  ****************************************/

FLINT_DLL void TEMPLATE(T, poly_get_coeff)(TEMPLATE(T, t) x, const TEMPLATE(T, poly_t) poly,
                            slong n, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_set_coeff)(TEMPLATE(T, poly_t) poly, slong n,
                                 const TEMPLATE(T, t) x,
                                 const TEMPLATE(T, ctx_t) ctx);

FQ_POLY_TEMPLATES_INLINE void 
TEMPLATE(T, poly_set_coeff_fmpz)(TEMPLATE(T, poly_t) poly, slong n,
                                 const fmpz_t x, const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_fit_length)(poly, n + 1, ctx);
    TEMPLATE(T, set_fmpz)(poly->coeffs + n, x, ctx);
    if (n + 1 > poly->length)
    {
        _TEMPLATE(T, poly_set_length)(poly, n + 1, ctx);
    }
    _TEMPLATE(T, poly_normalise)(poly, ctx);
}

FQ_POLY_TEMPLATES_INLINE int
TEMPLATE(T, poly_is_gen)(const TEMPLATE(T, poly_t) poly,
                         const TEMPLATE(T, ctx_t) ctx)
{
    return ((poly->length == 2) &&
            TEMPLATE(T, is_zero)(poly->coeffs, ctx) &&
            TEMPLATE(T, is_one)(poly->coeffs + 1, ctx));
}

/*  Comparison  **************************************************************/

FLINT_DLL int TEMPLATE(T, poly_equal)(const TEMPLATE(T, poly_t) poly1,
                            const TEMPLATE(T, poly_t) poly2,
                            const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, poly_equal_trunc)(const TEMPLATE(T, poly_t) poly1,
                            const TEMPLATE(T, poly_t) poly2,
                            slong n, const TEMPLATE(T, ctx_t) ctx);

FQ_POLY_TEMPLATES_INLINE int
TEMPLATE(T, poly_is_zero)(const TEMPLATE(T, poly_t) poly,
                          const TEMPLATE(T, ctx_t) ctx)
{
    return (poly->length == 0);
}

FQ_POLY_TEMPLATES_INLINE int
TEMPLATE(T, poly_is_one)(const TEMPLATE(T, poly_t) op,
                         const TEMPLATE(T, ctx_t) ctx)
{
    return (op->length == 1) && (TEMPLATE(T, is_one)(op->coeffs + 0, ctx));
}

FQ_POLY_TEMPLATES_INLINE int
TEMPLATE(T, poly_is_unit)(const TEMPLATE(T, poly_t) op,
                          const TEMPLATE(T, ctx_t) ctx)
{
    return (op->length == 1) && (!(TEMPLATE(T, is_zero)(op->coeffs + 0, ctx)));
}

FQ_POLY_TEMPLATES_INLINE int
TEMPLATE3(T, poly_equal, T)(const TEMPLATE(T, poly_t) poly,
                            const TEMPLATE(T, t) c,
                            const TEMPLATE(T, ctx_t) ctx)
{
    return ((poly->length == 0) && TEMPLATE(T, is_zero)(c, ctx)) ||
        ((poly->length == 1) && TEMPLATE(T, equal)(poly->coeffs, c, ctx));
}

/*  Addition and subtraction  ************************************************/

FLINT_DLL void _TEMPLATE(T, poly_add)(TEMPLATE(T, struct) *res,
                       const TEMPLATE(T, struct) *poly1, slong len1,
                       const TEMPLATE(T, struct) *poly2, slong len2,
                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_add)(TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op1,
                      const TEMPLATE(T, poly_t) op2,
                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_add_series)(TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op1,
                      const TEMPLATE(T, poly_t) op2,
                      slong n, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_sub)(TEMPLATE(T, struct) *res,
                       const TEMPLATE(T, struct) *poly1, slong len1,
                       const TEMPLATE(T, struct) *poly2, slong len2,
                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_sub)(TEMPLATE(T, poly_t) rop,
                      const TEMPLATE(T, poly_t) op1,
                      const TEMPLATE(T, poly_t) op2,
                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_sub_series)(TEMPLATE(T, poly_t) rop,
                      const TEMPLATE(T, poly_t) op1,
                      const TEMPLATE(T, poly_t) op2,
                      slong n, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_neg)(TEMPLATE(T, struct) *rop,
                       const TEMPLATE(T, struct) *op, slong len,
                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_neg)(TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op,
                      const TEMPLATE(T, ctx_t) ctx);

/*  Scalar multiplication and division  **************************************/

FLINT_DLL void _TEMPLATE3(T, poly_scalar_mul, T)(TEMPLATE(T, struct) *rop,
                                  const TEMPLATE(T, struct) *op, slong len,
                                  const TEMPLATE(T, t) x,
                                  const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE3(T, poly_scalar_mul, T)(TEMPLATE(T, poly_t) rop,
                                 const TEMPLATE(T, poly_t) op,
                                 const TEMPLATE(T, t) x,
                                 const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE3(T, poly_scalar_addmul, T)(TEMPLATE(T, struct) *rop,
                                     const TEMPLATE(T, struct) *op, slong len,
                                     const TEMPLATE(T, t) x,
                                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE3(T, poly_scalar_div, T)(TEMPLATE(T, struct) *rop,
                                  const TEMPLATE(T, struct) *op, slong len,
                                  const TEMPLATE(T, t) x,
                                  const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE3(T, poly_scalar_div, T)(TEMPLATE(T, poly_t) rop,
                                 const TEMPLATE(T, poly_t) op,
                                 const TEMPLATE(T, t) x,
                                 const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void  _TEMPLATE3(T, poly_scalar_addmul, T)(TEMPLATE(T, struct) *rop,
                                     const TEMPLATE(T, struct) *op, slong len,
                                     const TEMPLATE(T, t) x,
                                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE3(T, poly_scalar_addmul, T)(TEMPLATE(T, poly_t) rop,
                                    const TEMPLATE(T, poly_t) op,
                                    const TEMPLATE(T, t) x,
                                    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE3(T, poly_scalar_submul, T)(TEMPLATE(T, struct) *rop,
                                     const TEMPLATE(T, struct) *op, slong len,
                                     const TEMPLATE(T, t) x,
                                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE3(T, poly_scalar_submul, T)(TEMPLATE(T, poly_t) rop,
                                    const TEMPLATE(T, poly_t) op,
                                    const TEMPLATE(T, t) x,
                                    const TEMPLATE(T, ctx_t) ctx);

/*  Multiplication  **********************************************************/

FLINT_DLL void _TEMPLATE(T, poly_mul_classical)(TEMPLATE(T, struct) *rop,
                                 const TEMPLATE(T, struct) *op1, slong len1,
                                 const TEMPLATE(T, struct) *op2, slong len2,
                                 const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_mul_classical)(TEMPLATE(T, poly_t) rop,
                                const TEMPLATE(T, poly_t) op1,
                                const TEMPLATE(T, poly_t) op2,
                                const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_mul_reorder)(TEMPLATE(T, struct) *rop,
                               const TEMPLATE(T, struct) *op1, slong len1,
                               const TEMPLATE(T, struct) *op2, slong len2,
                               const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_mul_reorder)(TEMPLATE(T, poly_t) rop,
                              const TEMPLATE(T, poly_t) op1,
                              const TEMPLATE(T, poly_t) op2,
                              const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_mul_KS)(TEMPLATE(T, struct) *rop,
                          const TEMPLATE(T, struct) *op1, slong len1,
                          const TEMPLATE(T, struct) *op2, slong len2,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_mul_KS)(TEMPLATE(T, poly_t) rop,
                         const TEMPLATE(T, poly_t) op1,
                         const TEMPLATE(T, poly_t) op2,
                         const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_mul)(TEMPLATE(T, struct) *rop,
                       const TEMPLATE(T, struct) *op1, slong len1,
                       const TEMPLATE(T, struct) *op2, slong len2,
                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_mul)(TEMPLATE(T, poly_t) rop,
                      const TEMPLATE(T, poly_t) op1,
                      const TEMPLATE(T, poly_t) op2,
                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_mullow_classical)(TEMPLATE(T, struct) *rop,
                                    const TEMPLATE(T, struct) *op1, slong len1,
                                    const TEMPLATE(T, struct) *op2, slong len2,
                                    slong n,
                                    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_mullow_classical)(TEMPLATE(T, poly_t) rop,
                                   const TEMPLATE(T, poly_t) op1,
                                   const TEMPLATE(T, poly_t) op2,
                                   slong n, const
                                   TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_mullow_KS)(TEMPLATE(T, struct) *rop,
                             const TEMPLATE(T, struct) *op1, slong len1,
                             const TEMPLATE(T, struct) *op2, slong len2,
                             slong n,
                             const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_mullow_KS)(TEMPLATE(T, poly_t) rop,
                            const TEMPLATE(T, poly_t) op1,
                            const TEMPLATE(T, poly_t) op2,
                            slong n,
                            const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_mullow)(TEMPLATE(T, struct) *rop,
                          const TEMPLATE(T, struct) *op1, slong len1,
                          const TEMPLATE(T, struct) *op2, slong len2,
                          slong n,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_mullow)(TEMPLATE(T, poly_t) rop,
                         const TEMPLATE(T, poly_t) op1,
                         const TEMPLATE(T, poly_t) op2,
                         slong n,
                         const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_mulhigh_classical)(
    TEMPLATE(T, struct)* rop,
    const TEMPLATE(T, struct)* op1, slong len1,
    const TEMPLATE(T, struct)* op2, slong len2,
    slong start,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_mulhigh_classical)(TEMPLATE(T, poly_t) rop,
                                    const TEMPLATE(T, poly_t) op1,
                                    const TEMPLATE(T, poly_t) op2,
                                    slong start,
                                    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_mulhigh)(TEMPLATE(T, struct)* res,
                           const TEMPLATE(T, struct)* poly1, slong len1,
                           const TEMPLATE(T, struct)* poly2, slong len2,
                           slong n, TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_mulhigh)(TEMPLATE(T, poly_t) rop,
                          const TEMPLATE(T, poly_t) op1,
                          const TEMPLATE(T, poly_t) op2, slong start,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_mulmod)(TEMPLATE(T, struct) * res,
                          const TEMPLATE(T, struct) * poly1, slong len1,
                          const TEMPLATE(T, struct) * poly2, slong len2,
                          const TEMPLATE(T, struct) * f, slong lenf,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_mulmod)(TEMPLATE(T, poly_t) res,
                         const TEMPLATE(T, poly_t) poly1,
                         const TEMPLATE(T, poly_t) poly2,
                         const TEMPLATE(T, poly_t) f,
                         const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_mulmod_preinv)(TEMPLATE(T, struct) *res,
                                 const TEMPLATE(T, struct) *poly1, slong len1,
                                 const TEMPLATE(T, struct)  *poly2, slong len2,
                                 const TEMPLATE(T, struct) *f, slong lenf,
                                 const TEMPLATE(T, struct) *finv, slong lenfinv,
                                 const TEMPLATE(T, ctx_t) ctx);
FLINT_DLL void TEMPLATE(T, poly_mulmod_preinv)(TEMPLATE(T, poly_t) res,
                                const TEMPLATE(T, poly_t) poly1,
                                const TEMPLATE(T, poly_t) poly2,
                                const TEMPLATE(T, poly_t) f,
                                const TEMPLATE(T, poly_t) finv,
                                const TEMPLATE(T, ctx_t) ctx);

/* Squaring ******************************************************************/

FLINT_DLL void _TEMPLATE(T, poly_sqr_classical)(TEMPLATE(T, struct) *rop,
                                 const TEMPLATE(T, struct) *op, slong len,
                                 const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_sqr_classical)(TEMPLATE(T, poly_t) rop,
                                const TEMPLATE(T, poly_t) op,
                                const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_sqr_reorder)(TEMPLATE(T, struct) *rop,
                               const TEMPLATE(T, struct) *op, slong len,
                               const TEMPLATE(T, ctx_t) ctx);
                              
FLINT_DLL void TEMPLATE(T, poly_sqr_reorder)(TEMPLATE(T, poly_t) rop,
                              const TEMPLATE(T, poly_t) op,
                              const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_sqr_KS)(TEMPLATE(T, struct) *rop,
                          const TEMPLATE(T, struct) *op, slong len,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_sqr_KS)(TEMPLATE(T, poly_t) rop,
                         const TEMPLATE(T, poly_t) op,
                         const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_sqr)(TEMPLATE(T, struct) *rop,
                       const TEMPLATE(T, struct) *op, slong len,
                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_sqr)(TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op,
                      const TEMPLATE(T, ctx_t) ctx);

/*  Powering  ****************************************************************/

FLINT_DLL void _TEMPLATE(T, poly_pow)(TEMPLATE(T, struct) *rop,
                       const TEMPLATE(T, struct) *op, slong len,
                       ulong e,
                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_pow)(TEMPLATE(T, poly_t) rop,
                      const TEMPLATE(T, poly_t) op,
                      ulong e,
                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_powmod_fmpz_binexp)(TEMPLATE(T, struct) * res,
                                      const TEMPLATE(T, struct) * poly,
                                      const fmpz_t e,
                                      const TEMPLATE(T, struct) * f, slong lenf,
                                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_powmod_fmpz_binexp)(TEMPLATE(T, poly_t) res,
                                     const TEMPLATE(T, poly_t) poly,
                                     const fmpz_t e,
                                     const TEMPLATE(T, poly_t) f,
                                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_powmod_fmpz_binexp_preinv)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly,
    const fmpz_t e,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) *finv, slong lenfinv,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_powmod_fmpz_binexp_preinv)(
    TEMPLATE(T, poly_t) res,
    const TEMPLATE(T, poly_t) poly, const fmpz_t e,
    const TEMPLATE(T, poly_t) f, const TEMPLATE(T, poly_t) finv,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_powmod_ui_binexp)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly,
    ulong e,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_powmod_ui_binexp)(
    TEMPLATE(T, poly_t) res,
    const TEMPLATE(T, poly_t) poly, ulong e,
    const TEMPLATE(T, poly_t) f,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_powmod_ui_binexp_preinv)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly,
    ulong e,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * finv, slong lenfinv,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_powmod_ui_binexp_preinv)(
    TEMPLATE(T, poly_t) res,
    const TEMPLATE(T, poly_t) poly,
    ulong e,
    const TEMPLATE(T, poly_t) f,
    const TEMPLATE(T, poly_t) finv,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly,
    const fmpz_t e, ulong k,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * finv, slong lenfinv,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(
    TEMPLATE(T, poly_t) res,
    const TEMPLATE(T, poly_t) poly,
    const fmpz_t e, ulong k,
    const TEMPLATE(T, poly_t) f,
    const TEMPLATE(T, poly_t) finv,
    const TEMPLATE(T, ctx_t) ctx);


FLINT_DLL void _TEMPLATE(T, poly_powmod_x_fmpz_preinv)(
    TEMPLATE(T, struct) * res, const fmpz_t e,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct)* finv, slong lenfinv,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_powmod_x_fmpz_preinv)(
    TEMPLATE(T, poly_t) res, const fmpz_t e,
    const TEMPLATE(T, poly_t) f,
    const TEMPLATE(T, poly_t) finv,
    const TEMPLATE(T, ctx_t) ctx);

/*  Shifting  ****************************************************************/

FLINT_DLL void _TEMPLATE(T, poly_shift_left)(TEMPLATE(T, struct) *rop,
                              const TEMPLATE(T, struct) *op, slong len,
                              slong n,
                              const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_shift_left)(TEMPLATE(T, poly_t) rop,
                             const TEMPLATE(T, poly_t) op,
                             slong n,
                             const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_shift_right)(TEMPLATE(T, struct) *rop,
                               const TEMPLATE(T, struct) *op, slong len,
                               slong n,
                               const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_shift_right)(TEMPLATE(T, poly_t) rop,
                              const TEMPLATE(T, poly_t) op,
                              slong n,
                              const TEMPLATE(T, ctx_t) ctx);

/*  Norms  *******************************************************************/

FLINT_DLL slong _TEMPLATE(T, poly_hamming_weight)(const TEMPLATE(T, struct) *op, slong len,
                                  const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL slong TEMPLATE(T, poly_hamming_weight)(const TEMPLATE(T, poly_t) op,
                                 const TEMPLATE(T, ctx_t) ctx);

/*  Greatest common divisor  *************************************************/

FLINT_DLL void TEMPLATE(T, poly_gcd_euclidean)(TEMPLATE(T, poly_t) rop,
                                const TEMPLATE(T, poly_t) op1,
                                const TEMPLATE(T, poly_t) op2,
                                const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL slong _TEMPLATE(T, poly_gcd_euclidean)(TEMPLATE(T, struct)* G,
                                 const TEMPLATE(T, struct)* A, slong lenA,
                                 const TEMPLATE(T, struct)* B, slong lenB,
                                 const TEMPLATE(T, t) invB,
                                 const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL slong _TEMPLATE(T, poly_gcd)(TEMPLATE(T, struct)* G,
                       const TEMPLATE(T, struct)* A, slong lenA,
                       const TEMPLATE(T, struct)* B, slong lenB,
                       const TEMPLATE(T, t) invB,
                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_gcd)(TEMPLATE(T, poly_t) rop,
                      const TEMPLATE(T, poly_t) op1,
                      const TEMPLATE(T, poly_t) op2,
                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL slong _TEMPLATE(T, poly_gcd_euclidean_f)(TEMPLATE(T, t) f, TEMPLATE(T, struct)* G,
                                   const TEMPLATE(T, struct)* A, slong lenA,
                                   const TEMPLATE(T, struct)* B, slong lenB,
                                   const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_gcd_euclidean_f)(TEMPLATE(T, t) f, TEMPLATE(T, poly_t) G,
                                  const TEMPLATE(T, poly_t) A,
                                  const TEMPLATE(T, poly_t) B,
                                  const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL slong _TEMPLATE(T, poly_hgcd)(TEMPLATE(T, struct) **M, slong *lenM,
                        TEMPLATE(T, struct)* A, slong *lenA,
                        TEMPLATE(T, struct)* B, slong *lenB,
                        const TEMPLATE(T, struct)* a, slong lena,
                        const TEMPLATE(T, struct)* b, slong lenb,
                        const TEMPLATE(T, ctx_t) ctx);


FLINT_DLL slong _TEMPLATE(T, poly_gcd_hgcd)(TEMPLATE(T, struct)* G,
                            const TEMPLATE(T, struct)* A, slong lenA,
                            const TEMPLATE(T, struct)* B, slong lenB,
                            const TEMPLATE(T, t) invB,
                            const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_gcd_hgcd)(TEMPLATE(T, poly_t) G,
                                const TEMPLATE(T, poly_t) A,
                                const TEMPLATE(T, poly_t) B,
                                const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL slong _TEMPLATE(T, poly_xgcd_euclidean_f)(TEMPLATE(T, t) f, TEMPLATE(T, struct) *G,
                                    TEMPLATE(T, struct) *S,
                                    TEMPLATE(T, struct) *T,
                                    const TEMPLATE(T, struct) *A, slong lenA,
                                    const TEMPLATE(T, struct) *B, slong lenB,
                                    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_xgcd_euclidean_f)(TEMPLATE(T, t) f, TEMPLATE(T, poly_t) G,
                                   TEMPLATE(T, poly_t) S,
                                   TEMPLATE(T, poly_t) T,
                                   const TEMPLATE(T, poly_t) A,
                                   const TEMPLATE(T, poly_t) B,
                                   const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL slong _TEMPLATE(T, poly_xgcd_euclidean)(TEMPLATE(T, struct) *G,
                                  TEMPLATE(T, struct) *S,
                                  TEMPLATE(T, struct) *T,
                                  const TEMPLATE(T, struct) *A, slong lenA,
                                  const TEMPLATE(T, struct) *B, slong lenB,
                                  const TEMPLATE(T, t) invB,
                                  const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_xgcd_euclidean)(TEMPLATE(T, poly_t) G,
                                 TEMPLATE(T, poly_t) S, TEMPLATE(T, poly_t) T,
                                 const TEMPLATE(T, poly_t) A,
                                 const TEMPLATE(T, poly_t) B,
                                 const TEMPLATE(T, ctx_t) ctx);

FQ_POLY_TEMPLATES_INLINE slong
_TEMPLATE(T, poly_xgcd)(TEMPLATE(T, struct) *G,
                        TEMPLATE(T, struct) *S, TEMPLATE(T, struct) *T,
                        const TEMPLATE(T, struct) *A, slong lenA,
                        const TEMPLATE(T, struct) *B, slong lenB,
                        const TEMPLATE(T, t) invB,
                        const TEMPLATE(T, ctx_t) ctx)
{
    return _TEMPLATE(T, poly_xgcd_euclidean)(G, S, T, A, lenA,
                                             B, lenB, invB, ctx);
}


FQ_POLY_TEMPLATES_INLINE void
TEMPLATE(T, poly_xgcd)(TEMPLATE(T, poly_t) G,
                       TEMPLATE(T, poly_t) S, TEMPLATE(T, poly_t) T,
                       const TEMPLATE(T, poly_t) A,
                       const TEMPLATE(T, poly_t) B,
                       const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_xgcd_euclidean)(G, S, T, A, B, ctx);
}


/*  Euclidean division  ******************************************************/

FLINT_DLL ulong TEMPLATE(T, poly_remove)(TEMPLATE(T, poly_t) f,
                         const TEMPLATE(T, poly_t) g,
                         const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_div_basecase)(TEMPLATE(T, struct) *Q, TEMPLATE(T, struct) *R,
                               const TEMPLATE(T, struct) *A, slong lenA,
                               const TEMPLATE(T, struct) *B, slong lenB,
                               const TEMPLATE(T, t) invB,
                               const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_div_basecase)(TEMPLATE(T, poly_t) Q,
                               const TEMPLATE(T, poly_t) A,
                               const TEMPLATE(T, poly_t) B,
                               const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_divrem_basecase)(TEMPLATE(T, struct) *Q,
                                   TEMPLATE(T, struct) *R,
                                   const TEMPLATE(T, struct) *A, slong lenA,
                                   const TEMPLATE(T, struct) *B, slong lenB,
                                   const TEMPLATE(T, t) invB,
                                   const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_divrem_basecase)(TEMPLATE(T, poly_t) Q, TEMPLATE(T, poly_t) R,
                                  const TEMPLATE(T, poly_t) A,
                                  const TEMPLATE(T, poly_t) B,
                                  const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_divrem_divconquer_recursive)(
    TEMPLATE(T, struct) * Q, TEMPLATE(T, struct) * BQ,
    TEMPLATE(T, struct) * W,
    const TEMPLATE(T, struct) * A,
    const TEMPLATE(T, struct) * B, slong lenB,
    const TEMPLATE(T, t) invB,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_divrem_divconquer)(
    TEMPLATE(T, struct) *Q, TEMPLATE(T, struct) *R,
    const TEMPLATE(T, struct) *A, slong lenA,
    const TEMPLATE(T, struct) *B, slong lenB,
    const TEMPLATE(T, t) invB,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_divrem_divconquer)(TEMPLATE(T, poly_t) Q,
                                    TEMPLATE(T, poly_t) R,
                                    const TEMPLATE(T, poly_t) A,
                                    const TEMPLATE(T, poly_t) B,
                                    const TEMPLATE(T, ctx_t) ctx);

FQ_POLY_TEMPLATES_INLINE void
_TEMPLATE(T, poly_divrem)(TEMPLATE(T, struct) *Q, TEMPLATE(T, struct) *R,
                          const TEMPLATE(T, struct) *A, slong lenA,
                          const TEMPLATE(T, struct) *B, slong lenB,
                          const TEMPLATE(T, t) invB,
                          const TEMPLATE(T, ctx_t) ctx)
{
    _TEMPLATE(T, poly_divrem_divconquer)(Q, R, A, lenA, B, lenB, invB, ctx);
}

FQ_POLY_TEMPLATES_INLINE void
TEMPLATE(T, poly_divrem)(TEMPLATE(T, poly_t) Q, TEMPLATE(T, poly_t) R,
                         const TEMPLATE(T, poly_t) A,
                         const TEMPLATE(T, poly_t) B,
                         const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_divrem_divconquer)(Q, R, A, B, ctx);
}

FQ_POLY_TEMPLATES_INLINE void
_TEMPLATE(T, poly_rem)(TEMPLATE(T, struct) *R,
                       const TEMPLATE(T, struct) *A, slong lenA,
                       const TEMPLATE(T, struct) *B, slong lenB,
                       const TEMPLATE(T, t) invB,
                       const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, struct) *Q = _TEMPLATE(T, vec_init)(lenA - lenB + 1, ctx);

    if (lenA < lenB)
    {
        _TEMPLATE(T, vec_set)(R, A, lenA, ctx);
        _TEMPLATE(T, vec_zero)(R + lenA, lenB - 1 - lenA, ctx);
    }
    else
    {
        TEMPLATE(T, struct) *T = _TEMPLATE(T, vec_init)(lenA, ctx);
       _TEMPLATE(T, poly_divrem_divconquer)(Q, T, A, lenA, B, lenB, invB, ctx);
       _TEMPLATE(T, vec_set)(R, T, lenB - 1, ctx);
       _TEMPLATE(T, vec_clear)(T, lenA, ctx);
    }

    _TEMPLATE(T, vec_clear)(Q, lenA - lenB + 1, ctx);
}



FQ_POLY_TEMPLATES_INLINE void
TEMPLATE(T, poly_rem)(TEMPLATE(T, poly_t) R,
                      const TEMPLATE(T, poly_t) A, const TEMPLATE(T, poly_t) B,
                      const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_t) Q;
    TEMPLATE(T, poly_init)(Q, ctx);
    TEMPLATE(T, poly_divrem)(Q, R, A, B, ctx);
    TEMPLATE(T, poly_clear)(Q, ctx);
}

FLINT_DLL void _TEMPLATE(T, poly_inv_series_newton)(TEMPLATE(T, struct) * Qinv,
                                     const TEMPLATE(T, struct) * Q, slong n,
                                     const TEMPLATE(T, t) cinv,
                                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_inv_series_newton)(TEMPLATE(T, poly_t) Qinv,
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

FLINT_DLL void _TEMPLATE(T, poly_div_series) (TEMPLATE(T, struct) * Q, 
                          const TEMPLATE(T, struct) * A, slong Alen,
                          const TEMPLATE(T, struct) * B, slong Blen, 
                          slong n, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_div_series)(TEMPLATE(T, poly_t) Q, 
               const TEMPLATE(T, poly_t) A, const TEMPLATE(T, poly_t) B, 
                                        slong n, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_div_newton_n_preinv) (
    TEMPLATE(T, struct) *Q,
    const TEMPLATE(T, struct) *A, slong lenA,
    const TEMPLATE(T, struct)* B, slong lenB,
    const TEMPLATE(T, struct)* Binv, slong lenBinv,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_div_newton_n_preinv) (TEMPLATE(T, poly_t) Q,
                                       const TEMPLATE(T, poly_t) A,
                                       const TEMPLATE(T, poly_t) B,
                                       const TEMPLATE(T, poly_t) Binv,
                                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_divrem_newton_n_preinv) (
    TEMPLATE(T, struct)* Q, TEMPLATE(T, struct)* R,
    const TEMPLATE(T, struct)* A, slong lenA,
    const TEMPLATE(T, struct)* B, slong lenB,
    const TEMPLATE(T, struct)* Binv, slong lenBinv,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_divrem_newton_n_preinv)(
    TEMPLATE(T, poly_t) Q, TEMPLATE(T, poly_t) R,
    const TEMPLATE(T, poly_t) A, const TEMPLATE(T, poly_t) B,
    const TEMPLATE(T, poly_t) Binv, const TEMPLATE(T, ctx_t) ctx);


FLINT_DLL void _TEMPLATE(T, poly_divrem_f)(TEMPLATE(T, t) f,
                            TEMPLATE(T, struct)* Q, TEMPLATE(T, struct)* R,
                            const TEMPLATE(T, struct)* A, slong lenA,
                            const TEMPLATE(T, struct)* B, slong lenB,
                            const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_divrem_f)(TEMPLATE(T, t) f,
                           TEMPLATE(T, poly_t) Q, TEMPLATE(T, poly_t) R,
                           const TEMPLATE(T, poly_t) A,
                           const TEMPLATE(T, poly_t) B,
                           const TEMPLATE(T, ctx_t) ctx);

/*  Divisibility testing  ***************************************************/

FLINT_DLL int _TEMPLATE(T, poly_divides)(TEMPLATE(T, struct) *Q,
                           const TEMPLATE(T, struct) *A, slong lenA,
                           const TEMPLATE(T, struct) *B, slong lenB,
                           const TEMPLATE(T, t) invB,
                           const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, poly_divides)(TEMPLATE(T, poly_t) Q,
                          const TEMPLATE(T, poly_t) A, const TEMPLATE(T, poly_t) B,
                          const TEMPLATE(T, ctx_t) ctx);

/*  Derivative  **************************************************************/

FLINT_DLL void _TEMPLATE(T, poly_derivative)(TEMPLATE(T, struct) *rop,
                              const TEMPLATE(T, struct) *op, slong len,
                              const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_derivative)(TEMPLATE(T, poly_t) rop,
                             const TEMPLATE(T, poly_t) op,
                             const TEMPLATE(T, ctx_t) ctx);

/*  Evaluation  **************************************************************/

FLINT_DLL void _TEMPLATE3(T, poly_evaluate, T)(TEMPLATE(T, t) rop,
                                     const TEMPLATE(T, struct) *op, slong len,
                                     const TEMPLATE(T, t) a,
                                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE3(T, poly_evaluate, T)(TEMPLATE(T, t) res,
                                    const TEMPLATE(T, poly_t) f,
                                    const TEMPLATE(T, t) a,
                                    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE4(T, poly_evaluate, T, vec)(TEMPLATE(T, struct) * ys,
                                     const TEMPLATE(T, struct) *coeffs, slong len,
                                     const TEMPLATE(T, struct) *xs, slong n,
                                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE4(T, poly_evaluate, T, vec)(TEMPLATE(T, struct) * ys,
                                    const TEMPLATE(T, poly_t) poly,
                                    const TEMPLATE(T, struct) * xs, slong n,
                                    const TEMPLATE(T, ctx_t) ctx);
                                    
FLINT_DLL TEMPLATE(T, poly_struct) **
_TEMPLATE(T, poly_tree_alloc)(slong len, const TEMPLATE(T, ctx_t) ctx);
 
FLINT_DLL void _TEMPLATE(T, poly_tree_free)(TEMPLATE(T, poly_struct) ** tree, slong len,
                             const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_tree_build)(TEMPLATE(T, poly_struct) ** tree,
                              const TEMPLATE(T, struct) * roots,
                              slong len,
                              const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE4(T, poly_evaluate, T, vec_fast_precomp)
    (TEMPLATE(T, struct) * vs,
     const TEMPLATE(T, struct) * poly, slong plen,
     TEMPLATE(T, poly_struct) * const * tree, slong len,
     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE4(T, poly_evaluate, T, vec_fast)(TEMPLATE(T, struct) * ys,
                                          const TEMPLATE(T, struct) * poly, slong plen,
                                          const TEMPLATE(T, struct) * xs, slong n,
                                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE4(T, poly_evaluate, T, vec_fast)(TEMPLATE(T, struct) * ys,
                                         const TEMPLATE(T, poly_t) poly,
                                         const TEMPLATE(T, struct) *xs, slong n,
                                         const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE4(T, poly_evaluate, T, vec_iter)(TEMPLATE(T, struct) * ys,
                                          const TEMPLATE(T, struct) * coeffs, slong len,
                                          const TEMPLATE(T, struct) * xs, slong n,
                                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE4(T, poly_evaluate, T, vec_iter)(TEMPLATE(T, struct) * ys,
                                         const TEMPLATE(T, poly_t) poly,
                                         const TEMPLATE(T, struct) * xs, slong n,
                                         const TEMPLATE(T, ctx_t) ctx);

/*  Composition  *************************************************************/

FLINT_DLL void _TEMPLATE(T, poly_compose_divconquer)(
    TEMPLATE(T, struct) *rop,
    const TEMPLATE(T, struct) *op1, slong len1,
    const TEMPLATE(T, struct) *op2, slong len2,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_compose_divconquer)(TEMPLATE(T, poly_t) rop,
                                     const TEMPLATE(T, poly_t) op1,
                                     const TEMPLATE(T, poly_t) op2,
                                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_compose_horner)(TEMPLATE(T, struct) *rop,
                                  const TEMPLATE(T, struct) *op1, slong len1,
                                  const TEMPLATE(T, struct) *op2, slong len2,
                                  const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_compose_horner)(TEMPLATE(T, poly_t) rop,
                                 const TEMPLATE(T, poly_t) op1,
                                 const TEMPLATE(T, poly_t) op2,
                                 const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_compose)(TEMPLATE(T, struct) *rop,
                           const TEMPLATE(T, struct) *op1, slong len1,
                           const TEMPLATE(T, struct) *op2, slong len2,
                           const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_compose)(TEMPLATE(T, poly_t) rop,
                          const TEMPLATE(T, poly_t) op1,
                          const TEMPLATE(T, poly_t) op2,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_compose_mod)(TEMPLATE(T, struct) * res,
                     const TEMPLATE(T, struct) * f, slong lenf,
                     const TEMPLATE(T, struct) * g,
                     const TEMPLATE(T, struct) * h, slong lenh,
                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_compose_mod)(TEMPLATE(T, poly_t) res,
                              const TEMPLATE(T, poly_t) poly1,
                              const TEMPLATE(T, poly_t) poly2,
                              const TEMPLATE(T, poly_t) poly3,
                              const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_compose_mod_preinv)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * g,
    const TEMPLATE(T, struct) * h, slong lenh,
    const TEMPLATE(T, struct) * hinv, slong lenhinv,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_compose_mod_preinv)(TEMPLATE(T, poly_t) res,
                                     const TEMPLATE(T, poly_t) poly1,
                                     const TEMPLATE(T, poly_t) poly2,
                                     const TEMPLATE(T, poly_t) poly3,
                                     const TEMPLATE(T, poly_t) poly3inv,
                                     const TEMPLATE(T, ctx_t) ctx);


FLINT_DLL void _TEMPLATE(T, poly_compose_mod_horner)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * g,
    const TEMPLATE(T, struct) * h, slong lenh,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_compose_mod_horner)(TEMPLATE(T, poly_t) res,
                                     const TEMPLATE(T, poly_t) poly1,
                                     const TEMPLATE(T, poly_t) poly2,
                                     const TEMPLATE(T, poly_t) poly3,
                                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_compose_mod_horner_preinv)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * f, slong lenf,
    const TEMPLATE(T, struct) * g,
    const TEMPLATE(T, struct) * h, slong lenh,
    const TEMPLATE(T, struct) * hinv, slong lenhinv,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_compose_mod_horner_preinv)(TEMPLATE(T, poly_t) res,
                                            const TEMPLATE(T, poly_t) poly1,
                                            const TEMPLATE(T, poly_t) poly2,
                                            const TEMPLATE(T, poly_t) poly3,
                                            const TEMPLATE(T, poly_t) poly3inv,
                                            const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_compose_mod_brent_kung)(TEMPLATE(T, poly_t) res,
                                         const TEMPLATE(T, poly_t) poly1,
                                         const TEMPLATE(T, poly_t) poly2,
                                         const TEMPLATE(T, poly_t) poly3,
                                         const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_compose_mod_brent_kung)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly1, slong len1,
    const TEMPLATE(T, struct) * poly2,
    const TEMPLATE(T, struct) * poly3, slong len3,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_compose_mod_brent_kung_preinv)(
    TEMPLATE(T, struct) * res,
    const TEMPLATE(T, struct) * poly1, slong len1,
    const TEMPLATE(T, struct) * poly2,
    const TEMPLATE(T, struct) * poly3, slong len3,
    const TEMPLATE(T, struct) * poly3inv, slong len3inv,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_compose_mod_brent_kung_preinv)(
    TEMPLATE(T, poly_t) res,
    const TEMPLATE(T, poly_t) poly1,
    const TEMPLATE(T, poly_t) poly2,
    const TEMPLATE(T, poly_t) poly3,
    const TEMPLATE(T, poly_t) poly3inv,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_reduce_matrix_mod_poly) (TEMPLATE(T, mat_t) A,
                                           const TEMPLATE(T, mat_t) B,
                                           const TEMPLATE(T, poly_t) f,
                                           const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_precompute_matrix) (
    TEMPLATE(T, mat_t A),
    const TEMPLATE(T, struct)* poly1,
    const TEMPLATE(T, struct)* poly2, slong len2,
    const TEMPLATE(T, struct)* poly2inv, slong len2inv,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_precompute_matrix) (TEMPLATE(T, mat_t A),
                                     const TEMPLATE(T, poly_t) poly1,
                                     const TEMPLATE(T, poly_t) poly2,
                                     const TEMPLATE(T, poly_t) poly2inv,
                                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv)(
    TEMPLATE(T, struct)* res,
    const TEMPLATE(T, struct)* poly1, slong len1,
    const TEMPLATE(T, mat_t) A,
    const TEMPLATE(T, struct)* poly3, slong len3,
    const TEMPLATE(T, struct)* poly3inv, slong len3inv,
    const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv)(
    TEMPLATE(T, poly_t) res,
    const TEMPLATE(T, poly_t) poly1,
    const TEMPLATE(T, mat_t) A,
    const TEMPLATE(T, poly_t) poly3, 
    const TEMPLATE(T, poly_t) poly3inv,
    const TEMPLATE(T, ctx_t) ctx);

/*  Input and output  ********************************************************/

FLINT_DLL int _TEMPLATE(T, poly_fprint_pretty)(FILE *file,
                                 const TEMPLATE(T, struct) *poly, slong len,
                                 const char *x,
                                 const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, poly_fprint_pretty)(FILE * file,
                                const TEMPLATE(T, poly_t) poly,
                                const char *x,
                                const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int _TEMPLATE(T, poly_fprint)(FILE * file,
                          const TEMPLATE(T, struct) *poly, slong len,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL int TEMPLATE(T, poly_fprint)(FILE * file,
                         const TEMPLATE(T, poly_t) poly,
                         const TEMPLATE(T, ctx_t) ctx);

FQ_POLY_TEMPLATES_INLINE int
_TEMPLATE(T, poly_print)(const TEMPLATE(T, struct) *poly, slong len,
                         const TEMPLATE(T, ctx_t) ctx)
{
    return _TEMPLATE(T, poly_fprint)(stdout, poly, len, ctx);
}

FQ_POLY_TEMPLATES_INLINE int
TEMPLATE(T, poly_print)(const TEMPLATE(T, poly_t) poly,
                        const TEMPLATE(T, ctx_t) ctx)
{
    return TEMPLATE(T, poly_fprint)(stdout, poly, ctx);
}


FQ_POLY_TEMPLATES_INLINE int
_TEMPLATE(T, poly_print_pretty)(const TEMPLATE(T, struct) *poly, slong len,
                                const char *x, const TEMPLATE(T, ctx_t) ctx)
{
    return _TEMPLATE(T, poly_fprint_pretty)(stdout, poly, len, x, ctx);
}

FQ_POLY_TEMPLATES_INLINE int
TEMPLATE(T, poly_print_pretty)(const TEMPLATE(T, poly_t) poly,
                               const char *x,
                               const TEMPLATE(T, ctx_t) ctx)
{
    return TEMPLATE(T, poly_fprint_pretty)(stdout, poly, x, ctx);
}

FLINT_DLL char * _TEMPLATE(T, poly_get_str_pretty)(const TEMPLATE(T, struct) * poly, slong len,
                                  const char *x,
                                  const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL char * TEMPLATE(T, poly_get_str_pretty)(const TEMPLATE(T, poly_t) poly,
                                 const char *x,
                                 const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL char * _TEMPLATE(T, poly_get_str)(const TEMPLATE(T, struct) * poly, slong len,
                           const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL char * TEMPLATE(T, poly_get_str)(const TEMPLATE(T, poly_t) poly,
                          const TEMPLATE(T, ctx_t) ctx);

#ifdef __cplusplus
}
#endif

#endif
