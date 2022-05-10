/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#ifdef __cplusplus
extern "C" {
#endif

/* memory management **********************************************************/

FLINT_DLL void TEMPLATE(T, poly_init)(TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_init2)(TEMPLATE(T, poly_t) poly, slong alloc,
                        const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_clear)(TEMPLATE(T, poly_t) poly,
                        const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_fit_length)(TEMPLATE(T, poly_t) poly, slong len,
                             const TEMPLATE(T, ctx_t) ctx);

FQ_POLY_MINI_TEMPLATES_INLINE void
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

FLINT_DLL void _TEMPLATE(T, poly_normalise)(TEMPLATE(T, poly_t) poly,
                             const TEMPLATE(T, ctx_t) ctx);

/* assignments ****************************************************************/

FLINT_DLL void TEMPLATE(T, poly_set_coeff)(TEMPLATE(T, poly_t) poly, slong n,
                                 const TEMPLATE(T, t) x,
                                 const TEMPLATE(T, ctx_t) ctx);

FQ_POLY_MINI_TEMPLATES_INLINE void
_TEMPLATE(T, poly_zero)(TEMPLATE(T, struct) *rop, slong len,
                        const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    for (i = 0; i < len; i++)
        TEMPLATE(T, zero)(rop + i, ctx);
}

FQ_POLY_MINI_TEMPLATES_INLINE void
TEMPLATE(T, poly_zero)(TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx)
{
    _TEMPLATE(T, poly_set_length)(poly, 0, ctx);
}

FLINT_DLL void TEMPLATE(T, poly_one)(TEMPLATE(T, poly_t) poly, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_set)(TEMPLATE(T, struct) *rop, const TEMPLATE(T, struct) *op,
                       slong len, const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_set)(TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op,
                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE3(T, poly_set, T)(TEMPLATE(T, poly_t) poly, const TEMPLATE(T, t) c,
                          const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_swap)(TEMPLATE(T, poly_t) op1, TEMPLATE(T, poly_t) op2,
                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_make_monic)(TEMPLATE(T, struct) *rop,
                              const TEMPLATE(T, struct) *op, slong length,
                              const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_make_monic)(TEMPLATE(T, poly_t) rop,
                                  const TEMPLATE(T, poly_t) op,
                                  const TEMPLATE(T, ctx_t) ctx);

/* parameters *****************************************************************/

FQ_POLY_MINI_TEMPLATES_INLINE slong
TEMPLATE(T, poly_length)(const TEMPLATE(T, poly_t) poly,
                         const TEMPLATE(T, ctx_t) ctx)
{
    return poly->length;
}

FQ_POLY_MINI_TEMPLATES_INLINE slong
TEMPLATE(T, poly_degree)(const TEMPLATE(T, poly_t) poly,
                         const TEMPLATE(T, ctx_t) ctx)
{
    return poly->length - 1;
}

FQ_POLY_MINI_TEMPLATES_INLINE TEMPLATE(T, struct) *
TEMPLATE(T, poly_lead)(const TEMPLATE(T, poly_t) poly,
                       const TEMPLATE(T, ctx_t) ctx)
{
    return poly->length > 0 ? poly->coeffs + (poly->length - 1) : NULL;
}

/* comparisons ****************************************************************/

FQ_POLY_MINI_TEMPLATES_INLINE int
TEMPLATE(T, poly_is_zero)(const TEMPLATE(T, poly_t) poly,
                          const TEMPLATE(T, ctx_t) ctx)
{
    return (poly->length == 0);
}


/* arithmetic operations ******************************************************/

FLINT_DLL void _TEMPLATE(T, poly_neg)(TEMPLATE(T, struct) *rop,
                       const TEMPLATE(T, struct) *op, slong len,
                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_neg)(TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op,
                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_add)(TEMPLATE(T, struct) *res,
                       const TEMPLATE(T, struct) *poly1, slong len1,
                       const TEMPLATE(T, struct) *poly2, slong len2,
                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_add)(TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op1,
                      const TEMPLATE(T, poly_t) op2,
                      const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void _TEMPLATE(T, poly_sub)(TEMPLATE(T, struct) *res,
                       const TEMPLATE(T, struct) *poly1, slong len1,
                       const TEMPLATE(T, struct) *poly2, slong len2,
                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, poly_sub)(TEMPLATE(T, poly_t) rop,
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

FQ_POLY_MINI_TEMPLATES_INLINE void
_TEMPLATE(T, poly_divrem)(TEMPLATE(T, struct) *Q, TEMPLATE(T, struct) *R,
                          const TEMPLATE(T, struct) *A, slong lenA,
                          const TEMPLATE(T, struct) *B, slong lenB,
                          const TEMPLATE(T, t) invB,
                          const TEMPLATE(T, ctx_t) ctx)
{
    _TEMPLATE(T, poly_divrem_divconquer)(Q, R, A, lenA, B, lenB, invB, ctx);
}

FQ_POLY_MINI_TEMPLATES_INLINE void
TEMPLATE(T, poly_divrem)(TEMPLATE(T, poly_t) Q, TEMPLATE(T, poly_t) R,
                         const TEMPLATE(T, poly_t) A,
                         const TEMPLATE(T, poly_t) B,
                         const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_divrem_divconquer)(Q, R, A, B, ctx);
}

FLINT_DLL void _TEMPLATE(T, poly_rem)(TEMPLATE(T, struct) *R,
                                const TEMPLATE(T, struct) *A, slong lenA,
                                const TEMPLATE(T, struct) *B, slong lenB,
                                const TEMPLATE(T, t) invB,
                                const TEMPLATE(T, ctx_t) ctx);

FQ_POLY_MINI_TEMPLATES_INLINE void
TEMPLATE(T, poly_rem)(TEMPLATE(T, poly_t) R,
                      const TEMPLATE(T, poly_t) A, const TEMPLATE(T, poly_t) B,
                      const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_t) Q;
    TEMPLATE(T, poly_init)(Q, ctx);
    TEMPLATE(T, poly_divrem)(Q, R, A, B, ctx);
    TEMPLATE(T, poly_clear)(Q, ctx);
}

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

FLINT_DLL void _TEMPLATE3(T, poly_scalar_mul, T)(TEMPLATE(T, struct) *rop,
                                  const TEMPLATE(T, struct) *op, slong len,
                                  const TEMPLATE(T, t) x,
                                  const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE3(T, poly_scalar_mul, T)(TEMPLATE(T, poly_t) rop,
                                 const TEMPLATE(T, poly_t) op,
                                 const TEMPLATE(T, t) x,
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

/* evaluation *****************************************************************/

FLINT_DLL void _TEMPLATE3(T, poly_evaluate, T)(TEMPLATE(T, t) rop,
                                     const TEMPLATE(T, struct) *op, slong len,
                                     const TEMPLATE(T, t) a,
                                     const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE3(T, poly_evaluate, T)(TEMPLATE(T, t) res,
                                    const TEMPLATE(T, poly_t) f,
                                    const TEMPLATE(T, t) a,
                                    const TEMPLATE(T, ctx_t) ctx);

#ifdef __cplusplus
}
#endif

#endif
