/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T
#ifdef B

#include "templates.h"

#ifdef __cplusplus
extern "C" {
#endif

FLINT_DLL void TEMPLATE(T, embed_gens)(TEMPLATE(T, t) gen_sub,
                                       TEMPLATE(T, t) gen_sup,
                                       TEMPLATE(B, poly_t) minpoly,
                                       const TEMPLATE(T, ctx_t) sub_ctx,
                                       const TEMPLATE(T, ctx_t) sup_ctx);
FLINT_DLL void _TEMPLATE(T, embed_gens_naive)(TEMPLATE(T, t) gen_sub,
                                              TEMPLATE(T, t) gen_sup,
                                              TEMPLATE(B, poly_t) minpoly,
                                              const TEMPLATE(T, ctx_t) sub_ctx,
                                              const TEMPLATE(T, ctx_t) sup_ctx);
FLINT_DLL void _TEMPLATE(T, embed_gens_allombert)(TEMPLATE(T, t) gen_sub,
                                                  TEMPLATE(T, t) gen_sup,
                                                  TEMPLATE(B, poly_t) minpoly,
                                                  const TEMPLATE(T, ctx_t) sub_ctx,
                                                  const TEMPLATE(T, ctx_t) sup_ctx);

/* Given:

   - a context `sup_ctx` defining a field K of degree n,
   - a context `sub_ctx` defining a subfield k of degree m|n,
   - the n×m change-of-basis matrix from k to K,

   Return the m×n matrix of the trace from K to k.
   If m=n, res is the inverse of basis. */
FLINT_DLL void TEMPLATE(T, trace_matrix)(TEMPLATE(B, mat_t) res,
                                         const TEMPLATE(B, mat_t) basis,
                                         const TEMPLATE(T, ctx_t) sub_ctx,
                                         const TEMPLATE(T, ctx_t) sup_ctx);

/* Compute the matrix whose columns are (gen^0, gen^1, ..., gen^(trunc-1)) */
FLINT_DLL void TEMPLATE(T, composition_matrix_sub)(TEMPLATE(B, mat_t) matrix,
                                                   const TEMPLATE(T, t) gen,
                                                   const TEMPLATE(T, ctx_t) ctx,
                                                   slong trunc);
FQ_EMBED_TEMPLATES_INLINE
void TEMPLATE(T, composition_matrix)(TEMPLATE(B, mat_t) matrix,
                                     const TEMPLATE(T, t) gen,
                                     const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, composition_matrix_sub)(matrix, gen, ctx, 
                                        TEMPLATE(T, ctx_degree(ctx)));
}

FLINT_DLL void TEMPLATE(T, mul_matrix)(TEMPLATE(B, mat_t) matrix,
                                       const TEMPLATE(T, t) gen,
                                       const TEMPLATE(T, ctx_t) ctx);

FLINT_DLL void TEMPLATE(T, mono_to_dual_matrix)(TEMPLATE(B, mat_t) res,
                                                const TEMPLATE(T, ctx_t) ctx);
FLINT_DLL void TEMPLATE(T, dual_to_mono_matrix)(TEMPLATE(B, mat_t) res,
                                                const TEMPLATE(T, ctx_t) ctx);

FQ_EMBED_TEMPLATES_INLINE
void TEMPLATE(T, modulus_pow_series_inv)(TEMPLATE(B, poly_t) res,
                                         const TEMPLATE(T, ctx_t) ctx,
                                         slong trunc)
{
    TEMPLATE(B, poly_reverse)(res, 
                              TEMPLATE(T, ctx_modulus)(ctx), 
                              TEMPLATE(T, ctx_degree)(ctx) + 1);
    TEMPLATE(B, poly_inv_series)(res, res, trunc);
}

#ifdef __cplusplus
}
#endif

#endif
#endif
