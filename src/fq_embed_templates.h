/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T
#ifdef B

#include "templates.h"

#ifdef __cplusplus
extern "C" {
#endif

void TEMPLATE(T, embed_gens)(TEMPLATE(T, t) gen_sub,
                                       TEMPLATE(T, t) gen_sup,
                                       TEMPLATE(B, poly_t) minpoly,
                                       const TEMPLATE(T, ctx_t) sub_ctx,
                                       const TEMPLATE(T, ctx_t) sup_ctx);
void _TEMPLATE(T, embed_gens_naive)(TEMPLATE(T, t) gen_sub,
                                              TEMPLATE(T, t) gen_sup,
                                              TEMPLATE(B, poly_t) minpoly,
                                              const TEMPLATE(T, ctx_t) sub_ctx,
                                              const TEMPLATE(T, ctx_t) sup_ctx);
void _TEMPLATE(T, embed_gens_allombert)(TEMPLATE(T, t) gen_sub,
                                                  TEMPLATE(T, t) gen_sup,
                                                  TEMPLATE(B, poly_t) minpoly,
                                                  const TEMPLATE(T, ctx_t) sub_ctx,
                                                  const TEMPLATE(T, ctx_t) sup_ctx);

/* Convert to-from column vectors */
void TEMPLATE(T, embed_matrices)(TEMPLATE(B, mat_t) embed,
                                           TEMPLATE(B, mat_t) project,
                                           const TEMPLATE(T, t) gen_sub,
                                           const TEMPLATE(T, ctx_t) sub_ctx,
                                           const TEMPLATE(T, t) gen_sup,
                                           const TEMPLATE(T, ctx_t) sup_ctx,
                                           const TEMPLATE(B, poly_t) gen_minpoly);

/* Given:

   - a context `sup_ctx` defining a field K of degree n,
   - a context `sub_ctx` defining a subfield k of degree m|n,
   - the n×m change-of-basis matrix from k to K,

   Return the m×n matrix of the trace from K to k.
   If m=n, res is the inverse of basis.
*/
void TEMPLATE(T, embed_trace_matrix)(TEMPLATE(B, mat_t) res,
                                         const TEMPLATE(B, mat_t) basis,
                                         const TEMPLATE(T, ctx_t) sub_ctx,
                                         const TEMPLATE(T, ctx_t) sup_ctx);

/* Compute the matrix whose columns are (gen^0, gen^1, ..., gen^(trunc-1)) */
void TEMPLATE(T, embed_composition_matrix_sub)(TEMPLATE(B, mat_t) matrix,
                                                   const TEMPLATE(T, t) gen,
                                                   const TEMPLATE(T, ctx_t) ctx,
                                                   slong trunc);

void TEMPLATE(T, embed_composition_matrix)(TEMPLATE(B, mat_t) matrix,
                                     const TEMPLATE(T, t) gen,
                                     const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, embed_mul_matrix)(TEMPLATE(B, mat_t) matrix,
                                       const TEMPLATE(T, t) gen,
                                       const TEMPLATE(T, ctx_t) ctx);

void TEMPLATE(T, embed_mono_to_dual_matrix)(TEMPLATE(B, mat_t) res,
                                                const TEMPLATE(T, ctx_t) ctx);
void TEMPLATE(T, embed_dual_to_mono_matrix)(TEMPLATE(B, mat_t) res,
                                                const TEMPLATE(T, ctx_t) ctx);

#ifdef __cplusplus
}
#endif

#endif
#endif
