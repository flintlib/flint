/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_POLY_FACTOR_IMPL_H
#define GR_POLY_FACTOR_IMPL_H

#include "fmpz.h"
#include "gr.h"
#include "gr_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Information about a finite field F_q, q = p^d. Returns GR_UNABLE
   if ctx cannot be verified to be a finite field. */
GR_POLY_INLINE int
_gr_poly_factor_ff_info(fmpz_t q, fmpz_t p, slong * d, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (gr_ctx_is_field(ctx) != T_TRUE ||
        gr_ctx_is_finite_characteristic(ctx) != T_TRUE ||
        gr_ctx_is_finite(ctx) != T_TRUE)
        return GR_UNABLE;

    if (q != NULL)
        status |= gr_ctx_fq_order(q, ctx);
    if (p != NULL)
        status |= gr_ctx_fq_prime(p, ctx);
    if (d != NULL)
        status |= gr_ctx_fq_degree(d, ctx);

    return (status == GR_SUCCESS) ? GR_SUCCESS : GR_UNABLE;
}

/* Polynomials of length at most this are factored serially even when
   threads are available: below it, the cost of waking the threads
   exceeds the parallel work (measured with nmod, 8 threads: parallel
   sections help from degree ~300 on). The variable is provided so
   that the test code can exercise the parallel code paths with small
   inputs. */
#ifndef GR_POLY_FACTOR_THREADED_CUTOFF
# define GR_POLY_FACTOR_THREADED_CUTOFF 200
#endif

FLINT_DLL extern slong gr_poly_factor_threaded_cutoff;

/* Helper for functions computing a factorization of a monic polynomial:
   abort with GR_UNABLE if the status is not GR_SUCCESS or if the leading
   coefficient of pol cannot be verified to be nonzero. */
#define GR_POLY_FACTOR_CHECK(pol) \
    if (status != GR_SUCCESS || ((pol)->length > 0 && \
        gr_is_zero(gr_poly_coeff_srcptr((pol), (pol)->length - 1, ctx), ctx) != T_FALSE)) \
    { \
        status |= GR_UNABLE; \
        goto cleanup; \
    }

#define GR_POLY_FACTOR_CHECK_STATUS() \
    if (status != GR_SUCCESS) \
    { \
        status |= GR_UNABLE; \
        goto cleanup; \
    }

#ifdef __cplusplus
}
#endif

#endif
