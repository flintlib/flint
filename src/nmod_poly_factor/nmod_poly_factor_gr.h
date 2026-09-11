/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_FACTOR_GR_H
#define NMOD_POLY_FACTOR_GR_H

#include "nmod_poly.h"
#include "nmod_poly_factor.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    An nmod_poly_struct is a gr_poly_struct with the modulus appended, so
    the two cannot be cast to each other; but the coefficient data is
    laid out identically, so a polynomial can be viewed as the other type
    without copying the coefficients.

    NMOD_POLY_AS_GR views an nmod_poly as a gr_poly. The view shares the
    coefficient array: it must not be used for output (which could
    reallocate) unless the length is written back with
    NMOD_POLY_FROM_GR.
*/
#define NMOD_POLY_AS_GR(gpoly, npoly) \
    do { \
        (gpoly)->coeffs = (npoly)->coeffs; \
        (gpoly)->alloc = (npoly)->alloc; \
        (gpoly)->length = (npoly)->length; \
    } while (0)

#define NMOD_POLY_FROM_GR(npoly, gpoly) \
    do { \
        (npoly)->coeffs = (gpoly)->coeffs; \
        (npoly)->alloc = (gpoly)->alloc; \
        (npoly)->length = (gpoly)->length; \
    } while (0)

void _nmod_poly_factor_gr(nmod_poly_factor_t res, const nmod_poly_t f, int algorithm);

/* Move the polynomials in fac (with exponents exp) into res, appending to
   it; the entries of fac are left zero. */
NMOD_POLY_FACTOR_INLINE void
_nmod_poly_factor_set_gr(nmod_poly_factor_t res, gr_poly_vec_t fac,
    const fmpz_vec_t exp, nmod_t mod, gr_ctx_t gr_ctx)
{
    slong i, num = fac->length;

    nmod_poly_factor_fit_length(res, res->num + num);

    for (i = 0; i < num; i++)
    {
        nmod_poly_struct * p = res->p + res->num + i;

        nmod_poly_clear(p);
        NMOD_POLY_FROM_GR(p, fac->entries + i);
        p->mod = mod;

        gr_poly_init(fac->entries + i, gr_ctx);

        res->exp[res->num + i] = fmpz_get_si(exp->entries + i);
    }

    res->num += num;
}

#ifdef __cplusplus
}
#endif

#endif
