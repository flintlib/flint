/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_POLY_MINI_H
#define FMPZ_MOD_POLY_MINI_H

#ifdef FMPZ_MOD_POLY_MINI_INLINES_C
#define FMPZ_MOD_POLY_MINI_INLINE FLINT_DLL
#else
#define FMPZ_MOD_POLY_MINI_INLINE static __inline__
#endif

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

FMPZ_MOD_POLY_MINI_INLINE
void fmpz_mod_poly_init(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
{
    poly->coeffs = NULL;
    poly->alloc  = 0;
    poly->length = 0;
}

FLINT_DLL void fmpz_mod_poly_clear(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

FLINT_DLL void _fmpz_mod_poly_fit_length(fmpz_mod_poly_t poly, slong len);
FMPZ_MOD_POLY_MINI_INLINE
void fmpz_mod_poly_fit_length(fmpz_mod_poly_t poly, slong len, const fmpz_mod_ctx_t ctx)
{
    _fmpz_mod_poly_fit_length(poly, len);
}

FLINT_DLL void fmpz_mod_poly_reverse(fmpz_mod_poly_t res,
                const fmpz_mod_poly_t poly, slong n, const fmpz_mod_ctx_t ctx);

/* NOTE: The rest of the declarations are only here becuase of fq_embed.h. */

FLINT_DLL void fmpz_mod_poly_derivative(fmpz_mod_poly_t res,
                const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx);

FLINT_DLL void fmpz_mod_poly_inv_series_newton(fmpz_mod_poly_t Qinv, 
                   const fmpz_mod_poly_t Q, slong n, const fmpz_mod_ctx_t ctx);

FMPZ_MOD_POLY_MINI_INLINE void 
fmpz_mod_poly_inv_series(fmpz_mod_poly_t Qinv, const fmpz_mod_poly_t Q,
                                             slong n, const fmpz_mod_ctx_t ctx)
{
   fmpz_mod_poly_inv_series_newton(Qinv, Q, n, ctx);
}

#ifdef __cplusplus
}
#endif

#endif
