/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MOD_VEC_H
#define FMPZ_MOD_VEC_H

#ifdef FMPZ_MOD_VEC_INLINES_C
#define FMPZ_MOD_VEC_INLINE
#else
#define FMPZ_MOD_VEC_INLINE static inline \
    error fmpz_mod_vec/inlines.c does not exist
#endif

#include "fmpz_mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

void _fmpz_mod_vec_set_fmpz_vec(fmpz * A, const fmpz * B, slong len,
                                                     const fmpz_mod_ctx_t ctx);

void _fmpz_mod_vec_neg(fmpz * A, const fmpz * B, slong len,
                                                     const fmpz_mod_ctx_t ctx);

void _fmpz_mod_vec_add(fmpz * a, const fmpz * b, const fmpz * c,
                                            slong n, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_vec_sub(fmpz * a, const fmpz * b, const fmpz * c,
                                            slong n, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_vec_scalar_mul_fmpz_mod(fmpz * A, const fmpz * B,
                          slong len, const fmpz_t c, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_vec_scalar_addmul_fmpz_mod(fmpz * A, const fmpz * B,
                          slong len, const fmpz_t c, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_vec_mul(fmpz * A, const fmpz * B, const fmpz * C,
                                          slong len, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_vec_scalar_div_fmpz_mod(fmpz * A, const fmpz * B,
                          slong len, const fmpz_t c, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_vec_dot(fmpz_t d, const fmpz * A, const fmpz * B,
                                          slong len, const fmpz_mod_ctx_t ctx);

void _fmpz_mod_vec_dot_rev(fmpz_t r, const fmpz * a,
		          const fmpz * b, slong len, const fmpz_mod_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
