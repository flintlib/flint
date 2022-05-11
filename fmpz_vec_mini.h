/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_VEC_MINI_H
#define FMPZ_VEC_MINI_H

#ifdef FMPZ_VEC_INLINES_C
#define FMPZ_VEC_INLINE FLINT_DLL
#else
#define FMPZ_VEC_INLINE static __inline__
#endif

#include "flint.h"

#ifdef __cplusplus
extern "C" {
#endif

FLINT_DLL void _fmpz_vec_demote(fmpz * vec, slong len);

FLINT_DLL int _fmpz_vec_is_zero(const fmpz * vec, slong len);

FLINT_DLL void _fmpz_vec_height(fmpz_t height, const fmpz * vec, slong len);

FLINT_DLL slong _fmpz_vec_max_bits(const fmpz * vec, slong len);

FLINT_DLL mp_mock_size_t _fmpz_vec_max_limbs(const fmpz * vec, slong len);

FLINT_DLL void _fmpz_vec_content(fmpz_t res, const fmpz * vec, slong len);

#ifdef __cplusplus
}
#endif

#endif
