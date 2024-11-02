/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FMPZI_IMPL_H
#define FMPZI_IMPL_H

#include "fmpz_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#if FLINT_BITS == 64
# define GCD_MAX_D WORD(1125899906842623)
# define GCD_MIN_D WORD(-1125899906842623)
#else
# define GCD_MAX_D COEFF_MAX
# define GCD_MIN_D COEFF_MIN
#endif

void
_fmpzi_gcd_dddd(fmpzi_t res, double a, double b, double c, double d);

#ifdef __cplusplus
}
#endif

#endif /* FMPZI_IMPL_H */
