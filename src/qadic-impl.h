/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef QADIC_IMPL_H
#define QADIC_IMPL_H

#include "qadic.h"

#ifdef __cplusplus
extern "C" {
#endif

int
_artin_schreier_preimage(fmpz *rop, const fmpz *op, slong len,
                         const fmpz *a, const slong *j, slong lena);

#ifdef __cplusplus
}
#endif

#endif /* QADIC_IMPL_H */
