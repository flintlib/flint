/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef PADIC_POLY_IMPL_H
#define PADIC_POLY_IMPL_H

#include <stdio.h>
#include "fmpz_mod_poly.h"
#include "padic_poly.h"

/* Defined in compose.c */
void __padic_reduce(fmpz_t u, slong *v, slong N, const padic_ctx_t ctx);

#endif
