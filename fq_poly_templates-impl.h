/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FQ_POLY_TEMPLATES_IMPL_H
#define FQ_POLY_TEMPLATES_IMPL_H

#include "fmpz_poly.h"

/* defined in mul_reorder.c */
fmpz_poly_struct * __vec_init(slong len);
fmpz_poly_struct * __vec_init2(slong len, slong n);
void __vec_clear(fmpz_poly_struct * v, slong len);
void __scalar_addmul(fmpz_poly_struct * rop, const fmpz_poly_struct * op, slong len, const fmpz_poly_t x);
void __scalar_mul(fmpz_poly_struct * rop, const fmpz_poly_struct * op, slong len, const fmpz_poly_t x);

#endif
