/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_IMPL_H
#define ACB_IMPL_H

#include "acb_types.h"

void acb_gamma_stirling_eval(acb_t s, const acb_t z, slong nterms, int digamma, slong prec);

/* the elementary functions at one and two limbs by the register paths
   shared with acf (acb/elem_small.c): f(z) into r1 (sin into r1 and cos
   into r2, either NULL, for the pairs); returns 0 when the path does not
   apply (prec > 2 FLINT_BITS, a zero or special part, a midpoint of more
   than two limbs, a radius of 2^-8 or more, an argument the register
   paths decline) */
#define ACB_ELEM_EXP 0
#define ACB_ELEM_EXP_PI_I 1
#define ACB_ELEM_SIN_COS 2
#define ACB_ELEM_SIN_COS_PI 3
#define ACB_ELEM_LOG 4
int _acb_elem_small(acb_t r1, acb_t r2, const acb_t z, int f, slong prec);

#endif
