/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

int
unity_zp_equal(unity_zp f, unity_zp g)
{
    /*
        f and g can be reduced only by modylo x^{p^k} - 1,
        so reduce by cyclotomic polynomial
    */
    _unity_zp_reduce_cyclotomic(f);
    _unity_zp_reduce_cyclotomic(g);

    return fmpz_mod_poly_equal(f->poly, g->poly, f->ctx);
}

