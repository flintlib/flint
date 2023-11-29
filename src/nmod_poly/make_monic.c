/*
    Copyright (C) 2007, David Howden.
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void _nmod_poly_make_monic(mp_ptr output,
                            mp_srcptr input, slong len, nmod_t mod)
{
    mp_limb_t inv;

    inv = n_invmod(input[len - 1], mod.n);
    _nmod_vec_scalar_mul_nmod(output, input, len, inv, mod);
}

void nmod_poly_make_monic(nmod_poly_t output, const nmod_poly_t input)
{
    if (input->length == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_make_monic). Division by zero.\n");
    }

    nmod_poly_fit_length(output, input->length);
    _nmod_poly_make_monic(output->coeffs,
                            input->coeffs, input->length, input->mod);
    output->length = input->length;
}

