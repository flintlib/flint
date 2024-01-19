/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "nmod_poly_factor.h"

void
nmod_poly_factor_print(const nmod_poly_factor_t fac)
{
    slong i;

    for (i = 0; i < fac->num; i++)
    {
        nmod_poly_print(fac->p + i);
        flint_printf(" ^ %wd\n", fac->exp[i]);
    }
}

void nmod_poly_factor_print_pretty(const nmod_poly_factor_t fac, const char *var)
{
    slong i;
    for (i = 0; i < fac->num; i++)
    {
        nmod_poly_print_pretty(fac->p + i, var);
        flint_printf(" ^ %wd\n", fac->exp[i]);
    }
}
