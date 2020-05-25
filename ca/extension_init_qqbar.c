/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_extension_init_qqbar(ca_extension_t ext, const qqbar_t x)
{
    fmpq_poly_t t;

    ext->type = CA_EXT_QQBAR;
    qqbar_init(&ext->data.qqbar.x);
    qqbar_set(&ext->data.qqbar.x, x);

    /* nf_init wants an fmpq_poly_t, so mock up one */
    t->coeffs = QQBAR_POLY(x)->coeffs;
    t->den[0] = 1;
    t->length = QQBAR_POLY(x)->length;
    t->alloc = QQBAR_POLY(x)->alloc;

    nf_init(&ext->data.qqbar.nf, t);
}
