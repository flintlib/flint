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
ca_set_si(ca_t x, slong v, ca_ctx_t ctx)
{
    _ca_make_fmpq(x, ctx);
    fmpz_set_si(CA_FMPQ_NUMREF(x), v);
    fmpz_one(CA_FMPQ_DENREF(x));
}

