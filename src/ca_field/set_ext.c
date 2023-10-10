/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

/**

    Sets the extension number at position *i* (here indexed from 0) of *K*
    to the generator of the field with index *x_index* in *ctx*.
    (It is assumed that the generating field is a univariate field.)

    This only inserts a shallow reference: the field at index *x_index* must
    be kept alive until *K* has been cleared.
*/
void
ca_field_set_ext(ca_field_t K, slong i, ca_ext_srcptr x, ca_ctx_t ctx)
{
    CA_FIELD_EXT_ELEM(K, i) = (ca_ext_ptr) x;
    CA_FIELD_HASH(K) = CA_FIELD_HASH(K) * 100003 + CA_EXT_HASH(x);
}

