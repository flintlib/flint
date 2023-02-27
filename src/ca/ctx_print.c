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

void
ca_ctx_print(ca_ctx_t ctx)
{
    slong i;

    flint_printf("Calcium context with %wd cached fields:\n", CA_CTX_FIELD_CACHE(ctx)->length);
    for (i = 0; i < CA_CTX_FIELD_CACHE(ctx)->length; i++)
    {
        flint_printf("%wd   ", i);
        ca_field_print(CA_CTX_FIELD_CACHE(ctx)->items[i], ctx);
        flint_printf("\n");
    }
    flint_printf("\n");
}

