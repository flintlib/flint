/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

ca_field_ptr
_ca_ctx_get_field_const(ca_ctx_t ctx, calcium_func_code func)
{
    ca_ext_t ext;
    ca_ext_struct * ext_ptr[1];
    ca_field_ptr field;

    /* todo: shallow copy */
    ca_ext_init_const(ext, func, ctx);
    ext_ptr[0] = ca_ext_cache_insert(CA_CTX_EXT_CACHE(ctx), ext, ctx);
    field = ca_field_cache_insert_ext(CA_CTX_FIELD_CACHE(ctx), ext_ptr, 1, ctx);

    ca_ext_clear(ext, ctx);
    return field;
}
