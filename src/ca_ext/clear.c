/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_ext.h"
#include "ca_vec.h"

void
ca_ext_clear(ca_ext_t res, ca_ctx_t ctx)
{
    if (CA_EXT_HEAD(res) == CA_QQBar)
    {
        qqbar_clear(CA_EXT_QQBAR(res));
        nf_clear(CA_EXT_QQBAR_NF(res));
        flint_free(CA_EXT_QQBAR_NF(res));
    }
    else
    {
        if (CA_EXT_FUNC_NARGS(res) != 0)
            _ca_vec_clear(CA_EXT_FUNC_ARGS(res), CA_EXT_FUNC_NARGS(res), ctx);

        acb_clear(CA_EXT_FUNC_ENCLOSURE(res));

        if (res->data.func_data.qqbar != NULL)
        {
            qqbar_clear(res->data.func_data.qqbar);
            flint_free(res->data.func_data.qqbar);
        }
    }
}

