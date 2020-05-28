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
ca_field_init_mpoly_q(ca_field_t K, slong len)
{
    slong i;
    K->len = len;

    K->type = CA_FIELD_MPOLY_Q;
    fmpz_mpoly_ctx_init(&K->mctx, len, ORD_LEX);
    K->ext = flint_malloc(len * sizeof(ca_extension_struct *));

    for (i = 0; i < len; i++)
        K->ext[i] = NULL;

    K->ideal = NULL;
    K->ideal_len = 0;
}

void
ca_field_init_nf(ca_field_t K, ca_extension_struct * ext)
{
    K->len = 1;
    K->type = CA_FIELD_NF;
    K->nf_ext = ext;
    K->ideal = NULL;
    K->ideal_len = 0;
}

void
ca_field_init_qq(ca_field_t K)
{
    K->len = 0;
    K->type = CA_FIELD_QQ;
    K->ext = NULL;
    K->ideal = NULL;
    K->ideal_len = 0;
}

