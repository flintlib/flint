/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

ulong ca_hash_repr(const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_SPECIAL(x))
    {
        return 123;  /* not really interesting ... */
    }
    else if (CA_IS_QQ(x, ctx))
    {
        ulong a, b;

        a = calcium_fmpz_hash(CA_FMPQ_NUMREF(x));
        b = calcium_fmpz_hash(CA_FMPQ_DENREF(x));

        return a + 781237663 * b;
    }
    else if (CA_FIELD_IS_NF(CA_FIELD(x, ctx)))
    {
        const fmpz * num;
        const fmpz * den;
        slong i, len;
        ulong hash;

        if (CA_FIELD_NF(CA_FIELD(x, ctx))->flag & NF_LINEAR)
        {
            num = (fmpz *) LNF_ELEM_NUMREF(CA_NF_ELEM(x));
            den = LNF_ELEM_DENREF(CA_NF_ELEM(x));
            len = 1;
        }
        else if (CA_FIELD_NF(CA_FIELD(x, ctx))->flag & NF_QUADRATIC)
        {
            num = (fmpz *) QNF_ELEM_NUMREF(CA_NF_ELEM(x));
            den = QNF_ELEM_DENREF(CA_NF_ELEM(x));
            len = 2;
        }
        else
        {
            num = (fmpz *) NF_ELEM_NUMREF(CA_NF_ELEM(x));
            den = NF_ELEM_DENREF(CA_NF_ELEM(x));
            len = NF_ELEM(CA_NF_ELEM(x))->length;
        }

        hash = CA_EXT_HASH(CA_FIELD_EXT_ELEM(CA_FIELD(x, ctx), 0));

        hash = 1000003 * calcium_fmpz_hash(den) + hash;
        for (i = 0; i < len; i++)
            hash = calcium_fmpz_hash(num + i) * 1000003 + hash;

        return hash;
    }
    else
    {
        slong i;
        ulong hash;
        const fmpz_mpoly_struct *p, *q;

        hash = CA_FIELD_HASH(CA_FIELD(x, ctx));

        p = fmpz_mpoly_q_numref(CA_MPOLY_Q(x));
        q = fmpz_mpoly_q_numref(CA_MPOLY_Q(x));

        /* fixme: this only looks at coefficients -- accessing exponents is annoying */
        for (i = 0; i < p->length; i++)
            hash = calcium_fmpz_hash(p->coeffs + i) * 1000003 + hash;

        for (i = 0; i < q->length; i++)
            hash = calcium_fmpz_hash(q->coeffs + i) * 1000003 + hash;

        return hash;
    }
}
