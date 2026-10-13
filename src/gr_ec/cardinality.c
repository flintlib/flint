/*
    Copyright (C) 2026 Maël Hostettler
    Developed using Claude Opus 5 and Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "mpn_mod.h"
#include "gr.h"
#include "gr_ec.h"
#include "impl.h"

/* below this, counting naively is fastest */
#define GR_EC_CARDINALITY_NAIVE_CUTOFF WORD(10000)

/*
    Choose a point counting algorithm. SEA is used whenever it applies, as
    it beats BSGS from about 24 bits on; BSGS covers the long Weierstrass
    model in characteristic 2 and 3.
*/
int
_gr_ec_ctx_cardinality(fmpz_t res, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    fmpz_t q;
    int status;

    /* Z/NZ with N not known to be prime */
    if (gr_ctx_is_field(R) != T_TRUE)
        return gr_ec_ctx_cardinality_crt(res, ctx);

    fmpz_init(q);

    if (gr_ctx_cardinality_fmpz(q, R) != GR_SUCCESS)
    {
        fmpz_clear(q);
        return GR_UNABLE;
    }

    if (fmpz_cmp_si(q, GR_EC_CARDINALITY_NAIVE_CUTOFF) <= 0)
    {
        status = gr_ec_ctx_cardinality_naive(res, ctx);

        if (status == GR_SUCCESS)
        {
            fmpz_clear(q);
            return status;
        }
    }

    /* cheap checks for special curves first */
    if (gr_ec_ctx_cardinality_subfield(res, ctx) == GR_SUCCESS)
    {
        fmpz_clear(q);
        return GR_SUCCESS;
    }

    if (gr_ec_ctx_cardinality_cm(res, ctx) == GR_SUCCESS)
    {
        fmpz_clear(q);
        return GR_SUCCESS;
    }

    /* fmpz_mod is slow: count over nmod or mpn_mod instead */
    if (R->which_ring == GR_CTX_FMPZ_MOD && fmpz_size(q) <= MPN_MOD_MAX_LIMBS
            && _gr_ec_cardinality_mod_p(res, q, ctx) == GR_SUCCESS)
    {
        fmpz_clear(q);
        return GR_SUCCESS;
    }

    status = gr_ec_ctx_cardinality_sea(res, ctx);

    if (status != GR_SUCCESS)
        status = gr_ec_ctx_cardinality_bsgs(res, ctx);

    /* BSGS can fail when the group has a small exponent */
    if (status != GR_SUCCESS && fmpz_cmp_si(q, GR_EC_NAIVE_MAX_Q) <= 0)
        status = gr_ec_ctx_cardinality_naive(res, ctx);

    fmpz_clear(q);

    return status;
}

int
gr_ec_ctx_cardinality(fmpz_t res, gr_ec_ctx_t ctx)
{
    int status;

    if (!fmpz_is_zero(&GR_EC_CTX(ctx)->cardinality))
    {
        fmpz_set(res, &GR_EC_CTX(ctx)->cardinality);
        return GR_SUCCESS;
    }

    status = _gr_ec_ctx_cardinality(res, ctx);

    if (status == GR_SUCCESS)
        fmpz_set(&GR_EC_CTX(ctx)->cardinality, res);

    return status;
}
