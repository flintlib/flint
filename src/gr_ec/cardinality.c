/*
    Copyright (C) 2026 Maël Hostettler

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

/*
    Counting the points of E(F_q).

    Every algorithm returns #E(F_q) = q + 1 - t with |t| <= 2 sqrt(q), and
    each lives in a file of its own:

      cardinality_naive.c     O(q) field operations: walks over every x
                              and counts the roots in y. The reference, and
                              the only one independent of the group law.

      cardinality_bsgs.c      O(q^(1/4)) group operations, Shanks and
                              Mestre.

      cardinality_schoof.c    polynomial in log q: t mod l from Frobenius
                              on the l-torsion (torsion.c), Schoof's way
                              or, for SEA, Elkies' way (elkies.c, with the
                              modular polynomials of modular_poly.c).

      cardinality_cm.c        curves with complex multiplication by an
                              order of class number one, and supersingular
                              ones: the trace read off directly.

      cardinality_subfield.c  curves over F_{p^n} defined over F_p:
                              counted over F_p and lifted.

      cardinality_crt.c       curves over Z/NZ: counted modulo each prime
                              of N, lifted to p^k and multiplied.

    This file only chooses between them.
*/

/* ------------------------------------------------------------------ */
/* dispatch                                                           */
/* ------------------------------------------------------------------ */

/*
    Below this, walking the whole field is both quicker than setting up a
    baby-step table and free of any randomness.
*/
#define GR_EC_CARDINALITY_NAIVE_CUTOFF WORD(10000)

/*
    Past the naive range SEA wins at every size: measured against BSGS
    on random curves over prime fields, it is already ahead at 24 bits
    and twenty times faster at 64. BSGS remains the fallback for what SEA
    does not handle -- the long models of characteristic 2 and 3 -- and
    for a failure of SEA, which can only mean an unusable base ring.
*/

int
gr_ec_ctx_cardinality(fmpz_t res, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    fmpz_t q;
    int status;

    /*
        Over Z/NZ that is not known to be a field -- composite, or a prime
        nobody has vouched for -- count modulo each prime of N and put the
        counts together.
    */
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

    /*
        A curve over F_{p^n} whose coefficients happen to lie in F_p is a
        base change, and counting it over F_p and lifting the trace is
        enormously cheaper than counting it where it stands. Detecting
        that is five conversions, and over a prime field it declines
        immediately, so it is worth asking first.
    */
    if (gr_ec_ctx_cardinality_subfield(res, ctx) == GR_SUCCESS)
    {
        fmpz_clear(q);
        return GR_SUCCESS;
    }

    /*
        Complex multiplication next: it is a j-invariant comparison and a
        handful of scalar multiplications, and it gives up at once when the
        curve is not one of the special ones, so it costs almost nothing to
        try and saves everything when it applies.
    */
    if (gr_ec_ctx_cardinality_cm(res, ctx) == GR_SUCCESS)
    {
        fmpz_clear(q);
        return GR_SUCCESS;
    }

    /*
        What is left is a real count, and fmpz_mod is the slowest
        representation of F_p for it -- half the speed of mpn_mod at 256
        bits -- so the curve is copied over to nmod or mpn_mod when p fits
        one of them, and counted there.
    */
    if (R->which_ring == GR_CTX_FMPZ_MOD && fmpz_size(q) <= MPN_MOD_MAX_LIMBS)
    {
        status = _gr_ec_cardinality_mod_p(res, q, ctx);
        fmpz_clear(q);
        return status;
    }

    status = gr_ec_ctx_cardinality_sea(res, ctx);

    if (status != GR_SUCCESS)
        status = gr_ec_ctx_cardinality_bsgs(res, ctx);

    /* a curve whose group has small exponent can leave BSGS with more than
       one candidate; fall back while walking the field is still affordable */
    if (status != GR_SUCCESS && fmpz_cmp_si(q, GR_EC_NAIVE_MAX_Q) <= 0)
        status = gr_ec_ctx_cardinality_naive(res, ctx);

    fmpz_clear(q);

    return status;
}
