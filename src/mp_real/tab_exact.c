/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "arb.h"
#include "mp_real.h"
#include "impl.h"

/* Exact floors for the table entries.

   The tables guarantee that every entry equals the floor of the
   infinite-precision value at the full entry precision -- hence at
   every shorter truncation too, floors being nested.  The fast
   generation paths (native binary splitting and the fixed-point
   multi-summation) produce one-sided values at most a few ulps
   below the truth at the bottom (guard) limb, which coincide with
   the exact floor UNLESS the true value sits within that deficit
   above a limb boundary -- detectable, after the fact, as a guard
   limb within the error bound of 0xfff...: the deficit then
   borrowed through.  The ensure functions rescan for such guard
   limbs and recompute the (astronomically rare) offenders here, and
   route the small-index arb-combination entries through
   _mp_real_tab_store_floor with this as the retry fallback.  Nothing
   here may use the elementary functions of arb, which are built on
   these tables (the thread filling a table holds its lock). */

/* write e = floor(x 2^(64 nc)) if x's radius determines that floor
   uniquely; returns 0 otherwise.  prec must be at least the working
   precision of x: the floor itself is taken at that precision, since
   the value can sit arbitrarily close below an integer (the tail
   entries atan(2^-i) ~ 2^-i - 2^(-3i)/3 do, once the correction
   falls below the entry) and a fixed floor precision would never
   separate them from it. */
int
_mp_real_tab_store_floor(nn_ptr e, const arb_t x, slong nc, slong prec)
{
    arb_t y;
    fmpz_t f;
    int ok;

    arb_init(y);
    fmpz_init(f);
    arb_mul_2exp_si(y, x, FLINT_BITS * nc);
    arb_floor(y, y, FLINT_MAX(prec, FLINT_BITS * nc) + 64);
    ok = arb_get_unique_fmpz(f, y);
    if (ok)
    {
        FLINT_ASSERT(fmpz_sgn(f) >= 0
            && fmpz_bits(f) <= FLINT_BITS * nc);
        fmpz_get_ui_array(e, nc, f);
    }
    arb_clear(y);
    fmpz_clear(f);
    return ok;
}

/* e = floor(v 2^(64 nc)) exactly, v = log(1 + 2^-i) (which = 0) or
   atan(2^-i) (which = 1), retrying at increasing precision until
   the floor is unique.  The values come from the Machin-type binary
   splitting of mp_real (log(1 + 2^-i) = 2 atanh(1 / (2^(i+1) + 1)),
   atan(2^-i) = atan(1 / 2^i), and atan(1) = atan(1/2) + atan(1/3)),
   which uses no shared table: the elementary functions of arb are
   built on the tables being filled here. */
void
_mp_real_tab_entry_exact(nn_ptr e, int which, ulong i, slong nc)
{
    arb_t x;
    mp_real_t v, u;
    nn_ptr q;
    slong qn, wn;
    ulong p = 1;

    arb_init(x);
    mp_real_init(v);
    mp_real_init(u);

    /* q = 2^(i+1) + 1 or 2^i (an i of any size: the rescans of the
       tables call this for every entry) */
    qn = (slong) ((i + 1) / FLINT_BITS) + 1;
    q = flint_calloc(qn, sizeof(ulong));

    for (wn = nc + 1; ; wn *= 2)
    {
        if (which == 0)
        {
            flint_mpn_zero(q, qn);
            q[(i + 1) / FLINT_BITS] = UWORD(1) << ((i + 1) % FLINT_BITS);
            q[0] |= 1;
            _mp_real_atanh_frac_bsplit(v, &p, 1, q, qn, wn);
            mp_real_mul_2exp_si(v, v, 1);
        }
        else if (i == 0)
        {
            q[0] = 2;
            _mp_real_atan_frac_bsplit(v, &p, 1, q, 1, wn);
            q[0] = 3;
            _mp_real_atan_frac_bsplit(u, &p, 1, q, 1, wn);
            mp_real_add(v, v, u, wn);
        }
        else
        {
            flint_mpn_zero(q, qn);
            q[i / FLINT_BITS] = UWORD(1) << (i % FLINT_BITS);
            _mp_real_atan_frac_bsplit(v, &p, 1, q, qn, wn);
        }

        mp_real_get_arb(x, v);
        if (_mp_real_tab_store_floor(e, x, nc, FLINT_BITS * wn))
            break;
    }

    flint_free(q);
    arb_clear(x);
    mp_real_clear(v);
    mp_real_clear(u);
}
