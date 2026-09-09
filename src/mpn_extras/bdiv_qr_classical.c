/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

/*
    One balanced Hensel block (port of _radix_divmod_bn_block).

    Given a 2*bn-limb dividend window W, the divisor b (bn limbs) and its
    inverse binv = b^(-1) mod B^bn, compute the bn quotient limbs
    qb = W_low binv mod B^bn (so that qb b == W (mod B^bn)) and return the
    high part of the reduced window R = (W - qb b) / B^bn (low bn limbs)
    together with the borrow out of that subtraction.

    Both products use the assembly-optimised low and high multiplication:
    flint_mpn_mulhigh_n is a lower approximation of the top bn + 1 limbs of
    qb b whose deficit is below B (at most 2 bn units of the returned guard
    limb), and the true value of the guard limb, limb bn - 1 of qb b, is the
    known W[bn-1]; a borrow into the high part occurred iff the computed
    guard limb exceeds it. t is bn limbs of scratch. R may alias W.
*/
static mp_limb_t
_flint_mpn_bdiv_qr_block(mp_ptr qb, mp_ptr R, mp_srcptr W, mp_srcptr b,
    mp_srcptr binv, mp_ptr t, mp_size_t bn)
{
    mp_limb_t g;

    flint_mpn_mullow_n(qb, W, binv, bn);
    g = flint_mpn_mulhigh_n(t, qb, b, bn);
    if (g > W[bn - 1])
        mpn_add_1(t, t, bn, 1);
    return mpn_sub_n(R, W + bn, t, bn);
}

/*
    Classical (schoolbook) Hensel division, port of radix_divmod_bn_classical
    with B = 2^64: develops n limbs of q with q b == a (mod B^n), working from
    the least significant limb upwards in blocks of bn quotient limbs, each
    costing one bn x bn low product and one high product. The cost is
    O(n bn / bn * M(bn)), which beats the Karp-Markstein method when b is
    short relative to n.

    q may alias a. If r != NULL it receives (a - q b) / B^n mod B^bn (room
    for bn limbs), so that a == q b + B^n r (mod B^(n+bn)).
*/
/* the same with a precomputed inverse binv = b^(-1) mod B^bn (bn >= 2) */
void
_flint_mpn_bdiv_qr_classical_preinv(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an,
    mp_srcptr b, mp_size_t bn, mp_srcptr binv, mp_size_t n)
{
    mp_ptr W, t, qlast;
    mp_size_t blk, blocks, ext, pend;
    TMP_INIT;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(bn >= 2);
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(b[0] & 1);

    TMP_START;

    /* short quotient: a single low product */
    if (n <= bn && r == NULL)
    {
        if (an >= n)
        {
            flint_mpn_mullow_n(q, a, binv, n);
        }
        else
        {
            mp_ptr ap = TMP_ALLOC(n * sizeof(mp_limb_t));
            flint_mpn_copyi(ap, a, an);
            flint_mpn_zero(ap + an, n - an);
            flint_mpn_mullow_n(q, ap, binv, n);
        }
        TMP_END;
        return;
    }

    W = TMP_ALLOC((2 * bn + bn + bn) * sizeof(mp_limb_t));
    t = W + 2 * bn;             /* block product scratch */
    qlast = t + bn;             /* full final-block quotient */

    blocks = (n + bn - 1) / bn;
    ext = blocks * bn - n;      /* over-computed limbs in the last block */

    /* seed the window with the low (up to) 2*bn limbs of a */
    {
        mp_size_t avail = FLINT_MIN(an, 2 * bn);
        flint_mpn_copyi(W, a, avail);
        flint_mpn_zero(W + avail, 2 * bn - avail);
    }

    pend = 0;                   /* borrows pending above the window */

    for (blk = 0; blk < blocks; blk++)
    {
        mp_size_t qoff = blk * bn;

        if (blk < blocks - 1)
        {
            mp_limb_t bw = _flint_mpn_bdiv_qr_block(q + qoff, W, W, b, binv, t, bn);
            mp_size_t nextoff = (blk + 2) * bn;

            /* shift in the next bn limbs of a as the new high part */
            if (nextoff < an)
            {
                mp_size_t nextavail = FLINT_MIN(bn, an - nextoff);
                flint_mpn_copyi(W + bn, a + nextoff, nextavail);
                flint_mpn_zero(W + bn + nextavail, bn - nextavail);
            }
            else
            {
                flint_mpn_zero(W + bn, bn);
            }

            /* apply the discarded borrow (plus any carried from earlier
               blocks) to the freshly loaded high part */
            {
                mp_size_t tot = pend + (mp_size_t) bw;
                pend = 0;
                while (tot > 0)
                {
                    if (mpn_sub_1(W + bn, W + bn, bn, 1))
                        pend++;
                    tot--;
                }
            }
        }
        else
        {
            /* final block */
            mp_size_t qlen = n - qoff;      /* = bn - ext */

            if (r == NULL)
            {
                /* only the kept quotient limbs are needed */
                if (qlen == bn)
                    flint_mpn_mullow_n(q + qoff, W, binv, bn);
                else
                    flint_mpn_mulmid(q + qoff, W, bn, binv, bn, 0, qlen);
            }
            else
            {
                mp_ptr qb = (ext != 0) ? qlast : (q + qoff);
                (void) _flint_mpn_bdiv_qr_block(qb, W, W, b, binv, t, bn);
                if (ext != 0)
                    flint_mpn_copyi(q + qoff, qlast, qlen);
            }
        }
    }

    /*
        W[0, bn) now holds (a - q_full b) / B^(blocks bn) mod B^bn for the
        (blocks bn)-limb quotient q_full. When ext = 0 this is the remainder.
        Otherwise the discarded high ext limbs E = qlast[qlen, bn) of the
        last block give the remainder at precision n as
        E b + B^ext W (low bn limbs).
    */
    if (r != NULL)
    {
        if (ext == 0)
        {
            flint_mpn_copyi(r, W, bn);
        }
        else
        {
            mp_srcptr E = qlast + (bn - ext);
            flint_mpn_mulmid(r, b, bn, E, ext, 0, bn);
            mpn_add_n(r + ext, r + ext, W, bn - ext);
        }
    }

    TMP_END;
}

void
flint_mpn_bdiv_qr_classical(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an,
    mp_srcptr b, mp_size_t bn, mp_size_t n)
{
    mp_ptr binv;
    TMP_INIT;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(b[0] & 1);

    if (bn == 1)
    {
        flint_mpn_bdiv_qr_1(q, r, a, an, b[0], n);
        return;
    }

    TMP_START;
    binv = TMP_ALLOC(bn * sizeof(mp_limb_t));
    flint_mpn_binv(binv, b, bn, bn);
    _flint_mpn_bdiv_qr_classical_preinv(q, r, a, an, b, bn, binv, n);
    TMP_END;
}
