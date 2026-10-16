/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

/*
    One balanced Hensel block.

    Given a 2*bn-limb dividend window W, the divisor b (bn limbs) and its inverse
    binv = b^(-1) mod B^bn, compute the bn quotient limbs qb that cancel the low
    bn limbs of W (qb = W_low * binv mod B^bn, so that qb*b == W (mod B^bn)), and
    return the high part of the reduced window,

        R = (W - qb*b) / B^bn        (low bn limbs),

    which is W_high - (qb*b)_high. The borrow out of that subtraction is returned
    (it propagates into limbs above the window). t is 2*bn limbs of scratch.

    Only the high bn limbs of the 2*bn-limb product qb*b are needed, and its low
    bn limbs are known in advance: they equal W_low (since qb*b == W (mod B^bn)).
    When bn and the limb radix are both large enough this is exploited exactly as
    in radix_invmod_bn: the high part is formed by an approximate (carry-less)
    middle product with 3 guard limbs. Subtracting the known low limbs
    W_low[bn-3..bn) makes the guard limbs have true value zero, after which the
    high part is either exact or 1 ulp too small; the deficit occurs iff the
    carry-in dropped by the middle product was nonzero, which (because that carry
    is < B^2 when B >= bn) is signalled by a nonzero top guard limb. Otherwise
    the full 2*bn-limb product is formed.

    R may alias W_low (= W); the reads needed happen before R is written.
*/
static ulong
_radix_divmod_bn_block(nn_ptr qb, nn_ptr R, nn_srcptr W, nn_srcptr b,
    nn_srcptr binv, nn_ptr t, slong bn, const radix_t radix)
{
    radix_mulmid(qb, W, bn, binv, bn, 0, bn, radix);      /* qb = W_low * binv mod B^bn */

    if (bn > 3 && LIMB_RADIX(radix) >= (ulong) bn)
    {
        ulong one = 1;

        /* t[0..bn+3) = carry-less (qb*b)[bn-3 .. 2bn) (3 guard limbs at the bottom) */
        radix_mulmid(t, qb, bn, b, bn, bn - 3, 2 * bn, radix);
        /* Subtract the known low limbs so the guard limbs have true value zero. */
        radix_sub(t, t, bn + 3, W + (bn - 3), 3, radix);
        /* Nonzero top guard => high part is 1 ulp too small. */
        if (t[2] != 0)
            radix_add(t + 3, t + 3, bn, &one, 1, radix);
        /* t[3..bn+3) = (qb*b) >> bn; R = W_high - that. */
        return radix_sub(R, W + bn, bn, t + 3, bn, radix);
    }

    radix_mulmid(t, qb, bn, b, bn, 0, 2 * bn, radix);     /* t = qb * b (2*bn limbs) */
    /* t_low == W_low, so (W - t) is divisible by B^bn and (W - t)/B^bn = W_high - t_high. */
    return radix_sub(R, W + bn, bn, t + bn, bn, radix);
}

/*
    Classical (schoolbook) Hensel division.

    Given a of length an and b of length bn with b[0] invertible modulo the limb
    radix B, develop n limbs of the p-adic (Hensel) quotient q with

        q * b == a   (mod B^n),

    working from least significant to most significant limb. This is the p-adic
    analogue of the unbalanced Euclidean division in radix_divrem: there an
    an x bn division is run as a sequence of (2*bn) x bn block divisions from the
    top down; here the blocks run from the bottom up, each producing bn quotient
    limbs via one length-bn multiply and subtract (see _radix_divmod_bn_block).
    The bn == 1 case is exactly the Hensel analogue of radix_divrem_1.

    The running remainder window threads through the blocks just as in the
    Euclidean loop, so the remainder falls out of the final block (with an
    O(bn^2) adjustment when n is not a multiple of bn). The cost is O(n*bn) --
    linear in bn for fixed precision n -- hence cheaper than forming q = a*b^(-1)
    (two ~n-limb multiplications) when n is large and b is short.

    q receives n quotient limbs (caller provides room for n) and may alias a.
    Returns 1 on success and 0 if b[0] is not invertible modulo B (in which case
    nothing is written).

    If rem != NULL it receives the resume (Henselian) remainder (a - q*b)/B^n
    reduced to its low bn limbs (room for bn required), unnormalised; then
    a == q*b + B^n*rem (mod B^{n+bn}). The caller can test for an exact division
    with two zero tests:

        exact  <=>  flint_mpn_zero_p(rem, bn)
                    && flint_mpn_zero_p(a + (n + bn), FLINT_MAX(an - (n + bn), 0));

    (the second test being vacuous when an <= n + bn).

    Apart from q and rem the routine uses only O(bn) scratch, so it operates in
    place on q (q may alias a).
*/
int
radix_divmod_bn_classical(nn_ptr q, nn_ptr rem, nn_srcptr a, slong an,
    nn_srcptr b, slong bn, slong n, const radix_t radix)
{
    nn_ptr binv, W, t, qlast;
    slong blk, blocks, ext;
    slong pend;
    TMP_INIT;

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(n >= 1);

    /* The bn == 1 case has a dedicated, much leaner implementation. */
    if (bn == 1)
        return radix_divmod_bn_1(q, rem, a, an, b[0], n, radix);

    TMP_START;

    /* binv = b^(-1) mod B^bn; returns 0 iff b[0] is a non-unit. */
    binv = TMP_ALLOC(bn * sizeof(ulong));
    if (!radix_invmod_bn(binv, b, bn, bn, radix))
    {
        TMP_END;
        return 0;
    }

    W     = TMP_ALLOC((2 * bn) * sizeof(ulong));   /* running dividend window  */
    t     = TMP_ALLOC((2 * bn) * sizeof(ulong));   /* block product scratch    */
    qlast = TMP_ALLOC(bn * sizeof(ulong));         /* full final-block quotient */

    blocks = (n + bn - 1) / bn;
    ext = blocks * bn - n;                          /* over-computed limbs in last block, 0..bn-1 */

    /* Seed the window with the low (up to) 2*bn limbs of a. */
    {
        slong lo_avail = FLINT_MIN(an, bn);
        flint_mpn_copyi(W, a, lo_avail);
        flint_mpn_zero(W + lo_avail, 2 * bn - lo_avail);
        if (an > bn)
            flint_mpn_copyi(W + bn, a + bn, FLINT_MIN(an - bn, bn));
    }

    pend = 0;                                       /* borrow pending above the window */

    for (blk = 0; blk < blocks; blk++)
    {
        slong qoff = blk * bn;

        if (blk < blocks - 1)
        {
            /* Full block: bn quotient limbs straight into q, reduced window in W. */
            ulong bw = _radix_divmod_bn_block(q + qoff, W, W, b, binv, t, bn, radix);

            /* Shift in the next bn limbs of a as the new high part. */
            slong nextoff = (blk + 2) * bn;
            if (nextoff < an)
            {
                slong nextavail = FLINT_MIN(bn, an - nextoff);
                flint_mpn_copyi(W + bn, a + nextoff, nextavail);
                flint_mpn_zero(W + bn + nextavail, bn - nextavail);
            }
            else
            {
                flint_mpn_zero(W + bn, bn);
            }

            /* Apply the discarded borrow (plus any carried from earlier all-zero
               blocks) to the freshly loaded high part; a new borrow carries on. */
            {
                ulong one = 1;
                slong tot = pend + (slong) bw;
                pend = 0;
                while (tot > 0)
                {
                    if (radix_sub(W + bn, W + bn, bn, &one, 1, radix))
                        pend++;
                    tot--;
                }
            }
        }
        else
        {
            /* Final block. */
            slong qlen = n - qoff;                  /* = bn - ext */

            if (rem == NULL)
            {
                /* Only the kept quotient limbs are needed; skip the reduction. */
                radix_mulmid(q + qoff, W, bn, binv, bn, 0, qlen, radix);
            }
            else
            {
                /* Reduce the full block so the remainder can be read off W. */
                nn_ptr qb = (ext != 0) ? qlast : (q + qoff);
                (void) _radix_divmod_bn_block(qb, W, W, b, binv, t, bn, radix);
                if (ext != 0)
                    flint_mpn_copyi(q + qoff, qlast, qlen);
            }
        }
    }

    /*
        W[0..bn) now holds (a - q_full*b)/B^{blocks*bn} mod B^bn, where q_full is
        the (blocks*bn)-limb quotient. When n is a multiple of bn this is the
        resume remainder. Otherwise we kept only the low qlen = bn-ext limbs of
        the last block; with E = qlast[qlen..bn) the high ext discarded limbs, the
        resume remainder at precision n is

            Rem_n = E*b + B^ext * Rem_{blocks*bn}    (an O(bn^2) adjustment),

        of which we keep the low bn limbs.
    */
    if (rem != NULL)
    {
        if (ext == 0)
        {
            flint_mpn_copyi(rem, W, bn);
        }
        else
        {
            nn_srcptr E = qlast + (bn - ext);
            radix_mulmid(rem, b, bn, E, ext, 0, bn, radix);   /* low bn limbs of E*b */
            /* add W_low shifted up by ext digits (mod B^bn) */
            radix_add(rem + ext, rem + ext, bn - ext, W, bn - ext, radix);
        }
    }

    TMP_END;
    return 1;
}
