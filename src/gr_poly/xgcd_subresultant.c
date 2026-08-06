/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

/* Port of _fmpz_poly_pseudo_divrem_cohen to GR.
   Produces (Q, R) with lc(B)^e * A = Q*B + R, e = lenA - lenB + 1.
   R is computed in-place on the A buffer (R may alias A); Q must not alias A or B.
   On return, *lenR_out holds the normalised length of R. */
static int
_gr_poly_pseudo_divrem_cohen(gr_ptr Q, slong lenQ,
                                   gr_ptr R, slong * lenR_out,
                                   gr_srcptr A, slong lenA,
                                   gr_srcptr B, slong lenB,
                                   gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    gr_srcptr leadB = GR_ENTRY(B, lenB - 1, sz);
    slong lenR = lenA;
    slong e;

    if (R != A)
        status |= _gr_vec_set(R, A, lenA, ctx);

    status |= _gr_vec_zero(Q, lenQ, ctx);
    e = lenA - lenB;

    status |= gr_set(GR_ENTRY(Q, lenQ - 1, sz), GR_ENTRY(R, lenR - 1, sz), ctx);
    status |= _gr_vec_mul_scalar(R, R, lenR - 1, leadB, ctx);
    status |= _gr_vec_submul_scalar(GR_ENTRY(R, lenR - lenB, sz), B, lenB - 1,
                                    GR_ENTRY(R, lenR - 1, sz), ctx);
    status |= gr_zero(GR_ENTRY(R, lenR - 1, sz), ctx);
    lenR--;
    status |= _gr_vec_normalise(&lenR, R, lenR, ctx);
    if (status != GR_SUCCESS) goto done;

    while (lenR >= lenB)
    {
        status |= _gr_vec_mul_scalar(Q, Q, lenQ, leadB, ctx);
        status |= gr_add(GR_ENTRY(Q, lenR - lenB, sz),
                         GR_ENTRY(Q, lenR - lenB, sz),
                         GR_ENTRY(R, lenR - 1, sz), ctx);
        status |= _gr_vec_mul_scalar(R, R, lenR - 1, leadB, ctx);
        status |= _gr_vec_submul_scalar(GR_ENTRY(R, lenR - lenB, sz), B, lenB - 1,
                                        GR_ENTRY(R, lenR - 1, sz), ctx);
        status |= gr_zero(GR_ENTRY(R, lenR - 1, sz), ctx);
        lenR--;
        status |= _gr_vec_normalise(&lenR, R, lenR, ctx);
        if (status != GR_SUCCESS) goto done;
        e--;
    }

    if (e == 1)
    {
        status |= _gr_vec_mul_scalar(Q, Q, lenQ, leadB, ctx);
        status |= _gr_vec_mul_scalar(R, R, lenR, leadB, ctx);
    }
    else if (e > 1)
    {
        gr_ptr pow;
        GR_TMP_INIT(pow, ctx);
        status |= gr_pow_ui(pow, leadB, (ulong) e, ctx);
        status |= _gr_vec_mul_scalar(Q, Q, lenQ, pow, ctx);
        status |= _gr_vec_mul_scalar(R, R, lenR, pow, ctx);
        GR_TMP_CLEAR(pow, ctx);
    }

done:
    *lenR_out = lenR;
    return status;
}

static int
_gr_vec_content(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op gcd_op = GR_BINARY_OP(ctx, GCD);
    int status = GR_SUCCESS;
    slong i, sz = ctx->sizeof_elem;

    if (len == 0) return gr_zero(res, ctx);
    if (len == 1) return gr_set(res, vec, ctx);

    status |= gcd_op(res, GR_ENTRY(vec, 0, sz), GR_ENTRY(vec, 1, sz), ctx);
    for (i = 2; i < len; i++)
        status |= gcd_op(res, res, GR_ENTRY(vec, i, sz), ctx);
    return status;
}

/*
    Extended subresultant PRS for gr_poly.

    Given A, B in R[x] over a GCD domain R with lenA >= lenB >= 2, produces
    G, S, T satisfying  S*A + T*B = G  where G is an associate of gcd(A,B).

    Let a = cont(A), b = cont(B), d = gcd(a,b), a' = a/d, b' = b/d.
    Write A = a*A_p, B = b*B_p with A_p, B_p primitive.

    The core PRS runs on (A_p, B_p).  If it yields S_p, T_p, G_p with
        S_p * A_p + T_p * B_p = G_p,
    multiplying through by lcm(a,b) = ab/d gives:
        (b'*S_p) * A  +  (a'*T_p) * B  =  lcm(a,b) * G_p.
    So:  G = lcm(a,b)*G_p,  S = b'*S_p,  T = a'*T_p.

    T_p is recovered by exact polynomial division using the PRIMITIVE polynomials:
        T_p * B_p = G_p - S_p * A_p.
    Using the originals A, B would introduce non-integer denominators.
    Fixed copies Aprim and Bprim are stored before the loop alters curA/curB.

    lcm_ab = lcm(a,b) = a * b' is computed before the loop, since the scalar
    slots a and b are reused as scratch by the Brown cofactor update.

    Only the S-cofactor is tracked through the PRS; T_p is recovered at the end.
    At each step, pseudo-division gives
        lc(curB)^e * curA  =  Q * curB + R,   e = deg curA - deg curB + 1,
    and after dividing by the subresultant cofactor c = g_sc * h^delta:
        sR_new = (lc(curB)^e * sA - Q * sB) / c          [exact]
        R_new  = prem(curA, curB) / c                    [exact]
*/

int
_gr_poly_xgcd_subresultant(slong * lenG, gr_ptr G, gr_ptr S, gr_ptr T,
                            gr_srcptr A, slong lenA,
                            gr_srcptr B, slong lenB,
                            gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    FLINT_ASSERT(lenA >= lenB);
    FLINT_ASSERT(lenB >= 2);

    /*
    Offset            Size      Purpose
    ---------------   -------   -------
    0                 lenA      curA   - PRS polynomial  (starts as A_prim)
    lenA              lenB      curB   - PRS polynomial  (starts as B_prim)
    lenA+lenB         lenA      Aprim  - fixed copy of A_prim for T-recovery
    2*lenA+lenB       lenB      Bprim  - fixed copy of B_prim for T-recovery
    2*lenA+2*lenB     2*lenA    sA     - S-cofactor for curA  (starts [1])
    4*lenA+2*lenB     2*lenA    sB     - S-cofactor for curB  (starts [0])
    6*lenA+2*lenB     2*lenA    scr    - result buffer for sR each iteration
    8*lenA+2*lenB     2*lenA    qsb    - temporary for Q*sB product
    10*lenA+2*lenB    lenA      Q      - pseudo-quotient
    11*lenA+2*lenB    7         scalars: g_sc h a b ap bp lcm_ab

    Total: 11*lenA + 2*lenB + 7 elements.

    sA, sB, scr rotate as a group (3-way pointer swap each iteration);
    each holds a cofactor polynomial of degree < lenA by the Bézout bound,
    so 2*lenA capacity suffices.

    qsb is a separate fixed scratch buffer for the Q*sB product, which
    has length lenQ+lenSB-1 <= (lenA-lenB+1)+(lenA-1)-1+1 = 2*lenA-lenB
    <= 2*lenA, so 2*lenA capacity suffices.
     */
    slong Wsz = 11 * lenA + 2 * lenB + 7;
    gr_ptr W;
    gr_ptr curA, curB, Aprim, Bprim, sA, sB, scr, qsb, Q;
    gr_ptr g_sc, h, a, b, ap, bp, lcm_ab;

    GR_TMP_INIT_VEC(W, Wsz, ctx);
    curA   = GR_ENTRY(W,     0,            sz);
    curB   = GR_ENTRY(W,     lenA,         sz);
    Aprim  = GR_ENTRY(W,     lenA+lenB,    sz);
    Bprim  = GR_ENTRY(W,     2*lenA+lenB,  sz);
    sA     = GR_ENTRY(W,     2*lenA+2*lenB, sz);
    sB     = GR_ENTRY(W,     4*lenA+2*lenB, sz);
    scr    = GR_ENTRY(W,     6*lenA+2*lenB, sz);
    qsb    = GR_ENTRY(W,     8*lenA+2*lenB, sz);
    Q      = GR_ENTRY(W,     10*lenA+2*lenB, sz);
    g_sc   = GR_ENTRY(W,     11*lenA+2*lenB, sz);
    h      = GR_ENTRY(g_sc,  1, sz);
    a      = GR_ENTRY(h,     1, sz);
    b      = GR_ENTRY(a,     1, sz);
    ap     = GR_ENTRY(b,     1, sz);
    bp     = GR_ENTRY(ap,    1, sz);
    lcm_ab = GR_ENTRY(bp,    1, sz);

    slong lenCurA = lenA, lenCurB = lenB;
    slong lenSA, lenSB, lenSCR, lenQ, lenR;

    slong lenSmax = lenB - 1;
    slong lenTmax = lenA - 1;

    status |= _gr_vec_content(a, A, lenA, ctx);
    status |= _gr_vec_content(b, B, lenB, ctx);
    status |= _gr_vec_divexact_scalar(curA, A, lenA, a, ctx);
    status |= _gr_vec_divexact_scalar(curB, B, lenB, b, ctx);
    status |= _gr_vec_set(Aprim, curA, lenA, ctx);
    status |= _gr_vec_set(Bprim, curB, lenB, ctx);
    gr_ptr d = scr;
    status |= gr_gcd(d, a, b, ctx);
    status |= gr_divexact(ap, a, d, ctx);
    status |= gr_divexact(bp, b, d, ctx);
    status |= gr_mul(lcm_ab, a, bp, ctx);

    if (status != GR_SUCCESS)
        goto cleanup;

    status |= gr_one(sA, ctx);
    lenSA = 1;
    lenSB = 0;
    status |= gr_one(g_sc, ctx);
    status |= gr_one(h, ctx);

    while (lenCurB > 1)
    {
        slong delta = lenCurA - lenCurB;
        slong e     = lenCurA - lenCurB + 1;
        gr_srcptr lc_B = GR_ENTRY(curB, lenCurB - 1, sz);

        lenQ = lenCurA - lenCurB + 1;

        status |= _gr_poly_pseudo_divrem_cohen(Q, lenQ,
                curA, &lenR, curA, lenCurA, curB, lenCurB, ctx);

        if (status != GR_SUCCESS)
            goto cleanup;

        /* sR_new = lc(B)^e * sA  -  Q * sB  into scr */
        if (e == 1)
            status |= gr_set(a, lc_B, ctx);
        else
            status |= gr_pow_ui(a, lc_B, (ulong) e, ctx);

        status |= _gr_vec_mul_scalar(scr, sA, lenSA, a, ctx);
        lenSCR = lenSA;

        if (lenSB > 0)
        {
            slong lenQSB = lenQ + lenSB - 1;
            status |= _gr_poly_mul(qsb, Q, lenQ, sB, lenSB, ctx);
            status |= _gr_poly_sub(scr, scr, lenSCR, qsb, lenQSB, ctx);
            lenSCR = FLINT_MAX(lenSCR, lenQSB);
            status |= _gr_vec_normalise(&lenSCR, scr, lenSCR, ctx);
        }

        if (status != GR_SUCCESS)
            goto cleanup;

        FLINT_SWAP(gr_ptr, curA, curB);
        lenCurA = lenCurB;
        lenCurB = lenR;
        { gr_ptr _t = sA; sA = sB; sB = scr; scr = _t; }
        { slong  _l = lenSA; lenSA = lenSB; lenSB = lenSCR; lenSCR = _l; }

        if (lenCurB == 0)
            break;

        /* Divide curB and sB by the subresultant cofactor c = g_sc * h^delta. */
        if (delta == 0)
        {
            status |= _gr_vec_divexact_scalar(curB, curB, lenCurB, g_sc, ctx);
            status |= _gr_vec_divexact_scalar(sB,   sB,   lenSB,   g_sc, ctx);
        }
        else if (delta == 1)
        {
            status |= gr_mul(b, g_sc, h, ctx);
            status |= _gr_vec_divexact_scalar(curB, curB, lenCurB, b, ctx);
            status |= _gr_vec_divexact_scalar(sB,   sB,   lenSB,   b, ctx);
        }
        else
        {
            status |= gr_pow_ui(a, h, (ulong) delta, ctx);
            status |= gr_mul(b, g_sc, a, ctx);
            status |= _gr_vec_divexact_scalar(curB, curB, lenCurB, b, ctx);
            status |= _gr_vec_divexact_scalar(sB,   sB,   lenSB,   b, ctx);
        }

        status |= _gr_vec_normalise(&lenCurB, curB, lenCurB, ctx);

        /* Update Brown cofactors:
             g_sc_new = lc(curA),  h_new = h * lc(curA)^delta / h^delta. */
        if (delta == 0)
        {
            status |= gr_set(g_sc, GR_ENTRY(curA, lenCurA - 1, sz), ctx);
        }
        else if (delta == 1)
        {
            status |= gr_set(g_sc, GR_ENTRY(curA, lenCurA - 1, sz), ctx);
            status |= gr_set(h, g_sc, ctx);
        }
        else
        {
            status |= gr_pow_ui(a, h, (ulong) delta, ctx);
            status |= gr_pow_ui(b, GR_ENTRY(curA, lenCurA - 1, sz), (ulong) delta, ctx);
            status |= gr_mul(h, h, b, ctx);
            status |= gr_divexact(h, h, a, ctx);
            status |= gr_set(g_sc, GR_ENTRY(curA, lenCurA - 1, sz), ctx);
        }

        if (status != GR_SUCCESS)
            goto cleanup;
    }

    /* Termination.

        lenCurB == 0: G_prim = curA, S_p = sA.
        lenCurB >= 1: G_prim = curB, S_p = sB.

       Content scaling:
         G = lcm_ab * G_prim
         S = bp * S_p
         T = ap * T_p,  where T_p = (G_prim - S_p*Aprim) / Bprim  [exact] */
    {
        gr_srcptr G_prim;
        gr_srcptr Sp;
        slong lenG_prim, lenSp, lenTp;

        if (lenCurB == 0)
        {
            G_prim = curA;
            lenG_prim = lenCurA;
            Sp = sA;
            lenSp = lenSA;
        }
        else
        {
            G_prim = curB;
            lenG_prim = lenCurB;
            Sp = sB;
            lenSp = lenSB;
        }

        /* T_p = (G_prim - S_p * Aprim) / Bprim  using the primitive inputs. */
        if (lenSp == 0)
        {
            if (lenG_prim < lenB)
                lenTp = 0;
            else
            {
                lenTp = lenG_prim - lenB + 1;
                status |= _gr_poly_divexact(T, G_prim, lenG_prim, Bprim, lenB, ctx);
            }
        }
        else
        {
            slong lenSA_prod = lenSp + lenA - 1;
            status |= _gr_poly_mul(scr, Sp, lenSp, Aprim, lenA, ctx);
            status |= _gr_vec_normalise(&lenSA_prod, scr, lenSA_prod, ctx);
            status |= _gr_poly_sub(scr, G_prim, lenG_prim, scr, lenSA_prod, ctx);
            slong lenNum = FLINT_MAX(lenG_prim, lenSA_prod);
            status |= _gr_vec_normalise(&lenNum, scr, lenNum, ctx);
            if (lenNum < lenB)
                lenTp = 0;
            else
            {
                lenTp = lenNum - lenB + 1;
                status |= _gr_poly_divexact(T, scr, lenNum, Bprim, lenB, ctx);
            }
        }

        if (status != GR_SUCCESS)
            goto cleanup;

        /* Write G = lcm_ab * G_prim. */
        status |= _gr_vec_set(G, G_prim, lenG_prim, ctx);
        *lenG = lenG_prim;
        if (gr_is_one(lcm_ab, ctx) != T_TRUE)
            status |= _gr_vec_mul_scalar(G, G, lenG_prim, lcm_ab, ctx);

        /* Write S = bp * S_p. */
        status |= _gr_vec_set(S, Sp, lenSp, ctx);
        if (gr_is_one(bp, ctx) != T_TRUE && lenSp > 0)
            status |= _gr_vec_mul_scalar(S, S, lenSp, bp, ctx);

        /* Scale T = ap * T_p. */
        if (gr_is_one(ap, ctx) != T_TRUE && lenTp > 0)
            status |= _gr_vec_mul_scalar(T, T, lenTp, ap, ctx);

        /* The interface doesn't allow returning the lengths of S and T,
           so explicitly zero the high parts. */
        status |= _gr_vec_zero(GR_ENTRY(S, lenSp, sz), lenSmax - lenSp, ctx);
        status |= _gr_vec_zero(GR_ENTRY(T, lenTp, sz), lenTmax - lenTp, ctx);
    }

cleanup:
    GR_TMP_CLEAR_VEC(W, Wsz, ctx);
    if (status != GR_SUCCESS)
    {
        *lenG = 0;
    }

    return status;
}

int
gr_poly_xgcd_subresultant(gr_poly_t G, gr_poly_t S, gr_poly_t T,
                           const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx)
{
    return gr_poly_xgcd_wrapper((gr_method_poly_xgcd_op) _gr_poly_xgcd_subresultant, G, S, T, A, B, ctx);
}
