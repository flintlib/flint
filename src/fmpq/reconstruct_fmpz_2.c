/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "mpn_extras.h"
#include "fmpq.h"

#define FMPQ_RECONSTRUCT_ARRAY_LIMIT 12

/*
    hgcd for two-limb input, individual quotients not written
*/
static slong _hgcd_uiui_no_write(
    mp_limb_t A1, mp_limb_t A0,
    mp_limb_t B1, mp_limb_t B0,
    _ui_mat22_t M)
{
    slong written = 0; /* number of quotients generated */
    mp_limb_t last_written = 0;
    mp_limb_t d0, d1, t0, t1, t2, r0, r1;
    int det = 1;
    mp_limb_t m11 = 1;
    mp_limb_t m12 = 0;
    mp_limb_t m21 = 0;
    mp_limb_t m22 = 1;
    mp_limb_t a1 = A1;
    mp_limb_t a0 = A0;
    mp_limb_t b1 = B1;
    mp_limb_t b0 = B0;
    mp_limb_t q;

    FLINT_ASSERT(a1 != 0);
    FLINT_ASSERT(b1 < a1 || (b1 == a1 && b0 <= a0));

    if (b1 == 0 || !(b1 < a1 || (b1 == a1 && b0 < a0)))
        goto done;

    while (1)
    {
        eudiv_qrrnndd(q, r1, r0, a1, a0, b1, b0, tmp1_);

        t1 = m12 + q*m11;
        t2 = m22 + q*m21;

        if (r1 == 0)
            break;

        a0 = b0;
        a1 = b1;
        b0 = r0;
        b1 = r1;

        m12 = m11;
        m22 = m21;
        m11 = t1;
        m21 = t2;
        det *= -1;

        last_written = q;
        written++;
    }

    FLINT_ASSERT(a1 != 0);
    FLINT_ASSERT(b1 <= a1);
    FLINT_ASSERT(b1 < a1 || (b1 == a1 && b0 < a0));

    sub_ddmmss(d1,d0, a1,a0, b1,b0);
    if (det > 0)
    {
        if (b1 == 0 && b0 <= m21)
            goto fix;
        add_ssaaaa(t1,t0, 0,m11, 0,m12);
    }
    else
    {
        if (b1 == 0 && b0 <= m11)
            goto fix;
        add_ssaaaa(t1,t0, 0,m21, 0,m22);
    }
    if (d1 < t1 || (d1 == t1 && d0 < t0))
        goto fix;

fixed:

#ifdef FLINT_WANT_ASSERT

    /* should be ok */

    sub_ddmmss(d1,d0, a1,a0, b1,b0);
    if (det > 0)
    {
        FLINT_ASSERT(!(b1 == 0 && b0 <= m21));
        add_ssaaaa(t1,t0, 0,m11, 0,m12);
    }
    else
    {
        FLINT_ASSERT(!(b1 == 0 && b0 <= m11));
        add_ssaaaa(t1,t0, 0,m21, 0,m22);
    }
    FLINT_ASSERT(!(d1 < t1 || (d1 == t1 && d0 < t0)));

    /* should have {A, B} == {{m11, m12}, {m21, m22}} . {a, b} */

    umul_ppmm(t1,t0, a0, m11); t1 += a1*m11;
    umul_ppmm(d1,d0, b0, m12); d1 += b1*m12;
    add_sssaaaaaa(t2,t1,t0, 0,t1,t0, 0,d1,d0);
    FLINT_ASSERT(t0 == A0 && t1 == A1 && t2 == 0);

    umul_ppmm(t1,t0, a0, m21); t1 += a1*m21;
    umul_ppmm(d1,d0, b0, m22); d1 += b1*m22;
    add_sssaaaaaa(t2,t1,t0, 0,t1,t0, 0,d1,d0);
    FLINT_ASSERT(t0 == B0 && t1 == B1 && t2 == 0);

    /* should have det = Det[{{m11, m12}, {m21, m22}}] */

    umul_ppmm(t1,t0, m11, m22);
    umul_ppmm(d1,d0, m12, m21);
    sub_ddmmss(t1,t0, t1,t0, d1,d0);
    FLINT_ASSERT(t1 == FLINT_SIGN_EXT(t0));
    FLINT_ASSERT(t0 == det);

#endif

done:

    M->_11 = m11;
    M->_12 = m12;
    M->_21 = m21;
    M->_22 = m22;
    M->det = det;

    return written;

fix:

    written--;
    FLINT_ASSERT(written >= 0);

    FLINT_ASSERT(last_written != 0);

    q = last_written;

    t1 = m11 - q*m12;
    t2 = m21 - q*m22;
    m11 = m12;
    m21 = m22;
    m12 = t1;
    m22 = t2;
    det *= -1;

#ifdef FLINT_WANT_ASSERT
    umul_ppmm(t1,t0, a0, q); t1 += a1*q;
    add_ssaaaa(t1,t0, t1,t0, b1,b0);
    b0 = a0;
    b1 = a1;
    a0 = t0;
    a1 = t1;
#endif

    goto fixed;
}


/* u is odd */
static int coprime_ui(mp_limb_t u, mp_limb_t v)
{
    FLINT_ASSERT(u > 0);
    FLINT_ASSERT(v > 0);

    while ((v & 1) == 0) {v = v/2;}

again:

    FLINT_ASSERT(u & 1);
    FLINT_ASSERT(v & 1);

    if (u > v)
    {
        u -= v;
        do {
            u = u >> 1;
        } while ((u & 1) == 0);
        goto again;
    }
    else if (u < v)
    {
        v -= u;
        do {
            v = v >> 1;
        } while ((v & 1) == 0);
        goto again;
    }
    else
    {
        return u == 1;
    }
}

/* u is odd */
static int coprime_uiui(mp_limb_t u1, mp_limb_t u0, mp_limb_t v1, mp_limb_t v0)
{
    FLINT_ASSERT(u1 > 0 || u0 > 0);
    FLINT_ASSERT(v1 > 0 || v0 > 0);

    while ((v0 & 1) == 0)
    {
        v0 = MPN_RIGHT_SHIFT_LOW(v1, v0, 1);
        v1 = v1 >> 1;
    }

again:

    FLINT_ASSERT(u0 & 1);
    FLINT_ASSERT(v0 & 1);

    if (u1 > v1)
    {
        sub_ddmmss(u1, u0, u1, u0, v1, v0);
        do {
            u0 = MPN_RIGHT_SHIFT_LOW(u1, u0, 1);
            u1 = u1 >> 1;
        } while ((u0 & 1) == 0);
        goto again;
    }
    else if (v1 > u1)
    {
        sub_ddmmss(v1, v0, v1, v0, u1, u0);
        do {
            v0 = MPN_RIGHT_SHIFT_LOW(v1, v0, 1);
            v1 = v1 >> 1;
        } while ((v0 & 1) == 0);
        goto again;
    }
    else if (u0 > v0)
    {
        return coprime_ui(v0, u0 - v0);
    }
    else if (v0 > u0)
    {
        return coprime_ui(u0, v0 - u0);
    }
    else
    {
        return u1 == 0 && u0 == 1;
    }
}


int _fmpq_reconstruct_fmpz_2_ui(fmpz_t n, fmpz_t d,
              const fmpz_t a, const fmpz_t m, const fmpz_t NN, const fmpz_t DD)
{
    mp_limb_t Q, R, A, B, N;
    mp_limb_t m11 = 1, m12 = 0, t;
    int mdet = 1;

    FLINT_ASSERT(fmpz_size(m) == 1);

    A = fmpz_get_ui(m);
    B = fmpz_get_ui(a);
    N = fmpz_get_ui(NN);

gauss:

    FLINT_ASSERT(A > B && B > N);

    eudiv_qrnd(Q, R, A, B, tmp_);
    mdet *= -1;
    t = m12 + m11*Q;
    m12 = m11;
    m11 = t;
    A = B;
    B = R;

    if (B > N)
        goto gauss;

    FLINT_ASSERT(A > N && N >= B);

    if (fmpz_cmp_ui(DD, m11) < 0)
        return 0;

    if (mdet > 0)
        fmpz_set_ui(n, B);
    else
        fmpz_neg_ui(n, B);

    fmpz_set_ui(d, m11);

    FLINT_ASSERT(m11 != 0);
    if (B == 0)
        return m11 == 1;

    if (B & 1)
        return coprime_ui(B, m11);
    else if (m11 & 1)
        return coprime_ui(m11, B);
    else
        return 0;
}

int _fmpq_reconstruct_fmpz_2_uiui(fmpz_t n, fmpz_t d,
              const fmpz_t a, const fmpz_t m, const fmpz_t NN, const fmpz_t DD)
{
    mp_limb_t extra;
    mp_limb_t Q1, Q0, R1, R0, A1, A0, B1, B0, N1, N0, D1, D0;
    mp_limb_t m11[2] = {1, 0}, m12[2] = {0, 0}, t[2];
    int mdet = 1;

    FLINT_ASSERT(fmpz_size(m) == 2);

    fmpz_get_uiui(&A1, &A0, m);
    fmpz_get_uiui(&B1, &B0, a);
    fmpz_get_uiui(&N1, &N0, NN);
    fmpz_get_uiui(&D1, &D0, DD);

gauss:

    FLINT_ASSERT(A1 > B1 || (A1 == B1 && A0 > B0));
    FLINT_ASSERT(B1 > N1 || (B1 == N1 && B0 > N0));

    if (A1 != 0 && B1 != 0)
    {
        eudiv_qrrnndd(Q0, R1, R0, A1, A0, B1, B0, tmp1_);
        extra = 0;
    }
    else if (A1 == 0 && B1 == 0)
    {
        eudiv_qrnd(Q0, R0, A0, B0, tmp2_);
        R1 = 0;
        extra = 0;
    }
    else
    {
        eudiv_qqrnnd(Q1, Q0, R0, A1, A0, B0, tmp3_);
        R1 = 0;
        extra = m11[0]*Q1;
    }

    umul_ppmm(t[1], t[0], m11[0], Q0);
    add_ssaaaa(t[1], t[0], t[1], t[0], m12[1], m12[0]);
    t[1] += m11[1]*Q0 + extra;

    mdet *= -1;
    m12[1] = m11[1]; m12[0] = m11[0];
    m11[1] = t[1]; m11[0] = t[0];
    A1 = B1; A0 = B0;
    B1 = R1; B0 = R0;

    if (B1 > N1 || (B1 == N1 && B0 > N0))
        goto gauss;

    if (D1 < m11[1] || (D1 == m11[1] && D0 < m11[0]))
        return 0;

    if (mdet > 0)
        fmpz_set_uiui(n, B1, B0);
    else
        fmpz_neg_uiui(n, B1, B0);

    fmpz_set_uiui(d, m11[1], m11[0]);

    if (B1 == 0 && B0 == 0)
        return m11[1] == 0 && m11[0] == 1;

    if (m11[0] & 1)
        return coprime_uiui(m11[1], m11[0], B1, B0);
    else if (B0 & 1)
        return coprime_uiui(B1, B0, m11[1], m11[0]);
    else
        return 0;
}

int _fmpq_reconstruct_fmpz_2_ui_array(fmpz_t n, fmpz_t d,
              const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)
{
    mp_limb_t ex0, ex1, ex2, ex3, A1, A0, B1, B0;
    unsigned int n_lzcnt, a_lzcnt;
    _ui_mat22_t h;
    slong written;
    const mp_limb_t * n_ptr, * d_ptr;
    slong n_len, d_len;
    mp_limb_t A[FMPQ_RECONSTRUCT_ARRAY_LIMIT + 1];
    mp_limb_t B[FMPQ_RECONSTRUCT_ARRAY_LIMIT + 1];
    mp_limb_t Q[FMPQ_RECONSTRUCT_ARRAY_LIMIT + 1];
    mp_limb_t R[FMPQ_RECONSTRUCT_ARRAY_LIMIT + 1];
    mp_limb_t m11[FMPQ_RECONSTRUCT_ARRAY_LIMIT + 2];
    mp_limb_t m12[FMPQ_RECONSTRUCT_ARRAY_LIMIT + 2];
    slong Alen, Blen, Qlen, Rlen, m_len;
    int mdet = 1;

    FLINT_ASSERT(fmpz_size(m) <= FMPQ_RECONSTRUCT_ARRAY_LIMIT);

    fmpz_get_ui_array(A, FMPQ_RECONSTRUCT_ARRAY_LIMIT, m);
    Alen = FMPQ_RECONSTRUCT_ARRAY_LIMIT;
    MPN_NORM(A, Alen);

    fmpz_get_ui_array(B, FMPQ_RECONSTRUCT_ARRAY_LIMIT, a);
    Blen = FMPQ_RECONSTRUCT_ARRAY_LIMIT;
    MPN_NORM(B, Blen);

    if (COEFF_IS_MPZ(*N))
    {
        n_ptr = COEFF_TO_PTR(*N)->_mp_d;
        n_len = COEFF_TO_PTR(*N)->_mp_size;
    }
    else
    {
        n_ptr = (mp_srcptr) N; /* haha, dirty but works */
        n_len = 1;
    }

    if (COEFF_IS_MPZ(*D))
    {
        d_ptr = COEFF_TO_PTR(*D)->_mp_d;
        d_len = COEFF_TO_PTR(*D)->_mp_size;
    }
    else
    {
        d_ptr = (mp_srcptr) D; /* haha, dirty but works */
        d_len = 1;
    }

    /* m12 is supposed to be zero-extended to the length of m11 (= m_len) */
    flint_mpn_zero(m11, FMPQ_RECONSTRUCT_ARRAY_LIMIT);
    flint_mpn_zero(m12, FMPQ_RECONSTRUCT_ARRAY_LIMIT);
    /* m11 = 1, m12 = 0 */
    m_len = 1;
    m11[0] = 1;

    FLINT_ASSERT(n_len > 0);
    FLINT_ASSERT(d_len > 0);
    FLINT_ASSERT(n_ptr[n_len - 1] != 0);
    n_lzcnt = flint_clz(n_ptr[n_len - 1]);

again:

    FLINT_ASSERT(Alen > 0 && A[Alen - 1] > 0);
    FLINT_ASSERT(Blen > 0 && B[Blen - 1] > 0);
    FLINT_ASSERT(Alen > Blen || (Alen == Blen && mpn_cmp(A, B, Blen) > 0));
    FLINT_ASSERT(Blen > n_len || (Blen == n_len && mpn_cmp(B, n_ptr, n_len) > 0));
    FLINT_ASSERT(m_len > 0 && m11[m_len - 1] > 0);
    FLINT_ASSERT(mpn_cmp(m11, m12, m_len) >= 0);

    if (Alen < 3 || Blen <= n_len)
    {
        /* too small or too close to the end */
        goto gauss;
    }

    a_lzcnt = flint_clz(A[Alen - 1]);

    if (Alen - 1 > Blen)
    {
        /* large quotient */
        goto gauss;
    }

    if (Alen - 1 == n_len && n_lzcnt < a_lzcnt)
    {
        /* too small or too close to the end */
        goto gauss;
    }

    FLINT_ASSERT(Alen == Blen || Alen - 1 == Blen);

    /* zero-extend B to length of A in the case Alen - 1 == Blen */
    B[Blen] = 0;

    A1 = MPN_LEFT_SHIFT_HI(A[Alen - 1], A[Alen - 2], a_lzcnt);
    A0 = MPN_LEFT_SHIFT_HI(A[Alen - 2], A[Alen - 3], a_lzcnt);
    B1 = MPN_LEFT_SHIFT_HI(B[Alen - 1], B[Alen - 2], a_lzcnt);
    B0 = MPN_LEFT_SHIFT_HI(B[Alen - 2], B[Alen - 3], a_lzcnt);

    written = _hgcd_uiui_no_write(A1, A0, B1, B0, h);
    if (written <= 0)
    {
        /* difficult quotient */
        goto gauss;
    }

    /* (Q, R) will be the new values for (A, B) */
    if (h->det == 1)
    {
        Qlen = flint_mpn_fmms1(Q, h->_22, A, h->_12, B, Alen);
        Rlen = flint_mpn_fmms1(R, h->_11, B, h->_21, A, Alen);
    }
    else
    {
        Qlen = flint_mpn_fmms1(Q, h->_12, B, h->_22, A, Alen);
        Rlen = flint_mpn_fmms1(R, h->_21, A, h->_11, B, Alen);
    }

    FLINT_ASSERT(Qlen >= Rlen && Rlen >= 0);
    FLINT_ASSERT(Qlen > Rlen || Rlen == 0 || mpn_cmp(Q, R, Rlen) > 0);

    if (Qlen < n_len || (Qlen == n_len && mpn_cmp(Q, n_ptr, n_len) <= 0))
    {
        /* overshot with too many quotients. rare (impossible?) due to above
           lzcnt restriction. can trigger by using n_lzcnt + 1 < a_lzcnt. */
        goto gauss;
    }

    /* copy (Q, R) to (A, B) */
    Alen = Qlen;
    Blen = Rlen;
    flint_mpn_copyi(A, Q, FMPQ_RECONSTRUCT_ARRAY_LIMIT);
    flint_mpn_copyi(B, R, FMPQ_RECONSTRUCT_ARRAY_LIMIT);

    /* multiply first row of m by h, use R for temp */
    mdet *= h->det;
    ex0 = mpn_mul_1(R, m11, m_len, h->_11);
    ex1 = mpn_addmul_1(R, m12, m_len, h->_21);
    ex2 = mpn_mul_1(m12, m12, m_len, h->_22);
    ex3 = mpn_addmul_1(m12, m11, m_len, h->_12);
    add_ssaaaa(m12[m_len + 1], m12[m_len], 0, ex2, 0, ex3);
    flint_mpn_copyi(m11, R, m_len);
    add_ssaaaa(m11[m_len + 1], m11[m_len], 0, ex0, 0, ex1);
    m_len += (m11[m_len + 1] != 0) ? 2 : (m11[m_len] != 0);

    /* so A > N. see if further A > N >= B */
    if (Blen < n_len || (Blen == n_len && mpn_cmp(B, n_ptr, n_len) <= 0))
    {
        /* got lucky. can happen */
        goto done;
    }

    goto again;

gauss:

    FLINT_ASSERT(Alen > 0 && A[Alen - 1] > 0);
    FLINT_ASSERT(Blen > 0 && B[Blen - 1] > 0);
    FLINT_ASSERT(Alen > Blen || (Alen == Blen && mpn_cmp(A, B, Blen) > 0));
    FLINT_ASSERT(Blen > n_len || (Blen == n_len && mpn_cmp(B, n_ptr, n_len) > 0));
    FLINT_ASSERT(m_len > 0 && m11[m_len - 1] > 0);
    FLINT_ASSERT(mpn_cmp(m11, m12, m_len) >= 0);

    /* (A, B) = (B, A mod B) */
    mpn_tdiv_qr(Q, R, 0, A, Alen, B, Blen);
    Qlen = Alen - Blen + 1;
    MPN_NORM(Q, Qlen);
    Rlen = Blen;
    MPN_NORM(R, Rlen);
    Alen = Blen;
    Blen = Rlen;
    flint_mpn_copyi(A, B, FMPQ_RECONSTRUCT_ARRAY_LIMIT);
    flint_mpn_copyi(B, R, FMPQ_RECONSTRUCT_ARRAY_LIMIT);

    /* (m11, m12) = (m12 + Q*m11, m11), use R for temp */
    mdet *= -1;
    ex0 = (Qlen > m_len) ? flint_mpn_mul(R, Q, Qlen, m11, m_len)
                         : flint_mpn_mul(R, m11, m_len, Q, Qlen);
    m_len = Qlen + m_len - (ex0 == 0);
    ex0 = mpn_add_n(R, R, m12, m_len);
    R[m_len] = ex0;
    m_len += ex0;
    flint_mpn_copyi(m12, m11, m_len);
    flint_mpn_copyi(m11, R, m_len);

    /* see if further A > N >= B */
    if (Blen > n_len || (Blen == n_len && mpn_cmp(B, n_ptr, n_len) > 0))
        goto again;

done:

    FLINT_ASSERT(Alen > n_len || (Alen == n_len && mpn_cmp(A, n_ptr, n_len) > 0));
    FLINT_ASSERT(n_len > Blen || (n_len == Blen && mpn_cmp(n_ptr, B, Blen) >= 0));

    if (d_len < m_len || (d_len == m_len && mpn_cmp(d_ptr, m11, m_len) < 0))
        return 0;

    fmpz_set_ui_array(d, m11, m_len);

    if (Blen > 0)
    {
        fmpz_set_ui_array(n, B, Blen);
        if (mdet < 0)
            fmpz_neg(n, n);
    }
    else
    {
        fmpz_zero(n);
    }

    FLINT_ASSERT(m_len > 0);
    if (Blen == 0)
        return m_len == 1 && m11[0] == 1;

    if (((m11[0] | B[0]) & 1) == 0)
        return 0; /* gcd is even */

    if (m_len >= Blen)
    {
        if (mpn_gcd(R, m11, m_len, B, Blen) != 1)
            return 0;
    }
    else
    {
        if (mpn_gcd(R, B, Blen, m11, m_len) != 1)
            return 0;
    }

    return R[0] == 1;
}


/*
    A, B, and N come in with A > B > N > 0. Work on A/B and accumulate any
    quotients into M. The second row of M is ignored.
    S and T are temp working space.

    return -1: no need to use _lehmer again; leave with A > B > N > 0
    return 0: good idea to use _lehmer again; leave with A > B > N > 0
    return 1: done; leave with A > N >= B
*/
static int _lehmer(_fmpz_mat22_t M, fmpz_t A, fmpz_t B, const fmpz_t N,
                                                            fmpz_t S, fmpz_t T)
{
    int ret;
    slong written;
    mp_srcptr n_ptr;
    mpz_ptr a, b, s, t;
    mp_ptr a_ptr, b_ptr, s_ptr, t_ptr;
    mp_size_t a_len, b_len, n_len, s_len, t_len;
    _ui_mat22_t h;
    mp_limb_t A0, A1, B0, B1;
    unsigned int n_lzcnt, a_lzcnt;

    if (!COEFF_IS_MPZ(*A) || !COEFF_IS_MPZ(*B))
    {
        /* don't come back */
        return -1;
    }

    a = COEFF_TO_PTR(*A);
    b = COEFF_TO_PTR(*B);

    if (COEFF_IS_MPZ(*N))
    {
        n_ptr = COEFF_TO_PTR(*N)->_mp_d;
        n_len = COEFF_TO_PTR(*N)->_mp_size;
    }
    else
    {
        n_ptr = (mp_srcptr) N; /* haha, dirty but works */
        n_len = 1;
    }

    FLINT_ASSERT(n_ptr[n_len - 1] != 0);
    n_lzcnt = flint_clz(n_ptr[n_len - 1]);

    if (a->_mp_size < 3 || b->_mp_size <= n_len)
    {
        /* don't come back */
        return -1;
    }

    s = _fmpz_promote(S);
    t = _fmpz_promote(T);

    /* fit everything to length a->_mp_size */
    a_len = a->_mp_size;
    FLINT_MPZ_REALLOC(b, a_len);
    FLINT_MPZ_REALLOC(s, a_len);
    FLINT_MPZ_REALLOC(t, a_len);

again:

    a_ptr = a->_mp_d;
    b_ptr = b->_mp_d;
    s_ptr = s->_mp_d;
    t_ptr = t->_mp_d;

    a_len = a->_mp_size;
    b_len = b->_mp_size;

    /* supposed a > b > n > 0 */
    FLINT_ASSERT(a_len >= b_len);
    FLINT_ASSERT(b_len >= n_len);
    FLINT_ASSERT(n_len > 0);

    FLINT_ASSERT(a_ptr[a_len - 1] != 0);
    FLINT_ASSERT(b_ptr[b_len - 1] != 0);
    FLINT_ASSERT(n_ptr[n_len - 1] != 0);

    if (a_len < 3 || b_len <= n_len)
    {
        /* too small or too close to the end */
        ret = -1;
        goto cleanup;
    }

    if (a_len - 1 > b_len)
    {
        /* large quotient */
        ret = 0;
        goto cleanup;
    }

    FLINT_ASSERT(a_ptr[a_len - 1] != 0);
    a_lzcnt = flint_clz(a_ptr[a_len - 1]);

    if (a_len - 1 == n_len && n_lzcnt < a_lzcnt)
    {
        /* too small or too close to the end */
        ret = -1;
        goto cleanup;
    }

    FLINT_ASSERT(a_len == b_len || a_len - 1 == b_len);

    if (a_len - 1 == b_len)
        b_ptr[a_len - 1] = 0;

    A1 = MPN_LEFT_SHIFT_HI(a_ptr[a_len - 1], a_ptr[a_len - 2], a_lzcnt);
    A0 = MPN_LEFT_SHIFT_HI(a_ptr[a_len - 2], a_ptr[a_len - 3], a_lzcnt);
    B1 = MPN_LEFT_SHIFT_HI(b_ptr[a_len - 1], b_ptr[a_len - 2], a_lzcnt);
    B0 = MPN_LEFT_SHIFT_HI(b_ptr[a_len - 2], b_ptr[a_len - 3], a_lzcnt);

    written = _hgcd_uiui_no_write(A1, A0, B1, B0, h);
    if (written <= 0)
    {
        /* difficult quotient */
        ret = 0;
        goto cleanup;
    }

    if (h->det == 1)
    {
        s_len = flint_mpn_fmms1(s_ptr, h->_22, a_ptr, h->_12, b_ptr, a_len);
        t_len = flint_mpn_fmms1(t_ptr, h->_11, b_ptr, h->_21, a_ptr, a_len);
    }
    else
    {
        FLINT_ASSERT(h->det == -1);
        s_len = flint_mpn_fmms1(s_ptr, h->_12, b_ptr, h->_22, a_ptr, a_len);
        t_len = flint_mpn_fmms1(t_ptr, h->_21, a_ptr, h->_11, b_ptr, a_len);
    }

    /* s > t >= 0 */
    FLINT_ASSERT(s_len >= t_len && t_len >= 0);
    FLINT_ASSERT(s_len > t_len || t_len == 0 || mpn_cmp(s_ptr, t_ptr, t_len) > 0);

    if (s_len < n_len || (s_len == n_len && mpn_cmp(s_ptr, n_ptr, n_len) <= 0))
    {
        /* overshot with too many quotients */
        ret = 0;
        goto cleanup;
    }

    /* multiply first row of M by h, using second row as temp space */
    fmpz_mul_ui(M->_21, M->_11, h->_11);
    fmpz_addmul_ui(M->_21, M->_12, h->_21);
    fmpz_mul_ui(M->_12, M->_12, h->_22);
    fmpz_addmul_ui(M->_12, M->_11, h->_12);
    fmpz_swap(M->_11, M->_21);
    M->det *= h->det;

    /* a = s; b = t */
    s->_mp_size = s_len;
    t->_mp_size = t_len;
    FLINT_SWAP(mpz_ptr, a, s);
    FLINT_SWAP(mpz_ptr, b, t);

    /* so a > n. see if further a > n >= b. */
    if (t_len < n_len || (t_len == n_len && mpn_cmp(t_ptr, n_ptr, n_len) <= 0))
    {
        /* lucky finish */
        ret = 1;
        goto cleanup;
    }

    goto again;

cleanup:

    /* a/b are valid; make s/t valid */
    s->_mp_size = 0;
    t->_mp_size = 0;

    *A = PTR_TO_COEFF(a);
    *B = PTR_TO_COEFF(b);
    *S = PTR_TO_COEFF(s);
    *T = PTR_TO_COEFF(t);

    _fmpz_demote_val(A);
    _fmpz_demote_val(B);
    _fmpz_demote_val(S);
    _fmpz_demote_val(T);

    return ret;
}

/*
    return 0: not done, output still satisfies A > B > N
    return 1: we are done with A > N >= B (rare)
*/
static int _split(_fmpz_mat22_t M, fmpz_t A, fmpz_t B, const fmpz_t N)
{
    int ret;
    _fmpq_cfrac_list_t v;
    _fmpz_mat22_t H;
    fmpz_t As, Bs, Q, R;
    slong a, b, n = fmpz_size(N);
    slong s;    /* shift amount in words */

    fmpz_init(As);
    fmpz_init(Bs);
    fmpz_init(Q);
    fmpz_init(R);
    _fmpz_mat22_init(H);
    _fmpq_cfrac_list_init(v);

again:

    FLINT_ASSERT(fmpz_cmp(A, B) > 0 && fmpz_cmp(B, N) > 0);

    a = fmpz_size(A);
    b = fmpz_size(B);

    if (b - n < FMPQ_RECONSTRUCT_HGCD_CUTOFF)
    {
        /* relatively few remaining quotients */
        ret = 0;
        goto cleanup;
    }

    s = 1 + FLINT_MAX(0, 2*n - a);

    if (s >= b)
    {
gauss:
        /* we hit a hard quotient */
        fmpz_fdiv_qr(Q, A, A, B);
        fmpz_addmul(M->_12, M->_11, Q); fmpz_swap(M->_11, M->_12); M->det *= -1;
        fmpz_swap(A, B);
        if (fmpz_cmp(B, N) > 0)
            goto again;
        ret = 1;
        goto cleanup;
    }

    fmpz_fdiv_q_2exp(As, A, FLINT_BITS*s);
    fmpz_fdiv_q_2exp(Bs, B, FLINT_BITS*s);
    if (fmpz_cmp(As, Bs) <= 0)
        goto gauss;

    _fmpz_mat22_one(H);
    v->length = 0;
    _fmpq_hgcd(v, H, As, Bs);
    if (_fmpz_mat22_is_one(H))
        goto gauss;

    fmpz_fdiv_r_2exp(Q, A, FLINT_BITS*s);
    fmpz_fdiv_r_2exp(R, B, FLINT_BITS*s);
    fmpz_mul_2exp(A, As, FLINT_BITS*s);
    fmpz_mul_2exp(B, Bs, FLINT_BITS*s);
    if (H->det == 1)
    {
        fmpz_addmul(A, Q, H->_22);
        fmpz_submul(A, R, H->_12);
        fmpz_addmul(B, R, H->_11);
        fmpz_submul(B, Q, H->_21);
    }
    else
    {
        fmpz_addmul(A, R, H->_12);
        fmpz_submul(A, Q, H->_22);
        fmpz_addmul(B, Q, H->_21);
        fmpz_submul(B, R, H->_11);
    }

    /* multiply first row of M by H, using second row as temp space */
    fmpz_mul(M->_21, M->_11, H->_11);
    fmpz_addmul(M->_21, M->_12, H->_21);
    fmpz_mul(M->_12, M->_12, H->_22);
    fmpz_addmul(M->_12, M->_11, H->_12);
    fmpz_swap(M->_11, M->_21);
    M->det *= H->det;

    while (fmpz_cmp(A, N) <= 0)
    {
        /* unlikely (impossible?) with above choice of s. pop a quotient */
        FLINT_ASSERT(v->length > 0);
        v->length--;
        fmpz_addmul(B, A, v->array + v->length);
        fmpz_swap(A, B);
        fmpz_submul(M->_11, M->_12, v->array + v->length);
        fmpz_swap(M->_11, M->_12);
        M->det *= -1;
    }

    /* should have used at least one quotient! */
    FLINT_ASSERT(v->length > 0);

    if (fmpz_cmp(B, N) > 0)
        goto again;

    ret = 1;

cleanup:

    fmpz_clear(As);
    fmpz_clear(Bs);
    fmpz_clear(Q);
    fmpz_clear(R);
    _fmpz_mat22_clear(H);
    _fmpq_cfrac_list_clear(v);

    return ret;
}

int _fmpq_reconstruct_fmpz_2(fmpz_t n, fmpz_t d,
                const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)
{
    int ret, success;
    mp_size_t Asize, Nsize;
    fmpz_t Q, R, A, B;
    _fmpz_mat22_t M; /* only need first row of matrix M */
#ifdef FLINT_WANT_ASSERT
    int cqt_success;
    fmpz_t cqt_n, cqt_d;

    /* check answer against correct naive version */
    fmpz_init(cqt_n);
    fmpz_init(cqt_d);
    cqt_success = _fmpq_reconstruct_fmpz_2_naive(cqt_n, cqt_d, a, m, N, D);
#endif

    /* Quickly identify small integers */
    if (fmpz_cmp(a, N) <= 0)
    {
        fmpz_set(n, a);
        fmpz_one(d);
        success = 1;
        goto cleanup_assert;
    }

    fmpz_sub(n, a, m);
    if (fmpz_cmpabs(n, N) <= 0)
    {
        fmpz_one(d);
        success = 1;
        goto cleanup_assert;
    }

    Asize = fmpz_size(m);
    Nsize = fmpz_size(N);

    /* it is better to avoid the fmpz overhead & allocation at small sizes */
    if (Asize <= FMPQ_RECONSTRUCT_ARRAY_LIMIT)
    {
        if (Asize < 2)
            success = _fmpq_reconstruct_fmpz_2_ui(n, d, a, m, N, D);
        else if (Asize == 2)
            success = _fmpq_reconstruct_fmpz_2_uiui(n, d, a, m, N, D);
        else
            success = _fmpq_reconstruct_fmpz_2_ui_array(n, d, a, m, N, D);
        goto cleanup_assert;
    }

    _fmpz_mat22_init(M);
    _fmpz_mat22_one(M);

    fmpz_init_set(A, m);
    fmpz_init_set(B, a);
    fmpz_init(Q);
    fmpz_init(R);

    /* We have A > B > N > 0; accumulate quotients into M until A > N >= B */
    FLINT_ASSERT(fmpz_cmp(A, B) > 0 && fmpz_cmp(B, N) > 0 && fmpz_sgn(N) > 0);

    if (Asize - Nsize < 3)
        goto gauss;

    if (Asize - Nsize < FMPQ_RECONSTRUCT_HGCD_CUTOFF)
        goto lehmer;

    if (_split(M, A, B, N))
        goto write_answer;

lehmer:

    FLINT_ASSERT(fmpz_cmp(A, B) > 0 && fmpz_cmp(B, N) > 0);

    ret = _lehmer(M, A, B, N, Q, R);
    if (ret < 0)
        goto gauss;
    else if (ret > 0)
        goto write_answer;

    fmpz_fdiv_qr(Q, A, A, B);
    fmpz_addmul(M->_12, M->_11, Q); fmpz_swap(M->_11, M->_12); M->det *= -1;
    fmpz_swap(A, B);

    if (fmpz_cmp(B, N) > 0)
        goto lehmer;

    goto write_answer;

gauss:

    FLINT_ASSERT(fmpz_cmp(A, B) > 0 && fmpz_cmp(B, N) > 0);

    fmpz_fdiv_qr(Q, A, A, B);
    fmpz_addmul(M->_12, M->_11, Q); fmpz_swap(M->_11, M->_12); M->det *= -1;
    fmpz_swap(A, B);

    if (fmpz_cmp(B, N) > 0)
        goto gauss;

write_answer:

    FLINT_ASSERT(fmpz_cmp(A, N) > 0 && fmpz_cmp(N, B) >= 0);

    fmpz_swap(n, B);
    fmpz_swap(d, M->_11);
    if (M->det != 1)
    {
        FLINT_ASSERT(M->det == -1);
        fmpz_neg(n, n);
    }

    success = 0;

    FLINT_ASSERT(fmpz_sgn(d) > 0);
    if (fmpz_cmp(d, D) <= 0)
    {
        fmpz_gcd(R, n, d);
        success = fmpz_is_one(R);
    }

    fmpz_clear(Q);
    fmpz_clear(R);
    fmpz_clear(A);
    fmpz_clear(B);
    _fmpz_mat22_clear(M);

cleanup_assert:

#ifdef FLINT_WANT_ASSERT
    FLINT_ASSERT(success == cqt_success);
    if (success)
    {
        FLINT_ASSERT(fmpz_equal(n, cqt_n));
        FLINT_ASSERT(fmpz_equal(d, cqt_d));
    }
    fmpz_clear(cqt_n);
    fmpz_clear(cqt_d);
#endif

    return success;
}


int fmpq_reconstruct_fmpz_2(fmpq_t res, const fmpz_t a, const fmpz_t m,
                                                const fmpz_t N, const fmpz_t D)
{
    return _fmpq_reconstruct_fmpz_2(fmpq_numref(res),
                                    fmpq_denref(res), a, m, N, D);
}

