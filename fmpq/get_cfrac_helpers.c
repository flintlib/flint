/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "longlong.h"
#include "fmpq.h"
#include "fmpz_poly.h"

/***** gmp-impl.h ??? **************/

#define MPN_EXTRACT_NUMB(count, xh, xl)                \
  ((((xh) << ((count) - GMP_NAIL_BITS)) & GMP_NUMB_MASK) |    \
   ((xl) >> (GMP_LIMB_BITS - (count))))

/* Realloc for an mpz_t WHAT if it has less than NEEDED limbs.  */
#define MPZ_REALLOC(z,n) ((n) > ALLOC(z)            \
              ? (mp_ptr) _mpz_realloc(z,n)            \
              : PTR(z))
#define MPZ_NEWALLOC(z,n) ((n) > ALLOC(z)            \
               ? (mp_ptr) _mpz_newalloc(z,n)        \
               : PTR(z))

/* Field access macros.  */
#define SIZ(x) ((x)->_mp_size)
#define ABSIZ(x) ABS (SIZ (x))
#define PTR(x) ((x)->_mp_d)
#define EXP(x) ((x)->_mp_exp)
#define PREC(x) ((x)->_mp_prec)
#define ALLOC(x) ((x)->_mp_alloc)

#define MP_LIMB_T_SWAP(x, y)                        \
  do {                                    \
    mp_limb_t __mp_limb_t_swap__tmp = (x);                \
    (x) = (y);                                \
    (y) = __mp_limb_t_swap__tmp;                    \
  } while (0)
#define MP_SIZE_T_SWAP(x, y)                        \
  do {                                    \
    mp_size_t __mp_size_t_swap__tmp = (x);                \
    (x) = (y);                                \
    (y) = __mp_size_t_swap__tmp;                    \
  } while (0)


#define MP_SRCPTR_SWAP(x, y)                        \
  do {                                    \
    mp_srcptr __mp_srcptr_swap__tmp = (x);                \
    (x) = (y);                                \
    (y) = __mp_srcptr_swap__tmp;                    \
  } while (0)

#define MPN_PTR_SWAP(xp,xs, yp,ys)                    \
  do {                                    \
    MP_PTR_SWAP (xp, yp);                        \
    MP_SIZE_T_SWAP (xs, ys);                        \
  } while(0)
#define MPN_SRCPTR_SWAP(xp,xs, yp,ys)                    \
  do {                                    \
    MP_SRCPTR_SWAP (xp, yp);                        \
    MP_SIZE_T_SWAP (xs, ys);                        \
  } while(0)

#define MPZ_PTR_SWAP(x, y)                        \
  do {                                    \
    mpz_ptr __mpz_ptr_swap__tmp = (x);                    \
    (x) = (y);                                \
    (y) = __mpz_ptr_swap__tmp;                        \
  } while (0)
#define MPZ_SRCPTR_SWAP(x, y)                        \
  do {                                    \
    mpz_srcptr __mpz_srcptr_swap__tmp = (x);                \
    (x) = (y);                                \
    (y) = __mpz_srcptr_swap__tmp;                    \
  } while (0)

static void my_mpz_swap(mpz_ptr u, mpz_ptr v)
{
  MP_SIZE_T_SWAP(ALLOC(u), ALLOC(v));
  MP_SIZE_T_SWAP(SIZ(u), SIZ(v));
  MP_PTR_SWAP(PTR(v), PTR(u));
}


/* enable for debug printing of various types */
#if 0
static void _fmpz_mat22_print(const _fmpz_mat22_t M)
{
    printf("Matrix[\n   ");
    fmpz_print(M->_11); printf(",\n   ");
    fmpz_print(M->_12); printf(",\n   ");
    fmpz_print(M->_21); printf(",\n   ");
    fmpz_print(M->_22); printf("]\n");
}

static void _ui_mat22_print(const _ui_mat22_t M)
{
    flint_printf("Matrix[{{%wu, %wu}, {%wu, %wu}}]", M->_11, M->_12, M->_21, M->_22);
}

static void _fmpq_ball_print(const _fmpq_ball_t x)
{
    if (x->exact)
        printf("ExactBall[");
    else
        printf("Ball[");
    fmpz_print(x->left_num);
    flint_printf("/");
    fmpz_print(x->left_den);
    flint_printf(", ");
    fmpz_print(x->right_num);
    flint_printf("/");
    fmpz_print(x->right_den);
    flint_printf("]");
}
#endif


#if WANT_ASSERT

static int _fmpq_ball_is_ordered(const _fmpq_ball_t x)
{
    return x->exact || _fmpq_cmp(x->left_num, x->left_den,
                                 x->right_num, x->right_den) <= 0;
}

/* y = m^-1(x) */
static void _fmpq_ball_apply_mat22_inv(
    _fmpq_ball_t y,
    const _fmpz_mat22_t M,
    const _fmpq_ball_t x)
{
    y->exact = x->exact;

    if (M->det == 1)
    {
        if (!x->exact)
        {
            fmpz_mul(y->right_num, x->right_num, M->_22);
            fmpz_submul(y->right_num, x->right_den, M->_12);
            fmpz_mul(y->right_den, x->right_den, M->_11);
            fmpz_submul(y->right_den, x->right_num, M->_21);
        }

        fmpz_mul(y->left_num, x->left_num, M->_22);
        fmpz_submul(y->left_num, x->left_den, M->_12);
        fmpz_mul(y->left_den, x->left_den, M->_11);
        fmpz_submul(y->left_den, x->left_num, M->_21);
    }
    else
    {
        FLINT_ASSERT(M->det == -1);
        /* det = -1 swaps the endpoints */

        if (x->exact)
        {
            fmpz_mul(y->left_num, x->left_den, M->_12);
            fmpz_submul(y->left_num, x->left_num, M->_22);
            fmpz_mul(y->left_den, x->left_num, M->_21);
            fmpz_submul(y->left_den, x->left_den, M->_11);
        }
        else
        {
            fmpz_mul(y->right_num, x->left_den, M->_12);
            fmpz_submul(y->right_num, x->left_num, M->_22);
            fmpz_mul(y->right_den, x->left_num, M->_21);
            fmpz_submul(y->right_den, x->left_den, M->_11);

            fmpz_mul(y->left_num, x->right_den, M->_12);
            fmpz_submul(y->left_num, x->right_num, M->_22);
            fmpz_mul(y->left_den, x->right_num, M->_21);
            fmpz_submul(y->left_den, x->right_den, M->_11);
        }
    }
}

static int _fmpq_ball_equal(const _fmpq_ball_t x, _fmpq_ball_t y)
{
    if (x->exact != y->exact)
        return 0;

    if (!fmpz_equal(x->left_num, y->left_num) || !fmpz_equal(x->left_den, y->left_den))
        return 0;

    if (x->exact)
        return 1;

    if (!fmpz_equal(x->right_num, y->right_num) || !fmpz_equal(x->right_den, y->right_den))
        return 0;

    return 1;
}

#endif

/* y = [q 1; 1 0]^-1(x), also takes the left remainder */
static void _fmpq_ball_apply_mat22_inv_elem2(
    _fmpq_ball_t y,
    const fmpz_t q, fmpz_t r,
    const _fmpq_ball_t x)
{
    y->exact = x->exact;

    if (x->exact)
    {
        fmpz_set(y->left_num, x->left_den);
        fmpz_swap(y->left_den, r);
    }
    else
    {
        fmpz_set(y->right_num, x->left_den);
        fmpz_swap(y->right_den, r);

        fmpz_set(y->left_num, x->right_den);
        fmpz_set(y->left_den, x->right_num);
        fmpz_submul(y->left_den, x->right_den, q);
    }
}

/* is x canonical and bounded away from 1, i.e. 1 < x ? */
int _fmpq_ball_gt_one(const _fmpq_ball_t x)
{
    if (fmpz_sgn(x->left_num) <= 0)
        return 0;
    if (fmpz_sgn(x->left_den) <= 0)
        return 0;
    if (fmpz_cmp(x->left_den, x->left_num) >= 0)
        return 0;

    if (x->exact)
        return 1;

    if (fmpz_sgn(x->right_num) <= 0)
        return 0;
    if (fmpz_sgn(x->right_den) <= 0)
        return 0;
    if (fmpz_cmp(x->right_den, x->right_num) >= 0)
        return 0;

    return 1;
}

#define CFRAC_NEED_MATRIX 1
#define CFRAC_NEED_HGCD 2

/********************* hgcd *******************

    for a > b > 0,
        the hgcd returns terms for the open interval (a/(b+1), (a+1)/b)
        M is the product of the quotients {{q,1},{1,0}}

    if xa/xb = M(ya/yb), then

        det -1:
            M^-1 = {{m22, -m12}, {-m21, m11}}
              xa   xa+1       ya-m12  ya+m22
            (----, ----) = M((------, ------))
             xb+1   xb        yb+m11  yb-m21

        det -1:
            M^-1 = {{-m22, m12}, {m21, -m11}}
              xa   xa+1       ya-m22  ya+m12
            (----, ----) = M((------, ------))
             xb+1   xb        yb+m21  yb-m11

    we say that ya/yb is ok w.r.t. M if the above intervals are > 1, i.e.

        det +1:
            yb > m21, and
            ya - yb >= m11 + m12

        det -1:
            yb > m11, and
            ya - yb >= m21 + m22
*/
static int _hgcd_ok(const _fmpz_mat22_t M, const fmpz_t a, const fmpz_t b)
{
    int r = 1;
    fmpz_t t1, t2;
    fmpz_init(t1);
    fmpz_init(t2);
    FLINT_ASSERT(fmpz_sgn(M->_11) >= 0);
    FLINT_ASSERT(fmpz_sgn(M->_12) >= 0);
    FLINT_ASSERT(fmpz_sgn(M->_21) >= 0);
    FLINT_ASSERT(fmpz_sgn(M->_22) >= 0);
    r = r && (fmpz_cmp(a, b) > 0);
    r = r && (fmpz_sgn(b) > 0);
    if (M->det == 1)
    {
        r = r && (fmpz_cmp(a, M->_12) > 0);
        r = r && (fmpz_cmp(b, M->_21) > 0);
        fmpz_add(t2, M->_11, M->_12);
    }
    else
    {
        FLINT_ASSERT(M->det == -1);
        r = r && (fmpz_cmp(a, M->_22) > 0);
        r = r && (fmpz_cmp(b, M->_11) > 0);
        fmpz_add(t2, M->_21, M->_22);
    }
    fmpz_sub(t1, a, b);
    r = r && (fmpz_cmp(t1, t2) >= 0);
    fmpz_clear(t1);
    fmpz_clear(t2);
    return r;
}


static flint_bitcnt_t _hgcd_split(
    fmpz_t xa,
    fmpz_t xb,
    const fmpz_t ya,
    const fmpz_t yb,
    const _fmpz_mat22_t M,
    flint_bitcnt_t shift)
{
    flint_bitcnt_t r = 0;
    fmpz_t ta, tb;

    FLINT_ASSERT(_hgcd_ok(M, ya, yb));

    fmpz_init(ta);
    fmpz_init(tb);

    if (M->det == 1)
    {
        fmpz_sub(xa, ya, M->_12);
        fmpz_sub(xb, yb, M->_21);
        fmpz_add(ta, ya, M->_22);
        fmpz_add(tb, yb, M->_11);
    }
    else
    {
        FLINT_ASSERT(M->det == -1);
        fmpz_sub(xa, ya, M->_22);
        fmpz_sub(xb, yb, M->_11);
        fmpz_add(ta, ya, M->_12);
        fmpz_add(tb, yb, M->_21);
    }

    FLINT_ASSERT(fmpz_sgn(xa) > 0);
    FLINT_ASSERT(fmpz_sgn(xb) > 0);
    FLINT_ASSERT(fmpz_cmp(ta, xa) >= 0);
    FLINT_ASSERT(fmpz_cmp(tb, xb) >= 0);

    fmpz_fdiv_q_2exp(xa, xa, shift);
    fmpz_fdiv_q_2exp(ta, ta, shift);
    fmpz_fdiv_q_2exp(xb, xb, shift);
    fmpz_fdiv_q_2exp(tb, tb, shift);
    if (fmpz_sgn(xb) <= 0 || fmpz_cmp(xa, xb) <= 0)
        goto cleanup;

    while (!fmpz_equal(xa, ta) || !fmpz_equal(xb, tb))
    {
        shift++;
        fmpz_fdiv_q_2exp(xa, xa, 1);
        fmpz_fdiv_q_2exp(ta, ta, 1);
        fmpz_fdiv_q_2exp(xb, xb, 1);
        fmpz_fdiv_q_2exp(tb, tb, 1);
        if (fmpz_sgn(xb) <= 0 || fmpz_cmp(xa, xb) <= 0)
            goto cleanup;
    }

    r = shift;

cleanup:

    fmpz_clear(ta);
    fmpz_clear(tb);

    return r;
}

/* hgcd for two-limb input */
static slong _uiui_hgcd(
    mp_limb_t * s,
    mp_limb_t A1, mp_limb_t A0,
    mp_limb_t B1, mp_limb_t B0,
    _ui_mat22_t M)
{
    int i;
    slong written = 0;
    mp_limb_t d0, d1;
    mp_limb_t t0, t1, t2, r0, r1;
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
        FLINT_ASSERT(a1 != 0);
        FLINT_ASSERT(b1 != 0);
        FLINT_ASSERT(b1 <= a1);
        FLINT_ASSERT(b1 < a1 || (b1 == a1 && b0 < a0));
        q = 1;

        sub_ddmmss(r1,r0, a1,a0, b1,b0);

subtract:

        for (i = 2; i <= 4; i++)
        {
            sub_dddmmmsss(t2,t1,t0, 0,r1,r0, 0,b1,b0);
            if (t2 != 0)
                goto quotient_found;
            q += 1;
            r0 = t0;
            r1 = t1;
        }

        if (r1 != 0)
        {
            int ncnt, dcnt;
            mp_limb_t Q = 0;

            count_leading_zeros(ncnt, r1);
            count_leading_zeros(dcnt, b1);
            dcnt -= ncnt;
            if (dcnt <= 0)
                goto subtract;

            d1 = (b1 << dcnt) | (b0 >> (FLINT_BITS - dcnt));
            d0 = b0 << dcnt;

            while (dcnt >= 0)
            {
                sub_dddmmmsss(t2,t1,t0, 0,r1,r0, 0,d1,d0);
                Q = 2*Q + 1 + t2;
                r1 = t2 ? r1 : t1;
                r0 = t2 ? r0 : t0;
                d0 = (d1 << (FLINT_BITS - 1)) | (d0 >> 1);
                d1 = d1 >> 1;
                dcnt--;
            }

            q += Q;
        }

quotient_found:

        FLINT_ASSERT(r1 < b1 || (r1 == b1 && r0 < b0));

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

        s[written] = q;
        written++;
        FLINT_ASSERT(written <= 2*FLINT_BITS);
    }

    FLINT_ASSERT(a1 != 0);
    FLINT_ASSERT(b1 <= a1);
    FLINT_ASSERT(b1 < a1 || (b1 == a1 && b0 < a0));

    sub_ddmmss(d1,d0, a1,a0, b1,b0);
    if (det == 1)
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

#if WANT_ASSERT

    /* should be ok */

    sub_ddmmss(d1,d0, a1,a0, b1,b0);
    if (det == 1)
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

    q = s[written];

    t1 = m11 - q*m12;
    t2 = m21 - q*m22;
    m11 = m12;
    m21 = m22;
    m12 = t1;
    m22 = t2;
    det *= -1;

#if WANT_ASSERT
    umul_ppmm(t1,t0, a0, q); t1 += a1*q;
    add_ssaaaa(t1,t0, t1,t0, b1,b0);
    b0 = a0;
    b1 = a1;
    a0 = t0;
    a1 = t1;
#endif

    goto fixed;

}

/* y = a1*x1 - a2*x2 */
static mp_size_t _msub(mp_ptr y, mp_limb_t a1, mp_ptr x1,
                                 mp_limb_t a2, mp_ptr x2, mp_size_t n)
{
    mp_limb_t h0, h1;

    h0 =    mpn_mul_1(y, x1, n, a1);
    h1 = mpn_submul_1(y, x2, n, a2);

    if (h1 != h0)
        return -1;

    while (n > 0 && y[n - 1] == 0)
        n--;

    return n;
}


static void _lehmer_exact(
    _fmpz_vector_t s,
    _fmpz_mat22_t M,
    int flags,
    mpz_ptr xn,
    mpz_ptr xd,
    mpz_ptr yn,
    mpz_ptr yd)
{
    mp_limb_t s_temp[2*FLINT_BITS];
    slong written;
    mp_size_t xn_len;
    mp_size_t xd_len;
    mp_size_t yn_len;
    mp_size_t yd_len;
    mp_ptr xn_ptr;
    mp_ptr xd_ptr;
    mp_ptr yn_ptr;
    mp_ptr yd_ptr;
    _ui_mat22_t m;
    mp_limb_t A0, A1, B0, B1;
    mp_size_t n;

    /* fit everything to xn_len */
    n = SIZ(xn);

    MPZ_REALLOC(xd, n);
    MPZ_REALLOC(yn, n);
    MPZ_REALLOC(yd, n);

again:

    xn_len = SIZ(xn);
    xd_len = SIZ(xd);

    xn_ptr = PTR(xn);
    xd_ptr = PTR(xd);
    yn_ptr = PTR(yn);
    yd_ptr = PTR(yd);

    /* supposed xn > xd > 0 */
    FLINT_ASSERT(xn_len >= xd_len);
    FLINT_ASSERT(xd_len > 0);

    FLINT_ASSERT(xn_ptr[xn_len - 1] != 0);
    FLINT_ASSERT(xd_ptr[xd_len - 1] != 0);

    n = xn_len;
    if (n < 3)
        goto cleanup;

    if ((flags & CFRAC_NEED_HGCD) && xd_len <= 3 + _fmpz_mat22_bits(M)/FLINT_BITS)
    {
        goto cleanup;
    }

    if (n != xd_len && n != xd_len + 1)
        goto cleanup;

    if (n == xd_len + 1)
        xd_ptr[n - 1] = 0;

    if ((slong)(xn_ptr[n - 1]) < 0)
    {
        A1 = xn_ptr[n - 1];
        A0 = xn_ptr[n - 2];
        B1 = xd_ptr[n - 1];
        B0 = xd_ptr[n - 2];
    }
    else
    {
        int shift;
        count_leading_zeros(shift, xn_ptr[n - 1]);
        FLINT_ASSERT(shift > 0);
        A1 = MPN_EXTRACT_NUMB(shift, xn_ptr[n - 1], xn_ptr[n - 2]);
        A0 = MPN_EXTRACT_NUMB(shift, xn_ptr[n - 2], xn_ptr[n - 3]);
        B1 = MPN_EXTRACT_NUMB(shift, xd_ptr[n - 1], xd_ptr[n - 2]);
        B0 = MPN_EXTRACT_NUMB(shift, xd_ptr[n - 2], xd_ptr[n - 3]);
    }

    written = _uiui_hgcd(s_temp, A1, A0, B1, B0, m);
    if (written <= 0 || s->length + written > s->limit)
        goto cleanup;

    if (m->det == 1)
    {
        yn_len = _msub(yn_ptr, m->_22, xn_ptr, m->_12, xd_ptr, n);
        if (yn_len <= 0)
            goto cleanup;

        yd_len = _msub(yd_ptr, m->_11, xd_ptr, m->_21, xn_ptr, n);
        if (yd_len <= 0)
            goto cleanup;
    }
    else
    {
        yn_len = _msub(yn_ptr, m->_12, xd_ptr, m->_22, xn_ptr, n);
        if (yn_len <= 0)
            goto cleanup;

        yd_len = _msub(yd_ptr, m->_21, xn_ptr, m->_11, xd_ptr, n);
        if (yd_len <= 0)
            goto cleanup;
    }

    if (flags & CFRAC_NEED_HGCD)
    {
        /* over-strict but fast _hcgd_ok(M, yn, yd) */
        mp_size_t j;
        FLINT_ASSERT(yn_len >= yd_len);
        _fmpz_mat22_rmul_ui(M, m);
        for (j = 2 + _fmpz_mat22_bits(M)/FLINT_BITS; j < yn_len; j++)
        {
            mp_limb_t aa = yn_ptr[j];
            mp_limb_t bb = j < yd_len ? yd_ptr[j] : 0;
            if (aa > bb && aa - bb > 1)
                goto its_ok;
        }
        _fmpz_mat22_rmul_inv_ui(M, m);
        goto cleanup;
    }

    if (flags & CFRAC_NEED_MATRIX)
        _fmpz_mat22_rmul_ui(M, m);

its_ok:

    SIZ(yn) = yn_len;
    SIZ(yd) = yd_len;

    _fmpz_vector_append_ui(s, s_temp, written);

    my_mpz_swap(xn, yn);
    my_mpz_swap(xd, yd);
    goto again;

cleanup:

    /* never wrong */
    SIZ(yn) = 0;
    SIZ(yd) = 0;

    return;
}


static void _lehmer_inexact(
    _fmpz_vector_t s,
    _fmpz_mat22_t M,
    int needM,
    mpz_ptr xln,
    mpz_ptr xld,
    mpz_ptr xrn,
    mpz_ptr xrd,
    mpz_ptr yln,
    mpz_ptr yld,
    mpz_ptr yrn,
    mpz_ptr yrd)
{
    mp_limb_t s_temp[2*FLINT_BITS];
    slong written;
    mp_size_t xln_len;
    mp_size_t xld_len;
    mp_size_t xrn_len;
    mp_size_t xrd_len;
    mp_size_t yln_len;
    mp_size_t yld_len;
    mp_size_t yrn_len;
    mp_size_t yrd_len;
    mp_ptr xln_ptr;
    mp_ptr xld_ptr;
    mp_ptr xrn_ptr;
    mp_ptr xrd_ptr;
    mp_ptr yln_ptr;
    mp_ptr yld_ptr;
    mp_ptr yrn_ptr;
    mp_ptr yrd_ptr;
    _ui_mat22_t m;
    mp_limb_t A0, A1, B0, B1;
    mp_size_t n, nl, nr;

    /* fit everything to xln_len */
    nl = SIZ(xln);
    nr = SIZ(xrn);
    n = FLINT_MAX(nl, nr);

    MPZ_REALLOC(xln, n);
    MPZ_REALLOC(xld, n);
    MPZ_REALLOC(yln, n);
    MPZ_REALLOC(yld, n);
    MPZ_REALLOC(xrn, n);
    MPZ_REALLOC(xrd, n);
    MPZ_REALLOC(yrn, n);
    MPZ_REALLOC(yrd, n);

again:

    xln_len = SIZ(xln);
    xld_len = SIZ(xld);
    xrn_len = SIZ(xrn);
    xrd_len = SIZ(xrd);

    xln_ptr = PTR(xln);
    xld_ptr = PTR(xld);
    xrn_ptr = PTR(xrn);
    xrd_ptr = PTR(xrd);
    yln_ptr = PTR(yln);
    yld_ptr = PTR(yld);
    yrn_ptr = PTR(yrn);
    yrd_ptr = PTR(yrd);

    /* supposed xln > xld > 0 */
    FLINT_ASSERT(xln_len >= xld_len);
    FLINT_ASSERT(xld_len > 0);

    /* supposed xrn > xrd > 0 */
    FLINT_ASSERT(xrn_len >= xrd_len);
    FLINT_ASSERT(xrd_len > 0);

    FLINT_ASSERT(xln_ptr[xln_len-1] != 0);
    FLINT_ASSERT(xld_ptr[xld_len-1] != 0);
    FLINT_ASSERT(xrn_ptr[xrn_len-1] != 0);
    FLINT_ASSERT(xrd_ptr[xrd_len-1] != 0);

    nl = xln_len;
    nr = xrn_len;

    if (nl < 3 || nr < 3)
        goto cleanup;


    if (nl != xld_len && nl != xld_len + 1)
        goto cleanup;

    if (nl == xld_len + 1)
        xld_ptr[nl - 1] = 0;


    if (nr != xrd_len && nr != xrd_len + 1)
        goto cleanup;

    if (nr == xrd_len + 1)
        xrd_ptr[nr - 1] = 0;


    if ((slong)(xln_ptr[nl-1]) < 0)
    {
        A1 = xln_ptr[nl - 1];
        A0 = xln_ptr[nl - 2];
        B1 = xld_ptr[nl - 1];
        B0 = xld_ptr[nl - 2];
    }
    else
    {
        int shift;
        count_leading_zeros(shift, xln_ptr[nl - 1]);
        A1 = MPN_EXTRACT_NUMB(shift, xln_ptr[nl - 1], xln_ptr[nl - 2]);
        A0 = MPN_EXTRACT_NUMB(shift, xln_ptr[nl - 2], xln_ptr[nl - 3]);
        B1 = MPN_EXTRACT_NUMB(shift, xld_ptr[nl - 1], xld_ptr[nl - 2]);
        B0 = MPN_EXTRACT_NUMB(shift, xld_ptr[nl - 2], xld_ptr[nl - 3]);
    }

    written = _uiui_hgcd(s_temp, A1, A0, B1, B0, m);
    if (written <= 0 || s->length + written > s->limit)
        goto cleanup;

    if (m->det == 1)
    {
        yln_len = _msub(yln_ptr, m->_22, xln_ptr, m->_12, xld_ptr, nl);
        if (yln_len <= 0)
            goto cleanup;

        yld_len = _msub(yld_ptr, m->_11, xld_ptr, m->_21, xln_ptr, nl);
        if (yld_len <= 0)
            goto cleanup;

        yrn_len = _msub(yrn_ptr, m->_22, xrn_ptr, m->_12, xrd_ptr, nr);
        if (yrn_len <= 0)
            goto cleanup;

        yrd_len = _msub(yrd_ptr, m->_11, xrd_ptr, m->_21, xrn_ptr, nr);
        if (yrd_len <= 0)
            goto cleanup;
    }
    else
    {
        yrn_len = _msub(yrn_ptr, m->_12, xld_ptr, m->_22, xln_ptr, nl);
        if (yrn_len <= 0)
            goto cleanup;

        yrd_len = _msub(yrd_ptr, m->_21, xln_ptr, m->_11, xld_ptr, nl);
        if (yrd_len <= 0)
            goto cleanup;

        yln_len = _msub(yln_ptr, m->_12, xrd_ptr, m->_22, xrn_ptr, nr);
        if (yln_len <= 0)
            goto cleanup;

        yld_len = _msub(yld_ptr, m->_21, xrn_ptr, m->_11, xrd_ptr, nr);
        if (yld_len <= 0)
            goto cleanup;
    }

    /* check yl > 1 */
    if (yln_len <= yld_len)
    {
        if (yln_len != yld_len)
            goto cleanup;

        if (mpn_cmp(yln_ptr, yld_ptr, yln_len) <= 0)
            goto cleanup;
    }

    SIZ(yln) = yln_len;
    SIZ(yld) = yld_len;
    SIZ(yrn) = yrn_len;
    SIZ(yrd) = yrd_len;

    if (needM)
        _fmpz_mat22_rmul_ui(M, m);

    _fmpz_vector_append_ui(s, s_temp, written);

    my_mpz_swap(xln, yln);
    my_mpz_swap(xld, yld);
    my_mpz_swap(xrn, yrn);
    my_mpz_swap(xrd, yrd);
    goto again;

cleanup:

    /* never wrong */
    SIZ(yln) = 0;
    SIZ(yld) = 0;
    SIZ(yrn) = 0;
    SIZ(yrd) = 0;

    return;
}


static void _hgcd_step(
    _fmpz_mat22_t M,
    fmpz_t xa,
    fmpz_t xb,
    flint_bitcnt_t shift,
    _fmpz_mat22_t N,
    fmpz_t ya,
    fmpz_t yb)
{
    fmpz_fdiv_r_2exp(xa, xa, shift);
    fmpz_fdiv_r_2exp(xb, xb, shift);
    if (M->det == 1)
    {
        fmpz_sub(xa, xa, M->_12);
        fmpz_sub(xb, xb, M->_21);
        fmpz_fdiv_r_2exp(xa, xa, shift);
        fmpz_fdiv_r_2exp(xb, xb, shift);
        fmpz_add(xa, xa, M->_12);
        fmpz_add(xb, xb, M->_21);
    }
    else
    {
        FLINT_ASSERT(M->det == -1);
        fmpz_sub(xa, xa, M->_22);
        fmpz_sub(xb, xb, M->_11);
        fmpz_fdiv_r_2exp(xa, xa, shift);
        fmpz_fdiv_r_2exp(xb, xb, shift);
        fmpz_add(xa, xa, M->_22);
        fmpz_add(xb, xb, M->_11);
    }

    fmpz_mul_2exp(ya, ya, shift);
    fmpz_mul_2exp(yb, yb, shift);
    if (N->det == 1)
    {
        fmpz_addmul(ya, N->_22, xa);
        fmpz_submul(ya, N->_12, xb);
        fmpz_addmul(yb, N->_11, xb);
        fmpz_submul(yb, N->_21, xa);
    }
    else
    {
        FLINT_ASSERT(N->det == -1);
        fmpz_addmul(ya, N->_12, xb);
        fmpz_submul(ya, N->_22, xa);
        fmpz_addmul(yb, N->_21, xa);
        fmpz_submul(yb, N->_11, xb);
    }
    fmpz_swap(xa, ya);
    fmpz_swap(xb, yb);
    _fmpz_mat22_rmul(M, N);
}


void _fmpq_hgcd(
    _fmpz_vector_t s,
    _fmpz_mat22_t M,
    fmpz_t xa,
    fmpz_t xb)
{
    flint_bitcnt_t k, km, shift;
    fmpz_t ya, yb;
    _fmpz_mat22_t N;
#if WANT_ASSERT
    fmpz_t xa_org, xb_org;
    fmpz_init_set(xa_org, xa);
    fmpz_init_set(xb_org, xb);
#endif

    fmpz_init(ya);
    fmpz_init(yb);
    _fmpz_mat22_init(N);

    _fmpz_mat22_one(M);

again:

    FLINT_ASSERT(_hgcd_ok(M, xa, xb));

    if (s->length >= s->limit)
        goto cleanup;

    k = fmpz_bits(xa);
    km = _fmpz_mat22_bits(M);
    FLINT_ASSERT(k >= km);
    k -= km;

    if (k > 500*FLINT_BITS)
        goto split;
    if (k > 4*FLINT_BITS)
        goto lehmer;

gauss:

    FLINT_ASSERT(s->length <= s->limit);
    FLINT_ASSERT(_hgcd_ok(M, xa, xb));

    if (s->length >= s->limit)
        goto cleanup;

    fmpz_fdiv_qr(ya, yb, xa, xb);
    _fmpz_mat22_rmul_elem(M, ya);
    if (!_hgcd_ok(M, xb, yb))
    {
        _fmpz_mat22_rmul_inv_elem(M, ya);
        goto cleanup;
    }

    fmpz_swap(xa, xb);
    fmpz_swap(xb, yb);

    _fmpz_vector_push_back(s, ya);
    goto again;

lehmer:

    FLINT_ASSERT(s->length < s->limit);
    FLINT_ASSERT(_hgcd_ok(M, xa, xb));

    if (COEFF_IS_MPZ(*xa) && COEFF_IS_MPZ(*xb))
    {
        _lehmer_exact(s, M, CFRAC_NEED_MATRIX | CFRAC_NEED_HGCD,
                                         COEFF_TO_PTR(*xa), COEFF_TO_PTR(*xb),
                                         _fmpz_promote(ya), _fmpz_promote(yb));
        _fmpz_demote_val(ya);
        _fmpz_demote_val(yb);
        _fmpz_demote_val(xa);
        _fmpz_demote_val(xb);
    }

    goto gauss;

split:

    shift = _hgcd_split(ya, yb, xa, xb, M, km + k/2);
    if (shift == 0)
        goto gauss;

    _fmpq_hgcd(s, N, ya, yb);
    if (_fmpz_mat22_is_one(N))
        goto gauss;

    _hgcd_step(M, xa, xb, shift, N, ya, yb);
    FLINT_ASSERT(_hgcd_ok(M, xa, xb));

    km = _fmpz_mat22_bits(M);
    shift = _hgcd_split(ya, yb, xa, xb, M, km + 1);
    if (shift == 0)
        goto gauss;

    _fmpq_hgcd(s, N, ya, yb);
    if (_fmpz_mat22_is_one(N))
        goto gauss;

    _hgcd_step(M, xa, xb, shift, N, ya, yb);
    FLINT_ASSERT(_hgcd_ok(M, xa, xb));

    goto again;

cleanup:

#if WANT_ASSERT
    FLINT_ASSERT(_hgcd_ok(M, xa, xb));
    fmpz_mul(ya, M->_11, xa);
    fmpz_addmul(ya, M->_12, xb);
    fmpz_mul(yb, M->_21, xa);
    fmpz_addmul(yb, M->_22, xb);
    FLINT_ASSERT(fmpz_equal(xa_org, ya));
    FLINT_ASSERT(fmpz_equal(xb_org, yb));
    fmpz_clear(xa_org);
    fmpz_clear(xb_org);
#endif

    fmpz_clear(ya);
    fmpz_clear(yb);
    _fmpz_mat22_clear(N);

    FLINT_ASSERT(s->length <= s->limit);
    FLINT_ASSERT(_hgcd_ok(M, xa, xb));
    FLINT_ASSERT(M->det == 1 || M->det == -1);

    return;
}

/* given a >= b > 0, return the smallest k with floor(a/2^k) = floor(b/2^k) */
static flint_bitcnt_t _fmpz_tail_bits(const fmpz_t a, const fmpz_t b)
{
    flint_bitcnt_t k, j, max;
    max = k = fmpz_bits(a);
    FLINT_ASSERT(max >= fmpz_bits(b));

    for (j = 0; j < max; j++)
        if (fmpz_tstbit(a, j) != fmpz_tstbit(b, j))
            k = j + 1;

    return k;
}

void _fmpq_ball_get_cfrac(
    _fmpz_vector_t s,
    _fmpz_mat22_t M,
    int needM,
    _fmpq_ball_t x)
{
    flint_bitcnt_t k;
    fmpz_t q, r;
    _fmpq_ball_t y;
    _fmpz_mat22_t N;
#if WANT_ASSERT
    _fmpq_ball_t xorg;
    _fmpq_ball_init(xorg);
    xorg->exact = x->exact;
    fmpz_set(xorg->left_num, x->left_num);
    fmpz_set(xorg->left_den, x->left_den);
    fmpz_set(xorg->right_num, x->right_num);
    fmpz_set(xorg->right_den, x->right_den);
#endif

    fmpz_init(q);
    fmpz_init(r);
    _fmpq_ball_init(y);
    _fmpz_mat22_init(N);

    _fmpz_mat22_one(M);

    if (!x->exact)
    {
        if (fmpz_equal(x->left_num, x->right_num))
        {
            k = _fmpz_tail_bits(x->left_den, x->right_den);
            goto chop;
        }
        if (fmpz_equal(x->left_den, x->right_den))
        {
            k = _fmpz_tail_bits(x->right_num, x->left_num);
            goto chop;
        }
    }

again:

    FLINT_ASSERT(s->length <= s->limit);
    FLINT_ASSERT(_fmpq_ball_is_ordered(x));
    FLINT_ASSERT(_fmpq_ball_gt_one(x));
    FLINT_ASSERT(M->det == 1 || M->det == -1);

    if (s->length >= s->limit)
        goto cleanup;

    k = fmpz_bits(x->left_num);
    
    if (k > 500*FLINT_BITS)
        goto split;
    if (k > 4*FLINT_BITS)
        goto lehmer;

gauss:

    FLINT_ASSERT(s->length <= s->limit);
    FLINT_ASSERT(_fmpq_ball_is_ordered(x));
    FLINT_ASSERT(_fmpq_ball_gt_one(x));
    FLINT_ASSERT(M->det == 1 || M->det == -1);

    if (s->length >= s->limit)
        goto cleanup;

    fmpz_fdiv_qr(q, r, x->left_num, x->left_den);
    _fmpq_ball_apply_mat22_inv_elem2(y, q, r, x);
    if (!_fmpq_ball_gt_one(y))
        goto cleanup;

    _fmpq_ball_swap(x, y);

    if (needM)
        _fmpz_mat22_rmul_elem(M, q);

    _fmpz_vector_push_back(s, q);
    goto again;

lehmer:

    FLINT_ASSERT(s->length < s->limit);
    FLINT_ASSERT(_fmpq_ball_is_ordered(x));
    FLINT_ASSERT(_fmpq_ball_gt_one(x));
    FLINT_ASSERT(M->det == 1 || M->det == -1);

    if (x->exact)
    {
        if (COEFF_IS_MPZ(*x->left_num) && COEFF_IS_MPZ(*x->left_den))
        {
            _fmpz_mat22_one(N);
            _lehmer_exact(s, N, needM ? CFRAC_NEED_MATRIX : 0,
                       COEFF_TO_PTR(*x->left_num), COEFF_TO_PTR(*x->left_den),
                       _fmpz_promote(y->left_num), _fmpz_promote(y->left_den));
            _fmpz_demote_val(x->left_num);
            _fmpz_demote_val(x->left_den);
            _fmpz_demote_val(y->left_num);
            _fmpz_demote_val(y->left_den);

            if (needM)
                _fmpz_mat22_rmul(M, N);
        }
    }
    else
    {
        if (COEFF_IS_MPZ(*x->left_num) && COEFF_IS_MPZ(*x->left_den) &&
            COEFF_IS_MPZ(*x->right_num) && COEFF_IS_MPZ(*x->right_den))
        {
            _fmpz_mat22_one(N);
            _lehmer_inexact(s, N, needM,
                       COEFF_TO_PTR(*x->left_num), COEFF_TO_PTR(*x->left_den),
                     COEFF_TO_PTR(*x->right_num), COEFF_TO_PTR(*x->right_den),
                       _fmpz_promote(y->left_num), _fmpz_promote(y->left_den),
                     _fmpz_promote(y->right_num), _fmpz_promote(y->right_den));
            _fmpz_demote_val(x->left_num);
            _fmpz_demote_val(x->left_den);
            _fmpz_demote_val(x->right_num);
            _fmpz_demote_val(x->right_den);
            _fmpz_demote_val(y->left_num);
            _fmpz_demote_val(y->left_den);
            _fmpz_demote_val(y->right_num);
            _fmpz_demote_val(y->right_den);

            if (needM)
                _fmpz_mat22_rmul(M, N);
        }
    }

    goto gauss;

split:

    FLINT_ASSERT(s->length < s->limit);
    FLINT_ASSERT(_fmpq_ball_is_ordered(x));
    FLINT_ASSERT(_fmpq_ball_gt_one(x));
    FLINT_ASSERT(M->det == 1 || M->det == -1);

    k = k/2;

    if (x->exact)
    {
        fmpz_fdiv_q_2exp(y->left_num, x->left_num, k);
        fmpz_fdiv_q_2exp(y->left_den, x->left_den, k);
        if (fmpz_sgn(y->left_den) <= 0 || fmpz_cmp(y->left_num, y->left_den) <= 0)
            goto gauss;

        _fmpq_hgcd(s, N, y->left_num, y->left_den);
        if (_fmpz_mat22_is_one(N))
            goto gauss;

        fmpz_fdiv_r_2exp(q, x->left_num, k);
        fmpz_fdiv_r_2exp(r, x->left_den, k);
        fmpz_mul_2exp(x->left_num, y->left_num, k);
        fmpz_mul_2exp(x->left_den, y->left_den, k);
        if (N->det == 1)
        {
            fmpz_addmul(x->left_num, N->_22, q);
            fmpz_submul(x->left_num, N->_12, r);
            fmpz_addmul(x->left_den, N->_11, r);
            fmpz_submul(x->left_den, N->_21, q);
        }
        else
        {
            FLINT_ASSERT(N->det == -1);
            fmpz_addmul(x->left_num, N->_12, r);
            fmpz_submul(x->left_num, N->_22, q);
            fmpz_addmul(x->left_den, N->_21, q);
            fmpz_submul(x->left_den, N->_11, r);
        }
    }
    else
    {
        fmpz_fdiv_q_2exp(y->left_num, x->left_num, k);
        fmpz_fdiv_q_2exp(y->left_den, x->left_den, k);
        fmpz_add_ui(y->left_den, y->left_den, 1);
        fmpz_fdiv_q_2exp(y->right_num, x->right_num, k);
        fmpz_fdiv_q_2exp(y->right_den, x->right_den, k);
        fmpz_add_ui(y->right_num, y->right_num, 1);
        y->exact = 0;
        if (!_fmpq_ball_gt_one(y))
            goto gauss;

        _fmpq_ball_get_cfrac(s, N, 1, y);
        if (_fmpz_mat22_is_one(N))
            goto gauss;

        if (N->det == 1)
        {
            fmpz_one(r);
            fmpz_mul_2exp(r, r, k);
            fmpz_fdiv_r_2exp(q, x->left_den, k);
            fmpz_sub(r, r, q);
            fmpz_fdiv_r_2exp(q, x->left_num, k);

            fmpz_mul_2exp(x->left_num, y->left_num, k);
            fmpz_addmul(x->left_num, q, N->_22);
            fmpz_addmul(x->left_num, r, N->_12);
            fmpz_mul_2exp(x->left_den, y->left_den, k);
            fmpz_submul(x->left_den, q, N->_21);
            fmpz_submul(x->left_den, r, N->_11);

            fmpz_one(q);
            fmpz_mul_2exp(q, q, k);
            fmpz_fdiv_r_2exp(r, x->right_num, k);
            fmpz_sub(q, q, r);
            fmpz_fdiv_r_2exp(r, x->right_den, k);

            fmpz_mul_2exp(x->right_num, y->right_num, k);
            fmpz_submul(x->right_num, q, N->_22);
            fmpz_submul(x->right_num, r, N->_12);
            fmpz_mul_2exp(x->right_den, y->right_den, k);
            fmpz_addmul(x->right_den, q, N->_21);
            fmpz_addmul(x->right_den, r, N->_11);
        }
        else
        {
            FLINT_ASSERT(N->det == -1);

            fmpz_one(r);
            fmpz_mul_2exp(r, r, k);
            fmpz_fdiv_r_2exp(q, x->left_den, k);
            fmpz_sub(r, r, q);
            fmpz_fdiv_r_2exp(q, x->left_num, k);

            fmpz_mul_2exp(x->left_num, y->right_num, k);
            fmpz_submul(x->left_num, q, N->_22);
            fmpz_submul(x->left_num, r, N->_12);
            fmpz_mul_2exp(x->left_den, y->right_den, k);
            fmpz_addmul(x->left_den, q, N->_21);
            fmpz_addmul(x->left_den, r, N->_11);

            fmpz_one(q);
            fmpz_mul_2exp(q, q, k);
            fmpz_fdiv_r_2exp(r, x->right_num, k);
            fmpz_sub(q, q, r);
            fmpz_fdiv_r_2exp(r, x->right_den, k);

            fmpz_mul_2exp(x->right_num, y->left_num, k);
            fmpz_addmul(x->right_num, q, N->_22);
            fmpz_addmul(x->right_num, r, N->_12);
            fmpz_mul_2exp(x->right_den, y->left_den, k);
            fmpz_submul(x->right_den, q, N->_21);
            fmpz_submul(x->right_den, r, N->_11);

            fmpz_swap(x->right_num, x->left_num);
            fmpz_swap(x->right_den, x->left_den);
        }
    }

    FLINT_ASSERT(_fmpq_ball_gt_one(x));

    if (!needM)
        goto again;

    _fmpz_mat22_rmul(M, N);
    _fmpq_ball_get_cfrac(s, N, 1, x);
    _fmpz_mat22_rmul(M, N);

    goto cleanup;

chop:

    FLINT_ASSERT(s->length <= s->limit);
    FLINT_ASSERT(_fmpq_ball_is_ordered(x));
    FLINT_ASSERT(_fmpq_ball_gt_one(x));
    FLINT_ASSERT(_fmpz_mat22_is_one(M));

    fmpz_fdiv_q_2exp(q, x->left_num, k);
    fmpz_fdiv_q_2exp(r, x->left_den, k);
    if (fmpz_sgn(r) <= 0 || fmpz_cmp(q, r) <= 0)
        goto again;

    _fmpq_hgcd(s, M, q, r);
    if (_fmpz_mat22_is_one(M))
        goto again;

    fmpz_fdiv_r_2exp(y->left_num, x->left_num, k);
    fmpz_fdiv_r_2exp(y->left_den, x->left_den, k);
    fmpz_fdiv_r_2exp(y->right_num, x->right_num, k);
    fmpz_fdiv_r_2exp(y->right_den, x->right_den, k);

    fmpz_mul_2exp(x->left_num, q, k);
    fmpz_mul_2exp(x->left_den, r, k);
    fmpz_mul_2exp(x->right_num, q, k);
    fmpz_mul_2exp(x->right_den, r, k);

    if (M->det == 1)
    {
        fmpz_addmul(x->left_num, M->_22, y->left_num);
        fmpz_submul(x->left_num, M->_12, y->left_den);
        fmpz_addmul(x->left_den, M->_11, y->left_den);
        fmpz_submul(x->left_den, M->_21, y->left_num);
        fmpz_addmul(x->right_num, M->_22, y->right_num);
        fmpz_submul(x->right_num, M->_12, y->right_den);
        fmpz_addmul(x->right_den, M->_11, y->right_den);
        fmpz_submul(x->right_den, M->_21, y->right_num);
    }
    else
    {
        FLINT_ASSERT(M->det == -1);
        fmpz_addmul(x->left_num, M->_12, y->right_den);
        fmpz_submul(x->left_num, M->_22, y->right_num);
        fmpz_addmul(x->left_den, M->_21, y->right_num);
        fmpz_submul(x->left_den, M->_11, y->right_den);
        fmpz_addmul(x->right_num, M->_12, y->left_den);
        fmpz_submul(x->right_num, M->_22, y->left_num);
        fmpz_addmul(x->right_den, M->_21, y->left_num);
        fmpz_submul(x->right_den, M->_11, y->left_den);
    }

    goto gauss;

cleanup:

#if WANT_ASSERT
    FLINT_ASSERT(!needM || (_fmpq_ball_apply_mat22_inv(y, M, xorg), 
                                                       _fmpq_ball_equal(y, x)));
    _fmpq_ball_clear(xorg);
#endif

    fmpz_clear(q);
    fmpz_clear(r);
    _fmpq_ball_clear(y);
    _fmpz_mat22_clear(N);

    FLINT_ASSERT(s->length <= s->limit);
    FLINT_ASSERT(_fmpq_ball_gt_one(x));
    FLINT_ASSERT(M->det == 1 || M->det == -1);

    return;
}

