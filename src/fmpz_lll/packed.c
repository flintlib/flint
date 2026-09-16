/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Packed representation of the basis matrix used inside the floating-point
    LLL variants.

    Every entry is stored as a fixed number m of limbs in two's complement,
    the entries of a row being contiguous (n * m limbs per row). Rows are
    addressed through a pointer array so that row rotations are O(d).

    All row operations are computed slot-wise modulo 2^(m*FLINT_BITS), with
    no carries between entries. This is exact because a conservative upper
    bound on the bit length of the entries is maintained for every row: an
    operation is only performed if the bound on its result fits in the slots;
    otherwise the exact bit lengths are recomputed and, if necessary, the
    representation is grown to m + 1 limbs per entry. Growth is rare since
    size reduction shrinks the entries; conversely the slots are shrunk when
    all rows comfortably fit in m - 1 limbs.

    Working modulo 2^(m*FLINT_BITS) has the additional advantage that a row
    operation with a multi-limb scalar x only requires x mod 2^(m*FLINT_BITS).

    Exact dot products of rows (used to compute Gram matrix entries in the
    multiprecision LLL) are computed with the fixed-size signed mpn dot
    product kernels.
*/

#include <math.h>
#include "double_extras.h"
#include "d_vec.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"

/* entries with this many bits (magnitude) fit in slots of m limbs */
#define SLOT_MAX_BITS(m) ((m) * FLINT_BITS - 2)

/* headroom in bits when choosing the number of limbs */
#define PACKED_HEADROOM 16

slong
fmpz_lll_packed_limbs(const fmpz_mat_t B)
{
    slong bits = fmpz_mat_max_bits(B);
    bits = FLINT_ABS(bits);
    return (bits + 2 + PACKED_HEADROOM + FLINT_BITS - 1) / FLINT_BITS;
}

void
fmpz_lll_packed_init(fmpz_lll_packed_t P, slong d, slong n, slong m)
{
    slong i;

    P->d = d;
    P->n = n;
    P->m = m;
    P->entries = flint_malloc(sizeof(ulong) * FLINT_MAX(d * n * m, 1));
    P->rows = flint_malloc(sizeof(nn_ptr) * FLINT_MAX(d, 1));
    P->bits = flint_malloc(sizeof(slong) * FLINT_MAX(d, 1));

    for (i = 0; i < d; i++)
    {
        P->rows[i] = P->entries + i * n * m;
        P->bits[i] = 0;
    }
}

void
fmpz_lll_packed_clear(fmpz_lll_packed_t P)
{
    flint_free(P->entries);
    flint_free(P->rows);
    flint_free(P->bits);
}

void
fmpz_lll_packed_set_fmpz_mat(fmpz_lll_packed_t P, const fmpz_mat_t B)
{
    slong i, j, m = P->m;

    for (i = 0; i < P->d; i++)
    {
        slong bits = 0;

        for (j = 0; j < P->n; j++)
        {
            const fmpz * c = fmpz_mat_entry(B, i, j);
            bits = FLINT_MAX(bits, (slong) fmpz_bits(c));
            fmpz_get_signed_ui_array(P->rows[i] + j * m, m, c);
        }

        P->bits[i] = bits;
    }
}

void
fmpz_lll_packed_get_fmpz_mat(fmpz_mat_t B, const fmpz_lll_packed_t P)
{
    slong i, j, m = P->m;

    for (i = 0; i < P->d; i++)
        for (j = 0; j < P->n; j++)
            fmpz_set_signed_ui_array(fmpz_mat_entry(B, i, j), P->rows[i] + j * m, m);
}

/* bit length of the magnitude of an m-limb two's complement value */
static slong
_twos_bits(nn_srcptr x, slong m, nn_ptr tmp)
{
    slong size = m;

    if ((slong) x[m - 1] < 0)
    {
        mpn_neg(tmp, x, m);
        x = tmp;
    }

    while (size > 0 && x[size - 1] == 0)
        size--;

    if (size == 0)
        return 0;

    return (size - 1) * FLINT_BITS + FLINT_BIT_COUNT(x[size - 1]);
}

/* recompute the exact bit bound of row i */
void
fmpz_lll_packed_tighten(fmpz_lll_packed_t P, slong i)
{
    slong j, m = P->m, bits = 0;
    nn_ptr tmp;
    TMP_INIT;

    TMP_START;
    tmp = TMP_ALLOC(sizeof(ulong) * m);

    for (j = 0; j < P->n; j++)
        bits = FLINT_MAX(bits, _twos_bits(P->rows[i] + j * m, m, tmp));

    TMP_END;
    P->bits[i] = bits;
}

/* change the number of limbs per entry to m2 (values must fit) */
static void
_packed_resize(fmpz_lll_packed_t P, slong m2)
{
    slong i, j, d = P->d, n = P->n, m = P->m, mm = FLINT_MIN(m, m2);
    nn_ptr entries2 = flint_malloc(sizeof(ulong) * FLINT_MAX(d * n * m2, 1));

    for (i = 0; i < d; i++)
    {
        nn_srcptr src = P->rows[i];
        nn_ptr dst = entries2 + i * n * m2;

        for (j = 0; j < n; j++)
        {
            slong k;
            ulong ext;

            flint_mpn_copyi(dst + j * m2, src + j * m, mm);
            ext = (ulong) (((slong) src[j * m + m - 1]) >> (FLINT_BITS - 1));

            for (k = mm; k < m2; k++)
                dst[j * m2 + k] = ext;
        }
    }

    flint_free(P->entries);
    P->entries = entries2;
    P->m = m2;

    for (i = 0; i < d; i++)
        P->rows[i] = entries2 + i * n * m2;
}

void
fmpz_lll_packed_grow(fmpz_lll_packed_t P)
{
    _packed_resize(P, P->m + 1);
}

void
fmpz_lll_packed_maybe_shrink(fmpz_lll_packed_t P)
{
    slong i, maxbits = 0;

    if (P->m <= 1)
        return;

    for (i = 0; i < P->d; i++)
        maxbits = FLINT_MAX(maxbits, P->bits[i]);

    /* hysteresis: require half the initial headroom */
    if (maxbits + 2 + PACKED_HEADROOM / 2 <= SLOT_MAX_BITS(P->m - 1))
        _packed_resize(P, P->m - 1);
}

/* move row i to position j, shifting the rows in between */
void
fmpz_lll_packed_move_row(fmpz_lll_packed_t P, slong i, slong j)
{
    slong k;
    nn_ptr tmp = P->rows[i];
    slong tmpbits = P->bits[i];

    if (i > j)
    {
        for (k = i; k > j; k--)
        {
            P->rows[k] = P->rows[k - 1];
            P->bits[k] = P->bits[k - 1];
        }
    }
    else
    {
        for (k = i; k < j; k++)
        {
            P->rows[k] = P->rows[k + 1];
            P->bits[k] = P->bits[k + 1];
        }
    }

    P->rows[j] = tmp;
    P->bits[j] = tmpbits;
}

/*
    Make sure that the slots can hold the result of an operation on rows
    i and j whose entries are bounded by `bound` bits, where the bound was
    computed from the (possibly pessimistic) row bounds as
    max(bits[i], bits[j] + extra_j) + 1. First tighten the bounds, then
    grow if still necessary. Returns the recomputed bound.
*/
static slong
_ensure(fmpz_lll_packed_t P, slong i, slong j, slong extra_j, slong bound)
{
    if (bound <= SLOT_MAX_BITS(P->m))
        return bound;

    fmpz_lll_packed_tighten(P, i);
    fmpz_lll_packed_tighten(P, j);
    bound = FLINT_MAX(P->bits[i], P->bits[j] + extra_j) + 1;

    while (bound > SLOT_MAX_BITS(P->m))
        fmpz_lll_packed_grow(P);

    return bound;
}

/*
    Slot-wise kernels a_j <- a_j - u*b_j and a_j <- a_j + u*b_j modulo
    2^(m*FLINT_BITS), and a_j <- a_j -/+ b_j. The add and sub variants are
    generated as separate functions so that the inner loops are branch-free
    regardless of compiler optimizations.
*/

#define DEFINE_AORSMUL_1(name, AORS, add_ssaaaa_or_sub_ddmmss, add_sssaaaaaa_or_sub_dddmmmsss, mpn_aorsmul_1) \
static void \
name(nn_ptr a, nn_srcptr b, slong n, slong m, ulong u) \
{ \
    slong j; \
 \
    if (m == 1) \
    { \
        for (j = 0; j < n; j++) \
            a[j] = a[j] AORS b[j] * u; \
    } \
    else if (m == 2) \
    { \
        for (j = 0; j < n; j++) \
        { \
            ulong p1, p0, t1, t0; \
            umul_ppmm(p1, p0, b[2 * j], u); \
            p1 += b[2 * j + 1] * u; \
            add_ssaaaa_or_sub_ddmmss(t1, t0, a[2 * j + 1], a[2 * j], p1, p0); \
            a[2 * j] = t0; \
            a[2 * j + 1] = t1; \
        } \
    } \
    else if (m == 3) \
    { \
        for (j = 0; j < n; j++) \
        { \
            ulong p1, p0, q1, q0, t2, t1, t0; \
            umul_ppmm(p1, p0, b[3 * j], u); \
            umul_ppmm(q1, q0, b[3 * j + 1], u); \
            q1 += b[3 * j + 2] * u; \
            add_ssaaaa(q1, q0, q1, q0, 0, p1); \
            add_sssaaaaaa_or_sub_dddmmmsss(t2, t1, t0, a[3 * j + 2], a[3 * j + 1], a[3 * j], q1, q0, p0); \
            a[3 * j] = t0; \
            a[3 * j + 1] = t1; \
            a[3 * j + 2] = t2; \
        } \
    } \
    else \
    { \
        for (j = 0; j < n; j++) \
            mpn_aorsmul_1(a + j * m, b + j * m, m, u); \
    } \
}

DEFINE_AORSMUL_1(_submul_1, -, sub_ddmmss, sub_dddmmmsss, mpn_submul_1)
DEFINE_AORSMUL_1(_addmul_1, +, add_ssaaaa, add_sssaaaaaa, mpn_addmul_1)

#define DEFINE_AORS(name, AORS, add_ssaaaa_or_sub_ddmmss, mpn_aors_n) \
static void \
name(nn_ptr a, nn_srcptr b, slong n, slong m) \
{ \
    slong j; \
 \
    if (m == 1) \
    { \
        for (j = 0; j < n; j++) \
            a[j] = a[j] AORS b[j]; \
    } \
    else if (m == 2) \
    { \
        for (j = 0; j < n; j++) \
        { \
            ulong t1, t0; \
            add_ssaaaa_or_sub_ddmmss(t1, t0, a[2 * j + 1], a[2 * j], b[2 * j + 1], b[2 * j]); \
            a[2 * j] = t0; \
            a[2 * j + 1] = t1; \
        } \
    } \
    else \
    { \
        for (j = 0; j < n; j++) \
            mpn_aors_n(a + j * m, a + j * m, b + j * m, m); \
    } \
}

DEFINE_AORS(_sub, -, sub_ddmmss, mpn_sub_n)
DEFINE_AORS(_add, +, add_ssaaaa, mpn_add_n)

void
fmpz_lll_packed_row_sub(fmpz_lll_packed_t P, slong i, slong j)
{
    slong bound = FLINT_MAX(P->bits[i], P->bits[j]) + 1;
    bound = _ensure(P, i, j, 0, bound);
    _sub(P->rows[i], P->rows[j], P->n, P->m);
    P->bits[i] = bound;
}

void
fmpz_lll_packed_row_add(fmpz_lll_packed_t P, slong i, slong j)
{
    slong bound = FLINT_MAX(P->bits[i], P->bits[j]) + 1;
    bound = _ensure(P, i, j, 0, bound);
    _add(P->rows[i], P->rows[j], P->n, P->m);
    P->bits[i] = bound;
}

void
fmpz_lll_packed_row_submul_si(fmpz_lll_packed_t P, slong i, slong j, slong x)
{
    ulong u;
    slong bx, bound;

    if (x == 0)
        return;

    u = (x < 0) ? -(ulong) x : (ulong) x;
    bx = FLINT_BIT_COUNT(u);
    bound = FLINT_MAX(P->bits[i], P->bits[j] + bx) + 1;
    bound = _ensure(P, i, j, bx, bound);
    if (x > 0)
        _submul_1(P->rows[i], P->rows[j], P->n, P->m, u);
    else
        _addmul_1(P->rows[i], P->rows[j], P->n, P->m, u);
    P->bits[i] = bound;
}

void
fmpz_lll_packed_row_submul_fmpz(fmpz_lll_packed_t P, slong i, slong j, const fmpz_t x)
{
    slong bx, bound, m, n, k;
    nn_ptr xl, a;
    nn_srcptr b;
    TMP_INIT;

    if (fmpz_is_zero(x))
        return;

    if (fmpz_fits_si(x))
    {
        fmpz_lll_packed_row_submul_si(P, i, j, fmpz_get_si(x));
        return;
    }

    bx = fmpz_bits(x);
    bound = FLINT_MAX(P->bits[i], P->bits[j] + bx) + 1;
    bound = _ensure(P, i, j, bx, bound);
    m = P->m;
    n = P->n;

    TMP_START;
    xl = TMP_ALLOC(sizeof(ulong) * m);
    fmpz_get_signed_ui_array(xl, m, x);   /* x mod 2^(m*FLINT_BITS) */
    a = P->rows[i];
    b = P->rows[j];

    /* m == 1 cannot occur here: |x| >= 2^(FLINT_BITS-1) forces at least two limbs */
    if (m == 2)
    {
        for (k = 0; k < n; k++)
        {
            ulong p1, p0, t1, t0;
            umul_ppmm(p1, p0, b[2 * k], xl[0]);
            p1 += b[2 * k] * xl[1] + b[2 * k + 1] * xl[0];
            sub_ddmmss(t1, t0, a[2 * k + 1], a[2 * k], p1, p0);
            a[2 * k] = t0;
            a[2 * k + 1] = t1;
        }
    }
    else
    {
        nn_ptr t = TMP_ALLOC(sizeof(ulong) * 2 * m);
        for (k = 0; k < n; k++)
        {
            flint_mpn_mul_n(t, b + k * m, xl, m);   /* low m limbs are the product mod 2^(m*FLINT_BITS) */
            mpn_sub_n(a + k * m, a + k * m, t, m);
        }
    }

    TMP_END;
    P->bits[i] = bound;
}

/* a_j <- a_j -/+ u * 2^e * b_j modulo 2^(m*FLINT_BITS) */
void
fmpz_lll_packed_row_submul_si_2exp(fmpz_lll_packed_t P, slong i, slong j, slong x, ulong e)
{
    ulong u;
    slong bx, bound, k, m, n, s, r;
    nn_ptr tmp;
    mp_limb_t (* mpn_aorsmul_1)(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);
    TMP_INIT;

    if (x == 0)
        return;

    if (e == 0)
    {
        fmpz_lll_packed_row_submul_si(P, i, j, x);
        return;
    }

    u = (x < 0) ? -(ulong) x : (ulong) x;
    bx = FLINT_BIT_COUNT(u) + e;
    bound = FLINT_MAX(P->bits[i], P->bits[j] + bx) + 1;
    bound = _ensure(P, i, j, bx, bound);
    m = P->m;
    n = P->n;
    s = e / FLINT_BITS;
    r = e % FLINT_BITS;

    mpn_aorsmul_1 = (x < 0) ? mpn_addmul_1 : mpn_submul_1;

    TMP_START;
    tmp = TMP_ALLOC(sizeof(ulong) * m);

    for (k = 0; k < n; k++)
    {
        nn_srcptr bk = P->rows[j] + k * m;

        flint_mpn_zero(tmp, m);
        if (s < m)
        {
            if (r == 0)
                flint_mpn_copyi(tmp + s, bk, m - s);
            else
                mpn_lshift(tmp + s, bk, m - s, r);
        }

        mpn_aorsmul_1(P->rows[i] + k * m, tmp, m, u);
    }

    TMP_END;
    P->bits[i] = bound;
}

/* same semantics as _fmpz_vec_get_d_vec_2exp on the unpacked row; also sets the exact row bound */
slong
fmpz_lll_packed_get_d_vec_2exp(double * appv, fmpz_lll_packed_t P, slong i)
{
    slong j, n = P->n, m = P->m, maxexp = 0;
    nn_srcptr row = P->rows[i];
    slong * exp;
    nn_ptr tmp;
    TMP_INIT;

    TMP_START;
    exp = TMP_ALLOC(n * sizeof(slong));
    tmp = TMP_ALLOC(m * sizeof(ulong));

    for (j = 0; j < n; j++)
    {
        nn_srcptr e = row + j * m;
        slong size = m;
        int neg = ((slong) e[m - 1]) < 0;

        if (neg)
        {
            mpn_neg(tmp, e, m);
            e = tmp;
        }

        while (size > 0 && e[size - 1] == 0)
            size--;

        if (size == 0)
        {
            appv[j] = 0.0;
            exp[j] = 0;
        }
        else
        {
            exp[j] = (size - 1) * FLINT_BITS + FLINT_BIT_COUNT(e[size - 1]);
            appv[j] = flint_mpn_get_d(e, size, neg ? -1 : 1, -exp[j]);
            if (exp[j] > maxexp)
                maxexp = exp[j];
        }
    }

    for (j = 0; j < n; j++)
        appv[j] = d_mul_2exp(appv[j], exp[j] - maxexp);

    P->bits[i] = maxexp;

    TMP_END;
    return maxexp;
}

/* exact dot product of the first len entries of rows i and j */
void
fmpz_lll_packed_dot(fmpz_t res, const fmpz_lll_packed_t P, slong i, slong j, slong len)
{
    slong k, m = P->m;

    if (len == 0)
    {
        fmpz_zero(res);
        return;
    }

    /* signed fixed-width kernels; the output size is chosen from the row bounds */
    if (m <= FLINT_MPN_DOT_TAB_N)
    {
        ulong acc[2 * FLINT_MPN_DOT_TAB_N + 1];
        slong need = P->bits[i] + P->bits[j] + FLINT_BIT_COUNT(len) + 1;   /* signed result bits */
        slong sz = FLINT_MAX((need + FLINT_BITS - 1) / FLINT_BITS, 2 * m - 1);

        if (sz <= 2 * m + 1)
        {
            slong idx = sz - (2 * m - 1);

            if (m <= FLINT_MPN_DOT_DEDICATED_TAB_N && flint_mpn_dot_tab[1][m][m][idx] != NULL)
            {
                flint_mpn_dot_tab[1][m][m][idx](acc, P->rows[i], P->rows[j], len);
                fmpz_set_signed_ui_array(res, acc, sz);
                return;
            }
            else if (flint_mpn_dot_strided_tab[1][m][m][idx] != NULL)
            {
                flint_mpn_dot_strided_tab[1][m][m][idx](acc, P->rows[i], m, P->rows[j], m, len);
                fmpz_set_signed_ui_array(res, acc, sz);
                return;
            }
        }
    }

    /* generic fallback */
    {
        fmpz_t ta, tb;

        fmpz_init(ta);
        fmpz_init(tb);
        fmpz_zero(res);

        for (k = 0; k < len; k++)
        {
            fmpz_set_signed_ui_array(ta, P->rows[i] + k * m, m);
            fmpz_set_signed_ui_array(tb, P->rows[j] + k * m, m);
            fmpz_addmul(res, ta, tb);
        }

        fmpz_clear(ta);
        fmpz_clear(tb);
    }
}

/* same semantics as fmpz_lll_heuristic_dot */
double
fmpz_lll_packed_heuristic_dot(const double * vec1, const double * vec2, slong len2,
                              const fmpz_lll_packed_t P, slong k, slong j, slong exp_adj)
{
    double sum = _d_vec_dot(vec1, vec2, len2);
    double tmp = _d_vec_norm(vec1, len2);
    double tmp2 = _d_vec_norm(vec2, len2);

    tmp = tmp * tmp2 * ldexp(1.0, -70);
    tmp2 = sum * sum;

    if (tmp2 <= tmp)
    {
        slong exp;
        fmpz_t sp;
        fmpz_init(sp);
        fmpz_lll_packed_dot(sp, P, k, j, len2);
        sum = fmpz_get_d_2exp(&exp, sp);
        sum = ldexp(sum, exp - exp_adj);
        fmpz_clear(sp);
    }

    return sum;
}
