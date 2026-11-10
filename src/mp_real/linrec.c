/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "mag.h"
#include "thread_support.h"
#include "mp_real.h"
#include "impl.h"

/* Matrix binary splitting for linear recurrences with running sums and
   a scalar-geometric forcing (see impl.h for the recurrence).

   EXACT PRODUCTS.  The product over [a, b) is an integer matrix
   T = [A B; SA SB] ((d + s) x (d + g)) over a common denominator q,
   with the forcing factor phi (also over q).  Two halves combine as

       [A B; SA SB] = [A2 B2 0; SA2 SB2 q2 I] [A1 B1; 0 phi1 I; SA1 SB1],
       q = q1 q2,  phi = phi1 phi2,

   which for large entries is done as this single product, so that the
   multiplications by phi1 and q2 share the transforms of the matrix
   product (the transformed-integer product of fft_small transforms every
   entry once and skips the zeros, blocking the computation when the
   transforms exceed flint_fft_small_max_transformed_ring_size); for
   small entries the zeros of the (typically sparse) matrices are
   skipped by a classical product and the scalar parts are separate
   multiplications.

   EVALUATION.  The recurrence is applied to a state: [a, b) is cut into
   chunks whose exact products have about half the working precision in
   bits, and each chunk's matrix is applied to the vector (X, f v), so
   that above the chunks there are matrix-vector rather than
   matrix-matrix products.  The state is kept scaled by the product of
   the denominators so far (a floating-point scale, divided out once at
   the end), so a chunk costs no division.  For large entries the
   application is one truncated transformed product
   (fmpz_mat_mul_fft_small_trunc: each state entry is transformed once
   for all the rows, and the rows are accumulated in the transformed
   domain), each state entry in fixed point with its own exponent so
   that it keeps its relative accuracy (the entries can differ in size
   by many bits, as for the convergents of a continued fraction); for
   small entries it is done entry by entry in ball arithmetic.  Measured
   at 10^6 digits against a tree of truncated fixed-point matrix
   products over chunks of 2 to 4 times the working precision (the
   earlier design), this is faster (11% for the dual terms of
   Khinchin's constant before the changes above) and takes half the
   memory: one chunk's products and their transforms, as there is no
   fixed-point matrix per level. */

/* ball <-> fixed point <-> mag *********************************************/

void
_mp_real_get_rad_mag(mag_t res, const mp_real_t x)
{
    if (x->err == 0)
        mag_zero(res);
    else
        mag_set_ui_2exp_si(res, x->err,
            FLINT_BITS * ((x->size == 0) ? x->exp : x->exp - x->size));
}

/* the top bits of the mantissa: |mid| in [m, m + 1) 2^-s B^(exp - 1)
   with m normalized; kept to 31 bits for mag */
#define LEAD_SHIFT (FLINT_BITS - 31)

void
_mp_real_get_mag(mag_t res, const mp_real_t x)
{
    mag_t r;

    if (x->size == 0)
    {
        _mp_real_get_rad_mag(res, x);
        return;
    }

    {
        int s;
        ulong m = _mp_real_lead(x, &s);
        mag_set_ui_2exp_si(res, (m >> LEAD_SHIFT) + 1,
            LEAD_SHIFT - s + FLINT_BITS * (x->exp - 1));
    }
    mag_init(r);
    _mp_real_get_rad_mag(r, x);
    mag_add(res, res, r);
    mag_clear(r);
}

void
_mp_real_get_mag_lower(mag_t res, const mp_real_t x)
{
    mag_t r;

    if (x->size == 0)
    {
        mag_zero(res);
        return;
    }

    {
        int s;
        ulong m = _mp_real_lead(x, &s);
        mag_set_ui_lower(res, m >> LEAD_SHIFT);
        mag_mul_2exp_si(res, res, LEAD_SHIFT - s + FLINT_BITS * (x->exp - 1));
    }
    mag_init(r);
    _mp_real_get_rad_mag(r, x);
    mag_sub_lower(res, res, r);
    mag_clear(r);
}

void
_mp_real_get_fmpz_fixed(fmpz_t r, mag_t err, const mp_real_t x, slong e)
{
    mag_t t;

    if (x->size == 0)
    {
        fmpz_zero(r);
    }
    else
    {
        slong s = FLINT_BITS * (x->exp - x->size) - e;

        fmpz_set_ui_array(r, x->d, x->size);
        if (x->negative)
            fmpz_neg(r, r);
        if (s >= 0)
            fmpz_mul_2exp(r, r, s);
        else
            fmpz_fdiv_q_2exp(r, r, -s);
    }

    mag_init(t);
    _mp_real_get_rad_mag(t, x);
    mag_add(err, err, t);
    mag_set_ui_2exp_si(t, 1, e);
    mag_add(err, err, t);
    mag_clear(t);
}

void
_mp_real_add_error_mag(mp_real_t x, const mag_t m)
{
    slong e, q;
    double v;

    if (mag_is_zero(m))
        return;

    if (mag_is_inf(m) || COEFF_IS_MPZ(MAG_EXP(m)))
        flint_throw(FLINT_ERROR, "_mp_real_add_error_mag: radius out of range\n");

    /* m = man 2^(exp - MAG_BITS) = (man 2^r) B^q with 0 <= r < FLINT_BITS;
       man 2^r has at most MAG_BITS + FLINT_BITS - 1 bits, exact as a double */
    e = MAG_EXP(m) - MAG_BITS;
    q = e >> MP_REAL_LGB;
    v = ldexp((double) MAG_MAN(m), (int) (e - q * FLINT_BITS));
    _mp_real_add_error_ulps_at(x, v, q);
}

/* Below this many bits per entry the classical product skipping zeros is
   used; above, the transformed-integer product of fft_small, called
   directly since fmpz_mat_mul does not choose it for small matrices
   (measured for 3 x 2 by 2 x 2: 15 times slower at 10^3 bits, as fast at
   1.6 10^4 bits, and 30% faster than separate products at 10^5 to 10^6
   bits) */
#define LINREC_FFT_CUTOFF 12000

/* ... and only for augmented products with at least this inner
   dimension: for the 3 x 3 by 3 x 2 products of the k-terms of
   Khinchin's constant (d = 2, s = 1), with their zeros and scalar
   blocks, separate products were faster at all sizes measured (up to
   5 10^5 bits: 20 to 25%), for the 8 x 11 by 11 x 7 products of the
   dual terms the shared transforms 20% faster from about 10^5 bits */
#define LINREC_FFT_MIN_INNER 6

/* the evaluation applies exact products of about LINREC_CHUNK_FRAC times
   the working precision (measured best among 1/4, 1/2, 1, 2) */
#define LINREC_CHUNK_FRAC 0.5

/* fork the halves of exact products on the caller's thread budget only
   when they are large enough to amortize a wake-up */
#define LINREC_FORK_BITS 1000000

/* a matrix whose entries alias other fmpz values, for read-only use
   (released with flint_free on the entries alone) */
static void
_shallow_init(fmpz_mat_t M, slong r, slong c)
{
    M->entries = flint_calloc(r * c, sizeof(fmpz));
    M->r = r;
    M->c = c;
    M->stride = c;
}

#define SH(M, i, j) ((M)->entries[(i) * (M)->stride + (j)])

/* classical product skipping zero entries */
static void
_sparse_mul(fmpz_mat_t P, const fmpz_mat_t L, const fmpz_mat_t R)
{
    slong i, j, k, n = fmpz_mat_ncols(L);

    for (i = 0; i < fmpz_mat_nrows(L); i++)
        for (j = 0; j < fmpz_mat_ncols(R); j++)
        {
            fmpz * c = fmpz_mat_entry(P, i, j);
            fmpz_zero(c);
            for (k = 0; k < n; k++)
            {
                const fmpz * x = fmpz_mat_entry(L, i, k);
                const fmpz * y = fmpz_mat_entry(R, k, j);
                if (!fmpz_is_zero(x) && !fmpz_is_zero(y))
                    fmpz_addmul(c, x, y);
            }
        }
}

/* exact products ***********************************************************/


typedef struct
{
    fmpz_mat_t T;
    fmpz_t phi, q;
}
_exact_node_struct;

static void
_exact_init(_exact_node_struct * N, const _mp_real_linrec_struct * rec)
{
    fmpz_mat_init(N->T, rec->d + rec->s, rec->d + rec->g);
    fmpz_init(N->phi);
    fmpz_init(N->q);
}

static void
_exact_clear(_exact_node_struct * N)
{
    fmpz_mat_clear(N->T);
    fmpz_clear(N->phi);
    fmpz_clear(N->q);
}

static void
_exact_leaf(_exact_node_struct * N, const _mp_real_linrec_struct * rec, slong l)
{
    fmpz_mat_zero(N->T);
    fmpz_zero(N->phi);
    rec->leaf(N->T, N->phi, N->q, l, rec->ctx);
}

/* N1 = N1 * N2 (N1 left, N2 right) */
static void
_exact_combine(_exact_node_struct * N1, _exact_node_struct * N2,
    const _mp_real_linrec_struct * rec)
{
    slong d = rec->d, s = rec->s, g = rec->g, i, j;
    fmpz_mat_t P;

    fmpz_mat_init(P, d + s, d + g);

    if (FLINT_MIN(FLINT_ABS(fmpz_mat_max_bits(N1->T)),
            FLINT_ABS(fmpz_mat_max_bits(N2->T))) >= LINREC_FFT_CUTOFF
        && d + g + s >= LINREC_FFT_MIN_INNER)
    {
        fmpz_mat_t L, R;

        _shallow_init(L, d + s, d + g + s);
        _shallow_init(R, d + g + s, d + g);
        for (i = 0; i < d + s; i++)
            for (j = 0; j < d + g; j++)
                SH(L, i, j) = *fmpz_mat_entry(N2->T, i, j);
        for (i = 0; i < s; i++)
            SH(L, d + i, d + g + i) = *N2->q;
        for (i = 0; i < d; i++)
            for (j = 0; j < d + g; j++)
                SH(R, i, j) = *fmpz_mat_entry(N1->T, i, j);
        for (i = 0; i < g; i++)
            SH(R, d + i, d + i) = *N1->phi;
        for (i = 0; i < s; i++)
            for (j = 0; j < d + g; j++)
                SH(R, d + g + i, j) = *fmpz_mat_entry(N1->T, d + i, j);
        if (!fmpz_mat_mul_fft_small(P, L, R))
            fmpz_mat_mul(P, L, R);
        flint_free(L->entries);
        flint_free(R->entries);
    }
    else
    {
        fmpz_mat_t L, R;

        fmpz_mat_window_init(L, N2->T, 0, 0, d + s, d);
        fmpz_mat_window_init(R, N1->T, 0, 0, d, d + g);
        _sparse_mul(P, L, R);
        fmpz_mat_window_clear(L);
        fmpz_mat_window_clear(R);
        if (g > 0 && !fmpz_is_zero(N1->phi))
            for (i = 0; i < d + s; i++)
                for (j = d; j < d + g; j++)
                    fmpz_addmul(fmpz_mat_entry(P, i, j),
                        fmpz_mat_entry(N2->T, i, j), N1->phi);
        for (i = d; i < d + s; i++)
            for (j = 0; j < d + g; j++)
                fmpz_addmul(fmpz_mat_entry(P, i, j),
                    fmpz_mat_entry(N1->T, i, j), N2->q);
    }

    fmpz_mat_swap(N1->T, P);
    fmpz_mat_clear(P);
    if (g > 0)
        fmpz_mul(N1->phi, N1->phi, N2->phi);
    fmpz_mul(N1->q, N1->q, N2->q);
}

static void
_exact_bsplit(_exact_node_struct * R, const _mp_real_linrec_struct * rec,
    slong a, slong b)
{
    if (rec->bottom != NULL && b - a <= rec->bottom_len)
    {
        fmpz_mat_zero(R->T);
        fmpz_zero(R->phi);
        rec->bottom(R->T, R->phi, R->q, a, b, rec->ctx);
    }
    else if (b - a == 1)
    {
        _exact_leaf(R, rec, a);
    }
    else
    {
        _exact_node_struct N2;
        slong m = a + (b - a) / 2;

        _exact_init(&N2, rec);
        _exact_bsplit(R, rec, a, m);
        _exact_bsplit(&N2, rec, m, b);
        _exact_combine(R, &N2, rec);
        _exact_clear(&N2);
    }
}

/* len lb, saturating (a size estimate; the product overflows on 32-bit
   machines at high precision) */
static slong
_size_bits(slong len, slong lb)
{
    if (lb > 0 && len > WORD_MAX / lb)
        return WORD_MAX;
    return len * lb;
}


/* the exact bsplit, forking the halves of large nodes on the caller's
   thread budget */
typedef struct
{
    _exact_node_struct * R;
    const _mp_real_linrec_struct * rec;
    slong a, b, lb;
}
_exact_args_struct;

static void
_exact_bsplit_par(void * arg)
{
    _exact_args_struct * A = (_exact_args_struct *) arg;
    slong a = A->a, b = A->b;

    if (_size_bits(b - a, A->lb) < LINREC_FORK_BITS || flint_get_num_threads() < 2
        || (A->rec->bottom != NULL && b - a <= A->rec->bottom_len))
    {
        _exact_bsplit(A->R, A->rec, a, b);
    }
    else
    {
        _exact_node_struct N2;
        _exact_args_struct A1 = *A, A2 = *A;
        slong m = a + (b - a) / 2;

        _exact_init(&N2, A->rec);
        A1.b = m;
        A2.a = m;
        A2.R = &N2;
        _mp_real_parallel_pair(_exact_bsplit_par, &A1, _exact_bsplit_par, &A2);
        _exact_combine(A->R, &N2, A->rec);
        _exact_clear(&N2);
    }
}

static slong
_leaf_bits(const _mp_real_linrec_struct * rec, slong l)
{
    _exact_node_struct N;
    slong lb;

    _exact_init(&N, rec);
    _exact_leaf(&N, rec, l);
    lb = FLINT_MAX((slong) fmpz_bits(N.q), FLINT_ABS(fmpz_mat_max_bits(N.T)));
    _exact_clear(&N);
    return FLINT_MAX(lb, 1);
}

/* Y = T [X; f v] for large entries, by truncated transformed products,
   one per column of the state: each entry z_k of [X; f v] (d + g rows)
   is written as an integer Z_k 2^(e_k) with W' = FLINT_BITS n + 32 bits
   (so that each keeps its relative accuracy, as in ball arithmetic: the
   entries may differ in size by many bits), the shifts 2^(e_k - e) are
   folded into the columns of T, and T' Z is computed by
   fmpz_mat_mul_fft_small_trunc above a limb lo leaving W' + 64 bits
   below the smallest output row's terms.  Entry i of the result carries
   sum_k |T_ik| err_k from the inputs (err_k the radius of z_k plus the
   floor, 2^(e_k)), plus two units of 2^(e + FLINT_BITS lo) for the
   truncation and the transform and two more as a margin (as for the
   fixed-point products of arb_mat_mul).  Every entry of the state and
   each forcing term is transformed once for all the rows, and the
   accumulation is in the transformed domain.  Returns 0 (then the
   caller computes all of Y entry by entry) if the entries are small,
   their sizes too spread out, or fft_small declines. */
static int
_apply_fft(mp_real_ptr Y, const _exact_node_struct * N,
    const _mp_real_linrec_struct * rec, mp_real_srcptr X, const mp_real_t f,
    mp_real_srcptr v, slong r, slong n)
{
    slong d = rec->d, s = rec->s, g = rec->g, K = d + g, i, j, k;
    slong Wf = FLINT_BITS * n + 32, tb, base, lo, emax;
    slong * ex;
    fmpz_mat_t Z, P, T2;
    mp_real_struct * fv = NULL;
    mag_struct * err;
    mag_t m, t;
    int ok = 1;

    tb = FLINT_ABS(fmpz_mat_max_bits(N->T));
    if (tb < LINREC_FFT_CUTOFF || Wf < LINREC_FFT_CUTOFF)
        return 0;

    if (g > 0)
    {
        fv = flint_malloc(sizeof(mp_real_struct) * g);
        for (i = 0; i < g; i++)
        {
            mp_real_init(fv + i);
            mp_real_mul(fv + i, f, v + i, n);
        }
    }

#define ZE(k, j) (((k) < d) ? X + (k) * r + (j) : fv + ((k) - d))

    ex = flint_malloc(sizeof(slong) * K);
    err = flint_malloc(sizeof(mag_struct) * K);
    for (k = 0; k < K; k++)
        mag_init(err + k);
    fmpz_mat_init(Z, K, 1);
    fmpz_mat_init(P, d + s, 1);
    fmpz_mat_init(T2, d + s, K);
    mag_init(m);
    mag_init(t);

    /* sizes spread by more than a few working precisions: the shifted
       T would be large, so entry by entry instead */
    for (j = 0; j < r && ok; j++)
    {
        base = WORD_MAX;
        emax = WORD_MIN;
        for (k = 0; k < K; k++)
            if (!mp_real_is_zero(ZE(k, j)))
            {
                slong e = mp_real_abs_bound_lt_2exp_si(ZE(k, j)) - Wf;
                base = FLINT_MIN(base, e);
                emax = FLINT_MAX(emax, e);
            }
        if (base != WORD_MAX && emax - base > 4 * Wf)
            ok = 0;
    }

    for (j = 0; j < r && ok; j++)
    {
        slong minbits = WORD_MAX;

        base = WORD_MAX;
        for (k = 0; k < K; k++)
        {
            fmpz_zero(fmpz_mat_entry(Z, k, 0));
            mag_zero(err + k);
            if (!mp_real_is_zero(ZE(k, j)))
            {
                ex[k] = mp_real_abs_bound_lt_2exp_si(ZE(k, j)) - Wf;
                base = FLINT_MIN(base, ex[k]);
                _mp_real_get_fmpz_fixed(fmpz_mat_entry(Z, k, 0), err + k, ZE(k, j), ex[k]);
            }
        }

        if (base == WORD_MAX)
        {
            for (i = 0; i < d + s; i++)
                mp_real_zero(Y + i * r + j);
            continue;
        }

        /* T' = T with column k shifted by e_k - base; the bits of the
           terms of each row, for the truncation */
        for (i = 0; i < d + s; i++)
        {
            slong rb = WORD_MIN;
            for (k = 0; k < K; k++)
            {
                fmpz * c = fmpz_mat_entry(T2, i, k);
                if (mp_real_is_zero(ZE(k, j)) || fmpz_is_zero(fmpz_mat_entry(N->T, i, k)))
                {
                    fmpz_zero(c);
                    continue;
                }
                fmpz_mul_2exp(c, fmpz_mat_entry(N->T, i, k), ex[k] - base);
                rb = FLINT_MAX(rb, (slong) fmpz_bits(c) + (slong) fmpz_bits(fmpz_mat_entry(Z, k, 0)));
            }
            if (rb != WORD_MIN)
                minbits = FLINT_MIN(minbits, rb);
        }
        if (minbits == WORD_MAX)
            minbits = 0;
        lo = (minbits - Wf - 64) / FLINT_BITS;
        lo = FLINT_MAX(lo, 0);

        if (!fmpz_mat_mul_fft_small_trunc(P, T2, Z, lo))
        {
            ok = 0;
            break;
        }

        for (i = 0; i < d + s; i++)
        {
            mp_real_struct * y = Y + i * r + j;

            /* sum_k |T_ik| err_k + 4 units */
            mag_set_ui_2exp_si(m, 1, base + FLINT_BITS * lo + 2);
            for (k = 0; k < K; k++)
            {
                if (fmpz_is_zero(fmpz_mat_entry(N->T, i, k)) || mag_is_zero(err + k))
                    continue;
                mag_set_ui_2exp_si(t, 1, fmpz_bits(fmpz_mat_entry(N->T, i, k)));
                mag_addmul(m, t, err + k);
            }
            mp_real_set_fmpz(y, fmpz_mat_entry(P, i, 0));
            mp_real_mul_2exp_si(y, y, base + FLINT_BITS * lo);
            _mp_real_add_error_mag(y, m);
        }
    }

    fmpz_mat_clear(Z);
    fmpz_mat_clear(P);
    fmpz_mat_clear(T2);
    for (k = 0; k < K; k++)
        mag_clear(err + k);
    flint_free(err);
    flint_free(ex);
    mag_clear(m);
    mag_clear(t);
    if (g > 0)
    {
        for (i = 0; i < g; i++)
            mp_real_clear(fv + i);
        flint_free(fv);
    }
#undef ZE
    return ok;
}

void
_mp_real_linrec_eval(mp_real_ptr X, mp_real_ptr S, mp_real_t f,
    const _mp_real_linrec_struct * rec, slong a, slong b, mp_real_srcptr v,
    slong r, slong n)
{
    slong d = rec->d, s = rec->s, g = rec->g, i, j, k, l;
    slong chunk_bits = (slong) (LINREC_CHUNK_FRAC * FLINT_BITS * n);
    _exact_node_struct N;
    mp_real_struct * Y, * Bv;
    mp_real_t Q, t, u;
    int scaled = 0;

    if (g > 0 && r != 1)
        flint_throw(FLINT_ERROR, "_mp_real_linrec_eval: forcing needs one column\n");

    _exact_init(&N, rec);
    Y = flint_malloc(sizeof(mp_real_struct) * (d + s) * r);
    Bv = flint_malloc(sizeof(mp_real_struct) * (d + s));
    for (i = 0; i < (d + s) * r; i++)
        mp_real_init(Y + i);
    for (i = 0; i < d + s; i++)
        mp_real_init(Bv + i);
    mp_real_init(Q);
    mp_real_init(t);
    mp_real_init(u);
    mp_real_set_ui(Q, 1);

    for (l = a; l < b; )
    {
        _exact_args_struct A;
        slong len = chunk_bits / _leaf_bits(rec, l);

        len = FLINT_MAX(len, 1);
        len = FLINT_MIN(len, b - l);

        A.R = &N;
        A.rec = rec;
        A.a = l;
        A.b = l + len;
        A.lb = _leaf_bits(rec, l + len - 1);
        _exact_bsplit_par(&A);

        if (_apply_fft(Y, &N, rec, X, f, v, r, n))
            goto update;

        /* otherwise entry by entry: the forcing columns collapsed with v,
           times f */
        if (g > 0)
        {
            for (i = 0; i < d + s; i++)
            {
                mp_real_zero(Bv + i);
                for (k = 0; k < g; k++)
                {
                    if (fmpz_is_zero(fmpz_mat_entry(N.T, i, d + k)))
                        continue;
                    mp_real_set_fmpz(t, fmpz_mat_entry(N.T, i, d + k));
                    mp_real_mul(u, t, v + k, n);
                    mp_real_add(Bv + i, Bv + i, u, n);
                }
                mp_real_mul(Bv + i, Bv + i, f, n);
            }
        }

        /* Y = T [X; f v] (rows of X, then the sum increments) */
        for (i = 0; i < d + s; i++)
            for (j = 0; j < r; j++)
            {
                mp_real_struct * y = Y + i * r + j;

                if (g > 0)
                    mp_real_set(y, Bv + i);
                else
                    mp_real_zero(y);
                for (k = 0; k < d; k++)
                {
                    if (fmpz_is_zero(fmpz_mat_entry(N.T, i, k)) || mp_real_is_zero(X + k * r + j))
                        continue;
                    mp_real_set_fmpz(t, fmpz_mat_entry(N.T, i, k));
                    mp_real_mul(u, t, X + k * r + j, n);
                    mp_real_add(y, y, u, n);
                }
            }

update:
        /* the state scaled by the product Q of the denominators so far:
           X Q -> T [X Q; f Q v], S Q -> q S Q + ..., f Q -> phi f Q */
        for (i = 0; i < d * r; i++)
            mp_real_swap(X + i, Y + i);
        if (!fmpz_is_one(N.q))
        {
            mp_real_set_fmpz(t, N.q);
            for (i = 0; i < s * r; i++)
                mp_real_mul(S + i, S + i, t, n);
            mp_real_mul(Q, Q, t, n);
            scaled = 1;
        }
        for (i = 0; i < s * r; i++)
            mp_real_add(S + i, S + i, Y + d * r + i, n);
        if (g > 0)
        {
            mp_real_set_fmpz(t, N.phi);
            mp_real_mul(f, f, t, n);
        }

        l += len;
    }

    if (scaled)
    {
        mp_real_set_ui(t, 1);
        mp_real_div(t, t, Q, n);
        for (i = 0; i < d * r; i++)
            mp_real_mul(X + i, X + i, t, n);
        for (i = 0; i < s * r; i++)
            mp_real_mul(S + i, S + i, t, n);
        if (g > 0)
            mp_real_mul(f, f, t, n);
    }

    for (i = 0; i < (d + s) * r; i++)
        mp_real_clear(Y + i);
    for (i = 0; i < d + s; i++)
        mp_real_clear(Bv + i);
    flint_free(Y);
    flint_free(Bv);
    mp_real_clear(Q);
    mp_real_clear(t);
    mp_real_clear(u);
    _exact_clear(&N);
}
