/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_mat/impl.h"
#include "nmod_mat/delayed_impl.h"

/*
    Triangular solves by substitution with delayed reductions. With an r x r
    triangular matrix, each entry of the solution is b - sum_t c_t y_t, with b
    and the y_t reduced and at most r - 1 terms: the products are accumulated
    unreduced, by vectors (acc += (n - c) v with the routines
    _nmod_vec_nored_*scalar_addmul*), and each entry is reduced once. The
    accumulators have one, two or three limbs, from the bound (r - 1) (n - 1)^2
    + (n - 1); for n <= 2^(FLINT_BITS / 2) they have one limb, reduced every m
    products if this bound needs more (unless m = 1).
*/

typedef struct
{
    nmod_t mod;
    int nlimbs;
    int halflimb;   /* n <= 2^(FLINT_BITS / 2): products fit in one limb */
    int red;        /* variant of the reduction */
    slong m;        /* one limb: products between two reductions */
    ulong c0, c1;   /* precomputations of the reduction */
}
tri_red_struct;

static void
tri_red_init(tri_red_struct * R, slong r, nmod_t mod)
{
    ulong u2, u1, u0, t1, t0;

    R->mod = mod;
    R->halflimb = (mod.n <= (UWORD(1) << (FLINT_BITS / 2)));
    R->m = WORD_MAX;
    R->c0 = R->c1 = 0;

    if (R->halflimb)
    {
        ulong sq = (mod.n - 1) * (mod.n - 1);
        ulong m = (sq == 0) ? UWORD_MAX : (UWORD_MAX - (mod.n - 1)) / sq;

        /* with a reduction after every product (n close to 2^(FLINT_BITS
           / 2)), two limbs are faster unless r is small */
        if (m >= 2 || r < 8)
        {
            R->nlimbs = 1;
            if (m < (ulong) r - 1)
                R->m = m;
            else
                m = r - 1;

            /* Lemire's reduction for entries below 2^(FLINT_BITS / 2) */
            if (m * sq + (mod.n - 1) < (UWORD(1) << (FLINT_BITS / 2)))
            {
                R->red = 0;
                R->c0 = n_lemire_precomp(mod.n);
            }
            else
            {
                R->red = 1;
                R->c0 = n_barrett_precomp(mod.n);
            }
            return;
        }
    }

    /* (r - 1) (n - 1)^2 + (n - 1) */
    umul_ppmm(t1, t0, mod.n - 1, mod.n - 1);
    umul_ppmm(u1, u0, t0, r - 1);
    umul_ppmm(u2, t0, t1, r - 1);
    add_ssaaaa(u2, u1, u2, u1, 0, t0);
    add_sssaaaaaa(u2, u1, u0, u2, u1, u0, 0, 0, mod.n - 1);

    R->nlimbs = (u2 != 0) ? 3 : (u1 != 0) ? 2 : 1;

    if (R->nlimbs == 1)
    {
        R->red = 1;
        R->c0 = n_barrett_precomp(mod.n);
    }
    else if (R->nlimbs == 2)
    {
        if (u1 < mod.n)
            R->red = 0;
        else if (mod.norm > 0)
        {
            R->red = 1;
            n_ll_rem_l_precomp(&R->c1, &R->c0, mod.n);
        }
        else
            R->red = 2;
    }
    else
    {
        /* as in lu_classical_delayed.c: c1 = 2^FLINT_BITS,
           c0 = 2^(2 FLINT_BITS) mod n (not reduced when "limited") */
        if (mod.norm == 0 && mod.n < (UWORD(1) << (FLINT_BITS - 1))
                                     + (UWORD(1) << (FLINT_BITS / 2 - 2)))
        {
            R->red = 1;
            R->c1 = -mod.n;
            R->c0 = 4 * R->c1 * R->c1;
        }
        else
        {
            R->red = 0;
            R->c1 = nmod_set_ui(UWORD(1) << (FLINT_BITS - 1), mod);
            R->c1 = nmod_add(R->c1, R->c1, mod);
            R->c0 = nmod_mul(R->c1, R->c1, mod);
        }
    }
}

FLINT_FORCE_INLINE ulong
tri_red(nn_srcptr a, const tri_red_struct * R, const int nlimbs)
{
    ulong r;
    nmod_t mod = R->mod;

    if (nlimbs == 1)
    {
        if (R->red == 0)
            r = n_mod_lemire(a[0], mod.n, R->c0);
        else
            r = n_mod_barrett(a[0], mod.n, R->c0);
    }
    else if (nlimbs == 2)
    {
        if (R->red == 0)
            NMOD_RED2(r, a[1], a[0], mod);
        else if (R->red == 1)
            r = n_ll_rem_l_nonfullword(a[1], a[0], mod.n, R->c1, R->c0);
        else
            NMOD2_RED2(r, a[1], a[0], mod);
    }
    else
    {
        if (R->red == 1)
            r = n_lll_rem_l_fullword_limited(a[2], a[1], a[0], mod, R->c0, R->c1);
        else
            r = n_lll_rem_l(a[2], a[1], a[0], mod, R->c0, R->c1);
    }

    return r;
}

/* one limb: reduces len accumulators in place */
static void
tri_reduce_1(nn_ptr acc, slong len, const tri_red_struct * R)
{
    slong i;
    for (i = 0; i < len; i++)
        acc[i] = tri_red(acc + i, R, 1);
}

FLINT_FORCE_INLINE void
tri_set(nn_ptr acc, nn_srcptr b, slong len, const int nlimbs)
{
    slong i;

    if (nlimbs == 1)
    {
        for (i = 0; i < len; i++)
            acc[i] = b[i];
    }
    else if (nlimbs == 2)
    {
        for (i = 0; i < len; i++)
        {
            acc[2 * i] = b[i];
            acc[2 * i + 1] = 0;
        }
    }
    else
    {
        for (i = 0; i < len; i++)
        {
            acc[3 * i] = b[i];
            acc[3 * i + 1] = 0;
            acc[3 * i + 2] = 0;
        }
    }
}

/* acc[i] += c * v[i] for i < len, with v[i], c < n */
FLINT_FORCE_INLINE void
tri_addmul(nn_ptr acc, nn_srcptr v, slong len, ulong c,
           const tri_red_struct * R, const int nlimbs)
{
    if (nlimbs == 1)
        _nmod_vec_nored_scalar_addmul_halflimb(acc, v, len, c);
    else if (nlimbs == 2 && R->halflimb)
        _nmod_vec_nored_ll_scalar_addmul_halflimb(acc, v, len, c);
    else if (nlimbs == 2)
        _nmod_vec_nored_ll_scalar_addmul(acc, v, len, c);
    else
        _nmod_vec_nored_lll_scalar_addmul(acc, v, len, c);
}

/*
    Inverses of the diagonal entries of T (2r limbs): inv[j], and if
    n < 2^(FLINT_BITS - 1) their precomputations for n_mulmod_shoup in
    inv[r + j].
*/
static void
tri_inv_diag(nn_ptr inv, const nmod_mat_t T)
{
    slong j, r = T->r;
    nmod_t mod = T->mod;

    for (j = 0; j < r; j++)
        inv[j] = nmod_mat_entry(T, j, j);
    _nmod_vec_invert(inv, inv, r, mod);

    if (NMOD_CAN_USE_SHOUP(mod))
        for (j = 0; j < r; j++)
            inv[r + j] = _nmod_mulmod_precomp_shoup(inv[j], mod);
}

FLINT_FORCE_INLINE ulong
tri_mul_inv(ulong s, nn_srcptr inv, slong j, slong r, nmod_t mod)
{
    if (NMOD_CAN_USE_SHOUP(mod))
        return n_mulmod_shoup(inv[j], s, inv[r + j], mod.n);
    else
        return nmod_mul(s, inv[j], mod);
}

/*
    Y = T'^(-1) C by rows of Y, for the r x r triangular T' with entries
    T'[j][t] = tt[j * rs + t * cs], lower if lower != 0 and upper otherwise,
    and the r x w matrices Y and C (row strides ys and cstr; Y may alias C):
    Y[j] = (C[j] - sum_t T'[j][t] Y[t]) / T'[j][j], where t < j (lower) or
    t > j (upper). The accumulators acc hold w entries.
*/
FLINT_FORCE_INLINE void
_tri_rows(nn_ptr Y, slong ys, nn_srcptr C, slong cstr, slong w,
          nn_srcptr tt, slong rs, slong cs, slong r, int lower,
          nn_ptr acc, nn_srcptr inv, const tri_red_struct * R,
          const int nlimbs)
{
    slong jj, j, t, t0, t1, c, cnt;
    nmod_t mod = R->mod;
    ulong e;

    for (jj = 0; jj < r; jj++)
    {
        nn_ptr y;
        nn_srcptr tj;

        j = lower ? jj : r - 1 - jj;
        y = Y + j * ys;
        tj = tt + j * rs;

        if (jj == 0)
        {
            /* no terms */
            if (inv == NULL)
            {
                if (y != C + j * cstr)
                    for (c = 0; c < w; c++)
                        y[c] = C[j * cstr + c];
            }
            else
            {
                for (c = 0; c < w; c++)
                    y[c] = tri_mul_inv(C[j * cstr + c], inv, j, r, mod);
            }
            continue;
        }

        tri_set(acc, C + j * cstr, w, nlimbs);

        t0 = lower ? 0 : j + 1;
        t1 = lower ? j : r;
        cnt = 0;
        for (t = t0; t < t1; t++)
        {
            if ((e = tj[t * cs]) != 0)
            {
                if (nlimbs == 1 && cnt == R->m)
                {
                    tri_reduce_1(acc, w, R);
                    cnt = 0;
                }
                cnt++;
                tri_addmul(acc, Y + t * ys, w, mod.n - e, R, nlimbs);
            }
        }

        if (inv == NULL)
            for (c = 0; c < w; c++)
                y[c] = tri_red(acc + nlimbs * c, R, nlimbs);
        else
            for (c = 0; c < w; c++)
                y[c] = tri_mul_inv(tri_red(acc + nlimbs * c, R, nlimbs),
                                   inv, j, r, mod);
    }
}

static void
_tri_rows_any(nn_ptr Y, slong ys, nn_srcptr C, slong cstr, slong w,
              nn_srcptr tt, slong rs, slong cs, slong r, int lower,
              nn_ptr acc, nn_srcptr inv, const tri_red_struct * R)
{
    if (R->nlimbs == 1)
        _tri_rows(Y, ys, C, cstr, w, tt, rs, cs, r, lower, acc, inv, R, 1);
    else if (R->nlimbs == 2)
        _tri_rows(Y, ys, C, cstr, w, tt, rs, cs, r, lower, acc, inv, R, 2);
    else
        _tri_rows(Y, ys, C, cstr, w, tt, rs, cs, r, lower, acc, inv, R, 3);
}

/* columns of B processed at a time by the rows kernel */
#define TRI_ROWS_BLOCK 256

/* buffers of len limbs aligned on 64 bytes: the vectors of the
   accumulation routines then never cross a cache line or a page */
#define TRI_BUFFER(len) \
    ((nn_ptr) (((size_t) TMP_ALLOC((len) * sizeof(ulong) + 64) + 63) \
               & ~(size_t) 63))

#define TRI_ROUND8(w) (((w) + 7) & ~(slong) 7)

void
_nmod_mat_solve_tri_rows_delayed(nmod_mat_t X, const nmod_mat_t T,
                                 const nmod_mat_t B, int unit, int upper,
                                 int right)
{
    slong i, j, i0, w, ys, r = T->r, k = right ? B->r : B->c;
    nmod_t mod = T->mod;
    nn_ptr inv = NULL, acc, Y = NULL;
    tri_red_struct R[1];
    TMP_INIT;

    if (r == 0 || k == 0)
        return;

    if (mod.n == 1)
    {
        nmod_mat_zero(X);
        return;
    }

    TMP_START;
    tri_red_init(R, r, mod);
    acc = TRI_BUFFER(R->nlimbs * FLINT_MIN(k, TRI_ROWS_BLOCK));
    if (!unit)
    {
        inv = TMP_ALLOC(2 * r * sizeof(ulong));
        tri_inv_diag(inv, T);
    }

    if (!right)
    {
        /* T X = B, by rows of X */
        for (i0 = 0; i0 < k; i0 += TRI_ROWS_BLOCK)
        {
            w = FLINT_MIN(k - i0, TRI_ROWS_BLOCK);
            _tri_rows_any(X->entries + i0, X->stride, B->entries + i0,
                          B->stride, w, T->entries, T->stride, 1, r, !upper,
                          acc, inv, R);
        }
    }
    else
    {
        /* X T = B: T^T X^T = B^T by rows of X^T, for blocks of rows of B
           transposed in Y (rows of ys entries, aligned) */
        Y = TRI_BUFFER(r * TRI_ROUND8(FLINT_MIN(k, TRI_ROWS_BLOCK)));

        for (i0 = 0; i0 < k; i0 += TRI_ROWS_BLOCK)
        {
            w = FLINT_MIN(k - i0, TRI_ROWS_BLOCK);
            ys = TRI_ROUND8(w);

            for (i = 0; i < w; i++)
            {
                nn_srcptr b = nmod_mat_entry_ptr(B, i0 + i, 0);
                for (j = 0; j < r; j++)
                    Y[j * ys + i] = b[j];
            }

            _tri_rows_any(Y, ys, Y, ys, w, T->entries, 1, T->stride, r,
                          upper, acc, inv, R);

            for (i = 0; i < w; i++)
            {
                nn_ptr x = nmod_mat_entry_ptr(X, i0 + i, 0);
                for (j = 0; j < r; j++)
                    x[j] = Y[j * ys + i];
            }
        }
    }

    TMP_END;
}

/*
    One row x of X T = B (b the row of B, may alias x): for each j (from
    the last for T lower, from the first for T upper), x[j] is reduced, then
    the products with the row j of T are added to the accumulators of the
    entries not yet computed (at most r - 1 terms for each).
*/
FLINT_FORCE_INLINE void
_tri_right_axpy(nn_ptr x, nn_srcptr b, nn_ptr acc, const nmod_mat_t T,
                nn_srcptr inv, int upper, const tri_red_struct * R,
                const int nlimbs)
{
    slong jj, j, t, h, len, r = T->r, cnt = 0;
    nmod_t mod = R->mod;
    nn_srcptr tj;
    ulong s, c;

    tri_set(acc, b, r, nlimbs);

    for (jj = 0; jj < r; jj++)
    {
        j = upper ? jj : r - 1 - jj;

        s = (jj == 0) ? acc[nlimbs * j] : tri_red(acc + nlimbs * j, R, nlimbs);
        if (inv != NULL)
            s = tri_mul_inv(s, inv, j, r, mod);
        x[j] = s;

        if (s == 0 || jj == r - 1)
            continue;

        if (nlimbs == 1 && cnt == R->m)
        {
            if (upper)
                tri_reduce_1(acc + j + 1, r - j - 1, R);
            else
                tri_reduce_1(acc, j, R);
            cnt = 0;
        }
        cnt++;

        c = mod.n - s;
        if (upper)
        {
            /* entries j + 1, ..., r - 1; with one limb, the first ones
               separately so that the vectors of the routine start at the
               same positions for all j (modulo 4) */
            tj = nmod_mat_entry_ptr(T, j, j + 1);
            len = r - j - 1;
            h = (nlimbs == 1) ? (len & 3) : 0;
            for (t = 0; t < h; t++)
                acc[j + 1 + t] += c * tj[t];
            tri_addmul(acc + nlimbs * (j + 1 + h), tj + h, len - h, c,
                       R, nlimbs);
        }
        else
        {
            tri_addmul(acc, nmod_mat_entry_ptr(T, j, 0), j, c, R, nlimbs);
        }
    }
}

void
_nmod_mat_solve_tri_right_axpy_delayed(nmod_mat_t X, const nmod_mat_t T,
                                       const nmod_mat_t B, int unit,
                                       int upper)
{
    slong i, r = T->r, k = B->r;
    nmod_t mod = T->mod;
    nn_ptr inv = NULL, acc;
    tri_red_struct R[1];
    TMP_INIT;

    if (r == 0 || k == 0)
        return;

    if (mod.n == 1)
    {
        nmod_mat_zero(X);
        return;
    }

    TMP_START;
    tri_red_init(R, r, mod);
    acc = TRI_BUFFER(R->nlimbs * r);
    if (!unit)
    {
        inv = TMP_ALLOC(2 * r * sizeof(ulong));
        tri_inv_diag(inv, T);
    }

    if (R->nlimbs == 1)
        for (i = 0; i < k; i++)
            _tri_right_axpy(nmod_mat_entry_ptr(X, i, 0),
                nmod_mat_entry_ptr(B, i, 0), acc, T, inv, upper, R, 1);
    else if (R->nlimbs == 2)
        for (i = 0; i < k; i++)
            _tri_right_axpy(nmod_mat_entry_ptr(X, i, 0),
                nmod_mat_entry_ptr(B, i, 0), acc, T, inv, upper, R, 2);
    else
        for (i = 0; i < k; i++)
            _tri_right_axpy(nmod_mat_entry_ptr(X, i, 0),
                nmod_mat_entry_ptr(B, i, 0), acc, T, inv, upper, R, 3);

    TMP_END;
}
