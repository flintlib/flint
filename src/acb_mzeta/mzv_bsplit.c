/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "thread_support.h"
#include "fmpz_vec.h"
#include "impl.h"

/* Multiple zeta values by Hoelder convolution [BBBL2001] and binary
   splitting (quasilinear in the precision).

   Splitting the iterated integral of an admissible word w = e_1 ... e_k at
   1/2 (Chen's formula) and substituting t -> 1 - t in the part above 1/2:

     zeta(w) = sum_{j=0}^{k} I(w*[k-j..k)) I(w[j..k)),
     I(u) = int_{1/2 > t_1 > ... > t_L > 0} omega_{u_1}(t_1) ... omega_{u_L}(t_L),

   w* = the dual word (reversed, letters swapped; also admissible).  Every
   suffix of w (and of w*) ends with x1, so I(suffix) is a multiple
   polylogarithm at 1/2: for w = x0^(c_1 - 1) x1 ... x0^(c_d - 1) x1 and a
   suffix starting at offset o in block i (first entry c' = c_i - o),

     I = Li_{c', c_{i+1}, ..., c_d}(1/2) = sum_{n >= 1} z n^-c' H_{i+1}(n - 1) z^(n-1),

   z = 1/2, H_i(n) = sum_{n >= n_i > ... > n_d >= 1} n_i^-c_i ... n_d^-c_d.
   With G_i(n) = z^n H_i(n) (G_{d+1}(n) = z^n):

     G_i(n) = z G_i(n-1) + z n^-c_i G_{i+1}(n-1),
     A_j(n) = A_j(n-1) + z n^-c' G_{i+1}(n-1)        (one A_j per suffix j),

   a linear recurrence X(n) = M(n) X(n-1) of dimension d + k (the G_i with
   i >= 2 and the A_j), lower triangular in the order G_{d+1}, G_d, ...,
   G_2, A_0, ..., A_{k-1}.  With D(n) = 2 n^C, C = max c_i, the matrices
   D(n) M(n) have integer entries; their product over n = 1, ..., N is
   computed by binary splitting in ball arithmetic (exact at the leaves,
   rounded to the working precision higher up), and X(N) = P X(0) / prod D.

   Truncation: H_i(n) <= (1 + log n)^d, so each A_j omits at most
   sum_{n > N} 2^-n (1 + log n)^d <= 2^-N (1 + log(N+1))^d / (1 - q),
   q = (1 + 1/N)^d / 2. */

/* lower triangular square matrices of arb */
typedef struct
{
    arb_ptr e;
    slong n;
}
tmat_struct;

#define TE(M, i, j) ((M)->e + (i) * (M)->n + (j))

static void
tmat_init(tmat_struct * M, slong n)
{
    M->n = n;
    M->e = _arb_vec_init(n * n);
}

static void
tmat_clear(tmat_struct * M)
{
    _arb_vec_clear(M->e, M->n * M->n);
}

/* C = A B, lower triangular */
static void
tmat_mul(tmat_struct * C, const tmat_struct * A, const tmat_struct * B, slong prec)
{
    slong n = A->n, i, j, l;
    arb_t t;
    arb_init(t);
    for (i = 0; i < n; i++)
        for (j = 0; j <= i; j++)
        {
            arb_ptr c = TE(C, i, j);
            arb_zero(c);
            for (l = j; l <= i; l++)
            {
                if (arb_is_zero(TE(A, i, l)) || arb_is_zero(TE(B, l, j)))
                    continue;
                arb_addmul(c, TE(A, i, l), TE(B, l, j), prec);
            }
        }
    arb_clear(t);
}

typedef struct
{
    slong d, k, C, dim;
    const ulong * c;        /* the composition */
    slong * blk, * cp;      /* suffix j: block i and first entry c' */
}
holder_ctx;

/* the integer matrix D(n) M(n) */
static void
leaf(tmat_struct * Q, const holder_ctx * H, ulong n)
{
    slong d = H->d, k = H->k, C = H->C, i, j;
    fmpz_t t, u;

    fmpz_init(t);
    fmpz_init(u);
    _arb_vec_zero(Q->e, Q->n * Q->n);

    fmpz_set_ui(t, n);
    fmpz_pow_ui(t, t, C);           /* n^C = D(n) z */

    /* index of G_i (2 <= i <= d + 1): d + 1 - i;  A_j: d + j */
    arb_set_fmpz(TE(Q, 0, 0), t);   /* G_{d+1} = z^n */
    for (i = d; i >= 2; i--)
    {
        slong r = d + 1 - i;
        arb_set_fmpz(TE(Q, r, r), t);
        fmpz_set_ui(u, n);
        fmpz_pow_ui(u, u, C - H->c[i - 1]);
        arb_set_fmpz(TE(Q, r, r - 1), u);    /* G_{i+1} has index r - 1 */
    }
    for (j = 0; j < k; j++)
    {
        slong r = d + j, i1 = H->blk[j] + 1;     /* G_{i+1}, i = blk (1-based) */
        arb_set_fmpz(TE(Q, r, r), t);
        arb_mul_2exp_si(TE(Q, r, r), TE(Q, r, r), 1);   /* 2 n^C */
        fmpz_set_ui(u, n);
        fmpz_pow_ui(u, u, C - H->cp[j]);
        arb_set_fmpz(TE(Q, r, d + 1 - i1), u);
    }

    fmpz_clear(t);
    fmpz_clear(u);
}

/* exact version on fmpz (lower triangular, row-major), for the lower
   levels where the entries are small */
static void
leaf_fmpz(fmpz * Q, const holder_ctx * H, ulong n)
{
    slong d = H->d, k = H->k, C = H->C, dim = H->dim, i, j;
    fmpz * pw = _fmpz_vec_init(C + 1);     /* n^0, ..., n^C */
    fmpz_one(pw + 0);
    for (i = 1; i <= C; i++)
        fmpz_mul_ui(pw + i, pw + i - 1, n);
    _fmpz_vec_zero(Q, dim * dim);
    fmpz_set(Q + 0, pw + C);
    for (i = d; i >= 2; i--)
    {
        slong r = d + 1 - i;
        fmpz_set(Q + r * dim + r, pw + C);
        fmpz_set(Q + r * dim + r - 1, pw + C - H->c[i - 1]);
    }
    for (j = 0; j < k; j++)
    {
        slong r = d + j, i1 = H->blk[j] + 1;
        fmpz_mul_2exp(Q + r * dim + r, pw + C, 1);
        fmpz_set(Q + r * dim + (d + 1 - i1), pw + C - H->cp[j]);
    }
    _fmpz_vec_clear(pw, C + 1);
}

static void
bsplit_fmpz(fmpz * P, fmpz_t Den, const holder_ctx * H, ulong a, ulong b)
{
    slong dim = H->dim, i, j, l;
    if (b - a == 1)
    {
        leaf_fmpz(P, H, b);
        fmpz_set_ui(Den, b);
        fmpz_pow_ui(Den, Den, H->C);
        fmpz_mul_2exp(Den, Den, 1);
    }
    else
    {
        ulong m = a + (b - a) / 2;
        fmpz * P1 = _fmpz_vec_init(dim * dim), * P2 = _fmpz_vec_init(dim * dim);
        fmpz_t D2;
        fmpz_init(D2);
        bsplit_fmpz(P1, Den, H, a, m);
        bsplit_fmpz(P2, D2, H, m, b);
        for (i = 0; i < dim; i++)
            for (j = 0; j <= i; j++)
            {
                fmpz * c = P + i * dim + j;
                fmpz_zero(c);
                for (l = j; l <= i; l++)
                    if (!fmpz_is_zero(P2 + i * dim + l) && !fmpz_is_zero(P1 + l * dim + j))
                        fmpz_addmul(c, P2 + i * dim + l, P1 + l * dim + j);
            }
        fmpz_mul(Den, Den, D2);
        _fmpz_vec_clear(P1, dim * dim);
        _fmpz_vec_clear(P2, dim * dim);
        fmpz_clear(D2);
    }
}

/* P = Q(b) ... Q(a + 1), Den = prod_{a < n <= b} D(n) */
static void
bsplit(tmat_struct * P, arb_t Den, const holder_ctx * H, ulong a, ulong b, slong prec)
{
    /* the entries have at most ~(b - a) (C log2 b + 2) bits: exact
       arithmetic below the working precision */
    if ((double) (b - a) * (H->C * FLINT_BIT_COUNT(b) + 2) < prec)
    {
        slong dim = H->dim, i;
        fmpz * Pz = _fmpz_vec_init(dim * dim);
        fmpz_t Dz;
        fmpz_init(Dz);
        bsplit_fmpz(Pz, Dz, H, a, b);
        for (i = 0; i < dim * dim; i++)
            arb_set_fmpz(P->e + i, Pz + i);
        arb_set_fmpz(Den, Dz);
        _fmpz_vec_clear(Pz, dim * dim);
        fmpz_clear(Dz);
        return;
    }

    if (b - a == 1)
    {
        fmpz_t t;
        fmpz_init(t);
        leaf(P, H, b);
        fmpz_set_ui(t, b);
        fmpz_pow_ui(t, t, H->C);
        fmpz_mul_2exp(t, t, 1);
        arb_set_fmpz(Den, t);
        fmpz_clear(t);
    }
    else
    {
        ulong m = a + (b - a) / 2;
        tmat_struct P1, P2;
        arb_t D2;
        tmat_init(&P1, H->dim);
        tmat_init(&P2, H->dim);
        arb_init(D2);
        bsplit(&P1, Den, H, a, m, prec);
        bsplit(&P2, D2, H, m, b, prec);
        tmat_mul(P, &P2, &P1, prec);
        arb_mul(Den, Den, D2, prec);
        tmat_clear(&P1);
        tmat_clear(&P2);
        arb_clear(D2);
    }
}

/* v = P e_0 for P = Q(b) ... Q(a + 1): only the first column is needed at
   the root, so the lowest segment carries a vector (the matrix-matrix
   products are halved in size at each level) */
static void
bsplit_col(arb_ptr v, arb_t Den, const holder_ctx * H, ulong a, ulong b, slong prec)
{
    slong n = H->dim, i, l;

    if (b - a <= 4)
    {
        tmat_struct P;
        tmat_init(&P, n);
        bsplit(&P, Den, H, a, b, prec);
        for (i = 0; i < n; i++)
            arb_set(v + i, TE(&P, i, 0));
        tmat_clear(&P);
    }
    else
    {
        ulong m = a + (b - a) / 2;
        tmat_struct P2;
        arb_ptr v1;
        arb_t D2;
        tmat_init(&P2, n);
        v1 = _arb_vec_init(n);
        arb_init(D2);
        bsplit_col(v1, Den, H, a, m, prec);
        bsplit(&P2, D2, H, m, b, prec);
        for (i = 0; i < n; i++)
        {
            arb_zero(v + i);
            for (l = 0; l <= i; l++)
                if (!arb_is_zero(TE(&P2, i, l)) && !arb_is_zero(v1 + l))
                    arb_addmul(v + i, TE(&P2, i, l), v1 + l, prec);
        }
        arb_mul(Den, Den, D2, prec);
        tmat_clear(&P2);
        _arb_vec_clear(v1, n);
        arb_clear(D2);
    }
}

/* res[j] = I(w[j..k)) for j = 0, ..., k (res[k] = 1), w admissible */
static void
suffix_values(arb_ptr res, const unsigned char * w, slong k, slong prec)
{
    holder_ctx H;
    ulong c[512];
    slong d, i, j, o, N, wp;
    arb_t Den;
    mag_t tail;

    /* the composition and, per suffix, its block and first entry */
    d = 0;
    {
        slong z = 0;
        for (i = 0; i < k; i++)
        {
            if (w[i] == 0)
                z++;
            else
            {
                c[d++] = z + 1;
                z = 0;
            }
        }
    }
    H.d = d;
    H.k = k;
    H.c = c;
    H.C = 0;
    for (i = 0; i < d; i++)
        H.C = FLINT_MAX(H.C, (slong) c[i]);
    H.blk = flint_malloc(sizeof(slong) * k);
    H.cp = flint_malloc(sizeof(slong) * k);
    j = 0;
    for (i = 0; i < d; i++)
        for (o = 0; o < (slong) c[i]; o++)
        {
            H.blk[j] = i + 1;            /* 1-based block */
            H.cp[j] = c[i] - o;
            j++;
        }
    H.dim = d + k;

    /* N: 2^-N (1 + log(N+1))^d / (1 - q) below 2^-(prec + 4), rigorous
       bound below */
    N = prec + 8;
    while (N * 1.0 - d * log2(1.0 + log(N + 1.0)) < prec + 8)
        N += 8;

    mag_init(tail);
    {
        arb_t t, q;
        arb_init(t);
        arb_init(q);
        /* q = (1 + 1/N)^d / 2 */
        arb_set_ui(q, N + 1);
        arb_div_ui(q, q, N, 30);
        arb_pow_ui(q, q, d, 30);
        arb_mul_2exp_si(q, q, -1);
        arb_sub_ui(q, q, 1, 30);
        arb_neg(q, q);                      /* 1 - q */
        arb_set_ui(t, N + 1);
        arb_log(t, t, 30);
        arb_add_ui(t, t, 1, 30);
        arb_pow_ui(t, t, d, 30);
        arb_mul_2exp_si(t, t, -N);
        arb_div(t, t, q, 30);
        if (arb_is_positive(q))
            arb_get_mag(tail, t);
        else
            mag_inf(tail);
        arb_clear(t);
        arb_clear(q);
    }

    wp = prec + 2 * FLINT_BIT_COUNT(N) + 2 * FLINT_BIT_COUNT(H.dim) + 16;

    {
        arb_ptr v = _arb_vec_init(H.dim);
        arb_init(Den);
        bsplit_col(v, Den, &H, 0, N, wp);

        /* X(N) = P X(0) / Den, X(0) = e_0 (G_{d+1}(0) = 1) */
        for (j = 0; j < k; j++)
        {
            arb_div(res + j, v + d + j, Den, wp);
            arb_add_error_mag(res + j, tail);
        }
        arb_one(res + k);
        _arb_vec_clear(v, H.dim);
        arb_clear(Den);
    }
    mag_clear(tail);
    flint_free(H.blk);
    flint_free(H.cp);
}

typedef struct
{
    arb_ptr res[2];
    const unsigned char * w[2];
    slong k, prec;
}
suffix_job;

static void
suffix_worker(slong i, void * vj)
{
    suffix_job * J = vj;
    suffix_values(J->res[i], J->w[i], J->k, J->prec);
}

void
_acb_mzeta_mzv_bsplit(arb_t res, const unsigned char * w, slong k, slong prec)
{
    unsigned char * ws;
    arb_ptr A, B;
    slong i, wp;

    if (k == 0)
    {
        arb_one(res);
        return;
    }
    if (w[0] != 0 || w[k - 1] != 1 || k > 500)
    {
        arb_indeterminate(res);
        return;
    }

    ws = flint_malloc(k);
    for (i = 0; i < k; i++)
        ws[i] = 1 - w[k - 1 - i];

    A = _arb_vec_init(k + 1);
    B = _arb_vec_init(k + 1);
    wp = prec + FLINT_BIT_COUNT(k) + 8;
    {
        suffix_job J;
        J.res[0] = A; J.res[1] = B;
        J.w[0] = ws; J.w[1] = w;
        J.k = k;
        J.prec = wp;
        flint_parallel_do(suffix_worker, &J, 2, flint_get_num_threads(), FLINT_PARALLEL_UNIFORM);
    }

    arb_zero(res);
    for (i = 0; i <= k; i++)
        arb_addmul(res, A + k - i, B + i, wp);
    arb_set_round(res, res, prec);

    _arb_vec_clear(A, k + 1);
    _arb_vec_clear(B, k + 1);
    flint_free(ws);
}
