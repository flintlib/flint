/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include "ulong_extras.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_lll.h"
#include "nmod.h"
#include "nmod_poly.h"
#include "mpn_extras.h"
#include "impl.h"

/* Shell handler: products of theta series of low-rank blocks.

   Z^d is split into mutually F-orthogonal blocks L_b of rank 2 (and one
   of rank 1 when d is odd), with bases (rows of B) taken from an
   LLL-reduced basis and the integral kernels of the previous blocks, so
   that L_1 + ... + L_k has finite index in Z^d.  On a coset v + sum L_b,
   with y = S v - s, c_b = B_b F y and G_b = B_b F B_b^T, the points
   n = v + sum_b B_b^T u_b have

     V(n) = sum_b (S u_b + t_b)^T G_b (S u_b + t_b),   t_b = G_b^-1 c_b,

   and Den (S u + t)^T G (S u + t) = Den (S^2 u^T G u + 2 S u.c) + Den c^T G^-1 c
   is an integer for Den the common denominator of the c^T G^-1 c.  So the
   generating function sum_n x^(Den V(n)) zeta^(N ell.n) of the coset is
   the product of the series of the blocks, each of which is enumerated
   directly in linear time (a rank-2 block has about as many points as
   exponents).  The k - 1 products run in Z[x][zeta]/(zeta^N - 1) by
   Kronecker substitution (slots of width 2N - 1, folded after each
   product) modulo a prime just above a bound for every count; a last
   rank-1 block with few terms is multiplied in directly.  The exponents of
   a block on a coset lie in one class rho_b mod q (q the gcd over the
   blocks of the differences, a multiple of Den), so the series are taken
   in x^q.  Among the choices of the first block (pairs of reduced basis
   vectors) the one with the least total length of the series is used. */

#define MAX_DEN 1024
#define MAX_INDEX 4096
#define MAX_BITS 40

typedef struct
{
    slong d, nb;
    slong rank[ACB_MZETA_SHELLS_MAX_DIM], off[ACB_MZETA_SHELLS_MAX_DIM];
    fmpz_mat_t B;           /* the basis rows of the blocks, stacked */
    slong G[ACB_MZETA_SHELLS_MAX_DIM][3];      /* A (, Bc, C) of each block */
    ulong cl[ACB_MZETA_SHELLS_MAX_DIM][2];     /* N ell.(basis row) mod N */
    ulong hdiag[ACB_MZETA_SHELLS_MAX_DIM];
    slong index;
    ulong Den;
    double work;            /* sum over the cosets of Den / q (series length per value) */
    ulong * C;
    ulong len;
}
prod_struct;

static ulong
slot_width(ulong N)
{
    return N == 1 ? 1 : 2 * N - 1;
}

/* rows of K: a Z-basis of {n in Z^d : P F n = 0} (P has r rows of rank r),
   reduced with respect to F; returns the number of rows */
static slong
kernel_basis(fmpz_mat_t K, const fmpz_mat_t P, const fmpz_mat_t F)
{
    slong d = fmpz_mat_nrows(F), r = fmpz_mat_nrows(P), k = d - r, i, j;
    fmpz_mat_t A, H, U, Kr, G, V;
    fmpz_lll_t fl;

    fmpz_mat_init(A, r, d);
    fmpz_mat_mul(A, P, F);
    fmpz_mat_init(H, d, r);
    fmpz_mat_init(U, d, d);
    {
        fmpz_mat_t At;
        fmpz_mat_init(At, d, r);
        fmpz_mat_transpose(At, A);
        /* U At = H: the rows of U after the rank of H span the kernel */
        fmpz_mat_hnf_transform(H, U, At);
        fmpz_mat_clear(At);
    }

    fmpz_mat_init(Kr, k, d);
    for (i = 0; i < k; i++)
        for (j = 0; j < d; j++)
            fmpz_set(fmpz_mat_entry(Kr, i, j), fmpz_mat_entry(U, r + i, j));

    /* LLL with respect to F: the Gram matrix Kr F Kr^T */
    fmpz_mat_init(G, k, k);
    fmpz_mat_init(V, k, d);
    fmpz_mat_mul(V, Kr, F);
    {
        fmpz_mat_t Krt;
        fmpz_mat_init(Krt, d, k);
        fmpz_mat_transpose(Krt, Kr);
        fmpz_mat_mul(G, V, Krt);
        fmpz_mat_clear(Krt);
    }
    {
        fmpz_mat_t W;
        fmpz_mat_init(W, k, k);
        fmpz_mat_one(W);
        fmpz_lll_context_init(fl, 0.99, 0.51, GRAM, EXACT);
        fmpz_lll(G, W, fl);
        fmpz_mat_mul(K, W, Kr);
        fmpz_mat_clear(W);
    }

    fmpz_mat_clear(A);
    fmpz_mat_clear(H);
    fmpz_mat_clear(U);
    fmpz_mat_clear(Kr);
    fmpz_mat_clear(G);
    fmpz_mat_clear(V);
    return k;
}

/* y = S v - s for the coset k, the block vectors c_b = B_b F y and the
   character of v; returns the residue */
static ulong
coset_data(fmpz * c, slong k, const prod_struct * P, const _acb_mzeta_shells_struct * T)
{
    slong d = T->d, i, j, kk = k;
    fmpz * v = _fmpz_vec_init(d), * Fy = _fmpz_vec_init(d);
    fmpz_t t;
    ulong cv;

    fmpz_init(t);
    for (i = 0; i < d; i++)
    {
        fmpz_set_ui(v + i, kk % P->hdiag[i]);
        kk /= P->hdiag[i];
    }
    _fmpz_vec_dot_general(t, NULL, 0, T->ell, v, 0, d);
    cv = fmpz_fdiv_ui(t, T->N);
    for (i = 0; i < d; i++)
    {
        fmpz_mul_ui(v + i, v + i, T->S);
        fmpz_sub(v + i, v + i, T->s + i);
    }
    for (i = 0; i < d; i++)
    {
        fmpz_zero(Fy + i);
        for (j = 0; j < d; j++)
            fmpz_addmul(Fy + i, fmpz_mat_entry(T->F, i, j), v + j);
    }
    for (i = 0; i < d; i++)
        _fmpz_vec_dot_general(c + i, NULL, 0, fmpz_mat_row(P->B, i), Fy, 0, d);

    fmpz_clear(t);
    _fmpz_vec_clear(v, d);
    _fmpz_vec_clear(Fy, d);
    return cv;
}

/* c^T G^-1 c = num / den (reduced) for the block b */
static void
block_const(fmpz_t num, fmpz_t den, const fmpz * c, const slong * G, slong rank)
{
    fmpz_t t, g;
    fmpz_init(t);
    fmpz_init(g);
    if (rank == 1)
    {
        fmpz_mul(num, c, c);
        fmpz_set_si(den, G[0]);
    }
    else
    {
        /* adj G = [[C, -Bc], [-Bc, A]] */
        fmpz_mul(num, c, c);
        fmpz_mul_si(num, num, G[2]);
        fmpz_mul(t, c + 1, c + 1);
        fmpz_addmul_si(num, t, G[0]);
        fmpz_mul(t, c, c + 1);
        fmpz_mul_si(t, t, 2 * G[1]);
        fmpz_sub(num, num, t);
        fmpz_set_si(den, G[0]);
        fmpz_mul_si(den, den, G[2]);
        fmpz_set_si(t, G[1]);
        fmpz_submul(den, t, t);
    }
    fmpz_gcd(g, num, den);
    fmpz_divexact(num, num, g);
    fmpz_divexact(den, den, g);
    fmpz_clear(t);
    fmpz_clear(g);
}

/* for the block b on a coset: K0 = Den c^T G^-1 c (the exponent at u = 0)
   and q = the gcd of the differences e(u) - e(0) = Den (S^2 u^T G u + 2 S u.c),
   which is the gcd of the values at +-e_i and the cross terms 2 S^2 G_ij;
   0 if they do not fit */
static int
block_step(slong * K0, ulong * q, const fmpz * c, const slong * G, slong rank,
    ulong Den, ulong S)
{
    fmpz_t num, den, t, u, g;
    slong i;
    int ok;

    fmpz_init(num);
    fmpz_init(den);
    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(g);
    block_const(num, den, c, G, rank);
    fmpz_mul_ui(num, num, Den);
    fmpz_divexact(num, num, den);
    for (i = 0; i < rank; i++)
    {
        fmpz_set_si(t, G[2 * i]);           /* A or C */
        fmpz_mul_ui(t, t, S * S);
        fmpz_mul_ui(u, c + i, 2 * S);
        fmpz_add(den, t, u);
        fmpz_gcd(g, g, den);
        fmpz_sub(den, t, u);
        fmpz_gcd(g, g, den);
    }
    if (rank == 2)
    {
        fmpz_set_si(t, G[1]);
        fmpz_mul_ui(t, t, 2 * S * S);
        fmpz_gcd(g, g, t);
    }
    fmpz_mul_ui(g, g, Den);
    ok = (fmpz_bits(num) <= 61) && !fmpz_is_zero(g) && fmpz_bits(g) <= 61;
    if (ok)
    {
        *K0 = fmpz_get_si(num);
        *q = fmpz_get_ui(g);
    }
    fmpz_clear(num);
    fmpz_clear(den);
    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(g);
    return ok;
}

/* the blocks with first block the rows i, j of R (i = j: rank 1);
   0 if not usable */
static int
blocks_build(prod_struct * P, const fmpz_mat_t R, slong i, slong j, const _acb_mzeta_shells_struct * T)
{
    slong d = T->d, r = 0, b, k, l, m;
    fmpz_mat_t H;
    fmpz * c;
    fmpz_t num, den;
    ulong Den;
    int ok = 1;

    /* the bases */
    P->nb = 0;
    for (k = 0; k < d; k++)
        fmpz_set(fmpz_mat_entry(P->B, 0, k), fmpz_mat_entry(R, i, k));
    r = 1;
    if (j != i)
    {
        for (k = 0; k < d; k++)
            fmpz_set(fmpz_mat_entry(P->B, 1, k), fmpz_mat_entry(R, j, k));
        r = 2;
    }
    P->rank[0] = r;
    P->off[0] = 0;
    P->nb = 1;
    while (r < d)
    {
        fmpz_mat_t Pw, K;
        slong kr, take;
        fmpz_mat_window_init(Pw, P->B, 0, 0, r, d);
        fmpz_mat_init(K, d - r, d);
        kr = kernel_basis(K, Pw, T->F);
        take = (kr >= 2) ? 2 : 1;
        for (l = 0; l < take; l++)
            for (k = 0; k < d; k++)
                fmpz_set(fmpz_mat_entry(P->B, r + l, k), fmpz_mat_entry(K, l, k));
        P->off[P->nb] = r;
        P->rank[P->nb] = take;
        P->nb++;
        r += take;
        fmpz_mat_window_clear(Pw);
        fmpz_mat_clear(K);
    }

    /* Gram matrices and characters of the blocks */
    {
        fmpz_mat_t BF, Bt, G;
        fmpz_mat_init(BF, d, d);
        fmpz_mat_init(Bt, d, d);
        fmpz_mat_init(G, d, d);
        fmpz_mat_mul(BF, P->B, T->F);
        fmpz_mat_transpose(Bt, P->B);
        fmpz_mat_mul(G, BF, Bt);
        for (b = 0; b < P->nb && ok; b++)
        {
            slong o = P->off[b];
            for (l = 0; l < P->rank[b]; l++)
            {
                fmpz_t t;
                fmpz_init(t);
                _fmpz_vec_dot_general(t, NULL, 0, T->ell, fmpz_mat_row(P->B, o + l), 0, d);
                P->cl[b][l] = fmpz_fdiv_ui(t, T->N);
                fmpz_clear(t);
            }
            for (l = 0; l < 3 && ok; l++)
            {
                const fmpz * e;
                if (P->rank[b] == 1 && l > 0)
                    break;
                e = (l == 0) ? fmpz_mat_entry(G, o, o) : (l == 1) ? fmpz_mat_entry(G, o, o + 1)
                    : fmpz_mat_entry(G, o + 1, o + 1);
                if (fmpz_bits(e) > MAX_BITS)
                    ok = 0;
                else
                    P->G[b][l] = fmpz_get_si(e);
            }
            /* the blocks must be orthogonal (by construction) */
            for (k = 0; k < d && ok; k++)
                for (l = 0; l < P->rank[b]; l++)
                    if ((k < o || k >= o + P->rank[b]) && !fmpz_is_zero(fmpz_mat_entry(G, o + l, k)))
                        flint_throw(FLINT_ERROR, "shells_products: blocks not orthogonal\n");
        }
        fmpz_mat_clear(BF);
        fmpz_mat_clear(Bt);
        fmpz_mat_clear(G);
    }
    if (!ok)
        return 0;

    /* coset representatives: the box of the diagonal of the HNF */
    fmpz_mat_init(H, d, d);
    fmpz_mat_hnf(H, P->B);
    P->index = 1;
    for (k = 0; k < d && ok; k++)
    {
        if (!fmpz_abs_fits_ui(fmpz_mat_entry(H, k, k)) || fmpz_cmp_ui(fmpz_mat_entry(H, k, k), MAX_INDEX) > 0)
            ok = 0;
        else
        {
            P->hdiag[k] = fmpz_get_ui(fmpz_mat_entry(H, k, k));
            P->index *= P->hdiag[k];
            if (P->index > MAX_INDEX || P->index == 0)
                ok = 0;
        }
    }
    fmpz_mat_clear(H);
    if (!ok)
        return 0;

    /* Den over the cosets and blocks */
    c = _fmpz_vec_init(d);
    fmpz_init(num);
    fmpz_init(den);
    Den = 1;
    for (m = 0; m < P->index && ok; m++)
    {
        coset_data(c, m, P, T);
        for (b = 0; b < P->nb && ok; b++)
        {
            block_const(num, den, c + P->off[b], P->G[b], P->rank[b]);
            if (fmpz_cmp_ui(den, MAX_DEN) > 0)
                ok = 0;
            else
            {
                ulong x = fmpz_get_ui(den);
                Den = (Den / n_gcd(Den, x)) * x;
                if (Den > MAX_DEN)
                    ok = 0;
            }
        }
    }
    P->Den = Den;

    /* the length of the series per shell value, summed over the cosets */
    P->work = 0.0;
    for (m = 0; m < P->index && ok; m++)
    {
        ulong q = 0, qb;
        slong K0;
        coset_data(c, m, P, T);
        for (b = 0; b < P->nb && ok; b++)
        {
            ok = block_step(&K0, &qb, c + P->off[b], P->G[b], P->rank[b], Den, T->S);
            q = n_gcd(q, qb);
        }
        if (ok)
            P->work += (double) Den / q;
    }
    _fmpz_vec_clear(c, d);
    fmpz_clear(num);
    fmpz_clear(den);
    return ok;
}

static void
prod_clear(void * data)
{
    prod_struct * P = data;
    fmpz_mat_clear(P->B);
    flint_free(P->C);
    flint_free(P);
}

static void *
prod_init(const _acb_mzeta_shells_struct * T)
{
    slong d = T->d, i, j, bi = -1, bj = -1;
    double best = INFINITY;
    prod_struct * P;
    fmpz_mat_t R, G;
    fmpz_lll_t fl;

    /* F must be positive definite (leading principal minors, exactly) */
    for (i = 1; i <= d; i++)
    {
        fmpz_mat_t W;
        fmpz_t det;
        int pos;
        fmpz_mat_window_init(W, T->F, 0, 0, i, i);
        fmpz_init(det);
        fmpz_mat_det(det, W);
        pos = (fmpz_sgn(det) > 0);
        fmpz_clear(det);
        fmpz_mat_window_clear(W);
        if (!pos)
            return NULL;
    }

    P = flint_calloc(1, sizeof(prod_struct));
    P->d = d;
    fmpz_mat_init(P->B, d, d);

    /* an F-reduced basis R of Z^d */
    fmpz_mat_init(R, d, d);
    fmpz_mat_init_set(G, T->F);
    fmpz_mat_one(R);
    fmpz_lll_context_init(fl, 0.99, 0.51, GRAM, EXACT);
    fmpz_lll(G, R, fl);
    fmpz_mat_clear(G);

    /* the first block: a pair of rows of R (a single row for d = 1) */
    for (i = 0; i < d; i++)
    {
        for (j = (d == 1) ? i : i + 1; j < d || (d == 1 && j == i); j++)
        {
            if (blocks_build(P, R, i, j, T))
            {
                double score = P->work;
                if (score < best)
                {
                    best = score;
                    bi = i;
                    bj = j;
                }
            }
            if (d == 1)
                break;
        }
    }

    if (bi == -1 || !blocks_build(P, R, bi, bj, T))
    {
        fmpz_mat_clear(R);
        prod_clear(P);
        return NULL;
    }
    fmpz_mat_clear(R);
    return P;
}

static double
prod_cost(const void * data, const _acb_mzeta_shells_struct * T, ulong len)
{
    const prod_struct * P = data;
    /* the series length of an average coset */
    double EW = P->work / P->index * ((double) T->g * len + T->r0 + 1) * slot_width(T->N);

    /* memory: the block series and two products of EW coefficients */
    if (EW > 1e8 || (double) len * T->N > 2e8)
        return INFINITY;
    /* the enumeration of the blocks (at most Den times the points, which
       are not compressed) and nb - 1 multiplications per coset */
    return P->index * (P->nb * 3.0 * EW * P->Den + (P->nb - 1) * EW * (log2(EW + 1) + 4) * 2.0);
}

/* the points of the block b with exponent e <= Emax and k = (e - rho)/q
   <= Kmax (e = K0 + Den (S^2 u^T G u + 2 S u.c) = rho mod q): either added
   at k W + r to the dense array f (length (Kmax + 1) W, zeroed) or listed
   as terms (k, r); returns the number of points, or -1 if the exact
   values could overflow */
static slong
block_points(ulong * f, ulong * te, ulong * tr, slong cap, const prod_struct * P, slong b,
    const fmpz * c, ulong r0, const _acb_mzeta_shells_struct * T, ulong Emax,
    slong K0, ulong q, ulong rho, ulong Kmax)
{
    const slong * G = P->G[b];
    ulong N = T->N, W = slot_width(N), Den = P->Den, S = T->S;
    double Qmax = (double) Emax / Den, cd[2] = {0.0, 0.0}, td[2] = {0.0, 0.0}, R2, bnd;
    slong n = 0, u1, u2, lo1 = 0, hi1 = -1, lo2, hi2, cs[2] = {0, 0}, i;
    int rank = P->rank[b];

    for (i = 0; i < rank; i++)
    {
        if (fmpz_bits(c + i) > FLINT_BITS - 24)
            return -1;
        cs[i] = fmpz_get_si(c + i);
        cd[i] = (double) cs[i];
    }

    if (rank == 1)
    {
        double A = G[0], R, U1;
        td[0] = cd[0] / A;
        R = sqrt(Qmax / A) * (1 + 1e-9) + 1e-6;
        lo1 = (slong) floor((-R - td[0]) / S) - 1;
        hi1 = (slong) ceil((R - td[0]) / S) + 1;
        lo2 = hi2 = 0;
        U1 = FLINT_MAX(fabs((double) lo1), fabs((double) hi1));
        bnd = Den * (double) S * S * A * U1 * U1 + 2.0 * Den * S * fabs(cd[0]) * U1 + fabs((double) K0);
    }
    else
    {
        double A = G[0], Bc = G[1], Cc = G[2], D = A * Cc - Bc * Bc, U1, U2;
        /* t = G^-1 c */
        td[0] = (Cc * cd[0] - Bc * cd[1]) / D;
        td[1] = (A * cd[1] - Bc * cd[0]) / D;
        R2 = sqrt(Qmax * A / D) * (1 + 1e-9) + 1e-6;
        lo2 = (slong) floor((-R2 - td[1]) / S) - 1;
        hi2 = (slong) ceil((R2 - td[1]) / S) + 1;
        U2 = FLINT_MAX(fabs((double) lo2), fabs((double) hi2));
        /* |z1| <= |Bc| |z2| / A + sqrt(Qmax / A) */
        U1 = ((fabs(Bc) * (R2 + S) / A + sqrt(Qmax / A) + 1) + fabs(td[0])) / S + 2;
        bnd = Den * (double) S * S * (A * U1 * U1 + 2 * fabs(Bc) * U1 * U2 + Cc * U2 * U2)
            + 2.0 * Den * S * (fabs(cd[0]) * U1 + fabs(cd[1]) * U2) + fabs((double) K0);
    }
    /* the exact exponents (slong) stay below 2^(FLINT_BITS - 3) */
    if (!(bnd < ldexp(1.0, FLINT_BITS - 3)))
        return -1;

    for (u2 = lo2; u2 <= hi2; u2++)
    {
        slong base = K0, lin = 0;
        ulong rbase = 0;

        if (rank == 2)
        {
            double A = G[0], Bc = G[1], D = A * G[2] - Bc * Bc;
            double z2 = (double) S * u2 + td[1], cen, rr, q;
            q = Qmax - D * z2 * z2 / A;
            if (q < -1e-6 * (1 + Qmax))
                continue;
            rr = sqrt(FLINT_MAX(q, 0.0) / A) * (1 + 1e-9) + 1e-6 * (1 + fabs(z2));
            cen = -Bc * z2 / A;
            lo1 = (slong) floor((cen - rr - td[0]) / S) - 1;
            hi1 = (slong) ceil((cen + rr - td[0]) / S) + 1;
            /* the terms that do not depend on u1 */
            base = (slong) Den * ((slong) (S * S) * G[2] * u2 * u2 + 2 * (slong) S * cs[1] * u2) + K0;
            lin = (slong) Den * (2 * (slong) (S * S) * G[1] * u2);
            rbase = (N == 1) ? 0 : (ulong) (((u2 % (slong) N + (slong) N) * P->cl[b][1]) % N);
        }

        for (u1 = lo1; u1 <= hi1; u1++)
        {
            slong e = base + u1 * lin
                + (slong) Den * ((slong) (S * S) * G[0] * u1 * u1 + 2 * (slong) S * cs[0] * u1);
            ulong r, k;
            if (e < 0 || ((ulong) e - rho) % q != 0)
                flint_throw(FLINT_ERROR, "shells_products: bad exponent\n");
            if ((ulong) e > Emax)
                continue;
            k = ((ulong) e - rho) / q;
            if (k > Kmax)
                continue;
            r = (N == 1) ? 0 : (rbase + ((u1 % (slong) N + (slong) N) * P->cl[b][0]) + r0) % N;
            if (f != NULL)
                f[k * W + r]++;
            else
            {
                if (n >= cap)
                    return -1;
                te[n] = k;
                tr[n] = r;
            }
            n++;
        }
    }

    return n;
}

/* fold the slots r >= N of a product into r - N, truncate to Emax */
static void
fold(nmod_poly_t f, ulong N, ulong W, ulong Emax)
{
    slong e, r, len = f->length;
    slong E = FLINT_MIN((ulong) (len + W - 1) / W, Emax + 1);
    for (e = 0; e < E; e++)
    {
        for (r = N; r < (slong) W; r++)
        {
            slong k = e * W + r;
            if (k >= len)
                break;
            f->coeffs[k - N] = nmod_add(f->coeffs[k - N], f->coeffs[k], f->mod);
            f->coeffs[k] = 0;
        }
    }
    if (len > (slong) ((Emax + 1) * W))
        len = (Emax + 1) * W;
    f->length = len;
    _nmod_poly_normalise(f);
}

typedef struct
{
    prod_struct * P;
    const _acb_mzeta_shells_struct * T;
    ulong len, Emax;
    int ok;
}
work_struct;

static void
coset_worker(slong k, void * vw)
{
    work_struct * w = vw;
    prod_struct * P = w->P;
    const _acb_mzeta_shells_struct * T = w->T;
    slong d = T->d, nb = P->nb, b, i, nsp = 0;
    ulong N = T->N, W = slot_width(N), Emax = w->Emax, EW, cv, p, M, q = 0, rsum = 0, Kmax;
    ulong * dense[ACB_MZETA_SHELLS_MAX_DIM], * te = NULL, * tr = NULL;
    ulong qb[ACB_MZETA_SHELLS_MAX_DIM], rho[ACB_MZETA_SHELLS_MAX_DIM];
    slong K0[ACB_MZETA_SHELLS_MAX_DIM];
    double sum[ACB_MZETA_SHELLS_MAX_DIM], mx[ACB_MZETA_SHELLS_MAX_DIM], fs, fm, pm;
    fmpz * c = _fmpz_vec_init(d);
    nmod_poly_t f, g, h;
    int ok = 1;

    cv = coset_data(c, k, P, T);

    /* the exponents of the blocks are rho_b mod q: the series are in
       x^q, with the exponents (e - rho_b)/q up to Kmax */
    for (b = 0; b < nb && ok; b++)
    {
        ok = block_step(K0 + b, qb + b, c + P->off[b], P->G[b], P->rank[b], P->Den, T->S);
        q = n_gcd(q, qb[b]);
    }
    for (b = 0; b < nb && ok; b++)
    {
        rho[b] = (ulong) K0[b] % q;
        rsum += rho[b];
    }
    if (!ok)
    {
        w->ok = 0;
        _fmpz_vec_clear(c, d);
        return;
    }
    if (rsum > Emax)
    {
        /* no point of this coset has V < g len + r0 */
        _fmpz_vec_clear(c, d);
        return;
    }
    Kmax = (Emax - rsum) / q;
    EW = (Kmax + 1) * W;

    /* the series of the blocks; the last one, of rank 1, as a list of
       terms if it has few (decided after counting them) */
    for (b = 0; b < nb; b++)
        dense[b] = NULL;
    for (b = 0; b < nb && ok; b++)
    {
        slong n;
        ulong r0 = (b == 0) ? cv : 0;
        if (b == nb - 1 && b > 0 && P->rank[b] == 1)
        {
            slong cap = 2 * (slong) (sqrt((double) Emax) + 4);
            te = flint_malloc(sizeof(ulong) * cap);
            tr = flint_malloc(sizeof(ulong) * cap);
            n = block_points(NULL, te, tr, cap, P, b, c + P->off[b], 0, T, Emax,
                K0[b], q, rho[b], Kmax);
            if (n >= 0 && (double) n < 2.0 * (log2((double) EW + 1.0) + 4.0))
            {
                nsp = n;
                sum[b] = n;
                mx[b] = 1;
                continue;
            }
        }
        dense[b] = flint_calloc(EW, sizeof(ulong));
        n = block_points(dense[b], NULL, NULL, 0, P, b, c + P->off[b], r0, T, Emax,
            K0[b], q, rho[b], Kmax);
        if (n < 0)
            ok = 0;
        else
        {
            sum[b] = n;
            mx[b] = (double) _acb_mzeta_shells_table_max(dense[b], EW, 1);
        }
    }

    /* a bound for every coefficient of the partial products */
    fs = sum[0];
    fm = mx[0];
    pm = mx[0];
    for (b = 1; b < nb && ok; b++)
    {
        fm = FLINT_MIN(fs * mx[b], fm * sum[b]);
        fs = fs * sum[b];
        pm = FLINT_MAX(pm, FLINT_MAX(fm, mx[b]));
    }
    if (!(pm * P->index < ldexp(1.0, FLINT_BITS - 2)))
        ok = 0;

    if (ok)
    {
        p = n_nextprime(FLINT_MAX((ulong) pm, UWORD(1) << 20), 1);
        nmod_poly_init(f, p);
        nmod_poly_init(g, p);
        nmod_poly_init(h, p);

        nmod_poly_fit_length(f, EW);
        flint_mpn_copyi(f->coeffs, dense[0], EW);
        f->length = EW;
        _nmod_poly_normalise(f);

        for (b = 1; b < nb; b++)
        {
            if (dense[b] == NULL)
            {
                /* the sparse last block */
                slong cap = EW + W, len = 0;
                nmod_poly_fit_length(h, cap);
                flint_mpn_zero(h->coeffs, cap);
                for (i = 0; i < nsp; i++)
                {
                    slong sh = te[i] * W + tr[i], j, m = FLINT_MIN(f->length, cap - sh);
                    for (j = 0; j < m; j++)
                        h->coeffs[sh + j] = nmod_add(h->coeffs[sh + j], f->coeffs[j], f->mod);
                    if (m > 0)
                        len = FLINT_MAX(len, sh + m);
                }
                h->length = len;
                _nmod_poly_normalise(h);
            }
            else
            {
                nmod_poly_fit_length(g, EW);
                flint_mpn_copyi(g->coeffs, dense[b], EW);
                g->length = EW;
                _nmod_poly_normalise(g);
                nmod_poly_mullow(h, f, g, FLINT_MIN(f->length + g->length - 1, (slong) (EW + W)));
            }
            if (W > 1)
                fold(h, N, W, Kmax);
            else
                nmod_poly_truncate(h, Kmax + 1);
            nmod_poly_swap(f, h);
        }
        if (W > 1)
            fold(f, N, W, Kmax);

        /* the exponents Den (g M + r0) = rsum + q K */
        for (M = 0; M < w->len; M++)
        {
            ulong e = P->Den * (T->g * M + T->r0), r, K;
            if (e < rsum || (e - rsum) % q != 0)
                continue;
            K = (e - rsum) / q;
            for (r = 0; r < N; r++)
            {
                slong pos = K * W + r;
                if (pos < f->length)
                    P->C[M * N + r] += f->coeffs[pos];
            }
        }
        nmod_poly_clear(f);
        nmod_poly_clear(g);
        nmod_poly_clear(h);
    }
    else
        w->ok = 0;

    for (b = 0; b < nb; b++)
        flint_free(dense[b]);
    flint_free(te);
    flint_free(tr);
    _fmpz_vec_clear(c, d);
}

static int
prod_compute(void * data, const _acb_mzeta_shells_struct * T, ulong len)
{
    prod_struct * P = data;
    work_struct w;
    slong k;

    w.P = P;
    w.T = T;
    w.len = len;
    if (!((double) P->Den * ((double) T->g * (len - 1) + T->r0) < ldexp(1.0, FLINT_BITS - 3)))
        return 0;
    w.Emax = P->Den * (T->g * (len - 1) + T->r0);
    w.ok = 1;

    flint_free(P->C);
    P->C = flint_calloc(len * T->N, sizeof(ulong));
    P->len = len;

    /* the cosets one after the other: the large multiplications use the
       threads themselves (running the cosets in parallel instead was
       slower, the nested multiplications competing for the threads) */
    for (k = 0; k < P->index && w.ok; k++)
        coset_worker(k, &w);

    return w.ok;
}

static void
prod_block(ulong * out, const void * data, const _acb_mzeta_shells_struct * T, ulong M0, ulong n)
{
    _acb_mzeta_shells_table_block(out, ((const prod_struct *) data)->C, T->N, M0, n);
}

static ulong
prod_max_count(const void * data, const _acb_mzeta_shells_struct * T)
{
    const prod_struct * P = data;
    return _acb_mzeta_shells_table_max(P->C, P->len, T->N);
}

static void
prod_residues(char * mask, const void * data, const _acb_mzeta_shells_struct * T)
{
    const prod_struct * P = data;
    _acb_mzeta_shells_table_residues(mask, P->C, P->len, T->N);
}

const _acb_mzeta_shells_handler_struct _acb_mzeta_shells_products =
{
    "products", prod_init, prod_clear, prod_cost, prod_compute, prod_block,
    prod_max_count, prod_residues
};
