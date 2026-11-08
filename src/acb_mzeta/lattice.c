/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "arb_mat.h"
#include "impl.h"

/* (1 + sqrt(8/pi)) prod_j (1 + sqrt(2 pi)/C_jj), the constant of the theta
   tail bound of acb_theta_sum_radius */
void
_acb_mzeta_theta_tail_const(arb_t K, const arb_mat_t C, slong prec)
{
    slong j, d = arb_mat_nrows(C);
    arb_t s, t;

    arb_init(s);
    arb_init(t);
    arb_const_pi(s, prec);
    arb_mul_2exp_si(s, s, 1);
    arb_sqrt(s, s, prec);
    arb_set_ui(K, 4);
    arb_div(K, K, s, prec);
    arb_add_ui(K, K, 1, prec);
    for (j = 0; j < d; j++)
    {
        arb_div(t, s, arb_mat_entry(C, j, j), prec);
        arb_add_ui(t, t, 1, prec);
        arb_mul(K, K, t, prec);
    }
    arb_clear(s);
    arb_clear(t);
}

/* upper triangular C with C^T C = M */
int
_acb_mzeta_upper_cho(arb_mat_t C, const arb_mat_t M, slong prec)
{
    arb_mat_t L;
    int ok;

    arb_mat_init(L, arb_mat_nrows(M), arb_mat_nrows(M));
    ok = arb_mat_cho(L, M, prec);
    arb_mat_transpose(C, L);
    arb_mat_clear(L);
    return ok;
}

static void
lattice_init(_acb_mzeta_lattice_struct * L, const arb_mat_t A0, slong prec)
{
    slong d = arb_mat_nrows(A0);
    arb_mat_t G0, M;
    fmpz_mat_t U;
    arb_t t;
    int ok = 1;

    L->d = d;
    L->prec = prec;
    arb_mat_init(L->A0, d, d);
    arb_mat_init(L->A, d, d);
    arb_mat_init(L->Ainv, d, d);
    arb_mat_init(L->G, d, d);
    arb_mat_init(L->Ginv, d, d);
    arb_mat_init(L->Cd, d, d);
    arb_mat_init(L->Cr, d, d);
    arb_init(L->vol);
    arb_init(L->lam);
    arb_init(L->lam2);
    arb_init(L->Kd);
    arb_init(L->Kr);

    arb_mat_init(G0, d, d);
    arb_mat_init(M, d, d);
    fmpz_mat_init(U, d, d);
    arb_init(t);

    arb_mat_set(L->A0, A0);

    /* LLL on the Gram matrix: new basis A0 U^T, Gram matrix U G0 U^T */
    arb_mat_transpose(M, A0);
    arb_mat_mul(G0, M, A0, prec);
    arb_mat_spd_lll_reduce(U, G0, 0.99, 0.51, prec);
    arb_mat_set_fmpz_mat(M, U);
    arb_mat_transpose(M, M);
    arb_mat_mul(L->A, A0, M, prec);

    arb_mat_transpose(M, L->A);
    arb_mat_mul(L->G, M, L->A, prec);
    ok = ok && arb_mat_inv(L->Ainv, L->A, prec);
    ok = ok && arb_mat_spd_inv(L->Ginv, L->G, prec);
    arb_mat_det(L->vol, L->A, prec);
    arb_abs(L->vol, L->vol);
    ok = ok && arb_is_positive(L->vol);

    /* lam = V^(1/d) balances the two sums (exact when V = 1) */
    if (arb_is_one(L->vol))
        arb_one(L->lam);
    else
        arb_root_ui(L->lam, L->vol, d, prec);
    arb_sqr(L->lam2, L->lam, prec);

    arb_const_pi(t, prec);
    arb_div(t, t, L->lam2, prec);
    arb_mat_scalar_mul_arb(M, L->G, t, prec);
    ok = ok && _acb_mzeta_upper_cho(L->Cd, M, prec);
    arb_const_pi(t, prec);
    arb_mul(t, t, L->lam2, prec);
    arb_mat_scalar_mul_arb(M, L->Ginv, t, prec);
    ok = ok && _acb_mzeta_upper_cho(L->Cr, M, prec);

    if (ok)
    {
        _acb_mzeta_theta_tail_const(L->Kd, L->Cd, prec);
        _acb_mzeta_theta_tail_const(L->Kr, L->Cr, prec);
        ok = arb_is_finite(L->Kd) && arb_is_finite(L->Kr);
    }
    L->ok = ok;

    arb_mat_clear(G0);
    arb_mat_clear(M);
    fmpz_mat_clear(U);
    arb_clear(t);
}

static void
lattice_clear(void * obj)
{
    _acb_mzeta_lattice_struct * L = obj;
    arb_mat_clear(L->A0);
    arb_mat_clear(L->A);
    arb_mat_clear(L->Ainv);
    arb_mat_clear(L->G);
    arb_mat_clear(L->Ginv);
    arb_mat_clear(L->Cd);
    arb_mat_clear(L->Cr);
    arb_clear(L->vol);
    arb_clear(L->lam);
    arb_clear(L->lam2);
    arb_clear(L->Kd);
    arb_clear(L->Kr);
}

typedef struct
{
    const arb_mat_struct * A;
    slong prec;
}
lattice_key;

static int
lattice_match(const void * obj, const void * key)
{
    const _acb_mzeta_lattice_struct * L = obj;
    const lattice_key * k = key;
    return L->prec >= k->prec && arb_mat_nrows(L->A0) == arb_mat_nrows(k->A)
        && arb_mat_equal(L->A0, k->A);
}

/* The lattice data at a precision of at least prec + 64 bits (so that the
   derived quantities, which are inexact for most lattices, do not limit
   the accuracy); NULL if A is not a basis to that precision. */
const _acb_mzeta_lattice_struct *
_acb_mzeta_lattice_get(acb_mzeta_cache_t cache, const arb_mat_t A, slong prec)
{
    _acb_mzeta_lattice_struct * L;
    lattice_key key;

    key.A = A;
    key.prec = prec + 64;
    L = _acb_mzeta_cache_find(cache, _ACB_MZETA_CACHE_LATTICE, lattice_match, &key);
    if (L == NULL)
    {
        L = flint_malloc(sizeof(_acb_mzeta_lattice_struct));
        lattice_init(L, A, prec + 64);
        _acb_mzeta_cache_insert(cache, _ACB_MZETA_CACHE_LATTICE, L, lattice_clear);
    }
    return L->ok ? L : NULL;
}
