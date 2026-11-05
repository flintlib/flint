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
#include "fmpz_vec.h"
#include "fmpq_vec.h"
#include "mpn_extras.h"
#include "impl.h"

/* Shell counts of a shifted integral lattice with a linear character: the
   common part (S, N, g, r0) and the dispatch to the handlers.

   V(n) = F(S n - s) with s = S sigma takes the values g M + r0, M >= 0
   (g = the gcd of V(n) - V(0), computed from the values at the generators
   +-e_i and the cross terms 2 S^2 F_ij). */

#define MAX_MOD 64

static const _acb_mzeta_shells_handler_struct * const handlers[ACB_MZETA_SHELLS_NUM_HANDLERS] =
{
    &_acb_mzeta_shells_r3,
    &_acb_mzeta_shells_r3_odd,
    &_acb_mzeta_shells_products,
    &_acb_mzeta_shells_enum,
};

static ulong
fmpq_vec_den_lcm(ulong * res, const fmpq * x, slong d)
{
    fmpz_t L;
    int ok;
    fmpz_init(L);
    /* the denominators of an fmpq vector are at stride 2 */
    _fmpz_vec_lcm_strided(L, fmpq_denref(x), 2, d, NULL);
    ok = fmpz_abs_fits_ui(L) && fmpz_cmp_ui(L, MAX_MOD) <= 0;
    *res = ok ? fmpz_get_ui(L) : 0;
    fmpz_clear(L);
    return ok;
}

int
_acb_mzeta_shells_init(_acb_mzeta_shells_struct * T, const fmpz_mat_t F,
    const fmpq * sigma, const fmpq * ell)
{
    slong d = fmpz_mat_nrows(F), i, j;
    fmpz_t t, u, gg;
    fmpz * Fs;
    int ok = 1;

    if (d < 1 || d > ACB_MZETA_SHELLS_MAX_DIM)
    {
        T->d = 0;
        return 0;
    }

    T->d = d;
    T->method = 0;
    T->variant = 0;
    T->active = -1;
    T->len = 0;
    for (i = 0; i < ACB_MZETA_SHELLS_NUM_HANDLERS; i++)
        T->data[i] = NULL;
    fmpz_mat_init_set(T->F, F);
    T->s = _fmpz_vec_init(d);
    T->ell = _fmpz_vec_init(d);

    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(gg);
    Fs = _fmpz_vec_init(d);

    if (!fmpq_vec_den_lcm(&T->S, sigma, d) || !fmpq_vec_den_lcm(&T->N, ell, d))
    {
        ok = 0;
        goto cleanup;
    }

    for (i = 0; i < d; i++)
    {
        fmpz_mul_ui(t, fmpq_numref(sigma + i), T->S);
        fmpz_divexact(T->s + i, t, fmpq_denref(sigma + i));
        fmpz_mul_ui(t, fmpq_numref(ell + i), T->N);
        fmpz_divexact(t, t, fmpq_denref(ell + i));
        fmpz_mod_ui(T->ell + i, t, T->N);
    }

    /* g and r0 */
    for (i = 0; i < d; i++)
    {
        fmpz_zero(Fs + i);
        for (j = 0; j < d; j++)
            fmpz_addmul(Fs + i, fmpz_mat_entry(F, i, j), T->s + j);
    }
    fmpz_zero(gg);
    for (i = 0; i < d; i++)
    {
        /* f(+-e_i) = S^2 F_ii -+ 2 S (F s)_i */
        fmpz_mul_ui(t, fmpz_mat_entry(F, i, i), T->S * T->S);
        fmpz_mul_ui(u, Fs + i, 2 * T->S);
        fmpz_sub(u, t, u);
        fmpz_gcd(gg, gg, u);
        fmpz_mul_ui(u, Fs + i, 2 * T->S);
        fmpz_add(u, t, u);
        fmpz_gcd(gg, gg, u);
        for (j = i + 1; j < d; j++)
        {
            fmpz_mul_ui(t, fmpz_mat_entry(F, i, j), 2 * T->S * T->S);
            fmpz_gcd(gg, gg, t);
        }
    }
    /* V(0) = s^T F s */
    _fmpz_vec_dot_general(t, NULL, 0, T->s, Fs, 0, d);
    /* the shell values g M + r0 are computed in a ulong */
    if (fmpz_is_zero(gg) || !fmpz_abs_fits_ui(gg) || fmpz_bits(gg) > FLINT_BITS - 24
        || fmpz_bits(t) > FLINT_BITS - 4)
    {
        ok = 0;
        goto cleanup;
    }
    T->g = fmpz_get_ui(gg);
    T->r0 = fmpz_fdiv_ui(t, T->g);

    for (i = 0; i < ACB_MZETA_SHELLS_NUM_HANDLERS; i++)
        T->data[i] = handlers[i]->init(T);

cleanup:
    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(gg);
    _fmpz_vec_clear(Fs, d);

    if (!ok)
        _acb_mzeta_shells_clear(T);

    return ok;
}

void
_acb_mzeta_shells_clear(_acb_mzeta_shells_struct * T)
{
    slong i;
    if (T->d == 0)
        return;
    for (i = 0; i < ACB_MZETA_SHELLS_NUM_HANDLERS; i++)
        if (T->data[i] != NULL)
            handlers[i]->clear(T->data[i]);
    fmpz_mat_clear(T->F);
    _fmpz_vec_clear(T->s, T->d);
    _fmpz_vec_clear(T->ell, T->d);
    T->d = 0;
}

static double
handler_cost(const _acb_mzeta_shells_struct * T, slong i, ulong len)
{
    if (T->data[i] == NULL || (T->method != 0 && T->method != i + 1))
        return INFINITY;
    return handlers[i]->cost(T->data[i], T, len);
}

double
_acb_mzeta_shells_cost(const _acb_mzeta_shells_struct * T, ulong len)
{
    double c = INFINITY;
    slong i;
    for (i = 0; i < ACB_MZETA_SHELLS_NUM_HANDLERS; i++)
        c = FLINT_MIN(c, handler_cost(T, i, len));
    return c;
}

int
_acb_mzeta_shells_compute(_acb_mzeta_shells_struct * T, ulong len)
{
    double c[ACB_MZETA_SHELLS_NUM_HANDLERS];
    int tried[ACB_MZETA_SHELLS_NUM_HANDLERS] = {0};
    slong i, k;

    T->active = -1;
    T->len = 0;

    for (i = 0; i < ACB_MZETA_SHELLS_NUM_HANDLERS; i++)
        c[i] = handler_cost(T, i, len);

    /* the cheapest first, then the others in order of cost */
    for (k = 0; k < ACB_MZETA_SHELLS_NUM_HANDLERS; k++)
    {
        slong best = -1;
        for (i = 0; i < ACB_MZETA_SHELLS_NUM_HANDLERS; i++)
            if (!tried[i] && c[i] != INFINITY && (best == -1 || c[i] < c[best]))
                best = i;
        if (best == -1)
            break;
        tried[best] = 1;
        if (handlers[best]->compute(T->data[best], T, len))
        {
            T->active = best;
            T->len = len;
            return 1;
        }
    }

    return 0;
}

void
_acb_mzeta_shells_block(ulong * out, const _acb_mzeta_shells_struct * T, ulong M0, ulong n)
{
    handlers[T->active]->block(out, T->data[T->active], T, M0, n);
}

ulong
_acb_mzeta_shells_max_count(const _acb_mzeta_shells_struct * T)
{
    return handlers[T->active]->max_count(T->data[T->active], T);
}

void
_acb_mzeta_shells_residues(char * mask, const _acb_mzeta_shells_struct * T)
{
    if (handlers[T->active]->residues == NULL)
        memset(mask, 1, T->N);
    else
        handlers[T->active]->residues(mask, T->data[T->active], T);
}

const char *
_acb_mzeta_shells_name(const _acb_mzeta_shells_struct * T)
{
    return (T->active < 0) ? "none" : handlers[T->active]->name;
}

/* helpers for tables C[M N + r] */

void
_acb_mzeta_shells_table_block(ulong * out, const ulong * C, ulong N, ulong M0, ulong n)
{
    flint_mpn_copyi(out, C + M0 * N, n * N);
}

ulong
_acb_mzeta_shells_table_max(const ulong * C, ulong len, ulong N)
{
    ulong i, mx = 0;
    for (i = 0; i < len * N; i++)
        mx = FLINT_MAX(mx, C[i]);
    return mx;
}

void
_acb_mzeta_shells_table_residues(char * mask, const ulong * C, ulong len, ulong N)
{
    ulong M, r;
    memset(mask, 0, N);
    for (M = 0; M < len; M++)
        for (r = 0; r < N; r++)
            if (C[M * N + r] != 0)
                mask[r] = 1;
}
