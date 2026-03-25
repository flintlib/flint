/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mpoly_factor.h"
#include "fmpz_mod_mpoly_factor.h"

static void fmpz_mod_bpoly_set_fmpz_bpoly(
    fmpz_mod_bpoly_t A,
    const fmpz_bpoly_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mod_bpoly_fit_length(A, B->length, ctx);
    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_set_fmpz_poly(A->coeffs + i, B->coeffs + i, ctx);
        if (!fmpz_mod_poly_is_zero(A->coeffs + i, ctx))
            A->length = i + 1;
    }
}

static void fmpz_mod_bpoly_set_polyx(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mod_bpoly_fit_length(A, B->length, ctx);
    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_set_fmpz(A->coeffs + i, B->coeffs + i, ctx);
        if (!fmpz_mod_poly_is_zero(A->coeffs + i, ctx))
            A->length = i + 1;
    }
}

static void fmpz_mod_bpoly_set_polyy(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_bpoly_fit_length(A, 1, ctx);
    fmpz_mod_poly_set(A->coeffs + 0, B, ctx);
    A->length = !fmpz_mod_poly_is_zero(A->coeffs + 0, ctx);
}


static void fmpz_mod_bpoly_add_poly_shift(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_poly_t B,
    slong yshift,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_t c;

    FLINT_ASSERT(A->length > B->length);

    fmpz_init(c);

    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_get_coeff_fmpz(c, A->coeffs + i, yshift, ctx);
        FLINT_ASSERT(fmpz_is_zero(c));
        fmpz_mod_poly_set_coeff_fmpz(A->coeffs + i, yshift, B->coeffs + i, ctx);
    }

    fmpz_clear(c);
}



static void fmpz_bpoly_set(fmpz_bpoly_t A, const fmpz_bpoly_t B)
{
    slong i;

    FLINT_ASSERT(A != B);

    fmpz_bpoly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
        fmpz_poly_set(A->coeffs + i, B->coeffs + i);
}

static void fmpz_bpoly_make_primitive(fmpz_poly_t g, fmpz_bpoly_t A)
{
    slong Alen = A->length;
    slong i;
    fmpz_poly_t q;

    fmpz_poly_init(q);

    fmpz_poly_zero(g);
    for (i = 0; i < Alen; i++)
    {
        fmpz_poly_gcd(q, g, A->coeffs + i);
        fmpz_poly_swap(g, q);
    }

    if (Alen > 0 && fmpz_sgn(A->coeffs[Alen - 1].coeffs + A->coeffs[Alen - 1].length - 1) < 0)
        fmpz_poly_neg(g, g);

    for (i = 0; i < A->length; i++)
    {
        fmpz_poly_divexact(q, A->coeffs + i, g);
        fmpz_poly_swap(A->coeffs + i, q);
    }

    fmpz_poly_clear(q);
}

static int fmpz_bpoly_divides(fmpz_bpoly_t Q, fmpz_bpoly_t A, fmpz_bpoly_t B)
{
    slong i, qoff;
    int divides;
    fmpz_poly_t q, t;
    fmpz_bpoly_t R;

    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);
    FLINT_ASSERT(B->length > 0);

    fmpz_poly_init(q);
    fmpz_poly_init(t);
    fmpz_bpoly_init(R);
    fmpz_bpoly_set(R, A);

    Q->length = 0;

    while (R->length >= B->length)
    {
        divides = fmpz_poly_divides(q, R->coeffs + R->length - 1, B->coeffs + B->length - 1);
        if (!divides)
            goto cleanup;

        for (i = 0; i < B->length; i++)
        {
            fmpz_poly_mul(t, B->coeffs + i, q);
            fmpz_poly_sub(R->coeffs + i + R->length - B->length, R->coeffs + i + R->length - B->length, t);
        }

        qoff = R->length - B->length;

        FLINT_ASSERT(qoff >= 0);
        if (qoff >= Q->length)
        {
            fmpz_bpoly_fit_length(Q, qoff + 1);
            for (i = Q->length; i <= qoff; i++)
                fmpz_poly_zero(Q->coeffs + i);
            Q->length = qoff + 1;
        }

        fmpz_poly_set(Q->coeffs + qoff, q);

        while (R->length > 0 && fmpz_poly_is_zero(R->coeffs + R->length - 1))
            R->length--;
    }

    divides = (R->length == 0);

cleanup:

    fmpz_poly_clear(q);
    fmpz_poly_clear(t);
    fmpz_bpoly_clear(R);

    return divides;
}

static void fmpz_bpoly_set_fmpz_mod_bpoly(
    fmpz_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    fmpz_bpoly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
    {
        fmpz_poly_fit_length(A->coeffs + i, (B->coeffs + i)->length);
        (A->coeffs + i)->length = (B->coeffs + i)->length;
        _fmpz_vec_scalar_smod_fmpz((A->coeffs + i)->coeffs,
                         (B->coeffs + i)->coeffs, (B->coeffs + i)->length,
                                                    fmpz_mod_ctx_modulus(ctx));
    }
}

static void fmpz_bpoly_eval(fmpz_poly_t E, const fmpz_bpoly_t A, const fmpz_t alpha)
{
    slong i;
    fmpz_t t;

    fmpz_init(t);

    fmpz_poly_zero(E);
    for (i = A->length - 1; i >= 0; i--)
    {
        fmpz_poly_evaluate_fmpz(t, A->coeffs + i, alpha);
        fmpz_poly_set_coeff_fmpz(E, i, t);
    }

    fmpz_clear(t);
}

static void fmpz_bpoly_taylor_shift(fmpz_bpoly_t A, const fmpz_t alpha)
{
    slong i;
    for (i = A->length - 1; i >= 0; i--)
        _fmpz_poly_taylor_shift((A->coeffs + i)->coeffs, alpha, (A->coeffs + i)->length);
}

typedef struct {
    slong r; /* number of local factors */
    ulong k;
    slong lifting_prec;
    fmpz_t p;
    fmpz_t pk;
    fmpz_mod_ctx_t ctxp;
    fmpz_mod_ctx_t ctxpk;
    fmpz_mod_bpoly_t Btilde;                /* mod p^k */
    fmpz_mod_bpoly_struct * newBitilde;     /* mod p^k */
    fmpz_mod_poly_struct * P;               /* mod p^k */
    fmpz_mod_poly_struct * d;               /* mod p^k */
    fmpz_mod_poly_struct * Bitilde;         /* mod p^k */
    fmpz_mod_poly_struct * d1;              /* mod p */
    fmpz_mod_poly_struct * Bitilde1;        /* mod p */
} bpoly_info_struct;

typedef bpoly_info_struct bpoly_info_t[1];

static void bpoly_info_init(bpoly_info_t I, slong r, const fmpz_t p, ulong k)
{
    slong i;

    FLINT_ASSERT(r >= 2);

    I->r = r;

    I->lifting_prec = 0;

    I->k = k;
    fmpz_init_set(I->p, p);
    fmpz_init(I->pk);
    fmpz_pow_ui(I->pk, p, k);

    fmpz_mod_ctx_init(I->ctxp, I->p);
    fmpz_mod_ctx_init(I->ctxpk, I->pk);

    fmpz_mod_bpoly_init(I->Btilde, I->ctxpk);

    I->newBitilde = FLINT_ARRAY_ALLOC(I->r, fmpz_mod_bpoly_struct);
    I->P          = FLINT_ARRAY_ALLOC(I->r, fmpz_mod_poly_struct);
    I->d          = FLINT_ARRAY_ALLOC(I->r, fmpz_mod_poly_struct);
    I->Bitilde    = FLINT_ARRAY_ALLOC(I->r, fmpz_mod_poly_struct);
    I->d1         = FLINT_ARRAY_ALLOC(I->r, fmpz_mod_poly_struct);
    I->Bitilde1   = FLINT_ARRAY_ALLOC(I->r, fmpz_mod_poly_struct);

    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_bpoly_init(I->newBitilde + i, I->ctxpk);
        fmpz_mod_poly_init(I->P + i, I->ctxpk);
        fmpz_mod_poly_init(I->d + i, I->ctxpk);
        fmpz_mod_poly_init(I->Bitilde + i, I->ctxpk);
        fmpz_mod_poly_init(I->d1 + i, I->ctxp);
        fmpz_mod_poly_init(I->Bitilde1 + i, I->ctxp);
    }
}

static void bpoly_info_clear(bpoly_info_t I)
{
    slong i;

    fmpz_clear(I->p);
    fmpz_clear(I->pk);

    fmpz_mod_bpoly_clear(I->Btilde, I->ctxpk);

    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_bpoly_clear(I->newBitilde + i, I->ctxpk);
        fmpz_mod_poly_clear(I->P + i, I->ctxpk);
        fmpz_mod_poly_clear(I->d + i, I->ctxpk);
        fmpz_mod_poly_clear(I->Bitilde + i, I->ctxpk);
        fmpz_mod_poly_clear(I->d1 + i, I->ctxp);
        fmpz_mod_poly_clear(I->Bitilde1 + i, I->ctxp);
    }

    flint_free(I->newBitilde);
    flint_free(I->P);
    flint_free(I->d);
    flint_free(I->Bitilde);
    flint_free(I->d1);
    flint_free(I->Bitilde1);

    fmpz_mod_ctx_clear(I->ctxp);
    fmpz_mod_ctx_clear(I->ctxpk);
}

/*
    set out[i] so that
    1/(f[0]*f[1]*...*f[n-1]) = out[0]/f[0] + ... + out[n-1]/f[n-1]
*/
static int partial_fraction_coeffs(
    fmpz_mod_poly_struct * out,
    const fmpz_mod_poly_struct * f,
    slong n,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mod_poly_t num, den, a, b, g, t;

    FLINT_ASSERT(n > 1);

    fmpz_mod_poly_init(num, ctx);
    fmpz_mod_poly_init(den, ctx);
    fmpz_mod_poly_init(a, ctx);
    fmpz_mod_poly_init(b, ctx);
    fmpz_mod_poly_init(g, ctx);
    fmpz_mod_poly_init(t, ctx);

    fmpz_mod_poly_set_ui(num, 1, ctx);
    fmpz_mod_poly_mul(den, f + 0, f + 1, ctx);
    for (i = 2; i < n; i++)
        fmpz_mod_poly_mul(den, den, f + i, ctx);

    for (i = 0; i < n; i++)
    {
        fmpz_mod_poly_divrem(den, t, den, f + i, ctx);
        FLINT_ASSERT(fmpz_mod_poly_is_zero(t, ctx));
        fmpz_mod_poly_xgcd(g, a, b, f + i, den, ctx);
        if (fmpz_mod_poly_degree(g, ctx) != 0)
            return 0;
        FLINT_ASSERT(fmpz_is_one(g->coeffs + 0));
        fmpz_mod_poly_mul(t, b, num, ctx);
        fmpz_mod_poly_rem(out + i, t, f + i, ctx);
        fmpz_mod_poly_mul(t, a, num, ctx);
        fmpz_mod_poly_rem(num, t, den, ctx);
    }

    fmpz_mod_poly_clear(num, ctx);
    fmpz_mod_poly_clear(den, ctx);
    fmpz_mod_poly_clear(a, ctx);
    fmpz_mod_poly_clear(b, ctx);
    fmpz_mod_poly_clear(g, ctx);
    fmpz_mod_poly_clear(t, ctx);
    return 1;
}


static int bpoly_info_disolve(bpoly_info_t I)
{
    slong i, j;
    fmpz_t pj, t1;
    fmpz_mod_poly_t error, t, s, s1, s2;

    if (!partial_fraction_coeffs(I->d1, I->Bitilde1, I->r, I->ctxp))
        return 0;

    fmpz_init(pj);
    fmpz_init(t1);
    fmpz_mod_poly_init(error, I->ctxpk);
    fmpz_mod_poly_init(t, I->ctxpk);
    fmpz_mod_poly_init(s, I->ctxp);
    fmpz_mod_poly_init(s1, I->ctxp);
    fmpz_mod_poly_init(s2, I->ctxpk);

    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_poly_set_ui(I->P + i, 1, I->ctxpk);
        for (j = 0; j < I->r; j++)
        {
            if (i == j)
                continue;
            fmpz_mod_poly_mul(I->P + i, I->P + i, I->Bitilde + j, I->ctxpk);
        }
    }

    fmpz_mod_poly_set_ui(error, 1, I->ctxpk);
    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_poly_set(I->d + i, I->d1 + i, I->ctxpk); /* slight abuse because moduli are different */
        fmpz_mod_poly_mul(t, I->d + i, I->P + i, I->ctxpk);
        fmpz_mod_poly_sub(error, error, t, I->ctxpk);
    }

    fmpz_one(pj);
    for (j = 1; (ulong) j < I->k; j++)
    {
        fmpz_mul(pj, pj, I->p);
        fmpz_mod_poly_zero(s, I->ctxp);
        for (i = error->length - 1; i >= 0; i--)
        {
            FLINT_ASSERT(fmpz_divisible(error->coeffs + i, pj));
            fmpz_divexact(t1, error->coeffs + i, pj);
            fmpz_mod(t1, t1, I->p);
            fmpz_mod_poly_set_coeff_fmpz(s, i, t1, I->ctxp);
        }

        for (i = 0; i < I->r; i++)
        {
            fmpz_mod_poly_mul(s1, s, I->d1 + i, I->ctxp);
            fmpz_mod_poly_rem(s2, s1, I->Bitilde1 + i, I->ctxp);
            fmpz_mod_poly_scalar_mul_fmpz(s2, s2, pj, I->ctxpk);
            fmpz_mod_poly_add(I->d + i, I->d + i, s2, I->ctxpk);
        }

        fmpz_mod_poly_set_ui(error, 1, I->ctxpk);
        for (i = 0; i < I->r; i++)
        {
            fmpz_mod_poly_mul(t, I->d + i, I->P + i, I->ctxpk);
            fmpz_mod_poly_sub(error, error, t, I->ctxpk);
        }
    }

    FLINT_ASSERT(fmpz_mod_poly_is_zero(error, I->ctxpk));

    fmpz_clear(pj);
    fmpz_clear(t1);
    fmpz_mod_poly_clear(error, I->ctxpk);
    fmpz_mod_poly_clear(t, I->ctxpk);
    fmpz_mod_poly_clear(s, I->ctxp);
    fmpz_mod_poly_clear(s1, I->ctxp);
    fmpz_mod_poly_clear(s2, I->ctxpk);

    return 1;
}



static void _bivar_lift_quintic(bpoly_info_t I)
{
    slong i, j, k;
    fmpz_mod_bpoly_t tp, tp1, error;
    fmpz_mod_poly_t ss, tt;

    fmpz_mod_poly_init(ss, I->ctxpk);
    fmpz_mod_poly_init(tt, I->ctxpk);
    fmpz_mod_bpoly_init(tp, I->ctxpk);
    fmpz_mod_bpoly_init(tp1, I->ctxpk);
    fmpz_mod_bpoly_init(error, I->ctxpk);

    fmpz_mod_bpoly_mul_series(tp, I->newBitilde + 0, I->newBitilde + 1, I->lifting_prec, I->ctxpk);
    for (i = 2; i < I->r; i++)
    {
        fmpz_mod_bpoly_mul_series(tp1, tp, I->newBitilde + i, I->lifting_prec, I->ctxpk);
        fmpz_mod_bpoly_swap(tp1, tp, I->ctxpk);
    }
    fmpz_mod_bpoly_sub(error, I->Btilde, tp, I->ctxpk);

    for (j = 1; j < I->lifting_prec; j++)
    {
        fmpz_mod_poly_zero(ss, I->ctxpk);
        for (i = error->length - 1; i >= 0; i--)
        {
            fmpz_t ct;
            fmpz_init(ct);

            fmpz_mod_bpoly_get_coeff(ct, error, i, j, I->ctxpk);
            fmpz_mod_poly_set_coeff_fmpz(ss, i, ct, I->ctxpk);
            for (k = 0; k < j; k++)
            {
                fmpz_mod_bpoly_get_coeff(ct, error, i, k, I->ctxpk);
                FLINT_ASSERT(fmpz_is_zero(ct));
            }

            fmpz_clear(ct);
        }

        for (i = 0; i < I->r; i++)
        {
            fmpz_mod_poly_mul(tt, ss, I->d + i, I->ctxpk);
            fmpz_mod_poly_rem(tt, tt, I->Bitilde + i, I->ctxpk);
            fmpz_mod_bpoly_add_poly_shift(I->newBitilde + i, tt, j, I->ctxpk);
        }

        fmpz_mod_bpoly_mul_series(tp, I->newBitilde + 0, I->newBitilde + 1, I->lifting_prec, I->ctxpk);
        for (i = 2; i < I->r; i++)
        {
            fmpz_mod_bpoly_mul_series(tp1, tp, I->newBitilde + i, I->lifting_prec, I->ctxpk);
            fmpz_mod_bpoly_swap(tp1, tp, I->ctxpk);
        }
        fmpz_mod_bpoly_sub(error, I->Btilde, tp, I->ctxpk);
    }

    fmpz_mod_poly_clear(ss, I->ctxpk);
    fmpz_mod_poly_clear(tt, I->ctxpk);
    fmpz_mod_bpoly_clear(tp, I->ctxpk);
    fmpz_mod_bpoly_clear(tp1, I->ctxpk);
    fmpz_mod_bpoly_clear(error, I->ctxpk);
}


static void _bivar_lift_quartic2(bpoly_info_t I)
{
    slong i, j, k;
    fmpz_mod_poly_t t, t1;
    fmpz_mod_bpoly_t btilde;
    fmpz_mod_bpoly_struct newbitilde[2];

    FLINT_ASSERT(I->r == 2);

    fmpz_mod_poly_init(t, I->ctxpk);
    fmpz_mod_poly_init(t1, I->ctxpk);
    fmpz_mod_bpoly_init(btilde, I->ctxpk);
    fmpz_mod_bpoly_reverse_vars(btilde, I->Btilde, I->ctxpk);
    for (k = 0; k < I->r; k++)
    {
        fmpz_mod_bpoly_init(newbitilde + k, I->ctxpk);
        fmpz_mod_bpoly_reverse_vars(newbitilde + k, I->newBitilde + k, I->ctxpk);
        fmpz_mod_bpoly_fit_length(newbitilde + k, I->lifting_prec, I->ctxpk);
        FLINT_ASSERT((newbitilde + k)->length == 1);
    }

    for (j = 1; j < I->lifting_prec; j++)
    {
        if (j < btilde->length)
            fmpz_mod_poly_set(t, btilde->coeffs + j, I->ctxpk);
        else
            fmpz_mod_poly_zero(t, I->ctxpk);

        for (i = 1; i < j; i++)
        {
            fmpz_mod_poly_mul(t1, newbitilde[0].coeffs + i, newbitilde[1].coeffs + j - i, I->ctxpk);
            fmpz_mod_poly_sub(t, t, t1, I->ctxpk);
        }

        for (k = 0; k < I->r; k++)
        {
            fmpz_mod_poly_mul(t1, t, I->d + k, I->ctxpk);
            fmpz_mod_poly_rem(newbitilde[k].coeffs + j, t1, I->Bitilde + k, I->ctxpk);
            if (!fmpz_mod_poly_is_zero(newbitilde[k].coeffs + j, I->ctxpk))
                newbitilde[k].length = j + 1;
        }
    }

    for (k = 0; k < I->r; k++)
        fmpz_mod_bpoly_reverse_vars(I->newBitilde + k, newbitilde + k, I->ctxpk);

    fmpz_mod_poly_clear(t, I->ctxpk);
    fmpz_mod_poly_clear(t1, I->ctxpk);
    fmpz_mod_bpoly_clear(btilde, I->ctxpk);
    for (k = 0; k < I->r; k++)
        fmpz_mod_bpoly_clear(newbitilde + k, I->ctxpk);
}

static void _bivar_lift_quartic(bpoly_info_t I)
{
    slong i, j, k;
    fmpz_mod_poly_t t, t1;
    fmpz_mod_bpoly_t btilde;
    fmpz_mod_bpoly_struct * newbitilde, * U;

    FLINT_ASSERT(I->r > 2);

    newbitilde = (fmpz_mod_bpoly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_bpoly_struct));
    U = (fmpz_mod_bpoly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_bpoly_struct));

    fmpz_mod_poly_init(t, I->ctxpk);
    fmpz_mod_poly_init(t1, I->ctxpk);
    fmpz_mod_bpoly_init(btilde, I->ctxpk);
    fmpz_mod_bpoly_reverse_vars(btilde, I->Btilde, I->ctxpk);
    for (k = 0; k < I->r; k++)
    {
        fmpz_mod_bpoly_init(U + k, I->ctxpk);
        fmpz_mod_bpoly_fit_length(U + k, I->lifting_prec, I->ctxpk);
        for (i = 0; i < I->lifting_prec; i++)
        {
            fmpz_mod_poly_zero(U[k].coeffs + i, I->ctxpk);
        }

        fmpz_mod_bpoly_init(newbitilde + k, I->ctxpk);
        fmpz_mod_bpoly_reverse_vars(newbitilde + k, I->newBitilde + k, I->ctxpk);
        fmpz_mod_bpoly_fit_length(newbitilde + k, I->lifting_prec, I->ctxpk);
        FLINT_ASSERT(newbitilde[k].length == 1);
        for (i = 1; i < I->lifting_prec; i++)
        {
            fmpz_mod_poly_zero(newbitilde[k].coeffs + i, I->ctxpk);
        }
    }

    k = I->r - 2;
    fmpz_mod_poly_mul(U[k].coeffs + 0, newbitilde[k].coeffs + 0, newbitilde[k + 1].coeffs + 0, I->ctxpk);
    for (k--; k >= 1; k--)
        fmpz_mod_poly_mul(U[k].coeffs + 0, newbitilde[k].coeffs + 0, U[k + 1].coeffs + 0, I->ctxpk);

    for (j = 1; j < I->lifting_prec; j++)
    {
        k = I->r - 2;
        fmpz_mod_poly_zero(U[k].coeffs + j, I->ctxpk);
        for (i = 0; i <= j; i++)
        {
            fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + i, newbitilde[k + 1].coeffs + j - i, I->ctxpk);
            fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t1, I->ctxpk);
        }
        for (k--; k >= 1; k--)
        {
            fmpz_mod_poly_zero(U[k].coeffs + j, I->ctxpk);
            for (i = 0; i <= j; i++)
            {
                fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + i, U[k + 1].coeffs + j - i, I->ctxpk);
                fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t1, I->ctxpk);
            }
        }

        if (j < btilde->length)
            fmpz_mod_poly_set(t, btilde->coeffs + j, I->ctxpk);
        else
            fmpz_mod_poly_zero(t, I->ctxpk);

        for (i = 0; i <= j; i++)
        {
            fmpz_mod_poly_mul(t1, newbitilde[0].coeffs + i, U[1].coeffs + j - i, I->ctxpk);
            fmpz_mod_poly_sub(t, t, t1, I->ctxpk);
        }

        for (k = 0; k < I->r; k++)
        {
            fmpz_mod_poly_mul(t1, t, I->d + k, I->ctxpk);
            fmpz_mod_poly_rem(newbitilde[k].coeffs + j, t1, I->Bitilde + k, I->ctxpk);
            if (!fmpz_mod_poly_is_zero(newbitilde[k].coeffs + j, I->ctxpk))
                newbitilde[k].length = j + 1;
        }

        k = I->r - 2;
        fmpz_mod_poly_mul(t, newbitilde[k].coeffs + 0, newbitilde[k + 1].coeffs + j, I->ctxpk);
        fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + j, newbitilde[k + 1].coeffs + 0, I->ctxpk);
        fmpz_mod_poly_add(t, t, t1, I->ctxpk);
        fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t, I->ctxpk);
        for (k--; k >= 1; k--)
        {
            fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + 0, t, I->ctxpk);
            fmpz_mod_poly_swap(t, t1, I->ctxpk);
            fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + j, U[k + 1].coeffs + 0, I->ctxpk);
            fmpz_mod_poly_add(t, t, t1, I->ctxpk);
            fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t, I->ctxpk);
        }
    }

    for (k = 0; k < I->r; k++)
        fmpz_mod_bpoly_reverse_vars(I->newBitilde + k, newbitilde + k, I->ctxpk);

    fmpz_mod_poly_clear(t, I->ctxpk);
    fmpz_mod_poly_clear(t1, I->ctxpk);
    fmpz_mod_bpoly_clear(btilde, I->ctxpk);
    for (k = 0; k < I->r; k++)
    {
        fmpz_mod_bpoly_clear(U + k, I->ctxpk);
        fmpz_mod_bpoly_clear(newbitilde + k, I->ctxpk);
    }

    flint_free(newbitilde);
    flint_free(U);
}


static void _recombine_naive(
    fmpz_tpoly_t F,
    fmpz_bpoly_t B,
    fmpz_t alpha,
    bpoly_info_t I)
{
    fmpz_poly_t g;
    fmpz_bpoly_t Q, R, trymez;
    fmpz_mod_bpoly_t tryme, trymet;
    fmpz_mod_poly_t leadB;
    slong i, k, len;
    slong * subset;

    fmpz_poly_init(g);
    fmpz_bpoly_init(Q);
    fmpz_bpoly_init(R);
    fmpz_bpoly_init(trymez);
    fmpz_mod_bpoly_init(tryme, I->ctxpk);
    fmpz_mod_bpoly_init(trymet, I->ctxpk);
    fmpz_mod_poly_init(leadB, I->ctxpk);

    F->length = 0;

    FLINT_ASSERT(B->length > 0);
    fmpz_mod_poly_set_fmpz_poly(leadB, B->coeffs + B->length - 1, I->ctxpk);

    len = I->r;
    subset = (slong *) flint_malloc(len * sizeof(slong));
    for (k = 0; k < len; k++)
        subset[k] = k;

    for (k = 1; k <= len/2; k++)
    {
        zassenhaus_subset_first(subset, len, k);
        while (1)
        {
            fmpz_mod_bpoly_set_polyy(tryme, leadB, I->ctxpk);
            for (i = 0; i < len; i++)
            {
                if (subset[i] >= 0)
                {
                    fmpz_mod_bpoly_mul_series(trymet, tryme, I->newBitilde + subset[i],
                                                    I->lifting_prec, I->ctxpk);
                    fmpz_mod_bpoly_swap(trymet, tryme, I->ctxpk);
                }
            }
            fmpz_bpoly_set_fmpz_mod_bpoly(trymez, tryme, I->ctxpk);
            fmpz_bpoly_make_primitive(g, trymez);

            if (fmpz_bpoly_divides(Q, B, trymez))
            {
                fmpz_neg(alpha, alpha);
                fmpz_bpoly_taylor_shift(trymez, alpha);
                fmpz_neg(alpha, alpha);
                fmpz_tpoly_fit_length(F, F->length + 1);
                fmpz_bpoly_swap(F->coeffs + F->length, trymez);
                F->length++;
                fmpz_bpoly_swap(Q, B);
                FLINT_ASSERT(B->length > 0);
                fmpz_mod_poly_set_fmpz_poly(leadB, B->coeffs + B->length - 1, I->ctxpk);

                len -= k;
                if (!zassenhaus_subset_next_disjoint(subset, len + k))
                    break;
            }
            else
            {
                if (!zassenhaus_subset_next(subset, len))
                    break;
            }
        }
    }

    fmpz_neg(alpha, alpha);
    fmpz_bpoly_taylor_shift(B, alpha);
    fmpz_neg(alpha, alpha);
    if (B->length > 1)
    {
        fmpz_tpoly_fit_length(F, F->length + 1);
        fmpz_bpoly_swap(F->coeffs + F->length, B);
        F->length++;
    }
    else
    {
        FLINT_ASSERT(B->length == 1);
        FLINT_ASSERT(fmpz_poly_is_one(B->coeffs + 0));
    }

    fmpz_poly_clear(g);
    fmpz_bpoly_clear(Q);
    fmpz_bpoly_clear(R);
    fmpz_bpoly_clear(trymez);
    fmpz_mod_bpoly_clear(tryme, I->ctxpk);
    fmpz_mod_bpoly_clear(trymet, I->ctxpk);
    fmpz_mod_poly_clear(leadB, I->ctxpk);

    flint_free(subset);
}


void fmpz_bpoly_factor(fmpz_poly_t c, fmpz_tpoly_t F, fmpz_bpoly_t B)
{
    slong i;
    fmpz_t alpha;
    fmpz_poly_t Beval;
    fmpz_poly_factor_t Bevalfac;
    slong Blengthx, Blengthy;
    flint_bitcnt_t Bbits;
    slong sumabs, maxabs;
    ulong pkbits;
    ulong k;
    fmpz_t p;
    bpoly_info_t I;

    k = 1;
    fmpz_init_set_ui(p, UWORD(1) << (FLINT_BITS - 1));
    fmpz_init(alpha);
    fmpz_poly_init(Beval);
    fmpz_poly_factor_init(Bevalfac);
    bpoly_info_init(I, 2, p, k);

    Blengthx = B->length;
    FLINT_ASSERT(Blengthx > 1);

    fmpz_bpoly_make_primitive(c, B);


    /* New fast path for sparse 4-term bivariate cases */
    if (fmpz_bpoly_factor_try_sparse_4term_pairs(c, F, B))
        goto cleanup;

    fmpz_zero(alpha);
    goto got_alpha;

next_alpha:

    fmpz_neg(alpha, alpha);
    fmpz_add_ui(alpha, alpha, fmpz_sgn(alpha) >= 0);

got_alpha:

    fmpz_bpoly_eval(Beval, B, alpha);

    /* if killed leading coeff, get new alpha */
    if (Beval->length != Blengthx)
        goto next_alpha;

    fmpz_one(&Bevalfac->c);
    fmpz_poly_factor(Bevalfac, Beval);

    /* if multiple factors, get new alpha */
    for (i = 0; i < Bevalfac->num; i++)
    {
        if (Bevalfac->exp[i] != 1)
            goto next_alpha;
    }

    /* if one factor, B is irreducible */
    if (Bevalfac->num < 2)
    {
        fmpz_tpoly_fit_length(F, 1);
        F->length = 1;
        fmpz_bpoly_swap(F->coeffs + 0, B);
        goto cleanup;
    }

    fmpz_bpoly_taylor_shift(B, alpha);

    Blengthy = 0;
    Bbits = 0;
    for (i = 0; i < B->length; i++)
    {
        slong this_bits;
        Blengthy = FLINT_MAX(Blengthy, B->coeffs[i].length);
        this_bits = _fmpz_vec_max_bits(B->coeffs[i].coeffs, B->coeffs[i].length);
        Bbits = FLINT_MAX(Bbits, FLINT_ABS(this_bits));
    }

    pkbits = (FLINT_BIT_COUNT(Blengthx*Blengthy) + 1)/2;
    pkbits += Blengthx + Blengthy + Bbits - 3;

next_prime:

    fmpz_nextprime(p, p, 1);

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT((B->coeffs + B->length - 1)->length > 0);
    FLINT_ASSERT(!fmpz_is_zero((B->coeffs + B->length - 1)->coeffs + 0));

    if (fmpz_divisible((B->coeffs + B->length - 1)->coeffs + 0, p))
        goto next_prime;

    /*
        2^pkbits is strict upperbound on the coefficient of any factor of B.
        An upperbound on the coefficient of any factor of B*lc_x(B) is needed.
    */
    _fmpz_vec_sum_max_bits(&sumabs, &maxabs, B->coeffs[B->length - 1].coeffs,
                                             B->coeffs[B->length - 1].length);

    k = (pkbits + sumabs + fmpz_bits(p))/fmpz_bits(p);

    bpoly_info_clear(I);
    bpoly_info_init(I, Bevalfac->num, p, k);
    I->lifting_prec = Blengthy + (B->coeffs + B->length - 1)->length;

    fmpz_mod_bpoly_set_fmpz_bpoly(I->Btilde, B, I->ctxpk);
    fmpz_mod_bpoly_make_monic_series(I->Btilde, I->Btilde, I->lifting_prec, I->ctxpk);
    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_poly_set_fmpz_poly(I->Bitilde1 + i, Bevalfac->p + i, I->ctxp);
        fmpz_mod_poly_make_monic(I->Bitilde1 + i, I->Bitilde1 + i, I->ctxp);
        fmpz_mod_poly_set_fmpz_poly(I->Bitilde + i, Bevalfac->p + i, I->ctxpk);
        fmpz_mod_poly_make_monic(I->Bitilde + i, I->Bitilde + i, I->ctxpk);
        fmpz_mod_bpoly_set_polyx(I->newBitilde + i, I->Bitilde + i, I->ctxpk);
    }

    FLINT_ASSERT(I->r > 1);

    if (!bpoly_info_disolve(I))
        goto next_prime;

    if (I->r == 2)
        _bivar_lift_quartic2(I);
    else if (I->r < 20)
        _bivar_lift_quartic(I);
    else
        _bivar_lift_quintic(I);

    _recombine_naive(F, B, alpha, I);

cleanup:

    bpoly_info_clear(I);
    fmpz_poly_factor_clear(Bevalfac);
    fmpz_poly_clear(Beval);
    fmpz_clear(alpha);
    fmpz_clear(p);
}



/*
    return 1: ok
           0: lift failed
          -1: not primitive
*/
int fmpz_bpoly_factor_ordered(
    fmpz_poly_t c,
    fmpz_tpoly_t F,
    fmpz_bpoly_t B,
    const fmpz_t alpha,
    const fmpz_poly_factor_t Bevalf)
{
    int success;
    slong i;
    slong Blengthx, Blengthy;
    flint_bitcnt_t Bbits;
    ulong pkbits;
    ulong k;
    slong sumabs, maxabs;
    fmpz_t p, malpha;
    bpoly_info_t I;
    fmpz_bpoly_t Q, trymez;
    fmpz_mod_bpoly_t tryme, trymet;
    fmpz_mod_poly_t Blead;
    fmpz_poly_t g;

    k = 1;
    fmpz_init_set_ui(p, UWORD(1) << (FLINT_BITS - 1));
    bpoly_info_init(I, 2, p, k);

    fmpz_poly_init(g);
    fmpz_bpoly_init(Q);
    fmpz_bpoly_init(trymez);
    fmpz_mod_bpoly_init(tryme, I->ctxpk);
    fmpz_mod_bpoly_init(trymet, I->ctxpk);
    fmpz_mod_poly_init(Blead, I->ctxpk);

    Blengthx = B->length;
    FLINT_ASSERT(Blengthx > 1);

    fmpz_init(malpha);

    fmpz_bpoly_make_primitive(c, B);
    if (fmpz_poly_degree(c) > 0)
    {
        success = -1;
        goto cleanup;
    }

    fmpz_neg(malpha, alpha);
    fmpz_bpoly_taylor_shift(B, alpha);

    Blengthy = 0;
    Bbits = 0;
    for (i = 0; i < B->length; i++)
    {
        slong this_bits;
        Blengthy = FLINT_MAX(Blengthy, B->coeffs[i].length);
        this_bits = _fmpz_vec_max_bits(B->coeffs[i].coeffs, B->coeffs[i].length);
        Bbits = FLINT_MAX(Bbits, FLINT_ABS(this_bits));
    }

    pkbits = (FLINT_BIT_COUNT(Blengthx*Blengthy) + 1)/2;
    pkbits += Blengthx + Blengthy + Bbits - 3;

next_prime:

    fmpz_nextprime(p, p, 1);

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT((B->coeffs + B->length - 1)->length > 0);
    FLINT_ASSERT(!fmpz_is_zero((B->coeffs + B->length - 1)->coeffs + 0));

    if (fmpz_divisible((B->coeffs + B->length - 1)->coeffs + 0, p))
        goto next_prime;

    _fmpz_vec_sum_max_bits(&sumabs, &maxabs, B->coeffs[B->length - 1].coeffs,
                                             B->coeffs[B->length - 1].length);

    k = (pkbits + sumabs + fmpz_bits(p))/fmpz_bits(p);

    bpoly_info_clear(I);
    bpoly_info_init(I, Bevalf->num, p, k);
    I->lifting_prec = Blengthy + (B->coeffs + B->length - 1)->length;

    fmpz_mod_bpoly_set_fmpz_bpoly(I->Btilde, B, I->ctxpk);
    fmpz_mod_bpoly_make_monic_series(I->Btilde, I->Btilde, I->lifting_prec, I->ctxpk);
    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_poly_set_fmpz_poly(I->Bitilde1 + i, Bevalf->p + i, I->ctxp);
        fmpz_mod_poly_make_monic(I->Bitilde1 + i, I->Bitilde1 + i, I->ctxp);
        fmpz_mod_poly_set_fmpz_poly(I->Bitilde + i, Bevalf->p + i, I->ctxpk);
        fmpz_mod_poly_make_monic(I->Bitilde + i, I->Bitilde + i, I->ctxpk);
        fmpz_mod_bpoly_set_polyx(I->newBitilde + i, I->Bitilde + i, I->ctxpk);
    }

    FLINT_ASSERT(I->r > 1);

    if (!bpoly_info_disolve(I))
        goto next_prime;

    if (I->r == 2)
        _bivar_lift_quartic2(I);
    else if (I->r < 20)
        _bivar_lift_quartic(I);
    else
        _bivar_lift_quintic(I);

    fmpz_tpoly_fit_length(F, I->r);
    F->length = 0;
    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_poly_set_fmpz_poly(Blead, B->coeffs + B->length - 1, I->ctxpk);
        fmpz_mod_bpoly_set_polyy(tryme, Blead, I->ctxpk);
        fmpz_mod_bpoly_mul_series(trymet, tryme, I->newBitilde + i, I->lifting_prec, I->ctxpk);
        fmpz_mod_bpoly_swap(trymet, tryme, I->ctxpk);
        fmpz_bpoly_set_fmpz_mod_bpoly(trymez, tryme, I->ctxpk);
        fmpz_bpoly_make_primitive(g, trymez);
        if (!fmpz_bpoly_divides(Q, B, trymez))
        {
            success = 0;
            goto cleanup;
        }
        fmpz_bpoly_swap(B, Q);
        fmpz_bpoly_taylor_shift(trymez, malpha);
        fmpz_bpoly_swap(F->coeffs + F->length, trymez);
        F->length++;
    }

    success = 1;

cleanup:

    fmpz_poly_clear(g);
    fmpz_bpoly_clear(Q);
    fmpz_bpoly_clear(trymez);
    fmpz_mod_bpoly_clear(tryme, I->ctxpk);
    fmpz_mod_bpoly_clear(trymet, I->ctxpk);
    fmpz_mod_poly_clear(Blead, I->ctxpk);

    bpoly_info_clear(I);
    fmpz_clear(malpha);
    fmpz_clear(p);

    return success;
}

/* -------------------------------------------------------------------------- */
/*  Fast sparse special cases for bivariate factorization                     */
/* -------------------------------------------------------------------------- */

typedef struct
{
    fmpz_t c;
    slong ex;
    slong ey;
} sparse_term2_struct;

typedef sparse_term2_struct sparse_term2_t[1];

/* Pairings: (0,1)|(2,3), (0,2)|(1,3), (0,3)|(1,2) */
static const int fmpz_bpoly_factor_pairings_4term[3][4] =
{
    {0, 1, 2, 3},
    {0, 2, 1, 3},
    {0, 3, 1, 2}
};

static void fmpz_bpoly_sparse_term_init(sparse_term2_t t)
{
    t->ex = 0;
    t->ey = 0;
    fmpz_init(t->c);
}

static void fmpz_bpoly_sparse_term_clear(sparse_term2_t t)
{
    fmpz_clear(t->c);
}

static int fmpz_bpoly_equal_exact(const fmpz_bpoly_t A, const fmpz_bpoly_t B)
{
    slong i;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!fmpz_poly_equal(A->coeffs + i, B->coeffs + i))
            return 0;
    }

    return 1;
}

static void fmpz_bpoly_trim_exact(fmpz_bpoly_t A)
{
    while (A->length > 0 && fmpz_poly_is_zero(A->coeffs + A->length - 1))
        A->length--;
}

static void fmpz_bpoly_zero_exact(fmpz_bpoly_t A)
{
    slong i;

    for (i = 0; i < A->length; i++)
        fmpz_poly_zero(A->coeffs + i);

    A->length = 0;
}

static void fmpz_bpoly_add_monomial_fmpz(fmpz_bpoly_t A, const fmpz_t coeff, slong xexp, slong yexp)
{
    slong i;
    fmpz_t old;

    FLINT_ASSERT(xexp >= 0);
    FLINT_ASSERT(yexp >= 0);

    if (fmpz_is_zero(coeff))
        return;

    if (yexp >= A->length)
    {
        fmpz_bpoly_fit_length(A, yexp + 1);
        for (i = A->length; i <= yexp; i++)
            fmpz_poly_zero(A->coeffs + i);
        A->length = yexp + 1;
    }

    fmpz_init(old);
    fmpz_poly_get_coeff_fmpz(old, A->coeffs + yexp, xexp);
    fmpz_add(old, old, coeff);
    fmpz_poly_set_coeff_fmpz(A->coeffs + yexp, xexp, old);
    fmpz_clear(old);

    fmpz_bpoly_trim_exact(A);
}

static slong fmpz_bpoly_term_count_exact(const fmpz_bpoly_t A)
{
    slong i, total = 0;
    for (i = 0; i < A->length; i++)
        total += A->coeffs[i].length;
    return total;
}

static int fmpz_bpoly_extract_terms_4(sparse_term2_t out[4], const fmpz_bpoly_t A)
{
    slong i, j, k = 0;

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < A->coeffs[i].length; j++)
        {
            if (fmpz_is_zero(A->coeffs[i].coeffs + j))
                continue;

            if (k >= 4)
                return 0;

            out[k][0].ex = j;
            out[k][0].ey = i;
            fmpz_set(out[k][0].c, A->coeffs[i].coeffs + j);
            k++;
        }
    }

    return (k == 4);
}

static int fmpz_bpoly_extract_terms_2(sparse_term2_t out[2], const fmpz_bpoly_t A)
{
    slong i, j, k = 0;

    for (i = 0; i < A->length; i++)
    {
        for (j = 0; j < A->coeffs[i].length; j++)
        {
            if (fmpz_is_zero(A->coeffs[i].coeffs + j))
                continue;

            if (k >= 2)
                return 0;

            out[k][0].ex = j;
            out[k][0].ey = i;
            fmpz_set(out[k][0].c, A->coeffs[i].coeffs + j);
            k++;
        }
    }

    return (k == 2);
}

static int fmpz_bpoly_leading_coeff_positive(const fmpz_bpoly_t A)
{
    const fmpz_poly_struct * p;
    FLINT_ASSERT(A->length > 0);
    p = A->coeffs + (A->length - 1);
    FLINT_ASSERT(p->length > 0);
    return fmpz_sgn(p->coeffs + (p->length - 1)) > 0;
}

/*
Normalize a 2-term polynomial
    t0 + t1
into
    M * R
where M is a monomial with integer coefficient and R is primitive 2-term.
*/
static int fmpz_bpoly_normalize_pair(fmpz_bpoly_t M, fmpz_bpoly_t R, const sparse_term2_t t0, const sparse_term2_t t1)
{
    slong minx, miny;
    slong ax, ay, bx, by;
    fmpz_t g, a0, a1;
    slong i;

    minx = FLINT_MIN(t0->ex, t1->ex);
    miny = FLINT_MIN(t0->ey, t1->ey);

    ax = t0->ex - minx;
    ay = t0->ey - miny;
    bx = t1->ex - minx;
    by = t1->ey - miny;

    if (ax == bx && ay == by)
        return 0;

    fmpz_init(g);
    fmpz_init(a0);
    fmpz_init(a1);

    fmpz_gcd(g, t0->c, t1->c);
    if (fmpz_is_zero(g))
    {
        fmpz_clear(g);
        fmpz_clear(a0);
        fmpz_clear(a1);
        return 0;
    }

    fmpz_divexact(a0, t0->c, g);
    fmpz_divexact(a1, t1->c, g);

    fmpz_bpoly_zero_exact(M);
    fmpz_bpoly_zero_exact(R);

    fmpz_bpoly_add_monomial_fmpz(M, g, minx, miny);
    fmpz_bpoly_add_monomial_fmpz(R, a0, ax, ay);
    fmpz_bpoly_add_monomial_fmpz(R, a1, bx, by);

    if (!fmpz_bpoly_leading_coeff_positive(R))
    {
        fmpz_poly_neg(M->coeffs + 0, M->coeffs + 0);
        for (i = 0; i < R->length; i++)
            fmpz_poly_neg(R->coeffs + i, R->coeffs + i);
    }

    fmpz_clear(g);
    fmpz_clear(a0);
    fmpz_clear(a1);
    return 1;
}

static void fmpz_bpoly_append_factor(fmpz_tpoly_t F, fmpz_bpoly_t A)
{
    fmpz_tpoly_fit_length(F, F->length + 1);
    fmpz_bpoly_swap(F->coeffs + F->length, A);
    F->length++;
}

static void fmpz_bpoly_append_univar_factorization(fmpz_poly_t c, fmpz_tpoly_t F, fmpz_bpoly_t A)
{
    fmpz_poly_factor_t uf;
    slong i, j, e;

    FLINT_ASSERT(A->length == 1);

    fmpz_poly_factor_init(uf);
    fmpz_poly_factor(uf, A->coeffs + 0);

    fmpz_poly_mul(c, c, &uf->c);

    for (i = 0; i < uf->num; i++)
    {
        e = uf->exp[i];
        for (j = 0; j < e; j++)
        {
            fmpz_bpoly_t T;
            fmpz_bpoly_init(T);
            fmpz_bpoly_fit_length(T, 1);
            fmpz_poly_set(T->coeffs + 0, uf->p + i);
            T->length = 1;
            fmpz_bpoly_append_factor(F, T);
            fmpz_bpoly_clear(T);
        }
    }

    fmpz_poly_factor_clear(uf);
}

/*
Fast factorization of A + B*m^d where m = x^u y^v and one term is constant.
This handles factors like 1 + x^66 y^168.
*/
static int fmpz_bpoly_factor_try_binomial_monomial(fmpz_poly_t c, fmpz_tpoly_t F, fmpz_bpoly_t B)
{
    sparse_term2_t t[2];
    fmpz_poly_t U;
    fmpz_poly_factor_t Uf;
    fmpz_bpoly_t T;
    slong i, j, k;
    slong dx, dy, d;

    if (fmpz_bpoly_term_count_exact(B) != 2)
        return 0;

    fmpz_bpoly_sparse_term_init(t[0]);
    fmpz_bpoly_sparse_term_init(t[1]);

    if (!fmpz_bpoly_extract_terms_2(t, B))
        goto fail;

    if (!(t[0]->ex == 0 && t[0]->ey == 0) && (t[1]->ex == 0 && t[1]->ey == 0))
    {
        sparse_term2_struct tmp = t[0][0];
        t[0][0] = t[1][0];
        t[1][0] = tmp;
    }

    if (!(t[0]->ex == 0 && t[0]->ey == 0))
        goto fail;

    dx = t[1]->ex;
    dy = t[1]->ey;
    d = n_gcd(dx, dy);

    if (d <= 0)
        goto fail;

    fmpz_poly_init(U);
    fmpz_poly_factor_init(Uf);
    fmpz_bpoly_init(T);

    fmpz_poly_zero(U);
    fmpz_poly_set_coeff_fmpz(U, 0, t[0]->c);
    fmpz_poly_set_coeff_fmpz(U, d, t[1]->c);

    fmpz_poly_factor(Uf, U);
    fmpz_poly_mul(c, c, &Uf->c);

    for (i = 0; i < Uf->num; i++)
    {
        for (j = 0; j < Uf->exp[i]; j++)
        {
            fmpz_bpoly_zero_exact(T);

            for (k = 0; k < Uf->p[i].length; k++)
            {
                if (fmpz_is_zero(Uf->p[i].coeffs + k))
                    continue;

                fmpz_bpoly_add_monomial_fmpz(
                    T,
                    Uf->p[i].coeffs + k,
                    (dx / d) * k,
                    (dy / d) * k);
            }

            fmpz_bpoly_append_factor(F, T);
        }
    }

    fmpz_bpoly_clear(T);
    fmpz_poly_factor_clear(Uf);
    fmpz_poly_clear(U);

    fmpz_bpoly_sparse_term_clear(t[0]);
    fmpz_bpoly_sparse_term_clear(t[1]);
    return 1;

fail:
    fmpz_bpoly_sparse_term_clear(t[0]);
    fmpz_bpoly_sparse_term_clear(t[1]);
    return 0;
}

static void fmpz_bpoly_append_factorization(fmpz_poly_t c, fmpz_tpoly_t F, fmpz_bpoly_t A)
{
    if (A->length == 1)
    {
        fmpz_bpoly_append_univar_factorization(c, F, A);
        return;
    }

    if (fmpz_bpoly_factor_try_binomial_monomial(c, F, A))
        return;

    {
        fmpz_poly_t c2;
        fmpz_tpoly_t F2;
        slong i;

        fmpz_poly_init(c2);
        fmpz_tpoly_init(F2);

        fmpz_bpoly_factor(c2, F2, A);

        fmpz_poly_mul(c, c, c2);

        for (i = 0; i < F2->length; i++)
            fmpz_bpoly_append_factor(F, F2->coeffs + i);

        fmpz_poly_clear(c2);
        fmpz_tpoly_clear(F2);
    }
}

/*
Try the sparse 4-term "two pairs share a primitive binomial" pattern.
Assumes B has already been primitive-normalized by fmpz_bpoly_make_primitive(c, B).
Returns 1 if handled, 0 otherwise.
*/
static int fmpz_bpoly_factor_try_sparse_4term_pairs(fmpz_poly_t c, fmpz_tpoly_t F, fmpz_bpoly_t B)
{
    sparse_term2_t t[4];
    int i, pidx;

    if (fmpz_bpoly_term_count_exact(B) != 4)
        return 0;

    for (i = 0; i < 4; i++)
        fmpz_bpoly_sparse_term_init(t[i]);

    if (!fmpz_bpoly_extract_terms_4(t, B))
        goto fail;

    for (pidx = 0; pidx < 3; pidx++)
    {
        int a = fmpz_bpoly_factor_pairings_4term[pidx][0];
        int b = fmpz_bpoly_factor_pairings_4term[pidx][1];
        int cidx = fmpz_bpoly_factor_pairings_4term[pidx][2];
        int d = fmpz_bpoly_factor_pairings_4term[pidx][3];

        fmpz_bpoly_t M1, M2, G1, G2, G, H, Q;
        int ok = 0;

        fmpz_bpoly_init(M1);
        fmpz_bpoly_init(M2);
        fmpz_bpoly_init(G1);
        fmpz_bpoly_init(G2);
        fmpz_bpoly_init(G);
        fmpz_bpoly_init(H);
        fmpz_bpoly_init(Q);

        if (!fmpz_bpoly_normalize_pair(M1, G1, t[a], t[b]))
            goto cleanup_pair;

        if (!fmpz_bpoly_normalize_pair(M2, G2, t[cidx], t[d]))
            goto cleanup_pair;

        if (!fmpz_bpoly_equal_exact(G1, G2))
            goto cleanup_pair;

        fmpz_bpoly_set(G, G1);

        fmpz_bpoly_zero_exact(H);
        {
            slong y1 = 0, y2 = 0;
            slong x1, x2;

            while (y1 < M1->length && fmpz_poly_is_zero(M1->coeffs + y1))
                y1++;
            while (y2 < M2->length && fmpz_poly_is_zero(M2->coeffs + y2))
                y2++;

            FLINT_ASSERT(y1 < M1->length);
            FLINT_ASSERT(y2 < M2->length);

            x1 = M1->coeffs[y1].length - 1;
            x2 = M2->coeffs[y2].length - 1;

            fmpz_bpoly_add_monomial_fmpz(H, M1->coeffs[y1].coeffs + x1, x1, y1);
            fmpz_bpoly_add_monomial_fmpz(H, M2->coeffs[y2].coeffs + x2, x2, y2);
        }

        if (!fmpz_bpoly_divides(Q, B, G))
            goto cleanup_pair;

        if (!fmpz_bpoly_equal_exact(Q, H))
            goto cleanup_pair;

        F->length = 0;
        fmpz_bpoly_append_factorization(c, F, G);
        fmpz_bpoly_append_factorization(c, F, H);
        ok = 1;

cleanup_pair:
        fmpz_bpoly_clear(M1);
        fmpz_bpoly_clear(M2);
        fmpz_bpoly_clear(G1);
        fmpz_bpoly_clear(G2);
        fmpz_bpoly_clear(G);
        fmpz_bpoly_clear(H);
        fmpz_bpoly_clear(Q);

        if (ok)
        {
            for (i = 0; i < 4; i++)
                fmpz_bpoly_sparse_term_clear(t[i]);
            return 1;
        }
    }

fail:
    for (i = 0; i < 4; i++)
        fmpz_bpoly_sparse_term_clear(t[i]);

    return 0;
}
