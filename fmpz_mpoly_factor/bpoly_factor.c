/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "nmod_mpoly_factor.h"


typedef struct
{
    fmpz_mod_poly_struct * coeffs;
    slong alloc;
    slong length;
    fmpz_t modulus;
} fmpz_mod_bpoly_struct;

typedef fmpz_mod_bpoly_struct fmpz_mod_bpoly_t[1];


void fmpz_mod_bpoly_init(fmpz_mod_bpoly_t A, const fmpz_t modulus)
{
    A->coeffs = NULL;
    A->alloc = 0;
    A->length = 0;
    fmpz_init_set(A->modulus, modulus);
}

void fmpz_mod_bpoly_clear(fmpz_mod_bpoly_t A)
{
    fmpz_clear(A->modulus);
    if (A->coeffs)
    {
        slong i;
        for (i = 0; i < A->alloc; i++)
            fmpz_mod_poly_clear(A->coeffs + i);
        flint_free(A->coeffs);
    }
}

void fmpz_mod_bpoly_swap(fmpz_mod_bpoly_t A, fmpz_mod_bpoly_t B)
{
    fmpz_mod_bpoly_struct t = *A;
    *A = *B;
    *B = t;
}

void fmpz_mod_bpoly_print_pretty(fmpz_mod_bpoly_t A, const char * xvar, const char * yvar)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (fmpz_mod_poly_is_zero(A->coeffs + i))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        fmpz_mod_poly_print_pretty(A->coeffs + i, yvar);
        flint_printf(")*%s^%wd", xvar, i);
    }

    if (first)
        flint_printf("0");
}

void fmpz_mod_bpoly_fit_length(fmpz_mod_bpoly_t A, slong len)
{
    slong i;

    if (len <= A->alloc)
        return;

    if (len < 2 * A->alloc)
        len = 2 * A->alloc;

    if (A->alloc == 0)
        A->coeffs = (fmpz_mod_poly_struct *) flint_malloc(len * sizeof(fmpz_mod_poly_struct));
    else
        A->coeffs = (fmpz_mod_poly_struct *) flint_realloc(A->coeffs, len * sizeof(fmpz_mod_poly_struct));

    for (i = A->alloc; i < len; i++)
        fmpz_mod_poly_init(A->coeffs + i, A->modulus);

    A->alloc = len;
}

void fmpz_tpoly_print(fmpz_tpoly_t A, const char * xvar, const char * yvar, const char * zvar)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        fmpz_bpoly_print_pretty(A->coeffs + i, yvar, zvar);
        flint_printf(")*%s^%wd", xvar, i);
    }

    if (first)
        flint_printf("0");
}



void fmpz_mod_bpoly_set_fmpz_bpoly(fmpz_mod_bpoly_t A, const fmpz_bpoly_t B)
{
    slong i;
    fmpz_mod_bpoly_fit_length(A, B->length);
    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_set_fmpz_poly(A->coeffs + i, B->coeffs + i);
        if (!fmpz_mod_poly_is_zero(A->coeffs + i))
            A->length = i + 1;
    }
}

void fmpz_mod_bpoly_set_polyx(fmpz_mod_bpoly_t A, const fmpz_mod_poly_t B)
{
    slong i;
    fmpz_mod_bpoly_fit_length(A, B->length);
    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_set_fmpz(A->coeffs + i, B->coeffs + i);
        if (!fmpz_mod_poly_is_zero(A->coeffs + i))
            A->length = i + 1;
    }    
}

void fmpz_mod_bpoly_set_polyy(fmpz_mod_bpoly_t A, const fmpz_mod_poly_t B)
{
    fmpz_mod_bpoly_fit_length(A, 1);
    fmpz_mod_poly_set(A->coeffs + 0, B);
    A->length = !fmpz_mod_poly_is_zero(A->coeffs + 0);
}


void fmpz_mod_bpoly_get_coeff(fmpz_t c, const fmpz_mod_bpoly_t A, slong xi, slong yi)
{
    if (xi >= A->length)
        fmpz_zero(c);

    fmpz_mod_poly_get_coeff_fmpz(c, A->coeffs + xi, yi);
}

void fmpz_mod_bpoly_make_monic(fmpz_mod_bpoly_t A, slong order)
{
    slong i;
    fmpz_mod_poly_t t, lcinv;

    FLINT_ASSERT(A->length > 0);

    fmpz_mod_poly_init(t, A->modulus);
    fmpz_mod_poly_init(lcinv, A->modulus);
    fmpz_mod_poly_inv_series(lcinv, A->coeffs + A->length - 1, order);

    for (i = 0; i < A->length; i++)
    {
        fmpz_mod_poly_mullow(t, A->coeffs + i, lcinv, order);
        fmpz_mod_poly_swap(A->coeffs + i, t);
    }

    fmpz_mod_poly_clear(t);
    fmpz_mod_poly_clear(lcinv);
}

void fmpz_mod_bpoly_mul(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_bpoly_t C, slong order)
{
    slong i, j;
    fmpz_mod_poly_t t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    fmpz_mod_poly_init(t, B->modulus);

    fmpz_mod_bpoly_fit_length(A, B->length + C->length - 1);
    for (i = 0; i < B->length + C->length - 1; i++)
        fmpz_mod_poly_zero(A->coeffs + i);

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            fmpz_mod_poly_mullow(t, B->coeffs + i, C->coeffs + j, order);
            fmpz_mod_poly_add(A->coeffs + i + j, A->coeffs + i + j, t);
        }
    }

    A->length = B->length + C->length - 1;

    fmpz_mod_poly_clear(t);
}

void fmpz_mod_bpoly_add_poly_shift(fmpz_mod_bpoly_t A, const fmpz_mod_poly_t B, slong yshift)
{
    slong i;
    fmpz_t c;

    FLINT_ASSERT(A->length > B->length);

    fmpz_init(c);

    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_get_coeff_fmpz(c, A->coeffs + i, yshift);
        FLINT_ASSERT(fmpz_is_zero(c));
        fmpz_mod_poly_set_coeff_fmpz(A->coeffs + i, yshift, B->coeffs + i);
    }

    fmpz_clear(c);
}

void fmpz_mod_bpoly_sub(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_bpoly_t C)
{
    slong i;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    fmpz_mod_bpoly_fit_length(A, FLINT_MAX(B->length, C->length));

    for (i = 0; i < FLINT_MAX(B->length, C->length) - 1; i++)
    {
        if (i < B->length)
        {
            if (i < C->length)
            {
                fmpz_mod_poly_sub(A->coeffs + i, B->coeffs + i, C->coeffs + i);
            }
            else
            {
                fmpz_mod_poly_set(A->coeffs + i, B->coeffs + i);
            }
        }
        else
        {
            FLINT_ASSERT(i < C->length);
            fmpz_mod_poly_neg(A->coeffs + i, C->coeffs + i);
        }

        if (!fmpz_mod_poly_is_zero(A->coeffs + i))
            A->length = i + 1;
    }
}



void fmpz_bpoly_set(fmpz_bpoly_t A, const fmpz_bpoly_t B)
{
    slong i;

    FLINT_ASSERT(A != B);

    fmpz_bpoly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
        fmpz_poly_set(A->coeffs + i, B->coeffs + i);
}

void fmpz_bpoly_make_primitive(fmpz_poly_t g, fmpz_bpoly_t A)
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
        fmpz_poly_div(q, A->coeffs + i, g);
        fmpz_poly_swap(A->coeffs + i, q);
    }

    fmpz_poly_clear(q);
}

void fmpz_mod_bpoly_set_coeff(fmpz_mod_bpoly_t A, slong xi, slong yi, const fmpz_t c)
{
    slong i;

    if (xi >= A->length)
    {
        fmpz_mod_bpoly_fit_length(A, xi + 1);
        for (i = A->length; i <= xi; i++)
            fmpz_mod_poly_zero(A->coeffs + i);
        A->length = xi + 1;
    }

    fmpz_mod_poly_set_coeff_fmpz(A->coeffs + xi, yi, c);

    while (A->length > 0 && fmpz_poly_is_zero(A->coeffs + A->length - 1))
        A->length--;
}

void fmpz_mod_bpoly_zero(fmpz_mod_bpoly_t A)
{
    A->length = 0;
}

int fmpz_bpoly_divides(fmpz_bpoly_t Q, fmpz_bpoly_t A, fmpz_bpoly_t B)
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

void fmpz_bpoly_set_fmpz_mod_bpoly(fmpz_bpoly_t A, const fmpz_mod_bpoly_t B)
{
    slong i;

    fmpz_bpoly_fit_length(A, B->length);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
    {
        fmpz_poly_fit_length(A->coeffs + i, (B->coeffs + i)->length);
        (A->coeffs + i)->length = (B->coeffs + i)->length;
        _fmpz_vec_scalar_smod_fmpz((A->coeffs + i)->coeffs, (B->coeffs + i)->coeffs, (B->coeffs + i)->length, B->modulus);
    }
}

void fmpz_bpoly_eval(fmpz_poly_t E, const fmpz_bpoly_t A, const fmpz_t alpha)
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

void fmpz_bpoly_taylor_shift(fmpz_bpoly_t A, const fmpz_t alpha)
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
    fmpz_mod_bpoly_t Btilde;                /* mod p^k */
    fmpz_mod_bpoly_struct * newBitilde;     /* mod p^k */
    fmpz_mod_poly_struct * P;               /* mod p^k */
    fmpz_mod_poly_struct * d;               /* mod p^k */
    fmpz_mod_poly_struct * Bitilde;         /* mod p^k */
    fmpz_mod_poly_struct * d1;              /* mod p */
    fmpz_mod_poly_struct * Bitilde1;        /* mod p */
} bpoly_info_struct;

typedef bpoly_info_struct bpoly_info_t[1];

void bpoly_info_init(bpoly_info_t I, slong r, const fmpz_t p, ulong k)
{
    slong i;

    FLINT_ASSERT(r >= 2);

    I->r = r;

    I->lifting_prec = 0;

    I->k = k;
    fmpz_init_set(I->p, p);
    fmpz_init(I->pk);
    fmpz_pow_ui(I->pk, p, k);

    fmpz_mod_bpoly_init(I->Btilde, I->pk);

    I->newBitilde = (fmpz_mod_bpoly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_bpoly_struct));
    I->P          = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));
    I->d          = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));
    I->Bitilde    = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));
    I->d1         = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));
    I->Bitilde1   = (fmpz_mod_poly_struct *) flint_malloc(I->r*sizeof(fmpz_mod_poly_struct));

    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_bpoly_init(I->newBitilde + i, I->pk);
        fmpz_mod_poly_init(I->P + i, I->pk);
        fmpz_mod_poly_init(I->d + i, I->pk);
        fmpz_mod_poly_init(I->Bitilde + i, I->pk);
        fmpz_mod_poly_init(I->d1 + i, I->p);
        fmpz_mod_poly_init(I->Bitilde1 + i, I->p);
    }
}

void bpoly_info_clear(bpoly_info_t I)
{
    slong i;

    fmpz_clear(I->p);
    fmpz_clear(I->pk);

    fmpz_mod_bpoly_clear(I->Btilde);

    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_bpoly_clear(I->newBitilde + i);
        fmpz_mod_poly_clear(I->P + i);
        fmpz_mod_poly_clear(I->d + i);
        fmpz_mod_poly_clear(I->Bitilde + i);
        fmpz_mod_poly_clear(I->d1 + i);
        fmpz_mod_poly_clear(I->Bitilde1 + i);
    }

    flint_free(I->newBitilde);
    flint_free(I->P);
    flint_free(I->d);
    flint_free(I->Bitilde);
    flint_free(I->d1);
    flint_free(I->Bitilde1);
}

/*
    set out[i] so that
    1/(f[0]*f[1]*...*f[n-1]) = out[0]/f[0] + ... + out[n-1]/f[n-1]
*/
int partial_fraction_coeffs(fmpz_mod_poly_struct * out, const fmpz_mod_poly_struct * f, const fmpz_t p, slong n)
{
    slong i;
    fmpz_mod_poly_t num, den, a, b, g, t;

    FLINT_ASSERT(n > 1);

    fmpz_mod_poly_init(num, p);
    fmpz_mod_poly_init(den, p);
    fmpz_mod_poly_init(a, p);
    fmpz_mod_poly_init(b, p);
    fmpz_mod_poly_init(g, p);
    fmpz_mod_poly_init(t, p);

    fmpz_mod_poly_set_ui(num, 1);
    fmpz_mod_poly_mul(den, f + 0, f + 1);
    for (i = 2; i < n; i++)
        fmpz_mod_poly_mul(den, den, f + i);

    for (i = 0; i < n; i++)
    {
        fmpz_mod_poly_divrem(den, t, den, f + i);
        FLINT_ASSERT(fmpz_mod_poly_is_zero(t));
        fmpz_mod_poly_xgcd(g, a, b, f + i, den);
        if (fmpz_mod_poly_degree(g) != 0)
            return 0;
        FLINT_ASSERT(fmpz_is_one(g->coeffs + 0));
        fmpz_mod_poly_mul(t, b, num);
        fmpz_mod_poly_rem(out + i, t, f + i);
        fmpz_mod_poly_mul(t, a, num);
        fmpz_mod_poly_rem(num, t, den);
    }

    fmpz_mod_poly_clear(num);
    fmpz_mod_poly_clear(den);
    fmpz_mod_poly_clear(a);
    fmpz_mod_poly_clear(b);
    fmpz_mod_poly_clear(g);
    fmpz_mod_poly_clear(t);
    return 1;
}


int bpoly_info_disolve(bpoly_info_t I)
{
    slong i, j;
    fmpz_t pj, t1;
    fmpz_mod_poly_t error, t, s, s1, s2;

    if (!partial_fraction_coeffs(I->d1, I->Bitilde1, I->p, I->r))
        return 0;

    fmpz_init(pj);
    fmpz_init(t1);
    fmpz_mod_poly_init(error, I->pk);
    fmpz_mod_poly_init(t, I->pk);
    fmpz_mod_poly_init(s, I->p);
    fmpz_mod_poly_init(s1, I->p);
    fmpz_mod_poly_init(s2, I->pk);

    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_poly_set_ui(I->P + i, 1);
        for (j = 0; j < I->r; j++)
        {
            if (i == j)
                continue;
            fmpz_mod_poly_mul(I->P + i, I->P + i, I->Bitilde + j);
        }
    }

    fmpz_mod_poly_set_ui(error, 1);
    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_poly_set(I->d + i, I->d1 + i); /* slight abuse because moduli are different */
        fmpz_mod_poly_mul(t, I->d + i, I->P + i);
        fmpz_mod_poly_sub(error, error, t);
    }

    fmpz_one(pj);
    for (j = 1; j < I->k; j++)
    {
        fmpz_mul(pj, pj, I->p);
        fmpz_mod_poly_zero(s);
        for (i = error->length - 1; i >= 0; i--)
        {
            FLINT_ASSERT(fmpz_divisible(error->coeffs + i, pj));
            fmpz_divexact(t1, error->coeffs + i, pj);
            fmpz_mod(t1, t1, I->p);
            fmpz_mod_poly_set_coeff_fmpz(s, i, t1);
        }

        for (i = 0; i < I->r; i++)
        {
            fmpz_mod_poly_mul(s1, s, I->d1 + i);
            fmpz_set(&s2->p, I->p);
            fmpz_mod_poly_rem(s2, s1, I->Bitilde1 + i);
            fmpz_set(&s2->p, I->pk);
            fmpz_mod_poly_scalar_mul_fmpz(s2, s2, pj);
            fmpz_mod_poly_add(I->d + i, I->d + i, s2);
        }

        fmpz_mod_poly_set_ui(error, 1);
        for (i = 0; i < I->r; i++)
        {
            fmpz_mod_poly_mul(t, I->d + i, I->P + i);
            fmpz_mod_poly_sub(error, error, t);
        }
    }

    FLINT_ASSERT(fmpz_mod_poly_is_zero(error));

    fmpz_clear(pj);
    fmpz_clear(t1);
    fmpz_mod_poly_clear(error);
    fmpz_mod_poly_clear(t);
    fmpz_mod_poly_clear(s);
    fmpz_mod_poly_clear(s1);
    fmpz_mod_poly_clear(s2);

    return 1;
}



static void _bivar_lift_quintic(bpoly_info_t I)
{
    slong i, j, k;
    fmpz_mod_bpoly_t tp, tp1, error;
    fmpz_mod_poly_t ss, tt;

    fmpz_mod_poly_init(ss, I->pk);
    fmpz_mod_poly_init(tt, I->pk);
    fmpz_mod_bpoly_init(tp, I->pk);
    fmpz_mod_bpoly_init(tp1, I->pk);
    fmpz_mod_bpoly_init(error, I->pk);

    fmpz_mod_bpoly_mul(tp, I->newBitilde + 0, I->newBitilde + 1, I->lifting_prec);
    for (i = 2; i < I->r; i++)
    {
        fmpz_mod_bpoly_mul(tp1, tp, I->newBitilde + i, I->lifting_prec);
        fmpz_mod_bpoly_swap(tp1, tp);
    }
    fmpz_mod_bpoly_sub(error, I->Btilde, tp);

    for (j = 1; j < I->lifting_prec; j++)
    {
        fmpz_mod_poly_zero(ss);
        for (i = error->length - 1; i >= 0; i--)
        {
            fmpz_t ct;
            fmpz_init(ct);

            fmpz_mod_bpoly_get_coeff(ct, error, i, j);
            fmpz_mod_poly_set_coeff_fmpz(ss, i, ct);
            for (k = 0; k < j; k++)
            {
                fmpz_mod_bpoly_get_coeff(ct, error, i, k);
                FLINT_ASSERT(fmpz_is_zero(ct));
            }

            fmpz_clear(ct);
        }

        for (i = 0; i < I->r; i++)
        {
            fmpz_mod_poly_mul(tt, ss, I->d + i);
            fmpz_mod_poly_rem(tt, tt, I->Bitilde + i);
            fmpz_mod_bpoly_add_poly_shift(I->newBitilde + i, tt, j);
        }

        fmpz_mod_bpoly_mul(tp, I->newBitilde + 0, I->newBitilde + 1, I->lifting_prec);
        for (i = 2; i < I->r; i++)
        {
            fmpz_mod_bpoly_mul(tp1, tp, I->newBitilde + i, I->lifting_prec);
            fmpz_mod_bpoly_swap(tp1, tp);
        }
        fmpz_mod_bpoly_sub(error, I->Btilde, tp);
    }

    fmpz_mod_poly_clear(ss);
    fmpz_mod_poly_clear(tt);
    fmpz_mod_bpoly_clear(tp);
    fmpz_mod_bpoly_clear(tp1);
    fmpz_mod_bpoly_clear(error);
}

void fmpz_mod_bpoly_reverse_vars(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B)
{
    slong i, j;
    fmpz_mod_bpoly_zero(A);
    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_struct * Bi = B->coeffs + i;
        for (j = 0; j < Bi->length; j++)
        {
            if (!fmpz_is_zero(Bi->coeffs + j))
            {
                fmpz_mod_bpoly_set_coeff(A, j, i, Bi->coeffs + j);
            }
        }
    }
}

static void _bivar_lift_quartic2(bpoly_info_t I)
{
    slong i, j, k;
    fmpz_mod_poly_t t, t1;
    fmpz_mod_bpoly_t btilde;
    fmpz_mod_bpoly_struct newbitilde[2];

    FLINT_ASSERT(I->r == 2);

    fmpz_mod_poly_init(t, I->pk);
    fmpz_mod_poly_init(t1, I->pk);
    fmpz_mod_bpoly_init(btilde, I->pk);
    fmpz_mod_bpoly_reverse_vars(btilde, I->Btilde);
    for (k = 0; k < I->r; k++)
    {
        fmpz_mod_bpoly_init(newbitilde + k, I->pk);
        fmpz_mod_bpoly_reverse_vars(newbitilde + k, I->newBitilde + k);
        fmpz_mod_bpoly_fit_length(newbitilde + k, I->lifting_prec);
        FLINT_ASSERT((newbitilde + k)->length == 1);
    }

    for (j = 1; j < I->lifting_prec; j++)
    {
        if (j < btilde->length)
            fmpz_mod_poly_set(t, btilde->coeffs + j);
        else
            fmpz_mod_poly_zero(t);

        for (i = 1; i < j; i++)
        {
            fmpz_mod_poly_mul(t1, newbitilde[0].coeffs + i, newbitilde[1].coeffs + j - i);
            fmpz_mod_poly_sub(t, t, t1);
        }

        for (k = 0; k < I->r; k++)
        {
            fmpz_mod_poly_mul(t1, t, I->d + k);
            fmpz_mod_poly_rem(newbitilde[k].coeffs + j, t1, I->Bitilde + k);
            if (!fmpz_mod_poly_is_zero(newbitilde[k].coeffs + j))
                newbitilde[k].length = j + 1;
        }
    }

    for (k = 0; k < I->r; k++)
        fmpz_mod_bpoly_reverse_vars(I->newBitilde + k, newbitilde + k);

    fmpz_mod_poly_clear(t);
    fmpz_mod_poly_clear(t1);
    fmpz_mod_bpoly_clear(btilde);
    for (k = 0; k < I->r; k++)
    {
        fmpz_mod_bpoly_clear(newbitilde + k);
    }
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

    fmpz_mod_poly_init(t, I->pk);
    fmpz_mod_poly_init(t1, I->pk);
    fmpz_mod_bpoly_init(btilde, I->pk);
    fmpz_mod_bpoly_reverse_vars(btilde, I->Btilde);
    for (k = 0; k < I->r; k++)
    {
        fmpz_mod_bpoly_init(U + k, I->pk);
        fmpz_mod_bpoly_fit_length(U + k, I->lifting_prec);
        for (i = 0; i < I->lifting_prec; i++)
        {
            fmpz_mod_poly_zero(U[k].coeffs + i);
        }

        fmpz_mod_bpoly_init(newbitilde + k, I->pk);
        fmpz_mod_bpoly_reverse_vars(newbitilde + k, I->newBitilde + k);
        fmpz_mod_bpoly_fit_length(newbitilde + k, I->lifting_prec);
        FLINT_ASSERT(newbitilde[k].length == 1);
        for (i = 1; i < I->lifting_prec; i++)
        {
            fmpz_mod_poly_zero(newbitilde[k].coeffs + i);
        }
    }

    k = I->r - 2;
    fmpz_mod_poly_mul(U[k].coeffs + 0, newbitilde[k].coeffs + 0, newbitilde[k + 1].coeffs + 0);
    for (k--; k >= 1; k--)
        fmpz_mod_poly_mul(U[k].coeffs + 0, newbitilde[k].coeffs + 0, U[k + 1].coeffs + 0);

    for (j = 1; j < I->lifting_prec; j++)
    {
        k = I->r - 2;
        fmpz_mod_poly_zero(U[k].coeffs + j);
        for (i = 0; i <= j; i++)
        {
            fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + i, newbitilde[k + 1].coeffs + j - i);
            fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t1);
        }
        for (k--; k >= 1; k--)
        {
            fmpz_mod_poly_zero(U[k].coeffs + j);
            for (i = 0; i <= j; i++)
            {
                fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + i, U[k + 1].coeffs + j - i);
                fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t1);
            }
        }

        if (j < btilde->length)
            fmpz_mod_poly_set(t, btilde->coeffs + j);
        else
            fmpz_mod_poly_zero(t);

        for (i = 0; i <= j; i++)
        {
            fmpz_mod_poly_mul(t1, newbitilde[0].coeffs + i, U[1].coeffs + j - i);
            fmpz_mod_poly_sub(t, t, t1);
        }

        for (k = 0; k < I->r; k++)
        {
            fmpz_mod_poly_mul(t1, t, I->d + k);
            fmpz_mod_poly_rem(newbitilde[k].coeffs + j, t1, I->Bitilde + k);
            if (!fmpz_mod_poly_is_zero(newbitilde[k].coeffs + j))
                newbitilde[k].length = j + 1;
        }

        k = I->r - 2;
        fmpz_mod_poly_mul(t, newbitilde[k].coeffs + 0, newbitilde[k + 1].coeffs + j);
        fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + j, newbitilde[k + 1].coeffs + 0);
        fmpz_mod_poly_add(t, t, t1);
        fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t);
        for (k--; k >= 1; k--)
        {
            fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + 0, t);
            fmpz_mod_poly_swap(t, t1);
            fmpz_mod_poly_mul(t1, newbitilde[k].coeffs + j, U[k + 1].coeffs + 0);
            fmpz_mod_poly_add(t, t, t1);
            fmpz_mod_poly_add(U[k].coeffs + j, U[k].coeffs + j, t);
        }
    }

    for (k = 0; k < I->r; k++)
        fmpz_mod_bpoly_reverse_vars(I->newBitilde + k, newbitilde + k);

    fmpz_mod_poly_clear(t);
    fmpz_mod_poly_clear(t1);
    fmpz_mod_bpoly_clear(btilde);
    for (k = 0; k < I->r; k++)
    {
        fmpz_mod_bpoly_clear(U + k);
        fmpz_mod_bpoly_clear(newbitilde + k);
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
    fmpz_mod_bpoly_init(tryme, I->pk);
    fmpz_mod_bpoly_init(trymet, I->pk);
    fmpz_mod_poly_init(leadB, I->pk);

    F->length = 0;

    FLINT_ASSERT(B->length > 0);
    fmpz_mod_poly_set_fmpz_poly(leadB, B->coeffs + B->length - 1);

    len = I->r;
    subset = (slong *) flint_malloc(len * sizeof(slong));
    for (k = 0; k < len; k++)
        subset[k] = k;

    for (k = 1; k <= len/2; k++)
    {
        zassenhaus_subset_first(subset, len, k);
        while (1)
        {
            fmpz_mod_bpoly_set_polyy(tryme, leadB);
            for (i = 0; i < len; i++)
            {
                if (subset[i] >= 0)
                {
                    fmpz_mod_bpoly_mul(trymet, tryme,
                                   I->newBitilde + subset[i], I->lifting_prec);
                    fmpz_mod_bpoly_swap(trymet, tryme);
                }
            }
            fmpz_bpoly_set_fmpz_mod_bpoly(trymez, tryme);
            fmpz_bpoly_make_primitive(g, trymez);

            if (fmpz_bpoly_divides(Q, B, trymez))
            {
                fmpz_neg(alpha, alpha);
                fmpz_bpoly_taylor_shift(trymez, alpha);
                fmpz_neg(alpha, alpha);
                fmpz_tpoly_fit_length(F, F->length + 1);
                fmpz_bpoly_swap(F->coeffs + F->length, trymez);
                F->length++;
                fmpz_bpoly_swap(B, Q);
                FLINT_ASSERT(B->length > 0);
                fmpz_mod_poly_set_fmpz_poly(leadB, B->coeffs + B->length - 1);

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
    fmpz_mod_bpoly_clear(tryme);
    fmpz_mod_bpoly_clear(trymet);
    fmpz_mod_poly_clear(leadB);

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
    Bevalfac->num = 0;
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

    k = (pkbits + fmpz_bits(p))/fmpz_bits(p);

    bpoly_info_clear(I);
    bpoly_info_init(I, Bevalfac->num, p, k);
    I->lifting_prec = Blengthy + (B->coeffs + B->length - 1)->length;

    fmpz_mod_bpoly_set_fmpz_bpoly(I->Btilde, B);
    fmpz_mod_bpoly_make_monic(I->Btilde, I->lifting_prec);
    for (i = 0; i < I->r; i++)
    {
        fmpz_mod_poly_set_fmpz_poly(I->Bitilde1 + i, Bevalfac->p + i);
        fmpz_mod_poly_make_monic(I->Bitilde1 + i, I->Bitilde1 + i);
        fmpz_mod_poly_set_fmpz_poly(I->Bitilde + i, Bevalfac->p + i);
        fmpz_mod_poly_make_monic(I->Bitilde + i, I->Bitilde + i);
        fmpz_mod_bpoly_set_polyx(I->newBitilde + i, I->Bitilde + i);
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
