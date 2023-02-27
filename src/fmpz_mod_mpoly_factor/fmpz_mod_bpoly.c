/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_vec.h"
#include "fmpz_mod_mpoly_factor.h"


int fmpz_mod_bpoly_is_canonical(const fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
{
    slong i;

    if (A->length < 1)
        return A->length == 0;

    for (i = 0; i < A->length; i++)
    {
        if (!fmpz_mod_poly_is_canonical(A->coeffs + i, ctx))
            return 0;
        if (i + 1 == A->length && fmpz_mod_poly_is_zero(A->coeffs + i, ctx))
            return 0;
    }

    return 1;
}

void fmpz_mod_bpoly_clear(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_mod_poly_clear(A->coeffs + i, ctx);
    if (A->alloc > 0)
        flint_free(A->coeffs);
}

void fmpz_mod_bpoly_print_pretty(
    const fmpz_mod_bpoly_t A,
    const char * xvar,
    const char * yvar,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    int first;

    first = 1;
    for (i = A->length - 1; i >= 0; i--)
    {
        if (fmpz_mod_poly_is_zero(A->coeffs + i, ctx))
            continue;

        if (!first)
            flint_printf(" + ");

        first = 0;

        flint_printf("(");
        fmpz_mod_poly_print_pretty(A->coeffs + i, yvar, ctx);
        flint_printf(")*%s^%wd", xvar, i);
    }

    if (first)
        flint_printf("0");
}

void fmpz_mod_bpoly_fit_length(
    fmpz_mod_bpoly_t A,
    slong len,
    const fmpz_mod_ctx_t ctx)
{
    slong i = A->alloc;

    if (len <= i)
        return;

    len = FLINT_MAX(len, 2*i);

    A->coeffs = FLINT_ARRAY_REALLOC(A->coeffs, len, fmpz_mod_poly_struct);

    for ( ; i < len; i++)
        fmpz_mod_poly_init(A->coeffs + i, ctx);

    A->alloc = len;
}

void fmpz_mod_bpoly_zero(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
{
    A->length = 0;
}

void fmpz_mod_bpoly_set_coeff(
    fmpz_mod_bpoly_t A,
    slong xi, slong yi,
    const fmpz_t c,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    if (xi >= A->length)
    {
        fmpz_mod_bpoly_fit_length(A, xi + 1, ctx);
        for (i = A->length; i <= xi; i++)
            fmpz_mod_poly_zero(A->coeffs + i, ctx);
        A->length = xi + 1;
    }

    fmpz_mod_poly_set_coeff_fmpz(A->coeffs + xi, yi, c, ctx);
    while (A->length > 0 && fmpz_mod_poly_is_zero(A->coeffs + A->length - 1, ctx))
        A->length--;
}

void fmpz_mod_bpoly_set_poly_gen1(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_bpoly_fit_length(A, 1, ctx);
	fmpz_mod_poly_set(A->coeffs + 0, B, ctx);
	A->length = !fmpz_mod_poly_is_zero(A->coeffs + 0, ctx);
}


void fmpz_mod_bpoly_set_poly_gen0(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mod_bpoly_fit_length(A, B->length, ctx);
    A->length = B->length;
    for (i = 0; i < B->length; i++)
        fmpz_mod_poly_set_fmpz(A->coeffs + i, B->coeffs + i, ctx);
}

void fmpz_mod_bpoly_one(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_bpoly_fit_length(A, 1, ctx);
    A->length = 1;
    fmpz_mod_poly_one(A->coeffs + 0, ctx);
}

slong fmpz_mod_bpoly_degree1(
    const fmpz_mod_bpoly_t A,
    const fmpz_mod_ctx_t ctx)
{
    slong i, len = 0;
    for (i = 0; i < A->length; i++)
        len = FLINT_MAX(len, A->coeffs[i].length);
    return len - 1;    
}



int fmpz_mod_bpoly_equal(
    const fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < A->length; i++)
    {
        if (!fmpz_mod_poly_equal(A->coeffs + i, B->coeffs + i, ctx))
            return 0;
    }

    return 1;
}


void fmpz_mod_bpoly_set(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    if (A == B)
        return;

    fmpz_mod_bpoly_fit_length(A, B->length, ctx);
    A->length = B->length;

    for (i = 0; i < B->length; i++)
        fmpz_mod_poly_set(A->coeffs + i, B->coeffs + i, ctx);
}

void _fmpz_mod_poly_taylor_shift_horner(
    fmpz * a,
    const fmpz_t c,
    slong n,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;

    if (!fmpz_is_zero(c))
    {
        for (i = n - 2; i >= 0; i--)
            for (j = i; j < n - 1; j++)
                fmpz_mod_addmul(a + j, a + j, a + j + 1, c, ctx);
    }
}



void fmpz_mod_bpoly_taylor_shift_gen1(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_t c,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    fmpz_mod_bpoly_set(A, B, ctx);

    for (i = A->length - 1; i >= 0; i--)
        _fmpz_mod_poly_taylor_shift_horner(A->coeffs[i].coeffs, c,
                                                     A->coeffs[i].length, ctx);
}

void fmpz_mod_bpoly_sub(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_bpoly_t C,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    slong Alen = FLINT_MAX(B->length, C->length);

    fmpz_mod_bpoly_fit_length(A, Alen, ctx);

    for (i = 0; i < Alen; i++)
    {
        if (i < B->length)
        {
            if (i < C->length)
                fmpz_mod_poly_sub(A->coeffs + i, B->coeffs + i, C->coeffs + i, ctx);
            else
                fmpz_mod_poly_set(A->coeffs + i, B->coeffs + i, ctx);
        }
        else
        {
            FLINT_ASSERT(i < C->length);
            fmpz_mod_poly_neg(A->coeffs + i, C->coeffs + i, ctx);
        }
    }

    A->length = Alen;
    fmpz_mod_bpoly_normalise(A, ctx);
}


void fmpz_mod_bpoly_add(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_bpoly_t C,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    slong Alen = FLINT_MAX(B->length, C->length);

    fmpz_mod_bpoly_fit_length(A, Alen, ctx);

    for (i = 0; i < Alen; i++)
    {
        if (i < B->length)
        {
            if (i < C->length)
                fmpz_mod_poly_add(A->coeffs + i, B->coeffs + i, C->coeffs + i, ctx);
            else
                fmpz_mod_poly_set(A->coeffs + i, B->coeffs + i, ctx);
        }
        else
        {
            FLINT_ASSERT(i < C->length);
            fmpz_mod_poly_set(A->coeffs + i, C->coeffs + i, ctx);
        }
    }

    A->length = Alen;
    fmpz_mod_bpoly_normalise(A, ctx);
}


void fmpz_mod_bpoly_make_primitive(
    fmpz_mod_poly_t g,
    fmpz_mod_bpoly_t A,
    const fmpz_mod_ctx_t ctx)
{
    slong Alen = A->length;
    slong i;
    fmpz_mod_poly_t q, r;

    fmpz_mod_poly_init(q, ctx);
    fmpz_mod_poly_init(r, ctx);

    fmpz_mod_poly_zero(g, ctx);
    for (i = 0; i < Alen; i++)
    {
        fmpz_mod_poly_gcd(q, g, A->coeffs + i, ctx);
        fmpz_mod_poly_swap(g, q, ctx);
    }

    for (i = 0; i < Alen; i++)
    {
        fmpz_mod_poly_divrem(q, r, A->coeffs + i, g, ctx);
        FLINT_ASSERT(fmpz_mod_poly_is_zero(r, ctx));
        fmpz_mod_poly_swap(A->coeffs + i, q, ctx);
    }

    /* make lc_xy(A) one */
    if (Alen > 0)
    {
        fmpz * c = A->coeffs[Alen - 1].coeffs + A->coeffs[Alen - 1].length - 1;
        if (!fmpz_is_one(c))
        {
            fmpz_t cinv;
            fmpz_mod_poly_scalar_mul_fmpz(g, g, c, ctx);
            fmpz_init(cinv);
            fmpz_mod_inv(cinv, c, ctx);
            for (i = 0; i < Alen; i++)
                fmpz_mod_poly_scalar_mul_fmpz(A->coeffs + i,
                                              A->coeffs + i, cinv, ctx);
            fmpz_clear(cinv);
        }
    }

    fmpz_mod_poly_clear(q, ctx);
    fmpz_mod_poly_clear(r, ctx);
}

/* multiplication in F[y][x] */
void fmpz_mod_bpoly_mul(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_bpoly_t C,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    fmpz_mod_poly_struct * t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    if (B->length < 1 || C->length < 1)
    {
        A->length = 0;
        return;
    }

#if 0
    if (B->length > 2 && C->length > 2)
    {
        n_poly_t a, b, c;
        slong order = n_bpoly_degree1(B) + n_bpoly_degree1(C) + 1;

        n_poly_init(a);
        n_poly_init(b);
        n_poly_init(c);

        for (i = B->length - 1; i >= 0; i--)
        {
            n_poly_struct * Bi = B->coeffs + i;
            for (j = Bi->length - 1; j >= 0; j--)
                n_poly_set_coeff(b, order*i + j, Bi->coeffs[j]);
        }

        for (i = C->length - 1; i >= 0; i--)
        {
            n_poly_struct * Ci = C->coeffs + i;
            for (j = Ci->length - 1; j >= 0; j--)
                n_poly_set_coeff(c, order*i + j, Ci->coeffs[j]);
        }

        n_poly_mod_mul(a, b, c, ctx);

        A->length = 0;
        for (i = B->length + C->length - 1; i >= 0; i--)
        {
            for (j = order - 1; j >= 0; j--)
                n_bpoly_set_coeff(A, i, j, n_poly_get_coeff(a, order*i + j));
        }

        n_poly_clear(a);
        n_poly_clear(b);
        n_poly_clear(c);
        return;
    }
#endif

    fmpz_mod_bpoly_fit_length(A, B->length + C->length, ctx);
    for (i = 0; i < B->length + C->length - 1; i++)
        fmpz_mod_poly_zero(A->coeffs + i, ctx);

    t = A->coeffs + B->length + C->length - 1;

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            fmpz_mod_poly_mul(t, B->coeffs + i, C->coeffs + j, ctx);
            fmpz_mod_poly_add(A->coeffs + i + j, A->coeffs + i + j, t, ctx);
        }
    }

    A->length = B->length + C->length - 1;
    fmpz_mod_bpoly_normalise(A, ctx);
}


/*
    multiplication in (F[y]/y^order)[x]
        B, C need not be reduced mod y^order
        A should come out reduced mod y^order
*/
void fmpz_mod_bpoly_mul_series(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_bpoly_t C,
    slong order,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    fmpz_mod_poly_struct * t;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    if (B->length < 1 || C->length < 1)
    {
        A->length = 0;
        return;
    }
#if 0
    if (B->length > 2 && C->length > 2)
    {
        n_poly_t a, b, c;

        n_poly_init(a);
        n_poly_init(b);
        n_poly_init(c);

        for (i = B->length - 1; i >= 0; i--)
        {
            n_poly_struct * Bi = B->coeffs + i;
            for (j = FLINT_MIN(order, Bi->length) - 1; j >= 0; j--)
                n_poly_set_coeff(b, 2*order*i + j, Bi->coeffs[j]);
        }

        for (i = C->length - 1; i >= 0; i--)
        {
            n_poly_struct * Ci = C->coeffs + i;
            for (j = FLINT_MIN(order, Ci->length) - 1; j >= 0; j--)
                n_poly_set_coeff(c, 2*order*i + j, Ci->coeffs[j]);
        }

        n_poly_mod_mul(a, b, c, ctx);

        A->length = 0;
        for (i = B->length + C->length - 1; i >= 0; i--)
        {
            for (j = order - 1; j >= 0; j--)
                n_bpoly_set_coeff(A, i, j, n_poly_get_coeff(a, 2*order*i + j));
        }

        n_poly_clear(a);
        n_poly_clear(b);
        n_poly_clear(c);
        return;
    }
#endif

    fmpz_mod_bpoly_fit_length(A, B->length + C->length, ctx);
    for (i = 0; i < B->length + C->length - 1; i++)
        fmpz_mod_poly_zero(A->coeffs + i, ctx);

    t = A->coeffs + B->length + C->length - 1;

    for (i = 0; i < B->length; i++)
    {
        for (j = 0; j < C->length; j++)
        {
            fmpz_mod_poly_mullow(t, B->coeffs + i, C->coeffs + j, order, ctx);
            fmpz_mod_poly_add(A->coeffs + i + j, A->coeffs + i + j, t, ctx);
        }
    }

    A->length = B->length + C->length - 1;
    fmpz_mod_bpoly_normalise(A, ctx);
}


/*
    division in (F[y]/y^order)[x]
        A, B need not be reduced mod y^order
        Q, R should come out reduced mod y^order

    TODO: make this faster
*/
void fmpz_mod_bpoly_divrem_series(
    fmpz_mod_bpoly_t Q,
    fmpz_mod_bpoly_t R,
    const fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    slong order,
    const fmpz_mod_ctx_t ctx)
{
    slong i, qoff;
    fmpz_mod_poly_t q, t;

    FLINT_ASSERT(R != A);
    FLINT_ASSERT(R != B);
    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);
    FLINT_ASSERT(B->length > 0);

    fmpz_mod_poly_init(q, ctx);
    fmpz_mod_poly_init(t, ctx);

    fmpz_mod_bpoly_set(R, A, ctx);
    for (i = 0; i < R->length; i++)
        fmpz_mod_poly_truncate(R->coeffs + i, order, ctx);
    fmpz_mod_bpoly_normalise(R, ctx);

    Q->length = 0;

    while (R->length >= B->length)
    {
        fmpz_mod_poly_div_series(q, R->coeffs + R->length - 1,
                                    B->coeffs + B->length - 1, order, ctx);

        for (i = 0; i < B->length; i++)
        {
            fmpz_mod_poly_mullow(t, B->coeffs + i, q, order, ctx);
            fmpz_mod_poly_sub(R->coeffs + i + R->length - B->length,
                              R->coeffs + i + R->length - B->length, t, ctx);
        }

        qoff = R->length - B->length;

        FLINT_ASSERT(qoff >= 0);
        if (qoff >= Q->length)
        {
            fmpz_mod_bpoly_fit_length(Q, qoff + 1, ctx);
            for (i = Q->length; i <= qoff; i++)
                fmpz_mod_poly_zero(Q->coeffs + i, ctx);
            Q->length = qoff + 1;
        }

        fmpz_mod_poly_set(Q->coeffs + qoff, q, ctx);

        FLINT_ASSERT(fmpz_mod_poly_is_zero(R->coeffs + R->length - 1, ctx));

        fmpz_mod_bpoly_normalise(R, ctx);
    }

    fmpz_mod_poly_clear(q, ctx);
    fmpz_mod_poly_clear(t, ctx);
}

/*
    divisibility in F[y][x]
*/
int fmpz_mod_bpoly_divides(
    fmpz_mod_bpoly_t Q,
    const fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong i, qoff;
    int divides;
    fmpz_mod_poly_t q, t;
    fmpz_mod_bpoly_t R;

    FLINT_ASSERT(Q != A);
    FLINT_ASSERT(Q != B);
    FLINT_ASSERT(B->length > 0);
#if 0
    /* ksub not faster :( */
    if (0 && A->length > B->length && B->length > 2)
    {
        n_poly_t a, b, q, r;
        slong j;
        slong Adeg = n_bpoly_degree1(A);
        slong Bdeg = n_bpoly_degree1(B);
        slong Qdeg = Adeg - Bdeg;
        slong order = Adeg + 1;

        if (Qdeg < 0)
            return 0;

        n_poly_init(a);
        n_poly_init(b);
        n_poly_init(q);
        n_poly_init(r);

        for (i = B->length - 1; i >= 0; i--)
        {
            n_poly_struct * Bi = B->coeffs + i;
            for (j = Bi->length - 1; j >= 0; j--)
                n_poly_set_coeff(b, order*i + j, Bi->coeffs[j]);
        }

        for (i = A->length - 1; i >= 0; i--)
        {
            n_poly_struct * Ai = A->coeffs + i;
            for (j = Ai->length - 1; j >= 0; j--)
                n_poly_set_coeff(a, order*i + j, Ai->coeffs[j]);
        }

        n_poly_mod_divrem(q, r, a, b, ctx);

        if (r->length > 0 ||
            q->length - 1 < order*(A->length - B->length) ||
            q->length - 1 > Qdeg + order*(A->length - B->length))
        {
            divides = 0;
            goto cleanup_inner;
        }

        for (i = A->length - B->length; i >= 0; i--)
        {
            for (j = order - 1; j >= 0; j--)
            {
                mp_limb_t qc = n_poly_get_coeff(q, order*i + j);
                if (qc == 0)
                    continue;

                if (j > Qdeg)
                {
                    divides = 0;
                    goto cleanup_inner;
                }
                n_bpoly_set_coeff(Q, i, j, qc);
            }
        }

        divides = 1;

    cleanup_inner:

        n_poly_clear(a);
        n_poly_clear(b);
        n_poly_clear(q);
        n_poly_clear(r);
        return divides;
    }
#endif

    fmpz_mod_poly_init(q, ctx);
    fmpz_mod_poly_init(t, ctx);
    fmpz_mod_bpoly_init(R, ctx);

    fmpz_mod_bpoly_set(R, A, ctx);

    Q->length = 0;

    while (R->length >= B->length)
    {
        fmpz_mod_poly_divrem(q, t, R->coeffs + R->length - 1,
                                   B->coeffs + B->length - 1, ctx);
        if (!fmpz_mod_poly_is_zero(t, ctx))
        {
            divides = 0;
            goto cleanup;
        }

        for (i = 0; i < B->length; i++)
        {
            fmpz_mod_poly_mul(t, B->coeffs + i, q, ctx);
            fmpz_mod_poly_sub(R->coeffs + i + R->length - B->length,
                              R->coeffs + i + R->length - B->length, t, ctx);
        }

        qoff = R->length - B->length;

        FLINT_ASSERT(qoff >= 0);
        if (qoff >= Q->length)
        {
            fmpz_mod_bpoly_fit_length(Q, qoff + 1, ctx);
            for (i = Q->length; i <= qoff; i++)
                fmpz_mod_poly_zero(Q->coeffs + i, ctx);
            Q->length = qoff + 1;
        }

        fmpz_mod_poly_set(Q->coeffs + qoff, q, ctx);

        while (R->length > 0 && fmpz_mod_poly_is_zero(R->coeffs + R->length - 1, ctx))
            R->length--;
    }

    divides = (R->length == 0);

cleanup:

    fmpz_mod_poly_clear(q, ctx);
    fmpz_mod_poly_clear(t, ctx);
    fmpz_mod_bpoly_clear(R, ctx);

    return divides;
}

void fmpz_mod_bpoly_taylor_shift_gen0(
    fmpz_mod_bpoly_t A,
    const fmpz_t alpha,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    slong n = A->length;
    fmpz_mod_poly_struct * Acoeffs = A->coeffs;
    fmpz_t c, alpha_inv;

    FLINT_ASSERT(fmpz_mod_is_canonical(alpha, ctx));

    if (fmpz_is_zero(alpha))
        return;

    fmpz_init(c);
    fmpz_init(alpha_inv);
    fmpz_mod_inv(alpha_inv, alpha, ctx);

    fmpz_one(c);
    for (i = 1; i < n; i++)
    {
        fmpz_mod_mul(c, c, alpha, ctx);
        _fmpz_mod_vec_scalar_mul_fmpz_mod(Acoeffs[i].coeffs, Acoeffs[i].coeffs,
                                                    Acoeffs[i].length, c, ctx);
    }

    for (i = n - 2; i >= 0; i--)
    {
        for (j = i; j < n - 1; j++)
        {
            fmpz_mod_poly_add(Acoeffs + j, Acoeffs + j, Acoeffs + j + 1, ctx);
        }
    }

    fmpz_one(c);
    for (i = 1; i < n; i++)
    {
        fmpz_mod_mul(c, c, alpha_inv, ctx);
        _fmpz_mod_vec_scalar_mul_fmpz_mod(Acoeffs[i].coeffs, Acoeffs[i].coeffs,
                                                    Acoeffs[i].length, c, ctx);
    }

    fmpz_clear(c);
    fmpz_clear(alpha_inv);
}


void fmpz_mod_bpoly_derivative_gen0(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong i;

    FLINT_ASSERT(A != B);

    if (B->length < 2)
    {
        fmpz_mod_bpoly_zero(A, ctx);
        return;
    }

    fmpz_mod_bpoly_fit_length(A, B->length - 1, ctx);

    for (i = 1; i < B->length; i++)
        fmpz_mod_poly_scalar_mul_ui(A->coeffs + i - 1, B->coeffs + i, i, ctx);

    A->length = B->length - 1;
    fmpz_mod_bpoly_normalise(A, ctx);
}


void fmpz_mod_mpoly_get_fmpz_mod_bpoly(
    fmpz_mod_bpoly_t A,
    const fmpz_mod_mpoly_t B,
    slong varx,
    slong vary,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong j;
    slong NB;
    ulong Bexpx, Bexpy;
    slong Boffx, Bshiftx, Boffy, Bshifty;
    ulong mask;

    FLINT_ASSERT(B->bits <= FLINT_BITS);
    NB = mpoly_words_per_exp_sp(B->bits, ctx->minfo);

    mpoly_gen_offset_shift_sp(&Boffx, &Bshiftx, varx, B->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&Boffy, &Bshifty, vary, B->bits, ctx->minfo);
    mask = (-UWORD(1)) >> (FLINT_BITS - B->bits);

    fmpz_mod_bpoly_zero(A, ctx->ffinfo);
    for (j = 0; j < B->length; j++)
    {
        Bexpx = ((B->exps + NB*j)[Boffx] >> Bshiftx) & mask;
        Bexpy = ((B->exps + NB*j)[Boffy] >> Bshifty) & mask;
        fmpz_mod_bpoly_set_coeff(A, Bexpx, Bexpy, B->coeffs + j, ctx->ffinfo);
    }
}


void fmpz_mod_mpoly_set_fmpz_mod_bpoly(
    fmpz_mod_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mod_bpoly_t B,
    slong varx,
    slong vary,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong n = ctx->minfo->nvars;
    slong i, j;
    slong NA;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    ulong * Aexps;
    TMP_INIT;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);
    TMP_START;

    Aexps = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    for (i = 0; i < n; i++)
        Aexps[i] = 0;

    NA = mpoly_words_per_exp(Abits, ctx->minfo);
    fmpz_mod_mpoly_fit_length_reset_bits(A, 4, Abits, ctx);

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Alen = 0;
    for (i = 0; i < B->length; i++)
    {
        fmpz_mod_poly_struct * Bc = B->coeffs + i;
        _fmpz_mod_mpoly_fit_length(&Acoeff, &A->coeffs_alloc,
                               &Aexp, &A->exps_alloc, NA, Alen + Bc->length);

        for (j = 0; j < Bc->length; j++)
        {
            if (fmpz_is_zero(Bc->coeffs + j))
                continue;

            Aexps[varx] = i;
            Aexps[vary] = j;
            fmpz_set(Acoeff + Alen, Bc->coeffs + j);
            mpoly_set_monomial_ui(Aexp + NA*Alen, Aexps, Abits, ctx->minfo);
            Alen++;
        }
    }
    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->length = Alen;

    TMP_END;

    fmpz_mod_mpoly_sort_terms(A, ctx);
}
