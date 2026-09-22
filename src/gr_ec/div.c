/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Dividing a point by an integer: given P and n, the Q with n Q = P.

    The solutions of n Q = P, when there are any, form a coset of E[n](F_q),
    so there are either none of them or exactly #E[n](F_q). That is what
    makes this two operations rather than one:

      div            insists on a single answer, and reports GR_DOMAIN when
                     there is none or when there is more than one
      div_nonunique  returns one of them, and reports GR_DOMAIN only when
                     there is none

    Two routes get there. Whenever the context knows an annihilator m -- and
    here, unlike everywhere else in the module, it is worth counting the
    points to learn one, since the alternative is polynomial root finding --
    write n = n1 n2 with gcd(n1, m) = 1. Dividing by n1 is multiplying by
    its inverse modulo m, which is one scalar multiplication and always
    unique, because gcd(n1, m) = 1 forces E[n1](F_q) to be trivial. If n2
    comes out as 1, which is the overwhelmingly common case, that is the
    whole computation.

    What is left is dividing by n2, every prime of which divides m. There is
    no shortcut for that: the x-coordinates of the solutions are roots of a
    polynomial of degree n2^2 built from division polynomials, and they have
    to be found and then checked. That is affordable only for small n2, so
    above GR_EC_DIV_MAX_DEGREE this reports GR_UNABLE rather than building a
    polynomial nobody wants to wait for.
*/

#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_vec.h"
#include "gr_ec.h"
#include "impl.h"

/*
    The largest division polynomial this will build, as a degree. Dividing
    by n needs one of degree n^2, so this allows n up to 64. Root finding
    over F_q costs a q-power of a polynomial of that degree, which is where
    the real time goes.
*/
#define GR_EC_DIV_MAX_DEGREE WORD(4096)

/*
    F(x), whose roots include the x-coordinates of every Q with n Q = R.

    With psi the division polynomials and phi_n = x psi_n^2 - psi_{n-1}
    psi_{n+1}, the x-coordinate of n Q is phi_n / psi_n^2, so the equation
    is phi_n(x) - x_R psi_n(x)^2 = 0. In terms of the normalised Psi that
    gr_ec_ctx_division_poly returns -- psi_n for odd n, psi_n / psi_2 for
    even n -- and W = psi_2^2 this is

      n odd:   (x - x_R) Psi_n^2 - Psi_{n-1} Psi_{n+1} W
      n even:  (x - x_R) Psi_n^2 W - Psi_{n-1} Psi_{n+1}

    For R the point at infinity the equation degenerates to psi_n = 0, whose
    roots are the x-coordinates of E[n] other than infinity; for even n the
    2-torsion sits in the psi_2 that the normalisation divided out, so W
    joins the product. Roots are a superset of what is wanted either way,
    since a root of psi_n satisfies the first equation as well without
    solving the problem, which is why every candidate is checked afterwards.
*/
static int
_gr_ec_division_target(gr_poly_t F, ulong n, gr_srcptr xR, int R_is_inf,
        gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_poly_t W, A, B, C;
    int status = GR_SUCCESS;

    gr_poly_init(W, R);
    gr_poly_init(A, R);
    gr_poly_init(B, R);
    gr_poly_init(C, R);

    status |= gr_ec_ctx_psi2_sqr(W, ctx);
    status |= gr_ec_ctx_division_poly(A, n, ctx);

    if (R_is_inf)
    {
        if (n % 2 == 0)
            status |= gr_poly_mul(F, A, W, R);
        else
            status |= gr_poly_set(F, A, R);

        goto cleanup;
    }

    status |= gr_ec_ctx_division_poly(B, n - 1, ctx);
    status |= gr_ec_ctx_division_poly(C, n + 1, ctx);

    /* B <- Psi_{n-1} Psi_{n+1}, A <- Psi_n^2 */
    status |= gr_poly_mul(B, B, C, R);
    status |= gr_poly_mul(A, A, A, R);

    if (n % 2 == 0)
        status |= gr_poly_mul(A, A, W, R);
    else
        status |= gr_poly_mul(B, B, W, R);

    /* C <- x - x_R */
    status |= gr_poly_zero(C, R);
    status |= gr_poly_set_coeff_si(C, 1, 1, R);
    {
        gr_ptr t;
        GR_TMP_INIT(t, R);
        status |= gr_neg(t, xR, R);
        status |= gr_poly_set_coeff_scalar(C, 0, t, R);
        GR_TMP_CLEAR(t, R);
    }

    status |= gr_poly_mul(F, A, C, R);
    status |= gr_poly_sub(F, F, B, R);

cleanup:
    gr_poly_clear(W, R);
    gr_poly_clear(A, R);
    gr_poly_clear(B, R);
    gr_poly_clear(C, R);

    return status;
}

/*
    Solve n Q = P by root finding, for n small enough that the polynomial
    above is worth building.

    found is set to the number of distinct solutions seen, capped at 2,
    which is all the caller needs in order to tell "none", "exactly one" and
    "more than one" apart. res receives the first one found.
*/
static int
_gr_ec_div_by_roots(gr_ec_point_t res, int * found, const gr_ec_point_t P,
        ulong n, int stop_at_first, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_poly_t F;
    gr_vec_t roots;
    fmpz_vec_t mult;
    gr_ec_point_t Q, T, first;
    gr_ptr xP;
    slong i, sz = R->sizeof_elem;
    int P_is_inf, status = GR_SUCCESS;

    *found = 0;

    if (gr_ec_ctx_model(ctx) != GR_EC_SHORT_WEIERSTRASS)
        return GR_UNABLE;

    P_is_inf = (gr_ec_point_is_inf(P, ctx) == T_TRUE);

    gr_poly_init(F, R);
    gr_vec_init(roots, 0, R);
    fmpz_vec_init(mult, 0);
    gr_ec_point_init(Q, ctx);
    gr_ec_point_init(T, ctx);
    gr_ec_point_init(first, ctx);
    GR_TMP_INIT(xP, R);

    /* the point at infinity is always divisible by anything */
    if (P_is_inf)
    {
        status |= gr_ec_point_zero(first, ctx);
        *found = 1;

        if (stop_at_first)
            goto done;
    }
    else
    {
        gr_ptr x, y;
        GR_TMP_INIT2(x, y, R);
        status |= gr_ec_point_get_affine(x, y, P, ctx);
        status |= gr_set(xP, x, R);
        GR_TMP_CLEAR2(x, y, R);
    }

    if (status != GR_SUCCESS)
        goto done;

    status = _gr_ec_division_target(F, n, xP, P_is_inf, ctx);

    if (status != GR_SUCCESS)
        goto done;

    status = gr_poly_roots(roots, mult, F, 0, R);

    if (status != GR_SUCCESS)
        goto done;

    for (i = 0; i < roots->length && !(stop_at_first && *found > 0); i++)
    {
        int sign;

        if (gr_ec_point_lift_x(Q, GR_ENTRY(roots->entries, i, sz), ctx)
                != GR_SUCCESS)
            continue;

        /* both points over this x, which are Q and -Q */
        for (sign = 0; sign < 2; sign++)
        {
            if (sign == 1 && gr_ec_point_neg(Q, Q, ctx) != GR_SUCCESS)
                break;

            if (gr_ec_point_mul_ui(T, Q, n, ctx) != GR_SUCCESS)
                continue;

            if (gr_ec_point_equal(T, P, ctx) != T_TRUE)
                continue;

            if (*found == 0)
            {
                status |= gr_ec_point_set(first, Q, ctx);
                *found = 1;
            }
            else if (gr_ec_point_equal(Q, first, ctx) != T_TRUE)
            {
                *found = 2;
                goto done;
            }

            if (stop_at_first)
                goto done;
        }
    }

done:
    if (status == GR_SUCCESS && *found > 0)
        status = gr_ec_point_set(res, first, ctx);

    GR_TMP_CLEAR(xP, R);
    gr_ec_point_clear(first, ctx);
    gr_ec_point_clear(T, ctx);
    gr_ec_point_clear(Q, ctx);
    fmpz_vec_clear(mult);
    gr_vec_clear(roots, R);
    gr_poly_clear(F, R);

    return status;
}

/*
    n split as n1 n2 with gcd(n1, m) = 1 and every prime of n2 dividing m.
    Peeling gcds off n is enough: each step removes at least one prime of m
    from what is left.
*/
static void
_split_off_smooth_part(fmpz_t n1, fmpz_t n2, const fmpz_t n, const fmpz_t m)
{
    fmpz_t g;

    fmpz_init(g);
    fmpz_set(n1, n);
    fmpz_one(n2);

    while (1)
    {
        fmpz_gcd(g, n1, m);

        if (fmpz_is_one(g))
            break;

        fmpz_mul(n2, n2, g);
        fmpz_divexact(n1, n1, g);
    }

    fmpz_clear(g);
}

/* the common body; unique says whether more than one answer is a failure */
static int
_gr_ec_point_div_fmpz(gr_ec_point_t res, const gr_ec_point_t P,
        const fmpz_t n, int unique, gr_ec_ctx_t ctx)
{
    fmpz_t k, m, n1, n2;
    gr_ec_point_t T;
    gr_ec_order_kind_t kind;
    int found, status = GR_SUCCESS;

    if (fmpz_is_zero(n))
    {
        /* 0 Q = P has every Q as a solution when P is the identity, and
           none otherwise; either way it is not a single answer */
        if (unique || gr_ec_point_is_inf(P, ctx) != T_TRUE)
            return GR_DOMAIN;

        return gr_ec_point_zero(res, ctx);
    }

    fmpz_init(k);
    fmpz_init(m);
    fmpz_init(n1);
    fmpz_init(n2);
    gr_ec_point_init(T, ctx);

    fmpz_abs(k, n);

    /* here it is worth counting: the alternative is root finding */
    kind = gr_ec_ctx_get_cached_order(m, ctx);

    if (kind == GR_EC_ORDER_UNKNOWN && gr_ec_ctx_order(m, ctx) == GR_SUCCESS)
        kind = GR_EC_ORDER_EXACT;

    if (kind == GR_EC_ORDER_UNKNOWN)
    {
        /* nothing is known to be invertible, so all of n is the hard part */
        fmpz_one(n1);
        fmpz_set(n2, k);
    }
    else if (fmpz_is_one(m))
    {
        /* the group is trivial: P is the identity and so is the answer */
        status = gr_ec_point_zero(res, ctx);
        goto cleanup;
    }
    else
        _split_off_smooth_part(n1, n2, k, m);

    status = gr_ec_point_set(T, P, ctx);

    /* the easy factor: multiply by the inverse of n1 modulo m */
    if (status == GR_SUCCESS && !fmpz_is_one(n1))
    {
        fmpz_t inv;
        fmpz_init(inv);

        if (fmpz_invmod(inv, n1, m))
            status = gr_ec_point_mul_fmpz(T, T, inv, ctx);
        else
            status = GR_UNABLE;

        fmpz_clear(inv);
    }

    if (status != GR_SUCCESS)
        goto cleanup;

    if (fmpz_is_one(n2))
    {
        /* gcd(n, m) = 1 makes E[n](F_q) trivial, so this is the only answer */
        status = gr_ec_point_set(res, T, ctx);
        goto cleanup;
    }

    /* the hard factor, if it is small enough to be worth attempting */
    if (!fmpz_abs_fits_ui(n2))
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    {
        ulong d = fmpz_get_ui(n2);

        if (d > (ulong) GR_EC_DIV_MAX_DEGREE
                || d * d > (ulong) GR_EC_DIV_MAX_DEGREE)
        {
            status = GR_UNABLE;
            goto cleanup;
        }

        status = _gr_ec_div_by_roots(res, &found, T, d, !unique, ctx);

        if (status != GR_SUCCESS)
            goto cleanup;

        if (found == 0 || (unique && found > 1))
            status = GR_DOMAIN;
    }

cleanup:
    if (status == GR_SUCCESS && fmpz_sgn(n) < 0)
        status = gr_ec_point_neg(res, res, ctx);

    gr_ec_point_clear(T, ctx);
    fmpz_clear(k);
    fmpz_clear(m);
    fmpz_clear(n1);
    fmpz_clear(n2);

    return status;
}

int
gr_ec_point_div_fmpz(gr_ec_point_t res, const gr_ec_point_t P,
        const fmpz_t n, gr_ec_ctx_t ctx)
{
    return _gr_ec_point_div_fmpz(res, P, n, 1, ctx);
}

int
gr_ec_point_div_fmpz_nonunique(gr_ec_point_t res, const gr_ec_point_t P,
        const fmpz_t n, gr_ec_ctx_t ctx)
{
    return _gr_ec_point_div_fmpz(res, P, n, 0, ctx);
}

int
gr_ec_point_div_ui(gr_ec_point_t res, const gr_ec_point_t P, ulong n,
        gr_ec_ctx_t ctx)
{
    fmpz_t t;
    int status;

    fmpz_init_set_ui(t, n);
    status = gr_ec_point_div_fmpz(res, P, t, ctx);
    fmpz_clear(t);

    return status;
}

int
gr_ec_point_div_si(gr_ec_point_t res, const gr_ec_point_t P, slong n,
        gr_ec_ctx_t ctx)
{
    fmpz_t t;
    int status;

    fmpz_init(t);
    fmpz_set_si(t, n);
    status = gr_ec_point_div_fmpz(res, P, t, ctx);
    fmpz_clear(t);

    return status;
}

/*
    (a / b) P is the Q with b Q = a P. Multiplying first and dividing after
    gives the same set of answers as the other order, because gcd(a, b) = 1
    makes multiplication by a a bijection of E[b].
*/
static int
_gr_ec_point_mul_fmpq(gr_ec_point_t res, const gr_ec_point_t P,
        const fmpq_t c, int unique, gr_ec_ctx_t ctx)
{
    gr_ec_point_t T;
    int status;

    gr_ec_point_init(T, ctx);

    status = gr_ec_point_mul_fmpz(T, P, fmpq_numref(c), ctx);

    if (status == GR_SUCCESS)
        status = _gr_ec_point_div_fmpz(res, T, fmpq_denref(c), unique, ctx);

    gr_ec_point_clear(T, ctx);

    return status;
}

int
gr_ec_point_mul_fmpq(gr_ec_point_t res, const gr_ec_point_t P,
        const fmpq_t c, gr_ec_ctx_t ctx)
{
    return _gr_ec_point_mul_fmpq(res, P, c, 1, ctx);
}

int
gr_ec_point_mul_fmpq_nonunique(gr_ec_point_t res, const gr_ec_point_t P,
        const fmpq_t c, gr_ec_ctx_t ctx)
{
    return _gr_ec_point_mul_fmpq(res, P, c, 0, ctx);
}

int
gr_ec_point_div_fmpq(gr_ec_point_t res, const gr_ec_point_t P,
        const fmpq_t c, gr_ec_ctx_t ctx)
{
    fmpq_t d;
    int status;

    if (fmpq_is_zero(c))
        return GR_DOMAIN;

    fmpq_init(d);
    fmpq_inv(d, c);
    status = gr_ec_point_mul_fmpq(res, P, d, ctx);
    fmpq_clear(d);

    return status;
}

/*
    The other two representations, through the projective one. Division is
    root finding, not a ladder, so there is nothing to gain from doing it
    natively in either of them.
*/

#define GR_EC_DIV_VIA_PROJECTIVE(kind, short, op, argtype) \
int \
kind ## _ ## op(kind ## _t res, const kind ## _t P, argtype n, \
        gr_ec_ctx_t ctx) \
{ \
    gr_ec_point_t A; \
    int status; \
 \
    gr_ec_point_init(A, ctx); \
    status = gr_ec_point_set_ ## short(A, P, ctx); \
 \
    if (status == GR_SUCCESS) \
        status = gr_ec_point_ ## op(A, A, n, ctx); \
 \
    if (status == GR_SUCCESS) \
        status = kind ## _set_point(res, A, ctx); \
 \
    gr_ec_point_clear(A, ctx); \
 \
    return status; \
}

#define GR_EC_DIV_FAMILY(kind, short) \
    GR_EC_DIV_VIA_PROJECTIVE(kind, short, div_fmpz, const fmpz_t) \
    GR_EC_DIV_VIA_PROJECTIVE(kind, short, div_fmpz_nonunique, const fmpz_t) \
    GR_EC_DIV_VIA_PROJECTIVE(kind, short, div_ui, ulong) \
    GR_EC_DIV_VIA_PROJECTIVE(kind, short, div_si, slong) \
    GR_EC_DIV_VIA_PROJECTIVE(kind, short, mul_fmpq, const fmpq_t) \
    GR_EC_DIV_VIA_PROJECTIVE(kind, short, mul_fmpq_nonunique, const fmpq_t) \
    GR_EC_DIV_VIA_PROJECTIVE(kind, short, div_fmpq, const fmpq_t)

GR_EC_DIV_FAMILY(gr_ec_aff_point, aff_point)
GR_EC_DIV_FAMILY(gr_ec_jac_point, jac_point)
