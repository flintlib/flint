/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "qfb.h"
#include "ulong_extras.h"
#include "gr.h"
#include "gr_ec.h"
#include "impl.h"

/* ------------------------------------------------------------------ */
/* complex multiplication and supersingular curves                    */
/* ------------------------------------------------------------------ */

/*
    When E/F_p has complex multiplication by an order of small discriminant
    the trace is not something to search for. For the thirteen discriminants
    D of class number one the Hilbert class polynomial is linear, so
    j(E) = j_D is exactly the statement that E has that complex
    multiplication, and then 4p = t^2 + |D| v^2 has a solution, which
    qfb_cornacchia finds. Supersingular curves are the same story with
    t = 0, which is what happens when D is not a square modulo p.

    Cornacchia only gives |t|, and which of the twists with this
    j-invariant the curve actually is has to be decided separately. None of
    the thirteen cases needs a scalar multiplication for that.

    For j = 0 and j = 1728 the curve has sextic and quartic twists, and the
    sextic (respectively quartic) residue character of the coefficient
    picks one out, at the cost of one exponentiation. Those characters are
    classical evaluations of Jacobi sums -- see Ireland and Rosen, A
    Classical Introduction to Modern Number Theory, chapter 18.

    For the other eleven only the quadratic twist exists, and then two
    independent bits decide the sign, neither of which costs more than a
    Jacobi symbol:

      * a sign convention for |t| itself. For odd m = |D| the Jacobi symbol
        (t | m) does it, because -1 is a non-residue modulo every such m,
        so exactly one of +-t satisfies (t | m) = 1. Where m is a power of
        two (D = -8 and D = -16) a congruence on t and v takes its place.

      * which quadratic twist the curve is, which is the Legendre symbol
        (c_D a6 | p) for a constant c_D attached to the discriminant.

    That this is the shape of the answer is classical; see Rubin and
    Silverberg, Choosing the correct elliptic curve in the CM method,
    Math. Comp. 79 (2010). The constants c_D and the overall signs in
    cm_discs below were determined here by calibration against
    gr_ec_ctx_cardinality_bsgs over three thousand curves per discriminant,
    and each resulting rule was then confirmed on twenty thousand further
    curves.
*/

/* the representative of x modulo p in (-p/2, p/2] */
static void
_center(fmpz_t x, const fmpz_t p)
{
    fmpz_t h;
    fmpz_init(h);
    fmpz_tdiv_q_2exp(h, p, 1);

    if (fmpz_cmp(x, h) > 0)
        fmpz_sub(x, x, p);

    fmpz_clear(h);
}

/*
    The trace of y^2 = x^3 + a6 over F_p, which has j = 0.

    For p = 2 mod 3 the curve is supersingular. Otherwise 4p = A^2 + 27 B^2
    has a solution, normalised by A = 2 mod 3, and the sextic character of
    -108 a6 turns |A| into the trace of this particular sextic twist.
*/
static int
_cm_trace_j0(fmpz_t t, const fmpz_t a6, const fmpz_t p)
{
    fmpz_t D, s, A, B, d, e;
    int ok = 0;

    if (fmpz_fdiv_ui(p, 3) != 1)
        return fmpz_zero(t), 1;                 /* supersingular */

    fmpz_init(D); fmpz_init(s); fmpz_init(A);
    fmpz_init(B); fmpz_init(d); fmpz_init(e);

    fmpz_set_si(D, -27);
    fmpz_mod(D, D, p);

    if (fmpz_jacobi(D, p) == 1 && fmpz_sqrtmod(s, D, p)
            && qfb_cornacchia(A, B, p, -27, s))
    {
        if (fmpz_fdiv_ui(A, 3) == 1)
            fmpz_neg(A, A);

        fmpz_mul_si(d, a6, -108);
        fmpz_mod(d, d, p);

        fmpz_sub_ui(e, p, 1);
        fmpz_divexact_ui(e, e, 6);
        fmpz_powm(d, d, e, p);

        fmpz_mul(t, A, d);
        fmpz_mod(t, t, p);
        _center(t, p);
        ok = 1;
    }

    fmpz_clear(D); fmpz_clear(s); fmpz_clear(A);
    fmpz_clear(B); fmpz_clear(d); fmpz_clear(e);

    return ok;
}

/*
    The trace of y^2 = x^3 + a4 x over F_p, which has j = 1728.

    For p = 3 mod 4 the curve is supersingular. Otherwise 4p = A^2 + 4 B^2,
    normalised to the even member of the pair with A = 2 mod 8, and the
    quartic character of a4 selects the quartic twist.
*/
static int
_cm_trace_j1728(fmpz_t t, const fmpz_t a4, const fmpz_t p)
{
    fmpz_t D, s, A, B, e;
    int ok = 0;

    if (fmpz_fdiv_ui(p, 4) != 1)
        return fmpz_zero(t), 1;                 /* supersingular */

    fmpz_init(D); fmpz_init(s); fmpz_init(A); fmpz_init(B); fmpz_init(e);

    fmpz_set_si(D, -4);
    fmpz_mod(D, D, p);

    if (fmpz_jacobi(D, p) == 1 && fmpz_sqrtmod(s, D, p)
            && qfb_cornacchia(A, B, p, -4, s))
    {
        if (fmpz_fdiv_ui(A, 4) == 0)
            fmpz_set(A, B);

        if (fmpz_is_odd(A))
            fmpz_mul_2exp(A, A, 1);

        if (fmpz_fdiv_ui(A, 8) == 6)
            fmpz_neg(A, A);

        fmpz_sub_ui(e, p, 1);
        fmpz_tdiv_q_2exp(e, e, 2);
        fmpz_powm(e, a4, e, p);

        fmpz_mul(t, A, e);
        fmpz_mod(t, t, p);
        _center(t, p);
        ok = 1;
    }

    fmpz_clear(D); fmpz_clear(s); fmpz_clear(A); fmpz_clear(B); fmpz_clear(e);

    return ok;
}

/* how the sign convention for |t| is pinned down */
#define CM_SIGN_JACOBI 0        /* by (t | m), m odd */
#define CM_SIGN_D8     1        /* D = -8:  by t mod 16 and v mod 4 */
#define CM_SIGN_D16    2        /* D = -16: by t mod 8 */

/*
    The orders of discriminant D with class number one other than -3 and
    -4, the j-invariant of the corresponding curve, and the data of the
    rule above:

      Dc   the discriminant Cornacchia is run on. It is D itself except
           for -16 and -28, whose own forms do not represent every split
           p, so the maximal order is used and the normalisation picks the
           right member of the pair out.
      m    modulus of the Jacobi symbol on t, 1 when the kind says
           otherwise. For -27 the symbol modulo 27 is the symbol modulo 3.
      c    the constant in the Legendre symbol (c a6 | p).
      eps  the overall sign.
*/
typedef struct
{
    slong D;
    slong j;
    slong Dc;
    slong m;
    slong c;
    slong eps;
    int kind;
}
cm_disc_struct;

static const cm_disc_struct cm_discs[] = {
    {  -7, WORD(-3375),                  -7,   7, WORD(-2),      1, CM_SIGN_JACOBI },
    {  -8, WORD(8000),                   -8,   1, WORD(21),      1, CM_SIGN_D8     },
    { -11, WORD(-32768),                -11,  11, WORD(21),     -1, CM_SIGN_JACOBI },
    { -12, WORD(54000),                 -12,   3, WORD(22),     -1, CM_SIGN_JACOBI },
    { -16, WORD(287496),                 -4,   1, WORD(7),       1, CM_SIGN_D16    },
    { -19, WORD(-884736),               -19,  19, WORD(1),      -1, CM_SIGN_JACOBI },
    { -27, WORD(-12288000),             -27,   3, WORD(253),    -1, CM_SIGN_JACOBI },
    { -28, WORD(16581375),               -7,   7, WORD(-114),    1, CM_SIGN_JACOBI },
    { -43, WORD(-884736000),            -43,  43, WORD(21),     -1, CM_SIGN_JACOBI },
    { -67, WORD(-147197952000),         -67,  67, WORD(217),    -1, CM_SIGN_JACOBI },
    {-163, WORD(-262537412640768000),  -163, 163, WORD(185801), -1, CM_SIGN_JACOBI },
};

#define CM_NUM_DISCS (sizeof(cm_discs) / sizeof(cm_disc_struct))

/*
    The trace of a curve with j = j_D for one of the entries above, given
    its a6. Returns 0 if the trace could not be determined, which for a
    curve that really does have this j-invariant should not happen.
*/
static int
_cm_trace_disc(fmpz_t t, const cm_disc_struct * E, const fmpz_t a6,
        const fmpz_t p)
{
    fmpz_t D, s, a, b, u;
    slong sgn = E->eps;
    int ok = 0;

    fmpz_init(D); fmpz_init(s); fmpz_init(a); fmpz_init(b); fmpz_init(u);

    fmpz_set_si(D, E->Dc);
    fmpz_mod(D, D, p);

    if (fmpz_jacobi(D, p) != 1)
    {
        fmpz_zero(t);           /* inert: the reduction is supersingular */
        ok = 1;
        goto cleanup;
    }

    if (!fmpz_sqrtmod(s, D, p) || !qfb_cornacchia(a, b, p, E->Dc, s))
        goto cleanup;

    if (E->kind == CM_SIGN_JACOBI)
    {
        int e = n_jacobi((slong) fmpz_fdiv_ui(a, (ulong) E->m), (ulong) E->m);

        if (e == 0)
            goto cleanup;

        if (e < 0)
            sgn = -sgn;
    }
    else if (E->kind == CM_SIGN_D8)
    {
        /* here 4p = a^2 + 8 b^2 with a = 2 mod 4 */
        if (fmpz_fdiv_ui(a, 16) >= 8)
            sgn = -sgn;

        if (fmpz_fdiv_ui(b, 4) == 2)
            sgn = -sgn;
    }
    else
    {
        /* here 4p = a^2 + 4 b^2 and exactly one of a/2, b is odd; the
           trace is twice that one, normalised to 2 mod 4 */
        if (fmpz_fdiv_ui(a, 4) == 0)
            fmpz_swap(a, b);

        if (fmpz_is_odd(a))
            fmpz_mul_2exp(a, a, 1);

        if (fmpz_fdiv_ui(a, 8) == 6)
            sgn = -sgn;
    }

    fmpz_mul_si(u, a6, E->c);
    fmpz_mod(u, u, p);

    if (fmpz_is_zero(u))
        goto cleanup;

    if (fmpz_jacobi(u, p) < 0)
        sgn = -sgn;

    if (sgn > 0)
        fmpz_set(t, a);
    else
        fmpz_neg(t, a);

    ok = 1;

cleanup:
    fmpz_clear(D); fmpz_clear(s); fmpz_clear(a);
    fmpz_clear(b); fmpz_clear(u);

    return ok;
}

/*
    Is j(E) = J, without dividing? For y^2 = x^3 + a4 x + a6 the
    j-invariant is 6912 a4^3 / (4 a4^3 + 27 a6^2), so the question is
    whether 6912 a4^3 = J (4 a4^3 + 27 a6^2) in F_p. The numerator is the
    same for every J in the table, so the caller passes it reduced.
*/
static int
_j_equals(slong J, const fmpz_t num, const fmpz_t den, const fmpz_t p)
{
    fmpz_t v;
    int eq;

    fmpz_init(v);

    fmpz_mul_si(v, den, J);
    fmpz_mod(v, v, p);

    eq = fmpz_equal(num, v);

    fmpz_clear(v);

    return eq;
}

int
gr_ec_ctx_cardinality_cm(fmpz_t res, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    fmpz_t p, a4, a6, a4c, den, t;
    slong deg, i;
    int status = GR_SUCCESS, found = 0;

    if (gr_ctx_is_field(R) != T_TRUE)
        return GR_DOMAIN;

    /* the class number one theory used here is over a prime field, and the
       twist characters below are written for the short model */
    if (gr_ctx_fq_degree(&deg, R) != GR_SUCCESS || deg != 1
            || gr_ec_ctx_model(ctx) != GR_EC_SHORT_WEIERSTRASS)
        return GR_UNABLE;

    fmpz_init(p); fmpz_init(a4); fmpz_init(a6);
    fmpz_init(a4c); fmpz_init(den); fmpz_init(t);

    if (gr_ctx_fq_prime(p, R) != GR_SUCCESS || fmpz_cmp_ui(p, 3) <= 0
            || gr_get_fmpz(a4, GR_EC_A4(ctx), R) != GR_SUCCESS
            || gr_get_fmpz(a6, GR_EC_A6(ctx), R) != GR_SUCCESS)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    if (fmpz_is_zero(a4))
        found = _cm_trace_j0(t, a6, p);
    else if (fmpz_is_zero(a6))
        found = _cm_trace_j1728(t, a4, p);
    else
    {
        /* the numerator 6912 a4^3 and the denominator 4 a4^3 + 27 a6^2 of
           the j-invariant, both reduced once for the whole table */
        {
            fmpz_t w;
            fmpz_init(w);

            fmpz_powm_ui(a4c, a4, 3, p);
            fmpz_mul_ui(den, a4c, 4);
            fmpz_mul(w, a6, a6);
            fmpz_mul_ui(w, w, 27);
            fmpz_add(den, den, w);
            fmpz_mod(den, den, p);

            fmpz_mul_ui(a4c, a4c, 6912);
            fmpz_mod(a4c, a4c, p);

            fmpz_clear(w);
        }

        for (i = 0; i < (slong) CM_NUM_DISCS && !found; i++)
            if (_j_equals(cm_discs[i].j, a4c, den, p))
                found = _cm_trace_disc(t, cm_discs + i, a6, p);
    }

    if (found)
    {
        fmpz_add_ui(res, p, 1);
        fmpz_sub(res, res, t);
    }
    else
        status = GR_UNABLE;

cleanup:
    fmpz_clear(p); fmpz_clear(a4); fmpz_clear(a6);
    fmpz_clear(a4c); fmpz_clear(den); fmpz_clear(t);

    return status;
}
