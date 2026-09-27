/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "qfb.h"
#include "ulong_extras.h"
#include "nmod.h"
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

/* ------------------------------------------------------------------ */
/* class numbers 2 to 4                                               */
/* ------------------------------------------------------------------ */

/*
    For an order O_D of class number h > 1, j(E) is not rational but a root
    of the Hilbert class polynomial H_D, of degree h, and every root of H_D
    modulo p is the reduction of a j-invariant with complex multiplication
    by O_D. So if H_D(j(E)) = 0 in F_p, End(E) contains O_D, and then:

      * if p does not split in K = Q(sqrt(D0)), D0 the fundamental
        discriminant, E is supersingular and, as p > 3, t = 0 (Deuring);

      * if p splits, Frobenius is an element pi of norm p in O_K, so
        4p = t^2 + |D0| w^2, and every solution is a unit multiple of the
        one Cornacchia finds: t is +-a, and for D0 = -4 also +-2b, for
        D0 = -3 also +-(a + 3b)/2 and +-(a - 3b)/2.

    Which of those few candidates is the trace is settled with points: the
    true order kills every point, so it survives every test, and the answer
    is taken only when it is the only survivor. Unlike the class number one
    table above, this needs no calibration and no twist character, and it
    costs a couple of scalar multiplications -- nothing next to a count.

    Recognising the curve costs one evaluation of each H_D: some ten
    microseconds in all while p fits in a word, where the coefficients are
    reduced a word at a time, and about sixty above. Below
    GR_EC_CM_HILBERT_MIN_BITS a whole count costs little more than that,
    CM or not, so the check is not made.
*/

#include "cm_hilbert_table.h"

#define GR_EC_CM_HILBERT_MIN_BITS 40
#define GR_EC_CM_SELECT_POINTS 20

/* a coefficient of the table, reduced modulo p */
static void
_cm_hilbert_coeff(fmpz_t r, const cm_hilbert_coeff_struct * C, const fmpz_t p)
{
    slong i;

    fmpz_zero(r);

    for (i = C->len - 1; i >= 0; i--)
    {
        fmpz_mul_2exp(r, r, 32);
        fmpz_add_ui(r, r, cm_hilbert_words[C->off + i]);
    }

    if (C->sign < 0)
        fmpz_neg(r, r);

    fmpz_mod(r, r, p);
}

/*
    The order among the candidates q + 1 - t, decided by random points;
    0 if no single candidate is left standing.
*/
static int
_cm_select_order(fmpz_t res, const fmpz * t, slong n, const fmpz_t p,
        gr_ec_ctx_t ctx)
{
    gr_ec_point_t P, Q;
    flint_rand_t state;
    fmpz * N;
    char * alive;
    slong i, nalive = n, round;
    int ok = 0;

    N = _fmpz_vec_init(n);
    alive = flint_malloc(n);

    for (i = 0; i < n; i++)
    {
        fmpz_add_ui(N + i, p, 1);
        fmpz_sub(N + i, N + i, t + i);
        alive[i] = 1;
    }

    gr_ec_point_init(P, ctx);
    gr_ec_point_init(Q, ctx);
    flint_rand_init(state);

    for (round = 0; round < GR_EC_CM_SELECT_POINTS && nalive > 1; round++)
    {
        if (gr_ec_point_randtest(P, state, ctx) != GR_SUCCESS
                || gr_ec_point_is_inf(P, ctx) != T_FALSE)
            continue;

        for (i = 0; i < n; i++)
        {
            if (!alive[i])
                continue;

            if (gr_ec_point_mul_fmpz(Q, P, N + i, ctx) != GR_SUCCESS
                    || gr_ec_point_is_inf(Q, ctx) != T_TRUE)
            {
                alive[i] = 0;
                nalive--;
            }
        }
    }

    if (nalive == 1)
        for (i = 0; i < n; i++)
            if (alive[i])
            {
                fmpz_set(res, N + i);
                ok = 1;
            }

    flint_rand_clear(state);
    gr_ec_point_clear(Q, ctx);
    gr_ec_point_clear(P, ctx);
    flint_free(alive);
    _fmpz_vec_clear(N, n);

    return ok;
}

/* #E for a curve known to have CM by an order of fundamental disc. D0 */
static int
_cm_hilbert_count(fmpz_t res, slong D0, const fmpz_t p, gr_ec_ctx_t ctx)
{
    fmpz_t D, s, a, b;
    fmpz t[6];
    slong n = 0, i;
    int ok = 0;

    fmpz_init(D); fmpz_init(s); fmpz_init(a); fmpz_init(b);
    for (i = 0; i < 6; i++)
        fmpz_init(t + i);

    fmpz_set_si(D, D0);
    fmpz_mod(D, D, p);

    if (fmpz_jacobi(D, p) != 1)
    {
        /* inert or ramified: supersingular */
        fmpz_add_ui(res, p, 1);
        ok = 1;
        goto cleanup;
    }

    if (!fmpz_sqrtmod(s, D, p) || !qfb_cornacchia(a, b, p, D0, s))
        goto cleanup;

    /* 4p = a^2 + |D0| b^2 */
    fmpz_set(t + n++, a);

    if (D0 == -4)
        fmpz_mul_2exp(t + n++, b, 1);
    else if (D0 == -3)
    {
        fmpz_mul_ui(t + n, b, 3);
        fmpz_add(t + n, t + n, a);
        fmpz_fdiv_q_2exp(t + n, t + n, 1);
        n++;
        fmpz_mul_ui(t + n, b, 3);
        fmpz_sub(t + n, a, t + n);
        fmpz_fdiv_q_2exp(t + n, t + n, 1);
        n++;
    }

    for (i = 0, n *= 2; i < n / 2; i++)
        fmpz_neg(t + n / 2 + i, t + i);

    ok = _cm_select_order(res, t, n, p, ctx);

cleanup:
    for (i = 0; i < 6; i++)
        fmpz_clear(t + i);
    fmpz_clear(D); fmpz_clear(s); fmpz_clear(a); fmpz_clear(b);

    return ok;
}

/* the general version, for p of any size */
static int
_cm_hilbert_fmpz(fmpz_t res, const fmpz_t num, const fmpz_t den,
        const fmpz_t p, gr_ec_ctx_t ctx)
{
    fmpz np[5], dp[5];
    fmpz_t c, v, w;
    slong i, k;
    int found = 0;

    for (i = 0; i < 5; i++)
    {
        fmpz_init(np + i);
        fmpz_init(dp + i);
    }
    fmpz_init(c); fmpz_init(v); fmpz_init(w);

    fmpz_one(np + 0);
    fmpz_one(dp + 0);
    for (i = 1; i < 5; i++)
    {
        fmpz_mul(np + i, np + i - 1, num);
        fmpz_mod(np + i, np + i, p);
        fmpz_mul(dp + i, dp + i - 1, den);
        fmpz_mod(dp + i, dp + i, p);
    }

    for (k = 0; k < CM_HILBERT_NUM && !found; k++)
    {
        const cm_hilbert_struct * H = cm_hilbert + k;

        /* den^h H(num/den) = num^h + sum c_i num^i den^(h-i) */
        fmpz_set(v, np + H->h);

        for (i = 0; i < H->h; i++)
        {
            _cm_hilbert_coeff(c, cm_hilbert_coeffs + H->coeff + i, p);
            fmpz_mul(w, np + i, dp + H->h - i);
            fmpz_addmul(v, c, w);
        }

        fmpz_mod(v, v, p);

        if (fmpz_is_zero(v))
            found = _cm_hilbert_count(res, H->D0, p, ctx);
    }

    for (i = 0; i < 5; i++)
    {
        fmpz_clear(np + i);
        fmpz_clear(dp + i);
    }
    fmpz_clear(c); fmpz_clear(v); fmpz_clear(w);

    return found;
}

/* the same while p fits in a word: num, den and the result modulo p */
static slong
_cm_hilbert_find_ui(ulong num, ulong den, nmod_t mod)
{
    ulong np[5], dp[5], two32, c, v;
    slong i, k, w;

    two32 = nmod_set_ui(UWORD(1) << 16, mod);
    two32 = nmod_mul(two32, two32, mod);

    np[0] = dp[0] = nmod_set_ui(1, mod);
    for (i = 1; i < 5; i++)
    {
        np[i] = nmod_mul(np[i - 1], num, mod);
        dp[i] = nmod_mul(dp[i - 1], den, mod);
    }

    for (k = 0; k < CM_HILBERT_NUM; k++)
    {
        const cm_hilbert_struct * H = cm_hilbert + k;

        v = np[H->h];

        for (i = 0; i < H->h; i++)
        {
            const cm_hilbert_coeff_struct * C = cm_hilbert_coeffs + H->coeff + i;

            c = 0;
            for (w = C->len - 1; w >= 0; w--)
                c = nmod_add(nmod_mul(c, two32, mod),
                        nmod_set_ui(cm_hilbert_words[C->off + w], mod), mod);

            if (C->sign < 0)
                c = nmod_neg(c, mod);

            v = nmod_add(v, nmod_mul(c, nmod_mul(np[i], dp[H->h - i], mod), mod), mod);
        }

        if (v == 0)
            return k;
    }

    return -1;
}

/*
    Is j = num/den a root of some H_D in the table? Then count. Returns 1
    and sets res if so.
*/
static int
_cm_hilbert(fmpz_t res, const fmpz_t num, const fmpz_t den, const fmpz_t p,
        gr_ec_ctx_t ctx)
{
    if (fmpz_abs_fits_ui(p))
    {
        nmod_t mod;
        slong k;

        nmod_init(&mod, fmpz_get_ui(p));
        k = _cm_hilbert_find_ui(fmpz_get_ui(num), fmpz_get_ui(den), mod);

        if (k < 0)
            return 0;

        if (_cm_hilbert_count(res, cm_hilbert[k].D0, p, ctx))
            return 1;

        /* j is a root of H_D for more than one D modulo p (a rare
           collision); let the general version look at all of them */
    }

    return _cm_hilbert_fmpz(res, num, den, p, ctx);
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

        /* class numbers 2 to 4: these give the order directly */
        if (!found && fmpz_bits(p) >= GR_EC_CM_HILBERT_MIN_BITS
                && _cm_hilbert(res, a4c, den, p, ctx))
            goto cleanup;
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
