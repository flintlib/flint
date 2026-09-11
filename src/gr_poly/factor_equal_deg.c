/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "thread_support.h"
#include "fmpz.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "factor_impl.h"

/*
    Equal degree factorization of a monic squarefree polynomial f of degree
    n = r d which is a product of r irreducible factors of degree d over
    F_q, q = p^k.

    We use the Cantor-Zassenhaus / von zur Gathen-Shoup approach: for a
    random a in F_q[x]/(f), compute b = Tr(a) = a + a^q + ... + a^(q^(d-1)),
    which lies in the Berlekamp subalgebra (b is congruent to an element of
    F_q modulo each irreducible factor). Then a proper factor of f is found
    with probability >= 1/2 (roughly) by computing gcd(b^((q-1)/2) - 1, f)
    for odd q, or gcd(b + b^2 + ... + b^(2^(k-1)), f) for even q.

    The trace is computed using modular composition with x^q mod f
    (via a precomputed matrix), or with repeated exponentiation when q is
    small relative to the degree, or with the doubling algorithm of
    von zur Gathen and Shoup for large d.

    For d = 1 (root finding), a is chosen to be a random linear
    polynomial x + c so that powering can be done more cheaply.

    Compared to the fmpz_mod_poly implementation by Daniel Schultz, this
    version omits the Berlekamp-Massey based "deterministic simplification"
    which requires a generic Berlekamp-Massey implementation.
*/

/* Should we use modular composition with a precomputed matrix rather than
   exponentiation to apply the Frobenius? */
static int
_use_composition(const fmpz_t q, slong lenf)
{
    return fmpz_bits(q) > ((n_sqrt(lenf - 1) + 1) * 3) / 4;
}

typedef struct
{
    gr_poly_vec_t f;    /* polynomials to be factored */
    gr_poly_vec_t xp;   /* xp[i] = x^q mod f[i] (unused when d = 1) */
}
queue_struct;

typedef queue_struct queue_t[1];

static void
queue_init(queue_t Q, gr_ctx_t ctx)
{
    gr_poly_vec_init(Q->f, 0, ctx);
    gr_poly_vec_init(Q->xp, 0, ctx);
}

static void
queue_clear(queue_t Q, gr_ctx_t ctx)
{
    gr_poly_vec_clear(Q->f, ctx);
    gr_poly_vec_clear(Q->xp, ctx);
}

/* Push a monic divisor piece of the input to the queue if it still needs
   to be factored, or to the result if it is irreducible (degree d).
   The content of piece is clobbered. */
static int
_push_piece(gr_poly_vec_t res, queue_t Q, gr_poly_t piece, slong d,
        const gr_poly_t xp, gr_ctx_t ctx)
{
    slong Qlen = Q->f->length;
    int status = GR_SUCCESS;

    if (piece->length - 1 > d)
    {
        gr_poly_vec_fit_length(Q->f, Qlen + 1, ctx);
        gr_poly_vec_fit_length(Q->xp, Qlen + 1, ctx);
        gr_poly_swap(Q->f->entries + Qlen, piece, ctx);
        if (d > 1)
            status |= gr_poly_rem(Q->xp->entries + Qlen, xp, Q->f->entries + Qlen, ctx);
        Q->f->length = Q->xp->length = Qlen + 1;
    }
    else if (piece->length - 1 == d)
    {
        gr_poly_vec_append_swap(res, piece, ctx);
    }
    else
    {
        status = GR_UNABLE;
    }

    return status;
}

/* Given g | f, push g and f/g. Clobbers f and g. */
static int
_add_split(gr_poly_vec_t res, queue_t Q, gr_poly_t f, gr_poly_t g,
        slong d, const gr_poly_t xp, gr_poly_t tmp, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    status |= gr_poly_divexact(tmp, f, g, ctx);

    /* process the larger piece last (i.e. put it on top of the stack) */
    if (tmp->length < g->length)
        gr_poly_swap(tmp, g, ctx);

    status |= _push_piece(res, Q, g, d, xp, ctx);
    status |= _push_piece(res, Q, tmp, d, xp, ctx);

    return status;
}

/*
    Berlekamp-Massey algorithm over a field: given the sequence
    s_0, ..., s_{N-1}, computes the minimal polynomial
    h(y) = y^L + c_1 y^(L-1) + ... + c_L of the linear recurrence
    s_n + c_1 s_{n-1} + ... + c_L s_{n-L} = 0 (valid for L <= n < N),
    where L is minimal. Assumes N >= 2L so that h is uniquely determined.
*/
static int
_gr_poly_berlekamp_massey(gr_poly_t h, gr_srcptr s, slong N, gr_ctx_t ctx)
{
    gr_poly_t C, B, T;
    gr_ptr dd, b, u;
    slong n, i, L, m;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    gr_poly_init(C, ctx);
    gr_poly_init(B, ctx);
    gr_poly_init(T, ctx);
    GR_TMP_INIT3(dd, b, u, ctx);

    /* C(x) = 1 + c_1 x + ... + c_L x^L is the connection polynomial */
    status |= gr_poly_one(C, ctx);
    status |= gr_poly_one(B, ctx);
    status |= gr_one(b, ctx);
    L = 0;
    m = 1;

    for (n = 0; n < N && status == GR_SUCCESS; n++)
    {
        /* discrepancy dd = s_n + sum_{i=1}^L c_i s_{n-i} */
        status |= gr_set(dd, GR_ENTRY(s, n, sz), ctx);
        for (i = 1; i <= L && i < C->length; i++)
            status |= gr_addmul(dd, GR_ENTRY(C->coeffs, i, sz), GR_ENTRY(s, n - i, sz), ctx);

        if (gr_is_zero(dd, ctx) == T_TRUE)
        {
            m++;
        }
        else
        {
            /* T = C - (dd/b) x^m B */
            status |= gr_div(u, dd, b, ctx);
            status |= gr_neg(u, u, ctx);
            status |= gr_poly_shift_left(T, B, m, ctx);
            status |= gr_poly_mul_scalar(T, T, u, ctx);
            status |= gr_poly_add(T, T, C, ctx);

            if (2 * L <= n)
            {
                L = n + 1 - L;
                gr_poly_swap(B, C, ctx);
                status |= gr_set(b, dd, ctx);
                m = 1;
            }
            else
            {
                m++;
            }

            gr_poly_swap(C, T, ctx);
        }
    }

    /* h = reversal of C to length L + 1 */
    gr_poly_fit_length(h, L + 1, ctx);
    for (i = 0; i <= L; i++)
    {
        if (i < C->length)
            status |= gr_set(GR_ENTRY(h->coeffs, L - i, sz), GR_ENTRY(C->coeffs, i, sz), ctx);
        else
            status |= gr_zero(GR_ENTRY(h->coeffs, L - i, sz), ctx);
    }
    _gr_poly_set_length_normalise(h, L + 1, ctx);

    gr_poly_clear(C, ctx);
    gr_poly_clear(B, ctx);
    gr_poly_clear(T, ctx);
    GR_TMP_CLEAR3(dd, b, u, ctx);

    return status;
}

/*
    Given a monic polynomial f, an element b of F_q[x]/(f) which is
    congruent to a constant modulo each irreducible factor of f, and the
    list cs of the values which b takes (each value being attained for at
    least one factor), splits f into the pieces gcd(f, b - c), c in cs,
    which are pushed to the queue or the result. This is done by divide
    and conquer: with h_1 = prod_{c in cs_1} (y - c) for the first half
    cs_1 of the values, g_1 = gcd(f, h_1(b)) is the product of the factors
    where b takes a value in cs_1, computed with modular composition.
    Clobbers f and b.
*/
static int
_gr_poly_edf_split_by_values(gr_poly_vec_t res, queue_t Q,
    gr_poly_t f, gr_poly_t b, gr_srcptr cs, slong num, slong d,
    const gr_poly_t xp, gr_ctx_t ctx)
{
    gr_poly_t h, t, g1, g2, b2;
    slong mid;
    int status = GR_SUCCESS;

    if (num <= 1 || f->length - 1 <= d)
        return _push_piece(res, Q, f, d, xp, ctx);

    gr_poly_init(h, ctx);
    gr_poly_init(t, ctx);
    gr_poly_init(g1, ctx);
    gr_poly_init(g2, ctx);
    gr_poly_init(b2, ctx);

    mid = num / 2;

    /* h = prod (y - c) over the first half of the values */
    gr_poly_fit_length(h, mid + 1, ctx);
    status |= _gr_poly_product_roots(h->coeffs, cs, mid, ctx);
    _gr_poly_set_length(h, mid + 1, ctx);

    /* t = h(b) mod f */
    status |= gr_poly_rem(b, b, f, ctx);
    status |= gr_poly_compose_mod(t, h, b, f, ctx);
    status |= gr_poly_gcd(g1, t, f, ctx);

    if (status == GR_SUCCESS && g1->length > 1 && g1->length < f->length)
    {
        status |= gr_poly_divexact(g2, f, g1, ctx);
        status |= gr_poly_rem(b2, b, g2, ctx);
        status |= gr_poly_rem(b, b, g1, ctx);
        status |= _gr_poly_edf_split_by_values(res, Q, g1, b, cs, mid, d, xp, ctx);
        status |= _gr_poly_edf_split_by_values(res, Q, g2, b2, GR_ENTRY(cs, mid, ctx->sizeof_elem), num - mid, d, xp, ctx);
    }
    else
    {
        /* should not happen for valid input; let the main loop try again */
        status |= _push_piece(res, Q, f, d, xp, ctx);
    }

    gr_poly_clear(h, ctx);
    gr_poly_clear(t, ctx);
    gr_poly_clear(g1, ctx);
    gr_poly_clear(g2, ctx);
    gr_poly_clear(b2, ctx);

    return status;
}

/*
    Given b in the Berlekamp subalgebra of F_q[x]/(f) (i.e. b is congruent
    to a constant c_i modulo each irreducible factor f_i), split f into
    the pieces gcd(f, b - c) for the distinct values c. The values c are
    found as the roots of the minimal polynomial of the sequence
    s_i = [x^0] (b^i mod f), which divides prod_i (y - c_i) and is
    computed with the Berlekamp-Massey algorithm from 2r terms, r being
    an upper bound for the number of irreducible factors.
    (Algorithm 3 in David Marquis' thesis "Deterministic Factorization of
    Polynomials over Finite Fields", as implemented for fmpz_mod_poly by
    Daniel Schultz.) Pushes the resulting pieces to the queue or the
    result. Sets *found = 0 if no proper factor was found, in which case
    f is left unchanged. Otherwise f is clobbered.
*/
static int
_gr_poly_edf_split_minpoly(gr_poly_vec_t res, queue_t Q, int * found,
    gr_poly_t f, const gr_poly_preinv_t P, const gr_poly_t b, slong d, slong r,
    const gr_poly_t xp, const fmpz_t q, flint_rand_t state,
    gr_poly_t t, gr_poly_t g, gr_poly_t h, gr_ctx_t ctx)
{
    gr_ptr s, w;
    gr_poly_vec_t lin;
    slong i, n = f->length - 1, N = 2 * r;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    *found = 0;

    GR_TMP_INIT_VEC(s, N, ctx);
    GR_TMP_INIT_VEC(w, n, ctx);
    gr_poly_vec_init(lin, 0, ctx);

    /* s_i = <w, b^i mod f> for a random linear functional w: a fixed
       functional (for example the constant coefficient) can be
       degenerate on the powers of b, in which case the sequence
       satisfies a recurrence shorter than the minimal polynomial of b
       and no splitting is obtained. */
    for (i = 0; i < n; i++)
        status |= gr_randtest(GR_ENTRY(w, i, sz), state, ctx);

    status |= gr_poly_one(t, ctx);
    for (i = 0; i < N && status == GR_SUCCESS; i++)
    {
        if (i != 0)
            status |= gr_poly_preinv_mulmod(t, t, b, P, ctx);
        status |= _gr_vec_dot(GR_ENTRY(s, i, sz), NULL, 0, w, t->coeffs,
            FLINT_MIN(n, t->length), ctx);
    }

    status |= _gr_poly_berlekamp_massey(h, s, N, ctx);

    /* h has degree between 1 and r and (for valid input) is a product of
       distinct linear factors; find its roots (d = 1 equal degree
       factorization). To guarantee termination for invalid input, we
       first extract the product of the linear factors gcd(h, x^q - x). */
    if (status == GR_SUCCESS && h->length >= 2)
    {
        if (h->length == 2)
        {
            status |= gr_poly_vec_append(lin, h, ctx);
        }
        else
        {
            status |= gr_poly_gen(t, ctx);
            status |= gr_poly_powmod_fmpz_binexp(g, t, q, h, ctx);
            status |= gr_poly_sub(g, g, t, ctx);
            status |= gr_poly_gcd(g, g, h, ctx);

            if (status == GR_SUCCESS && g->length >= 2)
                status |= _gr_poly_factor_equal_deg_with_frob(lin, g, 1, g, state, ctx);
        }
    }

    if (status == GR_SUCCESS && lin->length >= 2)
    {
        /* Collect the roots c_i and split f recursively into the pieces
           gcd(f, b - c_i). */
        gr_ptr cs;
        GR_TMP_INIT_VEC(cs, lin->length, ctx);
        for (i = 0; i < lin->length; i++)
            status |= gr_neg(GR_ENTRY(cs, i, sz), lin->entries[i].coeffs, ctx);

        *found = 1;
        status |= gr_poly_set(t, b, ctx);
        status |= _gr_poly_edf_split_by_values(res, Q, f, t, cs, lin->length, d, xp, ctx);

        GR_TMP_CLEAR_VEC(cs, lin->length, ctx);
    }

    GR_TMP_CLEAR_VEC(s, N, ctx);
    GR_TMP_CLEAR_VEC(w, n, ctx);
    gr_poly_vec_clear(lin, ctx);

    return status;
}

/*
    Compute a = b + b^q + ... + b^(q^(d-1)) mod f, i.e. the trace of b
    into the Berlekamp subalgebra. Clobbers b, uses xi and t as scratch.
*/
static int
_gr_poly_trace(gr_poly_t a, gr_poly_t b, slong d,
    const gr_poly_t xq, const gr_poly_t f, const gr_poly_preinv_t P,
    const fmpz_t q, gr_poly_t xi, gr_poly_t t, gr_ctx_t ctx)
{
    slong i;
    int status = GR_SUCCESS;
    int composition = _use_composition(q, f->length);

    if (d < 2)
    {
        gr_poly_swap(a, b, ctx);
    }
    else if (!composition)
    {
        /* linear iteration with exponentiation */
        status |= gr_poly_preinv_powmod_fmpz_sliding(xi, b, q, 0, P, ctx);
        status |= gr_poly_add(a, b, xi, ctx);
        for (i = 2; i < d && status == GR_SUCCESS; i++)
        {
            status |= gr_poly_preinv_powmod_fmpz_sliding(t, xi, q, 0, P, ctx);
            gr_poly_swap(xi, t, ctx);
            status |= gr_poly_add(a, a, xi, ctx);
        }
    }
    else if (d < 16)
    {
        /* linear iteration with modular composition */
        gr_mat_t H;
        gr_mat_init(H, n_sqrt(f->length - 1) + 1, f->length - 1, ctx);
        status |= gr_poly_preinv_precompute_matrix(H, xq, P, ctx);
        status |= gr_poly_preinv_compose_mod_brent_kung_precomp(xi, b, H, P, ctx);
        status |= gr_poly_add(a, b, xi, ctx);
        for (i = 2; i < d && status == GR_SUCCESS; i++)
        {
            status |= gr_poly_preinv_compose_mod_brent_kung_precomp(t, xi, H, P, ctx);
            gr_poly_swap(xi, t, ctx);
            status |= gr_poly_add(a, a, xi, ctx);
        }
        gr_mat_clear(H, ctx);
    }
    else
    {
        /*
            Doubling algorithm from
                Computing Frobenius maps and factoring polynomials
                    von zur Gathen and Shoup
            Invariant: with j = 2^i the current chunk size, b = T_j(alpha),
            xi = x^(q^j) mod f, and a = T_m(alpha) where m is the number
            formed by the bits of d already processed (a = 0 if m = 0).
            Here T_m(alpha) = alpha + alpha^q + ... + alpha^(q^(m-1)).
        */
        gr_mat_t H;
        gr_mat_init(H, n_sqrt(f->length - 1) + 1, f->length - 1, ctx);

        status |= gr_poly_zero(a, ctx);
        status |= gr_poly_set(xi, xq, ctx);

        while (status == GR_SUCCESS)
        {
            status |= gr_poly_preinv_precompute_matrix(H, xi, P, ctx);

            if (d % 2 == 0)
            {
                status |= gr_poly_preinv_compose_mod_brent_kung_precomp(t, b, H, P, ctx);
                status |= gr_poly_add(b, b, t, ctx);
            }
            else if (a->length == 0)
            {
                gr_poly_swap(a, b, ctx);
                status |= gr_poly_preinv_compose_mod_brent_kung_precomp(b, a, H, P, ctx);
                status |= gr_poly_add(b, b, a, ctx);
            }
            else
            {
                status |= gr_poly_preinv_compose_mod_brent_kung_precomp(t, a, H, P, ctx);
                status |= gr_poly_add(a, b, t, ctx);
                status |= gr_poly_preinv_compose_mod_brent_kung_precomp(t, b, H, P, ctx);
                status |= gr_poly_add(b, b, t, ctx);
            }

            d = d / 2;

            if (d <= 1)
            {
                if (a->length == 0)
                {
                    gr_poly_swap(a, b, ctx);
                    break;
                }

                status |= gr_poly_preinv_compose_mod_brent_kung_precomp(t, xi, H, P, ctx);
                gr_poly_swap(xi, t, ctx);
                status |= gr_poly_preinv_precompute_matrix(H, xi, P, ctx);
                status |= gr_poly_preinv_compose_mod_brent_kung_precomp(t, a, H, P, ctx);
                status |= gr_poly_add(a, t, b, ctx);
                break;
            }

            status |= gr_poly_preinv_compose_mod_brent_kung_precomp(t, xi, H, P, ctx);
            gr_poly_swap(xi, t, ctx);
        }

        gr_mat_clear(H, ctx);
    }

    return status;
}

/*
    Attempt to find a proper factor g of the monic squarefree polynomial f,
    a product of irreducibles of degree d. If a factor is found, sets
    *found = 1, otherwise sets *found = 0.

    Requires the precomputed modulus P for f, and xq = x^q mod f
    (only when d > 1). Requires q = p^k, halfq = (q-1)/2 for odd q.
    Uses a, b, t as temporaries.
*/
/*
    Computes a random element b of the Berlekamp subalgebra of F_q[x]/(f):
    when d = 1, b is a random linear polynomial (all of F_q[x]/(f) is the
    Berlekamp subalgebra); otherwise b is the trace of a random polynomial.
    Sets *nonconst = 0 if b is a constant (a new element should be
    generated). Uses a, g, t as temporaries.
*/
static int
_gr_poly_edf_random_element(gr_poly_t b, int * nonconst, flint_rand_t state,
    const gr_poly_t f, const gr_poly_preinv_t P, const gr_poly_t xq, slong d,
    const fmpz_t q, gr_poly_t a, gr_poly_t g, gr_poly_t t, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    int char2 = fmpz_is_even(q);
    slong n = f->length - 1;

    *nonconst = 0;

    if (d == 1)
    {
        /* a = alpha x + c with random alpha != 0 (alpha = 1 in odd characteristic) */
        gr_poly_fit_length(b, 2, ctx);
        status |= gr_randtest(b->coeffs, state, ctx);
        if (char2)
            status |= gr_randtest_not_zero(GR_ENTRY(b->coeffs, 1, ctx->sizeof_elem), state, ctx);
        else
            status |= gr_one(GR_ENTRY(b->coeffs, 1, ctx->sizeof_elem), ctx);
        _gr_poly_set_length(b, 2, ctx);
        *nonconst = 1;
    }
    else
    {
        do
        {
            status |= gr_poly_randtest(a, state, n, ctx);
            if (status != GR_SUCCESS)
                return status;
        }
        while (a->length <= 1);

        status |= _gr_poly_trace(b, a, d, xq, f, P, q, g, t, ctx);
        if (status != GR_SUCCESS)
            return status;

        /* b must be nonconstant mod f to give a nontrivial split */
        *nonconst = (b->length > 1);
    }

    return status;
}

/*
    Attempt to find a proper factor g of the monic squarefree polynomial f,
    a product of irreducibles of degree d, from a nonconstant element b
    of the Berlekamp subalgebra. If a factor is found, sets *found = 1,
    otherwise sets *found = 0.

    Requires the precomputed modulus P for f. Requires q = p^k,
    halfq = (q-1)/2 for odd q. Uses a, t as temporaries.
*/
static int
_gr_poly_edf_split_power(gr_poly_t g, int * found,
    const gr_poly_t f, const gr_poly_preinv_t P, const gr_poly_t b, slong d,
    const fmpz_t q, const fmpz_t halfq, slong k,
    gr_poly_t a, gr_poly_t t, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i;
    int char2 = fmpz_is_even(q);
    slong n = f->length - 1;

    *found = 0;

    if (char2)
    {
        /* t = b + b^2 + ... + b^(2^(k-1)) mod f */
        status |= gr_poly_rem(t, b, f, ctx);
        status |= gr_poly_set(a, t, ctx);
        for (i = 1; i < k && status == GR_SUCCESS; i++)
        {
            status |= gr_poly_preinv_mulmod(a, a, a, P, ctx);
            status |= gr_poly_add(t, t, a, ctx);
        }
    }
    else if (d == 1 && fmpz_bits(q) > 16 && n > 2)
    {
        /* t = (x + c)^halfq mod f = R(x + c), R = x^halfq mod f(x - c) */
        gr_ptr c;
        GR_TMP_INIT(c, ctx);
        status |= gr_neg(c, b->coeffs, ctx);
        status |= gr_poly_taylor_shift(a, f, c, ctx);
        {
            gr_poly_preinv_t Pa;
            gr_poly_preinv_init(Pa, ctx);
            status |= gr_poly_preinv_set(Pa, a, ctx);
            status |= gr_poly_preinv_powmod_x_fmpz(t, halfq, Pa, ctx);
            gr_poly_preinv_clear(Pa, ctx);
        }
        status |= gr_poly_taylor_shift(t, t, b->coeffs, ctx);
        GR_TMP_CLEAR(c, ctx);
        status |= gr_poly_sub_ui(t, t, 1, ctx);
    }
    else
    {
        status |= gr_poly_preinv_powmod_fmpz_sliding(t, b, halfq, 0, P, ctx);
        status |= gr_poly_sub_ui(t, t, 1, ctx);
    }

    if (status != GR_SUCCESS)
        return status;

    status |= gr_poly_gcd(g, t, f, ctx);

    if (status != GR_SUCCESS)
        return status;

    if (g->length > 1 && g->length < f->length)
        *found = 1;

    return status;
}

/*
    Splits one piece f (monic, squarefree, all factors of degree d),
    with xq = x^q mod f, pushing the results to res (irreducible) or Q
    (to be split further). f and xq are clobbered.
*/
static int
_gr_poly_edf_piece(gr_poly_vec_t res, queue_t Q, gr_poly_t f, gr_poly_t xq,
    slong d, const fmpz_t q, const fmpz_t halfq, slong k, flint_rand_t state, gr_ctx_t ctx)
{
    gr_poly_t a, b, t, g, h;
    gr_poly_preinv_t P;
    slong n, r, attempts;
    int status = GR_SUCCESS;
    int found, use_minpoly;

    gr_poly_init(a, ctx);
    gr_poly_init(b, ctx);
    gr_poly_init(t, ctx);
    gr_poly_init(g, ctx);
    gr_poly_init(h, ctx);
    gr_poly_preinv_init(P, ctx);

    status |= gr_poly_preinv_set(P, f, ctx);

    n = f->length - 1;
    r = n / d;

    /* The minimal polynomial approach finds all factors at once from
       a single random element using 2r multiplications mod f, while
       the exponentiation approach needs log2(q) multiplications
       per split (and log2(r) rounds of splitting). */
    use_minpoly = (d > 1 && fmpz_bits(q) >= 8 &&
                   r <= 8 * (slong) FLINT_BIT_COUNT(d + 1) + (slong) fmpz_bits(q));

    found = 0;
    attempts = 0;
    while (status == GR_SUCCESS && !found)
    {
        /* Each attempt succeeds with probability >= 4/9 for valid
           inputs; guard against looping forever on invalid ones. */
        if (++attempts > 1000)
        {
            status = GR_UNABLE;
            break;
        }

        status |= _gr_poly_edf_random_element(b, &found, state, f, P, xq, d, q, a, g, t, ctx);
        if (status != GR_SUCCESS || !found)
            continue;

        if (use_minpoly)
        {
            status |= _gr_poly_edf_split_minpoly(res, Q, &found, f, P, b, d, r, xq, q, state, t, g, h, ctx);
        }
        else
        {
            status |= _gr_poly_edf_split_power(g, &found, f, P, b, d, q, halfq, k, a, t, ctx);
            if (status == GR_SUCCESS && found)
                status |= _add_split(res, Q, f, g, d, xq, t, ctx);
        }
    }

    gr_poly_clear(a, ctx);
    gr_poly_clear(b, ctx);
    gr_poly_clear(t, ctx);
    gr_poly_clear(g, ctx);
    gr_poly_clear(h, ctx);
    gr_poly_preinv_clear(P, ctx);

    return status;
}

/* Parallel processing of several queue pieces: each worker has its own
   piece, random state and outputs. */
typedef struct
{
    gr_poly_struct * pieces;
    gr_poly_struct * xps;
    gr_poly_vec_struct * res;
    queue_struct * Q;
    flint_rand_struct * states;
    slong d, k;
    const fmpz * q;
    const fmpz * halfq;
    gr_ctx_struct * ctx;
    int * status;
}
edf_piece_args_t;

static void
_edf_piece_worker(slong i, void * arg)
{
    edf_piece_args_t * a = arg;
    a->status[i] = _gr_poly_edf_piece(a->res + i, a->Q + i, a->pieces + i, a->xps + i,
        a->d, a->q, a->halfq, a->k, a->states + i, a->ctx);
}

int
_gr_poly_factor_equal_deg_with_frob(gr_poly_vec_t res, const gr_poly_t ff,
    slong d, const gr_poly_t frob, flint_rand_t state, gr_ctx_t ctx)
{
    fmpz_t q, halfq;
    slong k, n;
    gr_poly_t f, xq, g;
    queue_t Q;
    int status = GR_SUCCESS;

    n = ff->length - 1;

    if (d <= 0 || n % d != 0)
        return GR_DOMAIN;

    if (n == d)
        return gr_poly_vec_append(res, ff, ctx);

    fmpz_init(q);
    fmpz_init(halfq);

    status |= _gr_poly_factor_ff_info(q, NULL, &k, ctx);
    if (status != GR_SUCCESS)
    {
        fmpz_clear(q);
        fmpz_clear(halfq);
        return status;
    }

    if (fmpz_is_odd(q))
    {
        fmpz_sub_ui(halfq, q, 1);
        fmpz_fdiv_q_2exp(halfq, halfq, 1);
    }

    gr_poly_init(f, ctx);
    gr_poly_init(xq, ctx);
    gr_poly_init(g, ctx);
    queue_init(Q, ctx);

    status |= gr_poly_vec_append(Q->f, ff, ctx);
    if (d > 1)
        status |= gr_poly_vec_append(Q->xp, frob, ctx);
    else
        gr_poly_vec_set_length(Q->xp, 1, ctx);

    /* Zero roots can be split off directly. */
    if (d == 1 && gr_is_zero(ff->coeffs, ctx) == T_TRUE)
    {
        slong i = 1;
        while (i < ff->length && gr_is_zero(GR_ENTRY(ff->coeffs, i, ctx->sizeof_elem), ctx) == T_TRUE)
            i++;

        status |= gr_poly_gen(g, ctx);
        status |= gr_poly_vec_append(res, g, ctx);
        status |= gr_poly_shift_right(Q->f->entries, ff, i, ctx);

        if (Q->f->entries->length - 1 == d)
        {
            status |= gr_poly_vec_append(res, Q->f->entries, ctx);
            Q->f->length = 0;
            Q->xp->length = 0;
        }
    }

    while (status == GR_SUCCESS && Q->f->length > 0)
    {
        slong Qlen = Q->f->length;
        slong nthreads = flint_get_num_available_threads();
        slong np = FLINT_MIN(Qlen, nthreads);
        slong work = 0;
        slong i;

        /* total work in the pieces that would be processed in parallel;
           the pieces are generally of different sizes, so looking at the
           top piece alone is not meaningful */
        for (i = 0; i < np; i++)
            work += Q->f->entries[Qlen - 1 - i].length - 1;

        if (np > 1 && work > gr_poly_factor_threaded_cutoff && gr_ctx_is_threadsafe(ctx) == T_TRUE)
        {
            /* process the np pieces on top of the stack in parallel,
               each worker splitting its piece into its own output */
            edf_piece_args_t args;
            slong j;

            args.pieces = flint_malloc(np * sizeof(gr_poly_struct));
            args.xps = flint_malloc(np * sizeof(gr_poly_struct));
            args.res = flint_malloc(np * sizeof(gr_poly_vec_struct));
            args.Q = flint_malloc(np * sizeof(queue_struct));
            args.states = flint_malloc(np * sizeof(flint_rand_struct));
            args.status = flint_calloc(np, sizeof(int));
            args.d = d; args.q = q; args.halfq = halfq; args.k = k; args.ctx = ctx;

            for (i = 0; i < np; i++)
            {
                gr_poly_init(args.pieces + i, ctx);
                gr_poly_init(args.xps + i, ctx);
                gr_poly_swap(args.pieces + i, Q->f->entries + Qlen - 1 - i, ctx);
                gr_poly_swap(args.xps + i, Q->xp->entries + Qlen - 1 - i, ctx);
                gr_poly_vec_init(args.res + i, 0, ctx);
                queue_init(args.Q + i, ctx);
                flint_rand_init(args.states + i);
                flint_rand_set_seed(args.states + i, n_randlimb(state), n_randlimb(state));
            }
            Q->f->length = Q->xp->length = Qlen - np;

            flint_parallel_do(_edf_piece_worker, &args, np, -1, FLINT_PARALLEL_UNIFORM);

            for (i = 0; i < np; i++)
            {
                status |= args.status[i];
                for (j = 0; j < args.res[i].length; j++)
                    gr_poly_vec_append_swap(res, args.res[i].entries + j, ctx);
                for (j = 0; j < args.Q[i].f->length; j++)
                {
                    slong Ql = Q->f->length;
                    gr_poly_vec_fit_length(Q->f, Ql + 1, ctx);
                    gr_poly_vec_fit_length(Q->xp, Ql + 1, ctx);
                    gr_poly_swap(Q->f->entries + Ql, args.Q[i].f->entries + j, ctx);
                    gr_poly_swap(Q->xp->entries + Ql, args.Q[i].xp->entries + j, ctx);
                    Q->f->length = Q->xp->length = Ql + 1;
                }
                gr_poly_clear(args.pieces + i, ctx);
                gr_poly_clear(args.xps + i, ctx);
                gr_poly_vec_clear(args.res + i, ctx);
                queue_clear(args.Q + i, ctx);
                flint_rand_clear(args.states + i);
            }

            flint_free(args.pieces);
            flint_free(args.xps);
            flint_free(args.res);
            flint_free(args.Q);
            flint_free(args.states);
            flint_free(args.status);
        }
        else
        {
            gr_poly_swap(f, Q->f->entries + Qlen - 1, ctx);
            gr_poly_swap(xq, Q->xp->entries + Qlen - 1, ctx);
            Q->f->length = Q->xp->length = Qlen - 1;

            status |= _gr_poly_edf_piece(res, Q, f, xq, d, q, halfq, k, state, ctx);
        }
    }

    fmpz_clear(q);
    fmpz_clear(halfq);
    gr_poly_clear(f, ctx);
    gr_poly_clear(xq, ctx);
    gr_poly_clear(g, ctx);
    queue_clear(Q, ctx);

    return status;
}

int
gr_poly_factor_equal_deg_prob(gr_poly_t factor, flint_rand_t state,
    const gr_poly_t pol, slong d, gr_ctx_t ctx)
{
    fmpz_t q, halfq;
    slong k, n;
    gr_poly_t f, xq, a, b, t;
    gr_poly_preinv_t P;
    int status = GR_SUCCESS;
    int found = 0;

    n = pol->length - 1;

    if (d <= 0 || n < 2 || n % d != 0)
        return GR_DOMAIN;

    if (n == d)
        return gr_poly_one(factor, ctx);

    fmpz_init(q);
    fmpz_init(halfq);

    status |= _gr_poly_factor_ff_info(q, NULL, &k, ctx);
    if (status != GR_SUCCESS)
    {
        fmpz_clear(q);
        fmpz_clear(halfq);
        return status;
    }

    if (fmpz_is_odd(q))
    {
        fmpz_sub_ui(halfq, q, 1);
        fmpz_fdiv_q_2exp(halfq, halfq, 1);
    }

    gr_poly_init(f, ctx);
    gr_poly_preinv_init(P, ctx);
    gr_poly_init(xq, ctx);
    gr_poly_init(a, ctx);
    gr_poly_init(b, ctx);
    gr_poly_init(t, ctx);

    status |= gr_poly_make_monic(f, pol, ctx);
    status |= gr_poly_preinv_set(P, f, ctx);
    if (d > 1)
        status |= gr_poly_preinv_powmod_x_fmpz(xq, q, P, ctx);

    if (status == GR_SUCCESS)
        status |= _gr_poly_edf_random_element(b, &found, state, f, P, xq, d, q, a, factor, t, ctx);
    if (status == GR_SUCCESS && found)
        status |= _gr_poly_edf_split_power(factor, &found, f, P, b, d, q, halfq, k, a, t, ctx);

    if (status == GR_SUCCESS && !found)
        status |= gr_poly_one(factor, ctx);

    fmpz_clear(q);
    fmpz_clear(halfq);
    gr_poly_clear(f, ctx);
    gr_poly_preinv_clear(P, ctx);
    gr_poly_clear(xq, ctx);
    gr_poly_clear(a, ctx);
    gr_poly_clear(b, ctx);
    gr_poly_clear(t, ctx);

    return status;
}

int
gr_poly_factor_equal_deg(gr_poly_vec_t fac, const gr_poly_t pol, slong d, gr_ctx_t ctx)
{
    gr_poly_t f, frob;
    gr_poly_preinv_t P;
    fmpz_t q;
    flint_rand_t state;
    slong n;
    int status = GR_SUCCESS;

    gr_poly_vec_set_length(fac, 0, ctx);

    n = pol->length - 1;

    if (n < 1 || d <= 0 || n % d != 0)
        return GR_DOMAIN;

    fmpz_init(q);
    gr_poly_init(f, ctx);
    gr_poly_preinv_init(P, ctx);
    gr_poly_init(frob, ctx);

    status |= _gr_poly_factor_ff_info(q, NULL, NULL, ctx);
    if (status != GR_SUCCESS)
        goto cleanup;

    status |= gr_poly_make_monic(f, pol, ctx);
    GR_POLY_FACTOR_CHECK(f)

    if (n > d && d > 1)
    {
        status |= gr_poly_preinv_set(P, f, ctx);
        status |= gr_poly_preinv_powmod_x_fmpz(frob, q, P, ctx);
        GR_POLY_FACTOR_CHECK_STATUS()
    }

    flint_rand_init(state);
    status |= _gr_poly_factor_equal_deg_with_frob(fac, f, d, frob, state, ctx);
    flint_rand_clear(state);

cleanup:
    fmpz_clear(q);
    gr_poly_clear(f, ctx);
    gr_poly_preinv_clear(P, ctx);
    gr_poly_clear(frob, ctx);

    return status;
}
