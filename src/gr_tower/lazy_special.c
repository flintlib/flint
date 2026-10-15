/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Special functions of the lazy tower field.

    A special function value F(z) becomes a transcendental generator of
    the tower only for a canonical argument z: the known functional
    equations (shifts, reflection, parity), special values and
    transformations are applied first, so that values related by them
    are expressed through the same generator (and the zero test, which
    treats distinct generators as algebraically independent, does not
    have to rediscover the relations). Every rule used here is an
    identity; numerical evaluation only chooses representatives (with
    exact tests where enclosures do not decide).

    The canonical arguments are intrinsic (they depend on the value of
    z, not on its representation), so that different towers create the
    same generators, which merges then identify:

    - Gamma, digamma and polygamma functions: the argument is shifted
      by an integer to the strip 0 <= Re(z) < 1, and reflected (z ->
      1 - z) to 0 < Re(z) < 1/2, or Re(z) in {0, 1/2} with Im(z) > 0.
      Rational arguments use the closed forms at integers and
      half-integers, Gauss's digamma theorem, and normal forms in terms
      of basis values (the gamma lattice; the Hurwitz zeta lattice for
      psi^(m), m >= 1, with zeta(odd) and Catalan's constant).
    - erf: odd, so the argument has Re(z) > 0; on the imaginary axis,
      erf(i y) = i erfi(y) with a (real) generator erfi(y), y > 0;
      erfc(z) = 1 - erf(z) and erfi(z) = -i erf(i z).
    - zeta at integers: Bernoulli numbers (rational multiples of powers
      of pi at even integers); odd integers >= 3 give generators.
    - polylogarithms: rational functions for s <= 0, logarithms for s = 1,
      zeta values at z = 1 and z = -1; at other roots of unity, sums of
      Hurwitz zeta values (in their normal form).
    - complete elliptic integrals: through 2F1 (lazy_hypgeom.c), which
      gives K(0), E(0), E(1), K(1/2), E(1/2) and the imaginary-modulus
      transformation to |m - 1| <= 1; the singular values for r = 2, 3, 4
      and their complements here (Legendre's relation and Landen's
      transformation are found by the zero test, elliptic_relations.c).

    Functions of real arguments are real where the function is: the
    enclosures of such generators have exactly zero imaginary parts, so
    the real views (GR_TOWER_LAZY_REAL) accept them.
*/

#include <math.h>
#include "fmpq.h"
#include "fmpq_mat.h"
#include "fmpq_vec.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "arith.h"
#include "ulong_extras.h"
#include "bernoulli.h"
#include "acb.h"
#include "qqbar.h"
#include "gr.h"
#include "gr_special.h"
#include "gr_poly.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"
#include "gr_tower/impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

/* recurrences (shifts of the argument) over at most this many steps;
   beyond, the generator of the unshifted argument is used */
#define SHIFT_LIMIT 1000
/* exact values (factorials, Bernoulli numbers) up to this size */
#define EXACT_LIMIT 100000
/* Gamma at positive integers: exact factorials up to this size
   (Gamma(10^6) has 5.6 million digits: 0.2 s) */
#define FACTORIAL_LIMIT 2000000

#define FLAGS(ctx) gr_tower_lazy_ctx_field_flags(ctx)
#define ALG(ctx) ((FLAGS(ctx) & GR_TOWER_LAZY_ALGEBRAIC) && _gr_tower_lazy_outermost(ctx))
#define REAL(ctx) ((FLAGS(ctx) & GR_TOWER_LAZY_REAL) && _gr_tower_lazy_outermost(ctx))

#define CHECK_PREC 64

/* whether |x| > c */
static int
_abs_gt_ui(const fmpz_t x, ulong c)
{
    return !fmpz_fits_si(x) || (ulong) FLINT_ABS(fmpz_get_si(x)) > c;
}

/* -------------------------------------------------------------------- */
/* helpers                                                               */
/* -------------------------------------------------------------------- */

static int
_is_rational(fmpq_t c, gr_srcptr x, gr_ctx_t ctx)
{
    return _gr_tower_lazy_rational_repr_locked(c, x, ctx);
}

static int
_pi_times(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr p;
    int status;
    GR_TMP_INIT(p, ctx);
    status = gr_pi(p, ctx);
    status |= gr_mul(res, p, x, ctx);
    GR_TMP_CLEAR(p, ctx);
    return status;
}

/* sin(pi x) or cos(pi x) (which = 0, 1); real algebraic numbers at
   rational x (rather than expressions in roots of unity), in the tangent
   normal form for denominators up to the option GR_TOWER_OPT_TRIG_PI_LIMIT */

static int
_trig_pi(gr_ptr res, gr_srcptr x, int which, gr_ctx_t ctx)
{
    fmpq_t c;
    int status;

    fmpq_init(c);
    if (_is_rational(c, x, ctx) && fmpz_fits_si(fmpq_numref(c)) && fmpz_abs_fits_ui(fmpq_denref(c)) &&
        fmpz_cmp_ui(fmpq_denref(c), gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_TRIG_PI_LIMIT)) <= 0)
    {
        /* the tangent normal form: the values at one level share the
           generator tan(pi/M) */
        status = _gr_tower_lazy_trig_pi_real_locked(res, c, which, ctx);
    }
    else if (_is_rational(c, x, ctx) && fmpz_fits_si(fmpq_numref(c)) && fmpz_abs_fits_ui(fmpq_denref(c)) &&
        fmpz_cmp_si(fmpq_denref(c), gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_TRIG_ALGEBRAIC_LIMIT)) <= 0)
    {
        qqbar_t q;
        gr_ctx_t QQbar;
        qqbar_init(q);
        gr_ctx_init_real_qqbar(QQbar);
        if (which == 0)
            qqbar_sin_pi(q, fmpz_get_si(fmpq_numref(c)), fmpz_get_ui(fmpq_denref(c)));
        else
            qqbar_cos_pi(q, fmpz_get_si(fmpq_numref(c)), fmpz_get_ui(fmpq_denref(c)));
        status = gr_set_other(res, q, QQbar, ctx);
        gr_ctx_clear(QQbar);
        qqbar_clear(q);
    }
    else
    {
        status = _pi_times(res, x, ctx);
        status |= (which == 0) ? gr_sin(res, res, ctx) : gr_cos(res, res, ctx);
    }
    fmpq_clear(c);
    return status;
}

static int
_trig_pi_fmpq(gr_ptr res, const fmpq_t c, int which, gr_ctx_t ctx)
{
    int status = gr_set_fmpq(res, c, ctx);
    if (status == GR_SUCCESS)
        status = _trig_pi(res, res, which, ctx);
    return status;
}

/* n = floor(Re(x)) */
static int
_re_floor(slong * n, gr_srcptr x, gr_ctx_t ctx)
{
    acb_t z;
    arb_t f;
    fmpz_t a, b;
    int status;

    acb_init(z);
    arb_init(f);
    fmpz_init(a);
    fmpz_init(b);

    status = gr_tower_lazy_get_acb(z, x, CHECK_PREC, ctx);
    if (status == GR_SUCCESS)
    {
        /* a = floor of the lower bound, b = floor of the upper bound */
        arb_floor(f, acb_realref(z), CHECK_PREC);
        if (arb_is_exact(f) && arb_is_finite(f))
        {
            arf_get_fmpz(a, arb_midref(f), ARF_RND_FLOOR);
            fmpz_set(b, a);
        }
        else
            fmpz_one(a);
    }

    if (status == GR_SUCCESS && fmpz_equal(a, b) &&
        !arb_contains_int(acb_realref(z)) && fmpz_fits_si(a))
    {
        *n = fmpz_get_si(a);
    }
    else if (status == GR_SUCCESS)
    {
        /* exact */
        gr_ptr t;
        GR_TMP_INIT(t, ctx);
        status = gr_floor(t, x, ctx);     /* (of the real part) */
        if (status == GR_SUCCESS)
            status = gr_get_si(n, t, ctx);
        GR_TMP_CLEAR(t, ctx);
    }

    acb_clear(z);
    arb_clear(f);
    fmpz_clear(a);
    fmpz_clear(b);
    return status;
}

/*
    The canonical representative of z under integer shifts and the
    reflection z -> 1 - z: with z = z0 + n, 0 <= Re(z0) < 1, sets
    *reflect = 1 if z0 is to be replaced by 1 - z0 (Re(z0) > 1/2, or
    Re(z0) in {0, 1/2} with Im(z0) < 0). Sets *pole = 1 if z is a
    nonpositive integer. GR_UNABLE if the shift exceeds SHIFT_LIMIT.
*/
static int
_shift_reflect(gr_ptr z0, slong * n, int * reflect, int * pole, gr_srcptr z, gr_ctx_t ctx)
{
    fmpq_t half;
    int status, c, s;

    *reflect = 0;
    *pole = 0;
    *n = 0;

    status = _re_floor(n, z, ctx);
    if (status != GR_SUCCESS)
        return status;
    if (*n > SHIFT_LIMIT || *n < -SHIFT_LIMIT)
        return GR_UNABLE;

    status = gr_sub_si(z0, z, *n, ctx);
    if (status != GR_SUCCESS)
        return status;

    fmpq_init(half);
    fmpq_set_si(half, 1, 2);

    status = _gr_tower_lazy_re_cmp(&c, z0, half, ctx);
    if (status == GR_SUCCESS)
    {
        if (c > 0)
            *reflect = 1;
        else
        {
            int line = (c == 0) ? 1 : 0;     /* 1: Re(z0) = 1/2, 2: Re(z0) = 0 */

            if (!line)
            {
                fmpq_zero(half);
                status = _gr_tower_lazy_re_cmp(&c, z0, half, ctx);
                if (status == GR_SUCCESS && c == 0)
                    line = 2;
            }

            if (status == GR_SUCCESS && line)
            {
                status = _gr_tower_lazy_im_sign(&s, z0, ctx);
                if (status == GR_SUCCESS)
                {
                    if (s < 0)
                        *reflect = 1;
                    else if (s == 0)
                    {
                        /* z is an integer or a half-integer which is not
                           represented as a rational number (an unrefined
                           tower): only the poles are decided */
                        if (line == 2 && *n <= 0)
                            *pole = 1;
                        else
                            status = GR_UNABLE;
                    }
                }
            }
        }
    }

    fmpq_clear(half);
    return status;
}

/* -------------------------------------------------------------------- */
/* Gamma                                                                 */
/* -------------------------------------------------------------------- */

/*
    The column order of the unknowns k = 1, ..., q - 1 of the lattices at
    level q (gamma and Hurwitz zeta values at k/q): the non-candidates
    for the basis (k/q > 1/2) first; then the candidates by decreasing
    (denominator, numerator), so that a reduced row echelon form prefers
    small denominators in the basis.
*/
static void
_lattice_order(slong * order, slong q)
{
    slong m = 0, nc = 0, a, b, k;
    slong * cand = order + (q - 1) / 2;

    for (k = 1; k < q; k++)
        if (2 * k > q)
            order[m++] = k;
    /* (the candidates, k <= q/2, sorted in place after them) */
    for (k = 1; 2 * k <= q; k++)
        cand[nc++] = k;
    for (a = 1; a < nc; a++)
    {
        slong x = cand[a], gx = n_gcd(x, q);
        slong dx = q / gx, nx = x / gx;
        for (b = a - 1; b >= 0; b--)
        {
            slong y = cand[b], gy = n_gcd(y, q);
            slong dy = q / gy, ny = y / gy;
            if (dy > dx || (dy == dx && ny > nx))
                break;
            cand[b + 1] = cand[b];
        }
        cand[b + 1] = x;
    }
}

/*
    Gamma at rational arguments: a normal form modulo all the relations
    which follow from the reflection formula and Gauss's multiplication
    theorem (the universal distribution relations; by the Rohrlich-Lang
    conjecture, these are all the algebraic relations between the values
    Gamma(a), a rational, and pi).

    At level q, the unknowns are x_k = log Gamma(k/q), 0 < k < q, with the
    relations

        x_k + x_{q-k} = log pi - log sin(pi k/q),
        sum_{j<n} x_{k + j q/n} = (n-1)/2 log(2 pi) + (1/2 - n k/q) log n + x_{n k}
            (n | q, 0 < k < q/n),

    whose right-hand sides are logarithms of pi, of integers and of the
    (positive, algebraic) numbers sin(pi k/q). A reduced row echelon form
    with the columns ordered by preference (the arguments k/q > 1/2
    first, then those in (0, 1/2] from the least to the most preferred,
    the most preferred having the smallest denominator, then the smallest
    numerator) leaves as free columns the canonical basis at level q:
    greedily, the most preferred values independent of those preferred
    to them (one value Gamma(1/3) at level 6, Gamma(1/3), Gamma(1/4) at
    level 12; phi(q)/2 values in all). Every other value is a monomial in
    the basis values times rational powers of pi, integers and sines.
    The basis at level q contains the basis at each level d | q, since
    the values of denominator d are preferred to the others, so values
    computed at different levels are expressed consistently.

    Sets *is_basis if p/q is itself a basis value; otherwise sets
    *num_free and the arrays: the basis values (as numerators over q) and
    their integer exponents, and the rational exponents of the constants
    pi, n (2 <= n <= q) and sin(pi k/q) (1 <= k < q/2). Returns 0 if a
    basis exponent is not an integer (not expected).
*/
static int
_gamma_lattice(int * is_basis, slong * num_free, slong * free_k, fmpz * free_exp,
    fmpq * exp_pi, fmpq * exp_int, fmpq * exp_sin, slong p, slong q)
{
    slong nvar = q - 1, nsin = (q - 1) / 2;
    slong cpi = GR_TOWER_GAMMA_CPI(q), cint = GR_TOWER_GAMMA_CINT(q), csin = GR_TOWER_GAMMA_CSIN(q), ncols = GR_TOWER_GAMMA_NCOLS(q);
    slong * order, * col_of_k, * k_of_col;
    slong nrows = 0, i, j, r, rank;
    fmpq_mat_t A, B;
    int ok = 1;

    /* column order of the unknowns */
    order = flint_malloc(sizeof(slong) * nvar);
    col_of_k = flint_malloc(sizeof(slong) * q);
    k_of_col = flint_malloc(sizeof(slong) * nvar);
    _lattice_order(order, q);
    for (i = 0; i < nvar; i++)
    {
        col_of_k[order[i]] = i;
        k_of_col[i] = order[i];
    }

    fmpq_mat_init(A, _gr_tower_gamma_relations_rows(q), ncols);
    nrows = _gr_tower_gamma_relations(A, col_of_k, q);

    {
        fmpq_mat_t W;
        fmpq_mat_window_init(W, A, 0, 0, nrows, ncols);
        fmpq_mat_init(B, nrows, ncols);
        rank = fmpq_mat_rref(B, W);
        fmpq_mat_window_clear(W);
    }

    /* the target column */
    {
        slong ct = col_of_k[p], row = -1;
        for (r = 0; r < rank && row < 0; r++)
        {
            /* the pivot of row r is its first nonzero entry */
            for (j = 0; j < ncols; j++)
                if (!fmpq_is_zero(fmpq_mat_entry(B, r, j)))
                    break;
            if (j == ct)
                row = r;
        }

        if (row < 0)
        {
            *is_basis = 1;
        }
        else
        {
            /* x_t = -sum_{free} B[row, f] x_f - sum_c B[row, c] C_c */
            *is_basis = 0;
            *num_free = 0;
            for (j = 0; j < nvar && ok; j++)
            {
                if (j != ct && !fmpq_is_zero(fmpq_mat_entry(B, row, j)))
                {
                    if (!fmpz_is_one(fmpq_denref(fmpq_mat_entry(B, row, j))))
                        ok = 0;
                    else
                    {
                        free_k[*num_free] = k_of_col[j];
                        fmpz_neg(free_exp + *num_free, fmpq_numref(fmpq_mat_entry(B, row, j)));
                        (*num_free)++;
                    }
                }
            }
            fmpq_neg(exp_pi, fmpq_mat_entry(B, row, cpi));
            for (j = 0; j < q - 1; j++)
                fmpq_neg(exp_int + j, fmpq_mat_entry(B, row, cint + j));
            for (j = 0; j < nsin; j++)
                fmpq_neg(exp_sin + j, fmpq_mat_entry(B, row, csin + j));
        }
    }

    fmpq_mat_clear(A);
    fmpq_mat_clear(B);
    flint_free(order);
    flint_free(col_of_k);
    flint_free(k_of_col);
    return ok;
}

/* res *= x^e for a positive real x and a rational e (a principal root) */
static int
_mul_pow_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t e, gr_ctx_t ctx)
{
    gr_ptr t;
    int status = GR_SUCCESS;

    if (fmpq_is_zero(e))
        return GR_SUCCESS;
    if (!fmpz_abs_fits_ui(fmpq_denref(e)) || !fmpz_fits_si(fmpq_numref(e)))
        return GR_UNABLE;

    GR_TMP_INIT(t, ctx);
    if (fmpz_is_one(fmpq_denref(e)))
        status = gr_set(t, x, ctx);
    else
        status = gr_tower_lazy_root_ui(t, x, fmpz_get_ui(fmpq_denref(e)), ctx);
    if (status == GR_SUCCESS)
        status = gr_pow_si(t, t, fmpz_get_si(fmpq_numref(e)), ctx);
    status |= gr_mul(res, res, t, ctx);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

/* res = sin(pi k/q) = (z^k - z^(-k)) / (2i), z = exp(pi i/q), as an
   element of a cyclotomic field (a square root for q = 3, 4, 6); t is
   scratch space */
static int
_sin_pi_cyclo(gr_ptr res, gr_ptr t, slong k, slong q, gr_ctx_t ctx)
{
    qqbar_t z;
    gr_ctx_t QQbar;
    int status;

    if (q == 3 || q == 4 || q == 6)
    {
        /* quadratic irrationals (sqrt(3)/2, sqrt(2)/2, 1/2): written as
           square roots */
        fmpq_t c;
        fmpq_init(c);
        qqbar_init(z);
        qqbar_sin_pi(z, k, q);
        qqbar_sqr(z, z);
        qqbar_get_fmpq(c, z);
        status = gr_set_fmpq(res, c, ctx);
        status |= gr_sqrt(res, res, ctx);
        qqbar_clear(z);
        fmpq_clear(c);
        return status;
    }

    qqbar_init(z);
    gr_ctx_init_complex_qqbar(QQbar);
    qqbar_root_of_unity(z, k, 2 * q);
    status = gr_set_other(res, z, QQbar, ctx);
    status |= gr_inv(t, res, ctx);
    status |= gr_sub(res, res, t, ctx);
    qqbar_root_of_unity(z, 1, 4);
    status |= gr_set_other(t, z, QQbar, ctx);
    status |= gr_mul_ui(t, t, 2, ctx);
    status |= gr_div(res, res, t, ctx);
    gr_ctx_clear(QQbar);
    qqbar_clear(z);
    return status;
}

/* Gamma(p/q) for 0 < p < q, gcd(p, q) = 1, 3 <= q <= the option GR_TOWER_OPT_GAMMA_LATTICE_LIMIT:
   GR_UNABLE if the normal form is not available */
static int
_gamma_fmpq_lattice(gr_ptr res, slong p, slong q, gr_ctx_t ctx)
{
    int is_basis = 0, ok, status = GR_SUCCESS;
    slong num_free = 0, i;
    slong * free_k = flint_malloc(sizeof(slong) * q);
    fmpz * free_exp = _fmpz_vec_init(q);
    fmpq_t exp_pi;
    fmpq * exp_int = _fmpq_vec_init(q);
    fmpq * exp_sin = _fmpq_vec_init(q);

    fmpq_init(exp_pi);

    ok = _gamma_lattice(&is_basis, &num_free, free_k, free_exp, exp_pi, exp_int, exp_sin, p, q);

    if (!ok)
        status = GR_UNABLE;
    else if (is_basis)
    {
        fmpq_t c;
        fmpq_init(c);
        fmpq_set_si(c, p, q);
        status = gr_set_fmpq(res, c, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_special_gen_locked(res, res, GR_TOWER_GAMMA, 0, ctx);
        fmpq_clear(c);
    }
    else
    {
        gr_ptr t;
        fmpq_t c, e;
        GR_TMP_INIT(t, ctx);
        fmpq_init(c);
        fmpq_init(e);

        status = gr_one(res, ctx);

        for (i = 0; i < num_free && status == GR_SUCCESS; i++)
        {
            fmpq_set_si(c, free_k[i], q);
            status = gr_set_fmpq(t, c, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_special_gen_locked(t, t, GR_TOWER_GAMMA, 0, ctx);
            if (status == GR_SUCCESS)
                status = gr_pow_fmpz(t, t, free_exp + i, ctx);
            status |= gr_mul(res, res, t, ctx);
        }

        if (status == GR_SUCCESS && !fmpq_is_zero(exp_pi))
        {
            status = gr_pi(t, ctx);
            status |= _mul_pow_fmpq(res, t, exp_pi, ctx);
        }

        /* the algebraic factor prod n^(e_n) prod sin(pi k/q)^(e_k): for
           each denominator d of the exponents, one d-th root of a rational
           number (a structured radical) and one d-th root of a product of
           sines with integer exponents */
        if (status == GR_SUCCESS)
        {
            fmpz_t d, num, den;
            fmpq_t r, f;
            gr_ptr B;
            slong m, nd = 0;
            fmpz * dens = _fmpz_vec_init(2 * q);

            fmpz_init(d);
            fmpz_init(num);
            fmpz_init(den);
            fmpq_init(r);
            fmpq_init(f);
            GR_TMP_INIT(B, ctx);

            /* the distinct denominators */
            for (i = 0; i < q - 1 + (q - 1) / 2; i++)
            {
                const fmpq * x = (i < q - 1) ? exp_int + i : exp_sin + (i - (q - 1));
                if (fmpq_is_zero(x))
                    continue;
                for (m = 0; m < nd; m++)
                    if (fmpz_equal(dens + m, fmpq_denref(x)))
                        break;
                if (m == nd)
                    fmpz_set(dens + nd++, fmpq_denref(x));
            }

            for (m = 0; m < nd && status == GR_SUCCESS; m++)
            {
                ulong dd;
                fmpz_set(d, dens + m);
                if (!fmpz_abs_fits_ui(d))
                {
                    status = GR_UNABLE;
                    break;
                }
                dd = fmpz_get_ui(d);

                /* rational part: r = prod n^(e_n d) */
                fmpq_one(r);
                for (i = 0; i < q - 1; i++)
                {
                    if (!fmpq_is_zero(exp_int + i) && fmpz_equal(fmpq_denref(exp_int + i), d))
                    {
                        fmpq_set_si(f, i + 2, 1);
                        fmpq_pow_si(f, f, fmpz_get_si(fmpq_numref(exp_int + i)));
                        fmpq_mul(r, r, f);
                    }
                }
                if (!fmpq_is_one(r))
                {
                    status = gr_set_fmpq(B, r, ctx);
                    if (dd > 1 && status == GR_SUCCESS)
                        status = gr_tower_lazy_root_ui(B, B, dd, ctx);
                    status |= gr_mul(res, res, B, ctx);
                }

                /* sines: their product (with integer exponents), with
                   sin(pi k/q) = (z^k - z^(-k)) / (2i) in terms of the root of
                   unity z = exp(pi i/q), so that all values at one level
                   (and at levels dividing it) live in one cyclotomic field
                   (rather than in separate algebraic generators, whose
                   compositum the tower would have to discover) */
                {
                    int any = 0;
                    gr_ptr P, S, Z;
                    GR_TMP_INIT3(P, S, Z, ctx);
                    status |= gr_one(P, ctx);
                    for (i = 0; i < (q - 1) / 2 && status == GR_SUCCESS; i++)
                    {
                        if (!fmpq_is_zero(exp_sin + i) && fmpz_equal(fmpq_denref(exp_sin + i), d))
                        {
                            status = _sin_pi_cyclo(S, Z, i + 1, q, ctx);
                            status |= gr_pow_fmpz(S, S, fmpq_numref(exp_sin + i), ctx);
                            status |= gr_mul(P, P, S, ctx);
                            any = 1;
                        }
                    }
                    if (any && status == GR_SUCCESS && gr_is_one(P, ctx) != T_TRUE)
                    {
                        if (dd > 1)
                            status = gr_tower_lazy_root_ui(P, P, dd, ctx);
                        status |= gr_mul(res, res, P, ctx);
                    }
                    GR_TMP_CLEAR3(P, S, Z, ctx);
                }
            }

            _fmpz_vec_clear(dens, 2 * q);
            fmpz_clear(d);
            fmpz_clear(num);
            fmpz_clear(den);
            fmpq_clear(r);
            fmpq_clear(f);
            GR_TMP_CLEAR(B, ctx);
        }

        GR_TMP_CLEAR(t, ctx);
        fmpq_clear(c);
        fmpq_clear(e);
    }

    flint_free(free_k);
    _fmpz_vec_clear(free_exp, q);
    fmpq_clear(exp_pi);
    _fmpq_vec_clear(exp_int, q);
    _fmpq_vec_clear(exp_sin, q);
    return status;
}

static int
_gamma_fmpq(gr_ptr res, const fmpq_t c, gr_ctx_t ctx)
{
    fmpz_t n;
    fmpq_t c0, r, t;
    slong k, nn;
    int status = GR_SUCCESS;

    if (fmpz_is_one(fmpq_denref(c)))
    {
        if (fmpz_sgn(fmpq_numref(c)) <= 0)
            return GR_DOMAIN;
        if (fmpz_cmp_ui(fmpq_numref(c), FACTORIAL_LIMIT) > 0)
            return GR_UNABLE;
        fmpz_init(n);
        fmpz_fac_ui(n, fmpz_get_ui(fmpq_numref(c)) - 1);
        status = gr_set_fmpz(res, n, ctx);
        fmpz_clear(n);
        return status;
    }

    fmpz_init(n);
    fmpq_init(c0);
    fmpq_init(r);
    fmpq_init(t);

    /* c = c0 + n, 0 < c0 < 1 */
    fmpz_fdiv_q(n, fmpq_numref(c), fmpq_denref(c));
    fmpq_sub_fmpz(c0, c, n);

    if (_abs_gt_ui(n, EXACT_LIMIT))
    {
        status = GR_UNABLE;
        goto cleanup;
    }
    nn = fmpz_get_si(n);

    /* r = Gamma(c) / Gamma(c0), a rational number */
    fmpq_one(r);
    if (nn >= 0)
    {
        for (k = 0; k < nn; k++)
        {
            fmpq_add_si(t, c0, k);
            fmpq_mul(r, r, t);
        }
    }
    else
    {
        for (k = 1; k <= -nn; k++)
        {
            fmpq_sub_si(t, c0, k);
            fmpq_div(r, r, t);
        }
    }

    if (fmpz_equal_ui(fmpq_denref(c0), 2))
    {
        /* Gamma(1/2) = sqrt(pi) */
        status = gr_pi(res, ctx);
        status |= gr_sqrt(res, res, ctx);
    }
    else if (fmpz_cmp_ui(fmpq_denref(c0), gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_GAMMA_LATTICE_LIMIT)) <= 0 &&
        _gamma_fmpq_lattice(res, fmpz_get_si(fmpq_numref(c0)), fmpz_get_si(fmpq_denref(c0)), ctx) == GR_SUCCESS)
    {
        /* the normal form; otherwise (and beyond the limit) a generator,
           possibly dependent on others: the relation search of the
           zero test finds the distribution relations between gamma
           generators at rational arguments and eliminates them
           (certify.c) */
    }
    else
    {
        fmpq_t half;
        fmpq_init(half);
        fmpq_set_si(half, 1, 2);
        if (fmpq_cmp(c0, half) > 0)
        {
            /* Gamma(c0) = pi / (sin(pi c0) Gamma(1 - c0)) */
            gr_ptr s, g;
            GR_TMP_INIT2(s, g, ctx);
            fmpq_sub_si(t, c0, 1);
            fmpq_neg(t, t);
            /* (sin(pi c0) in terms of a root of unity, like the sines of
               the normal form: values at one level then share their
               cyclotomic field) */
            if (fmpz_fits_si(fmpq_numref(c0)) && fmpz_cmp_si(fmpq_denref(c0), gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_TRIG_ALGEBRAIC_LIMIT)) <= 0)
                status = _sin_pi_cyclo(s, g, fmpz_get_si(fmpq_numref(c0)), fmpz_get_si(fmpq_denref(c0)), ctx);
            else
                status = _trig_pi_fmpq(s, c0, 0, ctx);
            status |= gr_set_fmpq(g, t, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_special_gen_locked(g, g, GR_TOWER_GAMMA, 0, ctx);
            status |= gr_mul(s, s, g, ctx);
            status |= gr_pi(res, ctx);
            if (status == GR_SUCCESS)
                status = gr_div(res, res, s, ctx);
            GR_TMP_CLEAR2(s, g, ctx);
        }
        else
        {
            status = gr_set_fmpq(res, c0, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_special_gen_locked(res, res, GR_TOWER_GAMMA, 0, ctx);
        }
        fmpq_clear(half);
    }

    status |= gr_mul_fmpq(res, res, r, ctx);

cleanup:
    fmpz_clear(n);
    fmpq_clear(c0);
    fmpq_clear(r);
    fmpq_clear(t);
    return status;
}

static int _gamma(gr_ptr res, gr_srcptr z, gr_ctx_t ctx);

static int
_gamma(gr_ptr res, gr_srcptr z_in, gr_ctx_t ctx)
{
    fmpq_t c;
    gr_ptr z, z0, t, u;
    slong n;
    int reflect, pole, status;

    fmpq_init(c);
    if (_is_rational(c, z_in, ctx))
    {
        status = _gamma_fmpq(res, c, ctx);
        fmpq_clear(c);
        return status;
    }
    fmpq_clear(c);

    GR_TMP_INIT4(z, z0, t, u, ctx);
    status = gr_set(z, z_in, ctx);     /* (res may alias the argument) */
    if (status != GR_SUCCESS)
        goto cleanup;

    status = _shift_reflect(z0, &n, &reflect, &pole, z, ctx);

    if (status == GR_UNABLE && (n > SHIFT_LIMIT || n < -SHIFT_LIMIT))
    {
        /* (a far argument: its own generator) */
        status = _gr_tower_lazy_special_gen_locked(res, z, GR_TOWER_GAMMA, 0, ctx);
        goto cleanup;
    }

    if (status == GR_SUCCESS && pole)
        status = GR_DOMAIN;

    if (status != GR_SUCCESS)
        goto cleanup;

    if (reflect)
    {
        /* Gamma(z0) = pi / (sin(pi z0) Gamma(1 - z0)), with Gamma(1 - z0)
           canonical, or (when Re(z0) = 0) Gamma(1 - z0) = w Gamma(w)
           for the canonical w = -z0 */
        int c0;
        fmpq_t zero;
        fmpq_init(zero);
        status = _gr_tower_lazy_re_cmp(&c0, z0, zero, ctx);
        fmpq_clear(zero);

        status |= gr_neg(t, z0, ctx);
        if (status == GR_SUCCESS && c0 == 0)
        {
            status = _gr_tower_lazy_special_gen_locked(u, t, GR_TOWER_GAMMA, 0, ctx);
            status |= gr_mul(u, u, t, ctx);
        }
        else if (status == GR_SUCCESS)
        {
            status = gr_add_si(t, t, 1, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_special_gen_locked(u, t, GR_TOWER_GAMMA, 0, ctx);
        }

        status |= _trig_pi(t, z0, 0, ctx);
        status |= gr_mul(u, u, t, ctx);
        status |= gr_pi(t, ctx);
        if (status == GR_SUCCESS)
            status = gr_div(res, t, u, ctx);
    }
    else
    {
        status = _gr_tower_lazy_special_gen_locked(res, z0, GR_TOWER_GAMMA, 0, ctx);
    }

    /* Gamma(z0 + n) = Gamma(z0) (z0)_n, Gamma(z0 - n) = Gamma(z0) / (z0 - n)_n */
    if (status == GR_SUCCESS && n > 0)
    {
        status = gr_rising_ui(t, z0, n, ctx);
        status |= gr_mul(res, res, t, ctx);
    }
    else if (status == GR_SUCCESS && n < 0)
    {
        status = gr_rising_ui(t, z, -n, ctx);
        if (status == GR_SUCCESS)
            status = gr_div(res, res, t, ctx);
    }

cleanup:
    GR_TMP_CLEAR4(z, z0, t, u, ctx);
    return status;
}

int
gr_tower_lazy_gamma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _gr_tower_lazy_view_finish(_gamma(res, x, ctx), res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_rgamma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _gamma(res, x, ctx);
    if (status == GR_DOMAIN)
        status = gr_zero(res, ctx);    /* at the poles */
    else if (status == GR_SUCCESS)
        status = gr_inv(res, res, ctx);
    status = _gr_tower_lazy_view_finish(status, res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/*
    log Gamma(x) (the analytic continuation from the positive reals, with
    the branch cut on the negative reals, as acb_lgamma) = log(Gamma(x))
    + 2 pi i k, the integer k from enclosures
*/
static int
_lgamma(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr g, t;
    acb_t a, b, c;
    slong prec;
    int status, found = 0;
    fmpz_t k;

    GR_TMP_INIT2(g, t, ctx);
    acb_init(a);
    acb_init(b);
    acb_init(c);
    fmpz_init(k);

    status = _gamma(g, x, ctx);
    if (status == GR_SUCCESS)
        status = gr_log(g, g, ctx);

    for (prec = 64; status == GR_SUCCESS && !found && prec <= 4096; prec *= 2)
    {
        status = gr_tower_lazy_get_acb(a, x, prec, ctx);
        status |= gr_tower_lazy_get_acb(b, g, prec, ctx);
        if (status != GR_SUCCESS)
            break;
        acb_lgamma(a, a, prec);
        acb_sub(a, a, b, prec);
        acb_const_pi(c, prec);
        acb_mul_2exp_si(c, c, 1);
        acb_div_onei(a, a);
        acb_div(a, a, c, prec);
        /* (a = k, an integer, with a real part enclosure free of other integers) */
        if (arb_contains_zero(acb_imagref(a)) && mag_cmp_2exp_si(arb_radref(acb_realref(a)), -2) < 0)
        {
            arf_get_fmpz(k, arb_midref(acb_realref(a)), ARF_RND_NEAR);
            found = arb_contains_fmpz(acb_realref(a), k);
            if (!found)
                status = GR_UNABLE;
        }
    }
    if (status == GR_SUCCESS && !found)
        status = GR_UNABLE;

    if (status == GR_SUCCESS && !fmpz_is_zero(k))
    {
        /* + 2 pi i k */
        status = gr_pi(t, ctx);
        status |= gr_mul_fmpz(t, t, k, ctx);
        status |= gr_mul_2exp_si(t, t, 1, ctx);
        {
            gr_ptr u;
            GR_TMP_INIT(u, ctx);
            status |= gr_i(u, ctx);
            status |= gr_mul(t, t, u, ctx);
            GR_TMP_CLEAR(u, ctx);
        }
        status |= gr_add(g, g, t, ctx);
    }
    if (status == GR_SUCCESS)
        status = gr_set(res, g, ctx);

    GR_TMP_CLEAR2(g, t, ctx);
    acb_clear(a);
    acb_clear(b);
    acb_clear(c);
    fmpz_clear(k);
    return status;
}

int
gr_tower_lazy_lgamma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _gr_tower_lazy_view_finish(_lgamma(res, x, ctx), res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_beta(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
{
    gr_ptr a, b, c;
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    GR_TMP_INIT3(a, b, c, ctx);
    status = _gamma(a, x, ctx);
    if (status == GR_SUCCESS)
        status = _gamma(b, y, ctx);
    if (status == GR_SUCCESS)
        status = gr_add(c, x, y, ctx);
    if (status == GR_SUCCESS)
    {
        status = _gamma(c, c, ctx);
        if (status == GR_DOMAIN)
            status = gr_zero(res, ctx);     /* B(x, y) = 0 when x + y is a pole */
        else if (status == GR_SUCCESS)
        {
            status = gr_mul(a, a, b, ctx);
            if (status == GR_SUCCESS)
                status = gr_div(res, a, c, ctx);
        }
    }
    GR_TMP_CLEAR3(a, b, c, ctx);
    status = _gr_tower_lazy_view_finish(status, res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* error functions                                                       */
/* -------------------------------------------------------------------- */

static int
_erf(gr_ptr res, gr_srcptr z, gr_ctx_t ctx)
{
    fmpq_t zero;
    truth_t t;
    int c, status;

    t = gr_is_zero(z, ctx);
    if (t == T_TRUE)
        return gr_zero(res, ctx);
    if (t == T_UNKNOWN)
        return GR_UNABLE;

    /* odd: a canonical argument has Re(z) > 0, or Re(z) = 0 and Im(z) > 0 */
    fmpq_init(zero);
    status = _gr_tower_lazy_re_cmp(&c, z, zero, ctx);
    fmpq_clear(zero);
    if (status == GR_SUCCESS && c == 0)
        status = _gr_tower_lazy_im_sign(&c, z, ctx);
    if (status != GR_SUCCESS)
        return status;

    {
        gr_ptr w, i;
        int imag;
        fmpq_t zero;

        GR_TMP_INIT2(w, i, ctx);
        status = (c < 0) ? gr_neg(w, z, ctx) : gr_set(w, z, ctx);

        /* w = i y with y real (Re(w) = 0): erf(w) = i erfi(y), with
           the real generator erfi(y) */
        fmpq_init(zero);
        imag = 0;
        if (status == GR_SUCCESS)
        {
            int r;
            status = _gr_tower_lazy_re_cmp(&r, w, zero, ctx);
            imag = (status == GR_SUCCESS && r == 0);
        }
        fmpq_clear(zero);

        if (status == GR_SUCCESS && imag)
        {
            status = gr_i(i, ctx);
            status |= gr_mul(w, w, i, ctx);
            status |= gr_neg(w, w, ctx);      /* y = -i w > 0 */
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_special_gen_locked(res, w, GR_TOWER_ERF, 1, ctx);
            status |= gr_mul(res, res, i, ctx);
        }
        else if (status == GR_SUCCESS)
            status = _gr_tower_lazy_special_gen_locked(res, w, GR_TOWER_ERF, 0, ctx);

        if (c < 0)
            status |= gr_neg(res, res, ctx);
        GR_TMP_CLEAR2(w, i, ctx);
        return status;
    }
}

int
gr_tower_lazy_erf(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _gr_tower_lazy_view_finish(_erf(res, x, ctx), res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_erfc(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _erf(res, x, ctx);
    if (status == GR_SUCCESS)
    {
        status = gr_neg(res, res, ctx);
        status |= gr_add_si(res, res, 1, ctx);
    }
    status = _gr_tower_lazy_view_finish(status, res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_erfi(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_ptr i, t;
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    GR_TMP_INIT2(i, t, ctx);
    /* erfi(x) = -i erf(i x) */
    status = gr_i(i, ctx);
    status |= gr_mul(t, x, i, ctx);
    if (status == GR_SUCCESS)
        status = _erf(t, t, ctx);
    status |= gr_mul(res, t, i, ctx);
    status |= gr_neg(res, res, ctx);
    GR_TMP_CLEAR2(i, t, ctx);
    status = _gr_tower_lazy_view_finish(status, res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* Lambert W                                                             */
/* -------------------------------------------------------------------- */

static int
_lambertw(gr_ptr res, gr_srcptr z, slong k, gr_ctx_t ctx)
{
    truth_t t;
    acb_t a, b;
    int status = GR_SUCCESS, done = 0;

    t = gr_is_zero(z, ctx);
    if (t == T_TRUE)
        return (k == 0) ? gr_zero(res, ctx) : GR_DOMAIN;
    if (t == T_UNKNOWN)
        return GR_UNABLE;

    /* the branch point W_0(-1/e) = W_{-1}(-1/e) = -1 */
    if (k == 0 || k == -1)
    {
        acb_init(a);
        acb_init(b);
        if (gr_tower_lazy_get_acb(a, z, CHECK_PREC, ctx) == GR_SUCCESS)
        {
            arb_const_e(acb_realref(b), CHECK_PREC);
            arb_inv(acb_realref(b), acb_realref(b), CHECK_PREC);
            arb_neg(acb_realref(b), acb_realref(b));
            if (acb_overlaps(a, b))
            {
                gr_ptr e;
                GR_TMP_INIT(e, ctx);
                status = gr_set_si(e, -1, ctx);
                status |= gr_exp(e, e, ctx);
                status |= gr_add(e, e, z, ctx);
                if (status == GR_SUCCESS)
                {
                    t = gr_is_zero(e, ctx);
                    if (t == T_TRUE)
                    {
                        status = gr_set_si(res, -1, ctx);
                        done = 1;
                    }
                    else if (t == T_UNKNOWN)
                        status = GR_UNABLE;
                }
                GR_TMP_CLEAR(e, ctx);
            }
        }
        acb_clear(a);
        acb_clear(b);
    }

    if (status == GR_SUCCESS && !done)
        status = _gr_tower_lazy_special_gen_locked(res, z, GR_TOWER_LAMBERTW, k, ctx);

    return status;
}

int
gr_tower_lazy_lambertw_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t k, gr_ctx_t ctx)
{
    int status, real, alg;
    if (!fmpz_fits_si(k))
        return GR_UNABLE;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _gr_tower_lazy_view_finish(_lambertw(res, x, fmpz_get_si(k), ctx), res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_lambertw(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _gr_tower_lazy_view_finish(_lambertw(res, x, 0, ctx), res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* zeta at integers                                                      */
/* -------------------------------------------------------------------- */

static int
_zeta_si(gr_ptr res, slong s, gr_ctx_t ctx)
{
    fmpq_t b;
    int status;

    if (s == 1)
        return GR_DOMAIN;

    if (s == 0)
        return gr_set_si(res, -1, ctx) | gr_div_ui(res, res, 2, ctx);

    if (s >= 3 && (s % 2 == 1))
    {
        gr_ptr t;
        GR_TMP_INIT(t, ctx);
        status = gr_set_si(t, s, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_special_gen_locked(res, t, GR_TOWER_ZETA, 0, ctx);
        GR_TMP_CLEAR(t, ctx);
        return status;
    }

    if (s > EXACT_LIMIT || s < -EXACT_LIMIT)
        return GR_UNABLE;

    fmpq_init(b);

    if (s < 0)
    {
        /* zeta(-n) = -B_{n+1} / (n + 1) */
        bernoulli_fmpq_ui(b, 1 - s);
        {
            fmpz_t d;
            fmpz_init_set_si(d, 1 - s);
            fmpq_div_fmpz(b, b, d);
            fmpz_clear(d);
        }
        fmpq_neg(b, b);
        status = gr_set_fmpq(res, b, ctx);
    }
    else
    {
        /* zeta(2k) = (-1)^(k+1) B_{2k} (2 pi)^(2k) / (2 (2k)!) */
        gr_ptr p;
        GR_TMP_INIT(p, ctx);
        _gr_tower_zeta_even_over_pi(b, s);
        status = gr_pi(p, ctx);
        status |= gr_pow_ui(p, p, s, ctx);
        status |= gr_mul_fmpq(res, p, b, ctx);
        GR_TMP_CLEAR(p, ctx);
    }

    fmpq_clear(b);
    return status;
}

static int
_zeta(gr_ptr res, gr_srcptr s, gr_ctx_t ctx)
{
    fmpq_t c;
    int status;

    fmpq_init(c);
    if (_is_rational(c, s, ctx) && fmpz_is_one(fmpq_denref(c)))
    {
        if (fmpz_fits_si(fmpq_numref(c)))
            status = _zeta_si(res, fmpz_get_si(fmpq_numref(c)), ctx);
        else
            status = GR_UNABLE;
    }
    else
        status = _gr_tower_lazy_dirichlet_l_prim(res, s, 1, 1, ctx);   /* (the functional equation) */
    fmpq_clear(c);
    return status;
}

int
gr_tower_lazy_zeta(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _gr_tower_lazy_view_finish(_zeta(res, x, ctx), res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* digamma and polygamma functions                                       */
/* -------------------------------------------------------------------- */

/* res = cot(pi x) */
static int
_cot_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr s, c;
    int status;
    GR_TMP_INIT2(s, c, ctx);
    status = _trig_pi(c, x, 1, ctx);
    status |= _trig_pi(s, x, 0, ctx);
    if (status == GR_SUCCESS)
        status = gr_div(res, c, s, ctx);
    GR_TMP_CLEAR2(s, c, ctx);
    return status;
}

/* res = (d/dz)^m cot(pi z) = P_m(cot(pi z)), with P_0(c) = c and
   P_{j+1}(c) = -pi (1 + c^2) P_j'(c) */
static int
_cot_pi_derivative(gr_ptr res, gr_srcptr z, ulong m, gr_ctx_t ctx)
{
    gr_poly_t P, Q, R;
    gr_ptr c, p;
    ulong j;
    int status;

    gr_poly_init(P, ctx);
    gr_poly_init(Q, ctx);
    gr_poly_init(R, ctx);
    GR_TMP_INIT2(c, p, ctx);

    status = _cot_pi(c, z, ctx);
    status |= gr_pi(p, ctx);
    status |= gr_neg(p, p, ctx);

    status |= gr_poly_set_coeff_si(R, 0, 1, ctx);
    status |= gr_poly_set_coeff_si(R, 2, 1, ctx);   /* 1 + c^2 */
    status |= gr_poly_gen(P, ctx);
    for (j = 0; j < m && status == GR_SUCCESS; j++)
    {
        status |= gr_poly_derivative(Q, P, ctx);
        status |= gr_poly_mul(P, Q, R, ctx);
        status |= gr_poly_mul_scalar(P, P, p, ctx);
    }

    if (status == GR_SUCCESS)
        status = gr_poly_evaluate(res, P, c, ctx);

    gr_poly_clear(P, ctx);
    gr_poly_clear(Q, ctx);
    gr_poly_clear(R, ctx);
    GR_TMP_CLEAR2(c, p, ctx);
    return status;
}

/* Gauss's digamma theorem: psi(p/q) for 0 < p < q,
   -gamma - log(2q) - (pi/2) cot(pi p/q) + 2 sum_{k=1}^{floor((q-1)/2)} cos(2 pi k p/q) log sin(pi k/q) */
static int
_digamma_gauss(gr_ptr res, ulong p, ulong q, gr_ctx_t ctx)
{
    gr_ptr t, u, v;
    ulong k;
    int status;

    GR_TMP_INIT3(t, u, v, ctx);

    status = _gr_tower_lazy_special_gen_locked(res, NULL, GR_TOWER_CONSTANT, GR_TOWER_CONST_EULER, ctx);
    status |= gr_neg(res, res, ctx);

    status |= gr_set_ui(t, 2 * q, ctx);
    status |= gr_log(t, t, ctx);
    status |= gr_sub(res, res, t, ctx);

    if (2 * p != q)
    {
        fmpq_t r;
        fmpq_init(r);
        fmpq_set_ui(r, p, q);
        status |= gr_set_fmpq(t, r, ctx);
        status |= _cot_pi(t, t, ctx);
        status |= _pi_times(t, t, ctx);
        status |= gr_div_ui(t, t, 2, ctx);
        status |= gr_sub(res, res, t, ctx);
        fmpq_clear(r);
    }

    for (k = 1; k <= (q - 1) / 2 && status == GR_SUCCESS; k++)
    {
        fmpq_t r;
        fmpq_init(r);
        /* cos(2 pi k p / q) */
        fmpq_set_ui(r, 2 * k * p, q);
        status |= _trig_pi_fmpq(u, r, 1, ctx);
        /* log sin(pi k / q) */
        fmpq_set_ui(r, k, q);
        status |= _trig_pi_fmpq(v, r, 0, ctx);
        if (status == GR_SUCCESS)
            status = gr_log(v, v, ctx);
        status |= gr_mul(u, u, v, ctx);
        status |= gr_mul_ui(u, u, 2, ctx);
        status |= gr_add(res, res, u, ctx);
        fmpq_clear(r);
    }

    GR_TMP_CLEAR3(t, u, v, ctx);
    return status;
}

/*
    The normal form of zeta(s, p/q) = (-1)^s psi^(s-1)(p/q) / (s-1)!
    (s >= 2) at rational arguments, the additive analogue of the gamma
    lattice: by the reflection and distribution (multiplication)
    relations (special.c), which by a conjecture of the Chowla-Milnor
    type (the relations of Kubert's universal odd distribution) are all
    the linear relations between these values over the algebraic numbers
    and pi^s, the values at level q span a space with a basis of values
    zeta(s, k/q), k <= q/2, of the smallest denominators (and zeta(s) for
    odd s); every other value is a rational combination of the basis
    values plus an elementary part pi^s P(cot(pi k/q)). Catalan's
    constant replaces psi'(1/4) = pi^2 + 8G in the basis.

    Examples: psi'(3/4) = pi^2 - 8G, psi'(1/6) = 5 psi'(1/3) - 4 pi^2 ...,
    zeta(3, 1/2) = 7 zeta(3).
*/
/* beyond this level, the cotangents of the elementary parts are summed in
   a cyclotomic field (below, square roots appear for q = 3, 4, 6, and
   the small cyclotomic fields are cheap) */
#define HURWITZ_CYCLO_COT_MIN 12

/* res = cot(pi k/q) for 1 <= k <= q/2, in a cyclotomic field: i (z + 1) /
   (z - 1) with z = exp(2 pi i k/q) (square roots for q = 3, 4, 6) */
static int
_cot_pi_cyclo(gr_ptr res, gr_ptr t, slong k, slong q, gr_ctx_t ctx)
{
    qqbar_t z;
    gr_ctx_t QQbar;
    int status;
    slong g = n_gcd(k, q);

    k /= g;
    q /= g;

    if (q == 2)
        return gr_zero(res, ctx);

    if (q == 3 || q == 4 || q == 6)
    {
        /* cot(pi/3) = sqrt(1/3), cot(pi/4) = 1, cot(pi/6) = sqrt(3) */
        fmpq_t c;
        fmpq_init(c);
        fmpq_set_si(c, (q == 3) ? 1 : (q == 4) ? 1 : 3, (q == 3) ? 3 : 1);
        status = gr_set_fmpq(res, c, ctx);
        status |= gr_sqrt(res, res, ctx);
        fmpq_clear(c);
        return status;
    }

    qqbar_init(z);
    gr_ctx_init_complex_qqbar(QQbar);
    qqbar_root_of_unity(z, k, q);
    status = gr_set_other(t, z, QQbar, ctx);
    status |= gr_add_si(res, t, 1, ctx);
    status |= gr_sub_si(t, t, 1, ctx);
    if (status == GR_SUCCESS)
        status = gr_div(res, res, t, ctx);
    qqbar_root_of_unity(z, 1, 4);
    status |= gr_set_other(t, z, QQbar, ctx);
    status |= gr_mul(res, res, t, ctx);
    gr_ctx_clear(QQbar);
    qqbar_clear(z);
    return status;
}

/* res = (d/dz)^m cot(pi z) = pi^m P_m(cot(pi z)) at z = p/q (0 < p < q),
   with the cotangent in a cyclotomic field */
static int
_cot_pi_derivative_fmpq(gr_ptr res, slong p, slong q, ulong m, gr_ctx_t ctx)
{
    fmpz_poly_t P;
    gr_ptr c, t;
    slong i;
    int status;

    GR_TMP_INIT2(c, t, ctx);
    fmpz_poly_init(P);
    _gr_tower_cot_derivative_poly(P, m);

    if (2 * p > q)
    {
        status = _cot_pi_cyclo(c, t, q - p, q, ctx);
        status |= gr_neg(c, c, ctx);
    }
    else
        status = _cot_pi_cyclo(c, t, p, q, ctx);

    status |= gr_zero(res, ctx);
    for (i = fmpz_poly_degree(P); i >= 0 && status == GR_SUCCESS; i--)
    {
        status |= gr_mul(res, res, c, ctx);
        status |= gr_add_fmpz(res, res, P->coeffs + i, ctx);
    }
    status |= gr_pi(t, ctx);
    status |= gr_pow_ui(t, t, m, ctx);
    status |= gr_mul(res, res, t, ctx);

    fmpz_poly_clear(P);
    GR_TMP_CLEAR2(c, t, ctx);
    return status;
}

/* the reduced row echelon form of the relations at level q, with the
   column order of the unknowns: non-candidates (k/q > 1/2) first, then
   candidates by decreasing (denominator, numerator), zeta(s) last
   (gr_tower_hurwitz_lattice_struct: impl.h) */
static void
_hurwitz_lattice_init(gr_tower_hurwitz_lattice_struct * L, slong s, slong q)
{
    slong nvar = q, i, j, r, nrows;
    slong * order;
    fmpq_mat_t A;

    L->s = s;
    L->q = q;
    L->ncols = GR_TOWER_HURWITZ_NCOLS(q);
    L->col_of_k = flint_malloc(sizeof(slong) * (q + 1));
    L->k_of_col = flint_malloc(sizeof(slong) * nvar);
    L->row_of_col = flint_malloc(sizeof(slong) * L->ncols);
    order = flint_malloc(sizeof(slong) * nvar);

    _lattice_order(order, q);
    order[q - 1] = q;
    for (i = 0; i < nvar; i++)
    {
        L->col_of_k[order[i]] = i;
        L->k_of_col[i] = order[i];
    }
    flint_free(order);

    fmpq_mat_init(A, _gr_tower_hurwitz_relations_rows(q), L->ncols);
    nrows = _gr_tower_hurwitz_relations(A, L->col_of_k, s, q);
    {
        fmpq_mat_t W;
        fmpq_mat_window_init(W, A, 0, 0, nrows, L->ncols);
        fmpq_mat_init(L->B, nrows, L->ncols);
        L->rank = fmpq_mat_rref(L->B, W);
        fmpq_mat_window_clear(W);
    }
    fmpq_mat_clear(A);

    for (j = 0; j < L->ncols; j++)
        L->row_of_col[j] = -1;
    for (r = 0; r < L->rank; r++)
    {
        for (j = 0; j < L->ncols; j++)
            if (!fmpq_is_zero(fmpq_mat_entry(L->B, r, j)))
                break;
        if (j < L->ncols)
            L->row_of_col[j] = r;
    }
}

static void
_hurwitz_lattice_clear(gr_tower_hurwitz_lattice_struct * L)
{
    flint_free(L->col_of_k);
    flint_free(L->k_of_col);
    flint_free(L->row_of_col);
    fmpq_mat_clear(L->B);
}

/* the lattices are cached per thread by (s, q): the row echelon form at
   level q is the expensive part of a value in the normal form (and
   values at one level come in groups); an entry in use (refs > 0) is
   not replaced, since computing a value may need other levels */
#define HURWITZ_CACHE_SIZE 16

static FLINT_TLS_PREFIX gr_tower_hurwitz_lattice_struct * _hurwitz_cache[HURWITZ_CACHE_SIZE];
static FLINT_TLS_PREFIX slong _hurwitz_cache_refs[HURWITZ_CACHE_SIZE];
static FLINT_TLS_PREFIX slong _hurwitz_cache_next = 0;
static FLINT_TLS_PREFIX int _hurwitz_cache_registered = 0;

static void
_hurwitz_cache_cleanup(void)
{
    slong i;
    for (i = 0; i < HURWITZ_CACHE_SIZE; i++)
    {
        if (_hurwitz_cache[i] != NULL)
        {
            _hurwitz_lattice_clear(_hurwitz_cache[i]);
            flint_free(_hurwitz_cache[i]);
            _hurwitz_cache[i] = NULL;
        }
        _hurwitz_cache_refs[i] = 0;
    }
    _hurwitz_cache_next = 0;
    /* (flint_cleanup forgets the registration: register again when the
       cache is used afterwards) */
    _hurwitz_cache_registered = 0;
}

/* (to be released with _gr_tower_hurwitz_lattice_release) */
const gr_tower_hurwitz_lattice_struct *
_gr_tower_hurwitz_lattice(slong s, slong q)
{
    slong i, j;
    gr_tower_hurwitz_lattice_struct * L;

    for (i = 0; i < HURWITZ_CACHE_SIZE; i++)
    {
        if (_hurwitz_cache[i] != NULL && _hurwitz_cache[i]->s == s && _hurwitz_cache[i]->q == q)
        {
            _hurwitz_cache_refs[i]++;
            return _hurwitz_cache[i];
        }
    }

    if (!_hurwitz_cache_registered)
    {
        flint_register_cleanup_function(_hurwitz_cache_cleanup);
        _hurwitz_cache_registered = 1;
    }

    /* (round robin replacement among the entries not in use) */
    for (j = 0; j < HURWITZ_CACHE_SIZE; j++)
    {
        i = (_hurwitz_cache_next + j) % HURWITZ_CACHE_SIZE;
        if (_hurwitz_cache_refs[i] == 0)
            break;
    }

    if (j == HURWITZ_CACHE_SIZE)
    {
        /* all in use: not cached */
        L = flint_malloc(sizeof(gr_tower_hurwitz_lattice_struct));
        _hurwitz_lattice_init(L, s, q);
        return L;
    }

    _hurwitz_cache_next = (i + 1) % HURWITZ_CACHE_SIZE;
    if (_hurwitz_cache[i] != NULL)
        _hurwitz_lattice_clear(_hurwitz_cache[i]);
    else
        _hurwitz_cache[i] = flint_malloc(sizeof(gr_tower_hurwitz_lattice_struct));
    L = _hurwitz_cache[i];
    _hurwitz_lattice_init(L, s, q);
    _hurwitz_cache_refs[i] = 1;
    return L;
}

void
_gr_tower_hurwitz_lattice_release(const gr_tower_hurwitz_lattice_struct * L)
{
    slong i;
    for (i = 0; i < HURWITZ_CACHE_SIZE; i++)
    {
        if (_hurwitz_cache[i] == L)
        {
            _hurwitz_cache_refs[i]--;
            return;
        }
    }
    _hurwitz_lattice_clear((gr_tower_hurwitz_lattice_struct *) L);
    flint_free((gr_tower_hurwitz_lattice_struct *) L);
}

/* the value of the basis element zeta(s, k/q) (a generator, zeta(s), or
   pi^2 + 8 G for s = 2, k/q = 1/4) */
static int
_hurwitz_basis_value(gr_ptr res, slong s, slong k, slong q, gr_ctx_t ctx)
{
    fmpq_t c;
    fmpz_t f;
    int status;

    if (k == q)
        return _zeta_si(res, s, ctx);

    fmpq_init(c);
    fmpq_set_si(c, k, q);

    if (s == 2 && fmpz_is_one(fmpq_numref(c)) && fmpz_equal_ui(fmpq_denref(c), 4))
    {
        gr_ptr t;
        GR_TMP_INIT(t, ctx);
        status = _gr_tower_lazy_special_gen_locked(res, NULL, GR_TOWER_CONSTANT, GR_TOWER_CONST_CATALAN, ctx);
        status |= gr_mul_ui(res, res, 8, ctx);
        status |= gr_pi(t, ctx);
        status |= gr_sqr(t, t, ctx);
        status |= gr_add(res, res, t, ctx);
        GR_TMP_CLEAR(t, ctx);
        fmpq_clear(c);
        return status;
    }

    /* zeta(s, a) = (-1)^s psi^(s-1)(a) / (s-1)! */
    status = gr_set_fmpq(res, c, ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_special_gen_locked(res, res, GR_TOWER_POLYGAMMA, s - 1, ctx);
    fmpz_init(f);
    fmpz_fac_ui(f, s - 1);
    if (s % 2 == 1)
        fmpz_neg(f, f);
    status |= gr_div_fmpz(res, res, f, ctx);
    fmpz_clear(f);
    fmpq_clear(c);
    return status;
}

/* zeta(s, k/q) (0 < k <= q) in the normal form of the lattice */
static int
_hurwitz_nf(gr_ptr res, const gr_tower_hurwitz_lattice_struct * L, slong k, gr_ctx_t ctx)
{
    slong s = L->s, q = L->q, ct = L->col_of_k[k], row = L->row_of_col[ct], j;
    slong ce = GR_TOWER_HURWITZ_CE(q), cz = GR_TOWER_HURWITZ_CZ(q);
    gr_ptr t, u, E;
    int status = GR_SUCCESS;

    if (row < 0)
        return _hurwitz_basis_value(res, s, k, q, ctx);

    GR_TMP_INIT3(t, u, E, ctx);

    /* x_t = -sum_{j != t} B[row, j] x_j (basis values, E_k, Z) */
    status = gr_zero(res, ctx);
    for (j = 0; j < q && status == GR_SUCCESS; j++)
    {
        const fmpq * c = fmpq_mat_entry(L->B, row, j);
        if (j == ct || fmpq_is_zero(c))
            continue;
        status = _hurwitz_basis_value(t, s, L->k_of_col[j], q, ctx);
        status |= gr_mul_fmpq(t, t, c, ctx);
        status |= gr_sub(res, res, t, ctx);
    }

    /* elementary part: -sum_k B[row, E_k] E_k - B[row, Z] Z, with
       E_k = (-1)^(s-1) pi^s / (s-1)! P_{s-1}(cot(pi k/q)) */
    if (status == GR_SUCCESS)
    {
        fmpz_poly_t P;
        fmpz_t f;
        int any = 0;

        fmpz_poly_init(P);
        fmpz_init(f);
        _gr_tower_cot_derivative_poly(P, s - 1);
        status = gr_zero(E, ctx);
        if (q > HURWITZ_CYCLO_COT_MIN || REAL(ctx))
        {
            /* the cotangent sum computed in Q[x]/Phi_m, then evaluated at
               zeta_m once (rather than a division in the field per
               cotangent); in a real field, in the tangent normal form */
            ulong m = (q / n_gcd(q, 4)) * 4;
            fmpq * ck = _fmpq_vec_init(q / 2 + 1);
            fmpq_poly_t Ep;

            fmpq_poly_init(Ep);
            for (j = 1; 2 * j <= q; j++)
            {
                fmpq_set(ck + j, fmpq_mat_entry(L->B, row, ce + j - 1));
                if (!fmpq_is_zero(ck + j))
                    any = 1;
            }
            _gr_tower_cot_sum_cyclo(Ep, ck, q, P, m);
            if (q % 2 == 0 && !fmpq_is_zero(ck + q / 2))
            {
                /* cot(pi/2) = 0 */
                fmpq_t c0;
                fmpq_init(c0);
                fmpq_set_fmpz(c0, P->coeffs + 0);
                fmpq_mul(c0, c0, ck + q / 2);
                fmpq_poly_add_fmpq(Ep, Ep, c0);
                fmpq_clear(c0);
            }

            status = REAL(ctx) ? _gr_tower_lazy_real_cyclotomic_eval(E, Ep, m, ctx)
                               : _gr_tower_lazy_cyclotomic_eval(E, Ep, m, ctx);
            fmpq_poly_clear(Ep);
            _fmpq_vec_clear(ck, q / 2 + 1);
        }
        for (j = 1; 2 * j <= q && status == GR_SUCCESS && q <= HURWITZ_CYCLO_COT_MIN && !REAL(ctx); j++)
        {
            const fmpq * c = fmpq_mat_entry(L->B, row, ce + j - 1);
            if (fmpq_is_zero(c))
                continue;
            any = 1;
            status = _cot_pi_cyclo(u, t, j, q, ctx);
            /* t = P(u) by Horner */
            {
                slong i;
                status |= gr_zero(t, ctx);
                for (i = fmpz_poly_degree(P); i >= 0 && status == GR_SUCCESS; i--)
                {
                    status |= gr_mul(t, t, u, ctx);
                    status |= gr_add_fmpz(t, t, P->coeffs + i, ctx);
                }
            }
            status |= gr_mul_fmpq(t, t, c, ctx);
            status |= gr_add(E, E, t, ctx);
        }
        if (any && status == GR_SUCCESS)
        {
            /* times (-1)^s / (s-1)! (the minus sign of -B included) */
            fmpz_fac_ui(f, s - 1);
            if (s % 2 == 1)
                fmpz_neg(f, f);
            status = gr_div_fmpz(E, E, f, ctx);
            status |= gr_pi(t, ctx);
            status |= gr_pow_ui(t, t, s, ctx);
            status |= gr_mul(E, E, t, ctx);
            status |= gr_add(res, res, E, ctx);
        }
        if (status == GR_SUCCESS && !fmpq_is_zero(fmpq_mat_entry(L->B, row, cz)))
        {
            status = _zeta_si(t, s, ctx);
            status |= gr_mul_fmpq(t, t, fmpq_mat_entry(L->B, row, cz), ctx);
            status |= gr_sub(res, res, t, ctx);
        }
        fmpz_poly_clear(P);
        fmpz_clear(f);
    }

    GR_TMP_CLEAR3(t, u, E, ctx);
    return status;
}

/* zeta(s, p/q) for 0 < p <= q in the normal form (GR_UNABLE beyond the
   limits) */
static int
_hurwitz_fmpq_nf(gr_ptr res, slong s, slong p, slong q, gr_ctx_t ctx)
{
    if (s < 2 || s > gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_HURWITZ_WEIGHT_LIMIT) || q > gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_HURWITZ_LATTICE_LIMIT))
        return GR_UNABLE;

    {
        const gr_tower_hurwitz_lattice_struct * L = _gr_tower_hurwitz_lattice(s, q);
        int status = _hurwitz_nf(res, L, p, ctx);
        _gr_tower_hurwitz_lattice_release(L);
        return status;
    }
}

static int _polygamma_fmpq_unit(gr_ptr res, ulong m, const fmpq_t c0, gr_ctx_t ctx);

/* Li_s(exp(2 pi i a/N)) for 3 <= N <= HURWITZ_EXPAND_LIMIT, through the
   Hurwitz zeta values zeta(s, k/N): in the normal form for N <=
   the option GR_TOWER_OPT_HURWITZ_LATTICE_LIMIT, otherwise as generators (whose relations the
   zero test finds, hurwitz_relations.c) */
#define HURWITZ_EXPAND_LIMIT 240

static int
_polylog_root_of_unity(gr_ptr res, slong s, slong a, slong N, gr_ctx_t ctx)
{
    const gr_tower_hurwitz_lattice_struct * L = NULL;
    gr_ctx_t QQbar;
    qqbar_t w;
    gr_ptr t, u;
    fmpz_t f, fs;
    fmpq_t c;
    slong k;
    int lattice, status = GR_SUCCESS;

    if (s < 2 || s > gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_HURWITZ_WEIGHT_LIMIT) || N > HURWITZ_EXPAND_LIMIT)
        return GR_UNABLE;

    lattice = (N <= gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_HURWITZ_LATTICE_LIMIT));
    if (lattice)
        L = _gr_tower_hurwitz_lattice(s, N);
    gr_ctx_init_complex_qqbar(QQbar);
    qqbar_init(w);
    GR_TMP_INIT2(t, u, ctx);
    fmpz_init(f);
    fmpz_init(fs);
    fmpq_init(c);

    /* zeta(s, a) = (-1)^s psi^(s-1)(a) / (s-1)! */
    fmpz_fac_ui(fs, s - 1);
    if (s % 2 == 1)
        fmpz_neg(fs, fs);

    status = gr_zero(res, ctx);
    for (k = 1; k <= N && status == GR_SUCCESS; k++)
    {
        slong e = ((a * k) % N + N) % N;
        if (lattice)
            status = _hurwitz_nf(t, L, k, ctx);
        else if (k == N)
            status = _zeta_si(t, s, ctx);
        else
        {
            fmpq_set_si(c, k, N);
            status = _polygamma_fmpq_unit(t, s - 1, c, ctx);
            status |= gr_div_fmpz(t, t, fs, ctx);
        }
        if (e != 0 && status == GR_SUCCESS)
        {
            qqbar_root_of_unity(w, e, N);
            status = gr_set_other(u, w, QQbar, ctx);
            status |= gr_mul(t, t, u, ctx);
        }
        status |= gr_add(res, res, t, ctx);
    }

    fmpz_set_si(f, N);
    fmpz_pow_ui(f, f, s);
    status |= gr_div_fmpz(res, res, f, ctx);

    fmpz_clear(f);
    fmpz_clear(fs);
    fmpq_clear(c);
    GR_TMP_CLEAR2(t, u, ctx);
    qqbar_clear(w);
    gr_ctx_clear(QQbar);
    if (lattice)
        _gr_tower_hurwitz_lattice_release(L);
    return status;
}

/* psi^(m) at rational c0 in (0, 1] */
static int
_polygamma_fmpq_unit(gr_ptr res, ulong m, const fmpq_t c0, gr_ctx_t ctx)
{
    int status;

    if (m >= 1 && fmpz_cmp_ui(fmpq_denref(c0), 2) > 0 && fmpz_cmp_ui(fmpq_denref(c0), gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_HURWITZ_LATTICE_LIMIT)) <= 0 &&
        m + 1 <= gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_HURWITZ_WEIGHT_LIMIT))
    {
        /* psi^(m)(c0) = (-1)^(m+1) m! zeta(m+1, c0) */
        status = _hurwitz_fmpq_nf(res, m + 1, fmpz_get_si(fmpq_numref(c0)), fmpz_get_si(fmpq_denref(c0)), ctx);
        if (status == GR_SUCCESS)
        {
            fmpz_t f;
            fmpz_init(f);
            fmpz_fac_ui(f, m);
            if (m % 2 == 0)
                fmpz_neg(f, f);
            status = gr_mul_fmpz(res, res, f, ctx);
            fmpz_clear(f);
            return status;
        }
    }

    if (fmpz_is_one(fmpq_numref(c0)) && fmpz_cmp_ui(fmpq_denref(c0), 2) <= 0 && m >= 1)
    {
        /* psi^(m)(1) = (-1)^(m+1) m! zeta(m+1),
           psi^(m)(1/2) = (-1)^(m+1) m! (2^(m+1) - 1) zeta(m+1) */
        fmpz_t f, g;
        fmpz_init(f);
        fmpz_init(g);
        status = _zeta_si(res, m + 1, ctx);
        fmpz_fac_ui(f, m);
        if (m % 2 == 0)
            fmpz_neg(f, f);
        if (fmpz_equal_ui(fmpq_denref(c0), 2))
        {
            fmpz_one(g);
            fmpz_mul_2exp(g, g, m + 1);
            fmpz_sub_ui(g, g, 1);
            fmpz_mul(f, f, g);
        }
        status |= gr_mul_fmpz(res, res, f, ctx);
        fmpz_clear(f);
        fmpz_clear(g);
        return status;
    }

    if (m == 0 && fmpz_cmp_ui(fmpq_denref(c0), gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_GAUSS_DIGAMMA_LIMIT)) <= 0)
    {
        if (fmpz_is_one(fmpq_denref(c0)))
        {
            /* psi(1) = -gamma */
            status = _gr_tower_lazy_special_gen_locked(res, NULL, GR_TOWER_CONSTANT, GR_TOWER_CONST_EULER, ctx);
            return status | gr_neg(res, res, ctx);
        }
        return _digamma_gauss(res, fmpz_get_ui(fmpq_numref(c0)), fmpz_get_ui(fmpq_denref(c0)), ctx);
    }

    {
        fmpq_t half, w;
        fmpq_init(half);
        fmpq_init(w);
        fmpq_set_si(half, 1, 2);
        if (fmpq_cmp(c0, half) > 0)
        {
            /* psi^(m)(z) = (-1)^m psi^(m)(1 - z) - pi P_m(cot(pi z)) */
            gr_ptr t;
            GR_TMP_INIT(t, ctx);
            fmpq_sub_si(w, c0, 1);
            fmpq_neg(w, w);
            status = gr_set_fmpq(t, w, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_special_gen_locked(res, t, GR_TOWER_POLYGAMMA, m, ctx);
            if (m % 2 == 1)
                status |= gr_neg(res, res, ctx);
            /* (the cotangent in a cyclotomic field, like the normal form:
               values at one level then share their field) */
            if (fmpz_cmp_si(fmpq_denref(c0), gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_TRIG_ALGEBRAIC_LIMIT)) <= 0)
                status |= _cot_pi_derivative_fmpq(t, fmpz_get_si(fmpq_numref(c0)), fmpz_get_si(fmpq_denref(c0)), m, ctx);
            else
            {
                status |= gr_set_fmpq(t, c0, ctx);
                status |= _cot_pi_derivative(t, t, m, ctx);
            }
            status |= _pi_times(t, t, ctx);
            status |= gr_sub(res, res, t, ctx);
            GR_TMP_CLEAR(t, ctx);
        }
        else
        {
            status = gr_set_fmpq(res, c0, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_special_gen_locked(res, res, GR_TOWER_POLYGAMMA, m, ctx);
        }
        fmpq_clear(half);
        fmpq_clear(w);
    }

    return status;
}

/* res += sign * (-1)^m m! sum_{k} (x + k)^(-(m+1)), k = 0 .. n-1 */
static int
_polygamma_shift_sum(gr_ptr res, gr_srcptr x, slong n, ulong m, int sign, gr_ctx_t ctx)
{
    gr_ptr s, t;
    fmpz_t f;
    slong k;
    int status = GR_SUCCESS;

    GR_TMP_INIT2(s, t, ctx);
    fmpz_init(f);

    status |= gr_zero(s, ctx);
    for (k = 0; k < n && status == GR_SUCCESS; k++)
    {
        status |= gr_add_si(t, x, k, ctx);
        status |= gr_pow_ui(t, t, m + 1, ctx);
        if (status == GR_SUCCESS)
            status = gr_inv(t, t, ctx);
        status |= gr_add(s, s, t, ctx);
    }

    fmpz_fac_ui(f, m);
    if ((m % 2 == 1) != (sign < 0))
        fmpz_neg(f, f);
    status |= gr_mul_fmpz(s, s, f, ctx);
    status |= gr_add(res, res, s, ctx);

    fmpz_clear(f);
    GR_TMP_CLEAR2(s, t, ctx);
    return status;
}

static int
_polygamma(gr_ptr res, ulong m, gr_srcptr z_in, gr_ctx_t ctx)
{
    fmpq_t c;
    gr_ptr z, z0, t, u;
    slong n;
    int reflect, pole, status;

    fmpq_init(c);
    if (_is_rational(c, z_in, ctx))
    {
        fmpz_t nn;
        fmpq_t c0;

        if (fmpz_is_one(fmpq_denref(c)) && fmpz_sgn(fmpq_numref(c)) <= 0)
        {
            fmpq_clear(c);
            return GR_DOMAIN;
        }

        fmpz_init(nn);
        fmpq_init(c0);
        /* c = c0 + n with 0 < c0 <= 1 */
        fmpz_cdiv_q(nn, fmpq_numref(c), fmpq_denref(c));
        fmpz_sub_ui(nn, nn, 1);
        fmpq_sub_fmpz(c0, c, nn);

        if (_abs_gt_ui(nn, SHIFT_LIMIT))
            status = _gr_tower_lazy_special_gen_locked(res, z_in, GR_TOWER_POLYGAMMA, m, ctx);
        else
        {
            n = fmpz_get_si(nn);
            status = _polygamma_fmpq_unit(res, m, c0, ctx);
            GR_TMP_INIT(t, ctx);
            if (status == GR_SUCCESS && n > 0)
            {
                status = gr_set_fmpq(t, c0, ctx);
                status |= _polygamma_shift_sum(res, t, n, m, 1, ctx);
            }
            else if (status == GR_SUCCESS && n < 0)
            {
                status = gr_set_fmpq(t, c, ctx);
                status |= _polygamma_shift_sum(res, t, -n, m, -1, ctx);
            }
            GR_TMP_CLEAR(t, ctx);
        }

        fmpz_clear(nn);
        fmpq_clear(c0);
        fmpq_clear(c);
        return status;
    }
    fmpq_clear(c);

    GR_TMP_INIT4(z, z0, t, u, ctx);
    status = gr_set(z, z_in, ctx);     /* (res may alias the argument) */
    if (status != GR_SUCCESS)
        goto cleanup;

    status = _shift_reflect(z0, &n, &reflect, &pole, z, ctx);

    if (status == GR_UNABLE && (n > SHIFT_LIMIT || n < -SHIFT_LIMIT))
    {
        status = _gr_tower_lazy_special_gen_locked(res, z, GR_TOWER_POLYGAMMA, m, ctx);
        goto cleanup;
    }

    if (status == GR_SUCCESS && pole)
        status = GR_DOMAIN;
    if (status != GR_SUCCESS)
        goto cleanup;

    if (reflect)
    {
        /* psi^(m)(z0) = (-1)^m psi^(m)(1 - z0) - pi P_m(cot(pi z0)); when
           Re(z0) = 0, psi^(m)(1 + w) = psi^(m)(w) + (-1)^m m! w^(-m-1)
           with w = -z0 canonical */
        int c0;
        fmpq_t zero;
        fmpq_init(zero);
        status = _gr_tower_lazy_re_cmp(&c0, z0, zero, ctx);
        fmpq_clear(zero);

        status |= gr_neg(t, z0, ctx);
        if (status == GR_SUCCESS && c0 == 0)
        {
            status = _gr_tower_lazy_special_gen_locked(u, t, GR_TOWER_POLYGAMMA, m, ctx);
            status |= _polygamma_shift_sum(u, t, 1, m, 1, ctx);
        }
        else if (status == GR_SUCCESS)
        {
            status = gr_add_si(t, t, 1, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_special_gen_locked(u, t, GR_TOWER_POLYGAMMA, m, ctx);
        }

        if (m % 2 == 1)
            status |= gr_neg(u, u, ctx);
        status |= _cot_pi_derivative(t, z0, m, ctx);
        status |= _pi_times(t, t, ctx);
        status |= gr_sub(res, u, t, ctx);
    }
    else
    {
        status = _gr_tower_lazy_special_gen_locked(res, z0, GR_TOWER_POLYGAMMA, m, ctx);
    }

    /* psi^(m)(z0 + n) = psi^(m)(z0) + (-1)^m m! sum_{k<n} (z0 + k)^(-m-1);
       psi^(m)(z0 - n) = psi^(m)(z0) - (-1)^m m! sum_{k<n} (z + k)^(-m-1) */
    if (status == GR_SUCCESS && n > 0)
        status = _polygamma_shift_sum(res, z0, n, m, 1, ctx);
    else if (status == GR_SUCCESS && n < 0)
        status = _polygamma_shift_sum(res, z, -n, m, -1, ctx);

cleanup:
    GR_TMP_CLEAR4(z, z0, t, u, ctx);
    return status;
}

int
gr_tower_lazy_digamma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _gr_tower_lazy_view_finish(_polygamma(res, 0, x, ctx), res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* a nonnegative integer s (exactly) */
static int
_get_order(slong * m, gr_srcptr s, slong lo, gr_ctx_t ctx)
{
    fmpq_t c;
    int status = GR_UNABLE;
    fmpq_init(c);
    if (_is_rational(c, s, ctx) && fmpz_is_one(fmpq_denref(c)) && fmpz_fits_si(fmpq_numref(c)))
    {
        *m = fmpz_get_si(fmpq_numref(c));
        status = (*m >= lo && *m <= EXACT_LIMIT) ? GR_SUCCESS : GR_UNABLE;
    }
    fmpq_clear(c);
    return status;
}

int
gr_tower_lazy_polygamma(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t s, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    slong m;
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _get_order(&m, s, 0, ctx);
    if (status == GR_SUCCESS)
        status = _polygamma(res, m, x, ctx);
    status = _gr_tower_lazy_view_finish(status, res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* zeta(s, a) = (-1)^s psi^(s-1)(a) / (s-1)! for integer s >= 2 */
int
gr_tower_lazy_hurwitz_zeta(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t s, const gr_tower_lazy_elem_t a, gr_ctx_t ctx)
{
    slong m;
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _get_order(&m, s, 2, ctx);
    if (status != GR_SUCCESS)
    {
        /* other s: Bernoulli polynomials, character sums */
        status = _gr_tower_lazy_hurwitz_general(res, s, a, ctx);
        status = _gr_tower_lazy_view_finish(status, res, real, alg, ctx);
        _gr_tower_lazy_unlock(ctx);
        return status;
    }
    status = _gr_tower_lazy_hurwitz_zeta_int(res, m, a, ctx);
    status = _gr_tower_lazy_view_finish(status, res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* zeta(s, a) for an integer s >= 2 (GR_UNABLE above EXACT_LIMIT) */
int
_gr_tower_lazy_hurwitz_zeta_int(gr_ptr res, slong m, gr_srcptr a, gr_ctx_t ctx)
{
    int status;
    if (m < 2)
        return GR_DOMAIN;
    if (m > EXACT_LIMIT)
        return GR_UNABLE;
    status = _polygamma(res, m - 1, a, ctx);
    if (status == GR_SUCCESS)
    {
        fmpz_t f;
        fmpz_init(f);
        fmpz_fac_ui(f, m - 1);
        if (m % 2 == 1)
            fmpz_neg(f, f);
        status = gr_div_fmpz(res, res, f, ctx);
        fmpz_clear(f);
    }
    return status;
}

/* -------------------------------------------------------------------- */
/* polylogarithms                                                        */
/* -------------------------------------------------------------------- */

/*
    The dilogarithm under the anharmonic group {z, 1 - z, 1/z, 1/(1 - z),
    1 - 1/z, z/(z - 1)}: with principal branches (Li_2 continuous from
    below on [1, oo), as acb_polylog),

        s: Li_2(z) = -Li_2(1 - z) + pi^2/6 - log(z) log(1 - z)    (z not in (-oo, 0] or [1, oo)),
        t: Li_2(z) = -Li_2(1/z) - pi^2/6 - log(-z)^2/2            (z not in [0, oo)),
        r: Li_2(z) = -Li_2(z/(z - 1)) - log(1 - z)^2/2            (z not in [1, oo)),
        c: Li_2(z) = -Li_2(1/z) + pi^2/3 - log(z)^2/2 - i pi log(z)   (z > 1).

    The canonical argument is intrinsic: in [0, 1/2] for real z; for
    nonreal z, the orbit point in the closed region Re(w) <= 1/2,
    |w| <= 1, |w - 1| <= 1 (one of the six regions cut out by the lines
    and circles fixed by the involutions), with Im(w) > 0 on its boundary
    (where the other candidate is the conjugate). The orbit of a nonreal
    z avoids the real axis, so the steps s, t, r are valid there.
*/

/* one step: w = the orbit point, E = the elementary term (Li_2(z) = -Li_2(w) + E) */
static int
_dilog_step(gr_ptr w, gr_ptr E, int step, gr_srcptr z_in, gr_ctx_t ctx)
{
    gr_ptr a, b, p, z;
    int status;
    GR_TMP_INIT4(a, b, p, z, ctx);

    status = gr_set(z, z_in, ctx);     /* (w may alias z) */
    status |= gr_pi(p, ctx);
    status |= gr_sqr(p, p, ctx);

    if (step == 's')
    {
        /* w = 1 - z, E = pi^2/6 - log(z) log(1 - z) */
        status |= gr_sub_si(w, z, 1, ctx);
        status |= gr_neg(w, w, ctx);
        status |= gr_log(a, z, ctx);
        status |= gr_log(b, w, ctx);
        status |= gr_mul(a, a, b, ctx);
        status |= gr_div_ui(E, p, 6, ctx);
        status |= gr_sub(E, E, a, ctx);
    }
    else if (step == 't')
    {
        /* w = 1/z, E = -pi^2/6 - log(-z)^2/2 */
        status |= gr_inv(w, z, ctx);
        status |= gr_neg(a, z, ctx);
        status |= gr_log(a, a, ctx);
        status |= gr_sqr(a, a, ctx);
        status |= gr_div_ui(a, a, 2, ctx);
        status |= gr_div_ui(E, p, 6, ctx);
        status |= gr_add(E, E, a, ctx);
        status |= gr_neg(E, E, ctx);
    }
    else if (step == 'r')
    {
        /* w = z/(z - 1), E = -log(1 - z)^2/2 */
        status |= gr_sub_si(a, z, 1, ctx);
        status |= gr_div(w, z, a, ctx);
        status |= gr_neg(a, a, ctx);
        status |= gr_log(a, a, ctx);
        status |= gr_sqr(a, a, ctx);
        status |= gr_div_ui(a, a, 2, ctx);
        status |= gr_neg(E, a, ctx);
    }
    else
    {
        /* z > 1: w = 1/z, E = pi^2/3 - log(z)^2/2 - i pi log(z) */
        status |= gr_inv(w, z, ctx);
        status |= gr_log(a, z, ctx);
        status |= gr_sqr(b, a, ctx);
        status |= gr_div_ui(b, b, 2, ctx);
        status |= gr_div_ui(E, p, 3, ctx);
        status |= gr_sub(E, E, b, ctx);
        status |= gr_pi(b, ctx);
        status |= gr_mul(a, a, b, ctx);
        status |= gr_i(b, ctx);
        status |= gr_mul(a, a, b, ctx);
        status |= gr_sub(E, E, a, ctx);
    }

    GR_TMP_CLEAR4(a, b, p, z, ctx);
    return status;
}

/* the orbit point only (no elementary term) */
static int
_dilog_point(gr_ptr w, int step, gr_srcptr z, gr_ctx_t ctx)
{
    gr_ptr a;
    int status;
    if (step == 's')
    {
        status = gr_sub_si(w, z, 1, ctx);
        return status | gr_neg(w, w, ctx);
    }
    if (step == 't' || step == 'c')
        return gr_inv(w, z, ctx);
    GR_TMP_INIT(a, ctx);
    status = gr_sub_si(a, z, 1, ctx);
    status |= gr_div(w, z, a, ctx);
    GR_TMP_CLEAR(a, ctx);
    return status;
}

/* whether the nonreal w lies in the closed canonical region */
static int
_dilog_in_region(int * in, gr_srcptr w, gr_ctx_t ctx)
{
    gr_ptr a, b;
    fmpq_t half;
    int status, sg;

    GR_TMP_INIT2(a, b, ctx);
    fmpq_init(half);
    fmpq_set_si(half, 1, 2);

    *in = 0;
    status = _gr_tower_lazy_re_cmp(&sg, w, half, ctx);
    if (status == GR_SUCCESS && sg <= 0)
    {
        /* |w|^2 <= 1 */
        status = gr_conj(a, w, ctx);
        status |= gr_mul(a, a, w, ctx);
        status |= gr_sub_ui(a, a, 1, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_real_sign_fast(&sg, a, ctx);
        if (status == GR_SUCCESS && sg <= 0)
        {
            /* |w - 1|^2 <= 1 */
            status = gr_sub_ui(b, w, 1, ctx);
            status |= gr_conj(a, b, ctx);
            status |= gr_mul(a, a, b, ctx);
            status |= gr_sub_ui(a, a, 1, ctx);
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_real_sign_fast(&sg, a, ctx);
            if (status == GR_SUCCESS && sg <= 0)
                *in = 1;
        }
    }

    fmpq_clear(half);
    GR_TMP_CLEAR2(a, b, ctx);
    return status;
}

static int _polylog(gr_ptr res, slong s, gr_srcptr z, gr_ctx_t ctx);

/* Li_2(z) for z not 0, 1, not a root of unity: the canonical form */
static int
_dilog_reduce(gr_ptr res, gr_srcptr z_in, gr_ctx_t ctx)
{
    gr_ptr z, w, E, Etot;
    const char * word = "";
    int status, sign = 1, k;
    truth_t real;

    GR_TMP_INIT4(z, w, E, Etot, ctx);
    status = gr_set(z, z_in, ctx);      /* (res may alias z) */

    real = gr_tower_lazy_is_real(z, ctx);

    if (real == T_TRUE)
    {
        fmpq_t c;
        int s0, s1;
        fmpq_init(c);
        /* z < 0: r; z > 1: c; 1/2 < z < 1: s; then possibly s */
        status = _gr_tower_lazy_real_sign_fast(&s0, z, ctx);
        fmpq_one(c);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_re_cmp(&s1, z, c, ctx);
        if (status == GR_SUCCESS)
        {
            if (s0 < 0)
                word = "r";
            else if (s1 > 0)
                word = "c";
        }
        fmpq_clear(c);
    }
    else if (real == T_FALSE)
    {
        /* an orbit point which is a root of unity (where Li_2 is a sum of
           Hurwitz zeta values: Li_2(1 + i) through Li_2(-i)), Im > 0
           preferred; otherwise the orbit point in the canonical region
           (Im > 0 on its boundary) */
        static const char * words[6] = { "", "s", "t", "st", "ts", "r" };
        int found = -1, in, im, unity;
        for (k = 0; k < 6 && status == GR_SUCCESS; k++)
        {
            const char * p;
            fmpq_t r;
            fmpq_init(r);
            status = gr_set(w, z, ctx);
            for (p = words[k]; *p && status == GR_SUCCESS; p++)
                status = _dilog_point(w, *p, w, ctx);
            if (status == GR_SUCCESS && _gr_tower_lazy_root_of_unity_angle_locked(r, w, ctx) &&
                fmpz_cmp_ui(fmpq_denref(r), 2) > 0 && fmpz_cmp_ui(fmpq_denref(r), 240) <= 0)
            {
                if (found < 0 || fmpz_sgn(fmpq_numref(r)) > 0)
                    found = k;
            }
            fmpq_clear(r);
        }
        unity = (found >= 0);
        for (k = 0; k < 6 && status == GR_SUCCESS && !unity; k++)
        {
            /* the orbit point: apply the steps of the word (values only) */
            const char * p;
            status = gr_set(w, z, ctx);
            for (p = words[k]; *p && status == GR_SUCCESS; p++)
                status = _dilog_point(w, *p, w, ctx);
            if (status == GR_SUCCESS)
                status = _dilog_in_region(&in, w, ctx);
            if (status == GR_SUCCESS && in)
            {
                status = _gr_tower_lazy_im_sign(&im, w, ctx);
                if (status == GR_SUCCESS && im > 0)
                {
                    found = k;
                    break;
                }
                if (status == GR_SUCCESS && found < 0)
                    found = k;
            }
        }
        if (status == GR_SUCCESS && found < 0)
            status = GR_UNABLE;     /* (not expected) */
        if (status == GR_SUCCESS)
            word = words[found];
    }
    else
        status = GR_UNABLE;

    /* Li_2(z) = sign Li_2(w) + Etot along the word */
    if (status == GR_SUCCESS)
    {
        const char * p;
        status = gr_set(w, z, ctx);
        status |= gr_zero(Etot, ctx);
        for (p = word; *p && status == GR_SUCCESS; p++)
        {
            status = _dilog_step(w, E, *p, w, ctx);
            if (sign < 0)
                status |= gr_neg(E, E, ctx);
            status |= gr_add(Etot, Etot, E, ctx);
            sign = -sign;
        }

        /* real: from (0, 1) into (0, 1/2] */
        if (status == GR_SUCCESS && real == T_TRUE)
        {
            fmpq_t half;
            int sg;
            fmpq_init(half);
            fmpq_set_si(half, 1, 2);
            status = _gr_tower_lazy_re_cmp(&sg, w, half, ctx);
            if (status == GR_SUCCESS && sg > 0)
            {
                status = _dilog_step(w, E, 's', w, ctx);
                if (sign < 0)
                    status |= gr_neg(E, E, ctx);
                status |= gr_add(Etot, Etot, E, ctx);
                sign = -sign;
                word = "x";     /* (changed) */
            }
            fmpq_clear(half);
        }
    }

    if (status == GR_SUCCESS)
    {
        if (*word == '\0')
            status = _gr_tower_lazy_special_gen_locked(res, z, GR_TOWER_POLYLOG, 2, ctx);
        else
        {
            /* the canonical value (special values, roots of unity, or a
               generator) */
            status = _polylog(res, 2, w, ctx);
            if (sign < 0)
                status |= gr_neg(res, res, ctx);
            status |= gr_add(res, res, Etot, ctx);
        }
    }

    GR_TMP_CLEAR4(z, w, E, Etot, ctx);
    return status;
}

/*
    Li_s(z) for s >= 3 at z not 0, +-1, not a root of unity of the
    expanded orders: by the inversion formula

        Li_s(z) + (-1)^s Li_s(1/z) = -(2 pi i)^s / s! B_s(1/2 + log(-z) / (2 pi i))

    (z not in [0, 1]; on the cut z > 1, where Li_s is continuous from
    below, log(-z) = log(z) + pi i), the generators are kept at |z| < 1,
    or |z| = 1 with Im(z) > 0. Li_3(1/2) = 7 zeta(3)/8 - pi^2 log(2)/12 +
    log(2)^3/6 (Landen's identity at 1/2).
*/
static int _polylog(gr_ptr res, slong s, gr_srcptr z, gr_ctx_t ctx);

static int
_polylog_high(gr_ptr res, slong s, gr_srcptr z_in, gr_ctx_t ctx)
{
    gr_ptr z, w, u, v, L, E;
    fmpq_t c;
    int status, sg = 0, inv = 0, cut = 0;
    truth_t real;

    GR_TMP_INIT5(z, w, u, v, L, ctx);
    GR_TMP_INIT(E, ctx);
    fmpq_init(c);

    status = gr_set(z, z_in, ctx);      /* (res may alias z) */

    if (s == 3 && _is_rational(c, z, ctx) && fmpz_is_one(fmpq_numref(c)) && fmpz_equal_ui(fmpq_denref(c), 2))
    {
        status = _zeta_si(u, 3, ctx);
        status |= gr_mul_ui(u, u, 7, ctx);
        status |= gr_div_ui(u, u, 8, ctx);
        status |= gr_set_ui(L, 2, ctx);
        status |= gr_log(L, L, ctx);
        status |= gr_pi(v, ctx);
        status |= gr_sqr(v, v, ctx);
        status |= gr_mul(v, v, L, ctx);
        status |= gr_div_ui(v, v, 12, ctx);
        status |= gr_sub(u, u, v, ctx);
        status |= gr_pow_ui(v, L, 3, ctx);
        status |= gr_div_ui(v, v, 6, ctx);
        status |= gr_add(res, u, v, ctx);
        goto cleanup;
    }

    real = gr_tower_lazy_is_real(z, ctx);
    if (real == T_TRUE)
    {
        /* |z| > 1: z > 1 (the cut) or z < -1 */
        fmpq_one(c);
        status = _gr_tower_lazy_re_cmp(&sg, z, c, ctx);
        if (status == GR_SUCCESS && sg > 0)
            inv = cut = 1;
        else if (status == GR_SUCCESS)
        {
            fmpq_set_si(c, -1, 1);
            status = _gr_tower_lazy_re_cmp(&sg, z, c, ctx);
            inv = (status == GR_SUCCESS && sg < 0);
        }
    }
    else if (real == T_FALSE)
    {
        /* |z|^2 - 1 */
        status = gr_conj(u, z, ctx);
        status |= gr_mul(u, u, z, ctx);
        status |= gr_sub_ui(u, u, 1, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_real_sign_fast(&sg, u, ctx);
        if (status == GR_SUCCESS && sg > 0)
            inv = 1;
        else if (status == GR_SUCCESS && sg == 0)
        {
            status = _gr_tower_lazy_im_sign(&sg, z, ctx);
            inv = (status == GR_SUCCESS && sg < 0);
        }
    }
    else
        status = GR_UNABLE;

    if (status != GR_SUCCESS)
        goto cleanup;

    if (!inv)
    {
        status = _gr_tower_lazy_special_gen_locked(res, z, GR_TOWER_POLYLOG, s, ctx);
        goto cleanup;
    }

    /* L = log(-z), or log(z) + pi i on the cut */
    if (cut)
    {
        status = gr_log(L, z, ctx);
        status |= gr_pi(u, ctx);
        status |= gr_i(v, ctx);
        status |= gr_mul(u, u, v, ctx);
        status |= gr_add(L, L, u, ctx);
    }
    else
    {
        status = gr_neg(L, z, ctx);
        status |= gr_log(L, L, ctx);
    }

    /* E = -(2 pi i)^s / s! B_s(1/2 + L / (2 pi i)) */
    if (status == GR_SUCCESS)
    {
        fmpq_poly_t B;
        slong k;
        fmpz_t f;
        fmpq_poly_init(B);
        fmpz_init(f);
        arith_bernoulli_polynomial(B, s);
        /* v = 2 pi i, u = 1/2 + L / v */
        status = gr_pi(v, ctx);
        status |= gr_mul_2exp_si(v, v, 1, ctx);
        status |= gr_i(u, ctx);
        status |= gr_mul(v, v, u, ctx);
        status |= gr_div(u, L, v, ctx);
        fmpq_set_si(c, 1, 2);
        status |= gr_add_fmpq(u, u, c, ctx);
        status |= gr_zero(E, ctx);
        for (k = fmpq_poly_degree(B); k >= 0 && status == GR_SUCCESS; k--)
        {
            status = gr_mul(E, E, u, ctx);
            fmpq_poly_get_coeff_fmpq(c, B, k);
            status |= gr_add_fmpq(E, E, c, ctx);
        }
        status |= gr_pow_ui(v, v, s, ctx);
        status |= gr_mul(E, E, v, ctx);
        fmpz_fac_ui(f, s);
        status |= gr_div_fmpz(E, E, f, ctx);
        status |= gr_neg(E, E, ctx);
        fmpq_poly_clear(B);
        fmpz_clear(f);
    }

    /* Li_s(z) = E - (-1)^s Li_s(1/z) */
    if (status == GR_SUCCESS)
        status = gr_inv(w, z, ctx);
    if (status == GR_SUCCESS)
        status = _polylog(u, s, w, ctx);
    if (status == GR_SUCCESS)
    {
        if (s % 2 == 0)
            status = gr_sub(res, E, u, ctx);
        else
            status = gr_add(res, E, u, ctx);
    }

cleanup:
    GR_TMP_CLEAR5(z, w, u, v, L, ctx);
    GR_TMP_CLEAR(E, ctx);
    fmpq_clear(c);
    return status;
}

static int
_polylog(gr_ptr res, slong s, gr_srcptr z, gr_ctx_t ctx)
{
    truth_t t;
    gr_ptr u, v;
    int status = GR_SUCCESS;

    t = gr_is_zero(z, ctx);
    if (t == T_TRUE)
        return gr_zero(res, ctx);
    if (t == T_UNKNOWN)
        return GR_UNABLE;

    t = gr_is_one(z, ctx);
    if (t == T_UNKNOWN)
        return GR_UNABLE;
    if (t == T_TRUE)
        return (s <= 1) ? GR_DOMAIN : _zeta_si(res, s, ctx);

    GR_TMP_INIT2(u, v, ctx);

    if (s <= 0)
    {
        /* Li_{-n}(z) = sum_{k=0}^{n} k! S(n+1, k+1) (z / (1 - z))^(k+1) */
        slong n = -s, k;
        fmpz_t c, f;
        fmpz_init(c);
        fmpz_init(f);
        status = gr_sub_si(u, z, 1, ctx);
        status |= gr_neg(u, u, ctx);
        if (status == GR_SUCCESS)
            status = gr_div(u, z, u, ctx);
        status |= gr_zero(res, ctx);
        status |= gr_set(v, u, ctx);
        for (k = 0; k <= n && status == GR_SUCCESS; k++)
        {
            gr_ptr w;
            GR_TMP_INIT(w, ctx);
            arith_stirling_number_2(c, n + 1, k + 1);
            fmpz_fac_ui(f, k);
            fmpz_mul(c, c, f);
            status |= gr_mul_fmpz(w, v, c, ctx);
            status |= gr_add(res, res, w, ctx);
            status |= gr_mul(v, v, u, ctx);
            GR_TMP_CLEAR(w, ctx);
        }
        fmpz_clear(c);
        fmpz_clear(f);
    }
    else if (s == 1)
    {
        /* -log(1 - z) */
        status = gr_sub_si(u, z, 1, ctx);
        status |= gr_neg(u, u, ctx);
        if (status == GR_SUCCESS)
            status = gr_log(res, u, ctx);
        status |= gr_neg(res, res, ctx);
    }
    else
    {
        fmpq_t c;
        int done = 0;

        fmpq_init(c);
        if (_is_rational(c, z, ctx))
        {
            if (fmpz_equal_si(fmpq_numref(c), -1) && fmpz_is_one(fmpq_denref(c)))
            {
                /* Li_s(-1) = (2^(1-s) - 1) zeta(s) */
                fmpq_t r;
                fmpq_init(r);
                fmpq_one(r);
                fmpq_div_2exp(r, r, s - 1);
                fmpq_sub_si(r, r, 1);
                status = _zeta_si(res, s, ctx);
                status |= gr_mul_fmpq(res, res, r, ctx);
                fmpq_clear(r);
                done = 1;
            }
            else if (s == 2 && fmpz_is_one(fmpq_numref(c)) && fmpz_equal_ui(fmpq_denref(c), 2))
            {
                /* Li_2(1/2) = pi^2/12 - log(2)^2/2 */
                status = gr_pi(u, ctx);
                status |= gr_sqr(u, u, ctx);
                status |= gr_div_ui(u, u, 12, ctx);
                status |= gr_set_ui(v, 2, ctx);
                if (status == GR_SUCCESS)
                    status = gr_log(v, v, ctx);
                status |= gr_sqr(v, v, ctx);
                status |= gr_div_ui(v, v, 2, ctx);
                status |= gr_sub(res, u, v, ctx);
                done = 1;
            }
        }
        else if (s <= gr_tower_lazy_ctx_get_option(ctx, GR_TOWER_OPT_HURWITZ_WEIGHT_LIMIT) && _gr_tower_lazy_root_of_unity_angle_locked(c, z, ctx) &&
                 fmpz_cmp_ui(fmpq_denref(c), 2) > 0 && fmpz_cmp_ui(fmpq_denref(c), HURWITZ_EXPAND_LIMIT) <= 0)
        {
            /* Li_s(exp(2 pi i a/N)) = N^(-s) sum_{k=1}^{N} exp(2 pi i a k/N) zeta(s, k/N)
               in the normal form of the Hurwitz zeta values */
            status = _polylog_root_of_unity(res, s, fmpz_get_si(fmpq_numref(c)), fmpz_get_si(fmpq_denref(c)), ctx);
            done = (status == GR_SUCCESS);
        }
        fmpq_clear(c);

        if (!done && s == 2)
            status = _dilog_reduce(res, z, ctx);
        else if (!done)
            status = _polylog_high(res, s, z, ctx);
    }

    GR_TMP_CLEAR2(u, v, ctx);
    return status;
}

int
gr_tower_lazy_polylog(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t s, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    fmpq_t c;
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    fmpq_init(c);
    if (_is_rational(c, s, ctx) && fmpz_is_one(fmpq_denref(c)) && !_abs_gt_ui(fmpq_numref(c), EXACT_LIMIT))
        status = _polylog(res, fmpz_get_si(fmpq_numref(c)), x, ctx);
    else
    {
        /* other s: zeta(s) at 1, (2^(1-s) - 1) zeta(s) at -1, the Lerch
           transcendent z Phi(z, s, 1) at other roots of unity */
        fmpq_t r;
        gr_ptr t;
        fmpq_init(r);
        GR_TMP_INIT(t, ctx);
        if (gr_is_one(x, ctx) == T_TRUE)
            status = _zeta(res, s, ctx);
        else if (_gr_tower_lazy_root_of_unity_angle_locked(r, x, ctx))
        {
            status = gr_one(t, ctx);
            if (status == GR_SUCCESS)
                status = gr_lerch_phi(t, x, s, t, ctx);
            if (status == GR_SUCCESS)
                status = gr_mul(res, t, x, ctx);
        }
        else
            status = GR_UNABLE;
        GR_TMP_CLEAR(t, ctx);
        fmpq_clear(r);
    }
    fmpq_clear(c);
    status = _gr_tower_lazy_view_finish(status, res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_dilog(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _gr_tower_lazy_view_finish(_polylog(res, 2, x, ctx), res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* complete elliptic integrals                                           */
/* -------------------------------------------------------------------- */

static int _elliptic(gr_ptr res, gr_srcptr m, int kind, gr_ctx_t ctx);

/*
    Singular values (Chowla-Selberg): at m = k_r^2 with K(1-m)/K(m) =
    sqrt(r), K is a product of gamma values at rationals and pi, and
    E follows from the elliptic alpha function,
    E = K + (pi/(4K) - alpha(r) K) / sqrt(r); at 1 - m,
    K(1 - m) = sqrt(r) K(m) and E(1 - m) from Legendre's relation.
    (Each entry checked numerically to 90 digits.)

        r = 2: m = 3 - 2 sqrt(2),        K = (sqrt(2)+1)^(1/2) Gamma(1/8) Gamma(3/8) / (2^(13/4) sqrt(pi)),
               alpha = sqrt(2) - 1
        r = 3: m = (2 - sqrt(3))/4,      K = 3^(1/4) Gamma(1/3)^3 / (2^(7/3) pi),
               alpha = (sqrt(3) - 1)/2
        r = 4: m = 17 - 12 sqrt(2),      K = (sqrt(2)+1) Gamma(1/4)^2 / (2^(7/2) sqrt(pi)),
               alpha = 2 (sqrt(2) - 1)^2
*/
static int
_ell_sv_m(gr_ptr res, int r, gr_ctx_t ctx)
{
    const char * s = (r == 2) ? "3-2*sqrt(2)" : (r == 3) ? "(2-sqrt(3))/4" : "17-12*sqrt(2)";
    return gr_set_str(res, s, ctx);
}

static int
_ell_sv_K(gr_ptr res, int r, gr_ctx_t ctx)
{
    gr_ptr t;
    fmpq_t q;
    int status;
    GR_TMP_INIT(t, ctx);
    fmpq_init(q);
    if (r == 2)
    {
        status = gr_set_str(res, "sqrt(sqrt(2)+1)/(2^(13/4)*sqrt(pi))", ctx);
        fmpq_set_si(q, 1, 8);
        status |= _gamma_fmpq(t, q, ctx);
        status |= gr_mul(res, res, t, ctx);
        fmpq_set_si(q, 3, 8);
        status |= _gamma_fmpq(t, q, ctx);
        status |= gr_mul(res, res, t, ctx);
    }
    else if (r == 3)
    {
        status = gr_set_str(res, "3^(1/4)/(2^(7/3)*pi)", ctx);
        fmpq_set_si(q, 1, 3);
        status |= _gamma_fmpq(t, q, ctx);
        status |= gr_pow_ui(t, t, 3, ctx);
        status |= gr_mul(res, res, t, ctx);
    }
    else
    {
        status = gr_set_str(res, "(sqrt(2)+1)/(2^(7/2)*sqrt(pi))", ctx);
        fmpq_set_si(q, 1, 4);
        status |= _gamma_fmpq(t, q, ctx);
        status |= gr_sqr(t, t, ctx);
        status |= gr_mul(res, res, t, ctx);
    }
    fmpq_clear(q);
    GR_TMP_CLEAR(t, ctx);
    return status;
}

/* sets *done if m is a tabulated singular value or its complement */
static int
_elliptic_singular(gr_ptr res, int * done, gr_srcptr m, int kind, gr_ctx_t ctx)
{
    acb_t z;
    int r, comp, status = GR_SUCCESS;
    double mv;

    *done = 0;
    acb_init(z);
    /* (a real m written through nonreal generators has an imaginary
       part only containing zero; the exact comparison decides) */
    if (gr_tower_lazy_get_acb(z, m, 64, ctx) != GR_SUCCESS || !arb_contains_zero(acb_imagref(z)) ||
        !arb_is_finite(acb_realref(z)))
    {
        acb_clear(z);
        return GR_SUCCESS;
    }
    mv = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    acb_clear(z);

    for (r = 2; r <= 4 && !*done && status == GR_SUCCESS; r++)
    {
        /* m_2 = 0.1716, m_3 = 0.0670, m_4 = 0.0294 */
        double mr = (r == 2) ? 0.17157287525381 : (r == 3) ? 0.066987298107781 : 0.029437251522859;
        for (comp = 0; comp < 2 && !*done; comp++)
        {
            double target = comp ? 1.0 - mr : mr;
            gr_ptr t, K, E, Kc, u;
            if (fabs(mv - target) > 1e-8)
                continue;
            GR_TMP_INIT5(t, K, E, Kc, u, ctx);
            status = _ell_sv_m(t, r, ctx);
            if (comp)
            {
                status |= gr_sub_ui(t, t, 1, ctx);
                status |= gr_neg(t, t, ctx);
            }
            if (status == GR_SUCCESS && gr_equal(t, m, ctx) == T_TRUE)
            {
                /* K, E at m_r */
                status = _ell_sv_K(K, r, ctx);
                if (kind == GR_TOWER_ELLIPTIC_E || comp)
                {
                    const char * al = (r == 2) ? "sqrt(2)-1" : (r == 3) ? "(sqrt(3)-1)/2" : "2*(sqrt(2)-1)^2";
                    status |= gr_set_str(u, al, ctx);
                    status |= gr_mul(u, u, K, ctx);            /* alpha K */
                    status |= gr_pi(E, ctx);
                    status |= gr_div_ui(E, E, 4, ctx);
                    status |= gr_div(E, E, K, ctx);
                    status |= gr_sub(E, E, u, ctx);
                    status |= gr_set_ui(u, r, ctx);
                    status |= gr_sqrt(u, u, ctx);
                    status |= gr_div(E, E, u, ctx);
                    status |= gr_add(E, E, K, ctx);
                }
                if (!comp)
                    status |= gr_set(res, (kind == GR_TOWER_ELLIPTIC_K) ? K : E, ctx);
                else
                {
                    /* K' = sqrt(r) K, E' = (pi/2 - E K' + K K') / K */
                    status |= gr_set_ui(u, r, ctx);
                    status |= gr_sqrt(u, u, ctx);
                    status |= gr_mul(Kc, K, u, ctx);
                    if (kind == GR_TOWER_ELLIPTIC_K)
                        status |= gr_set(res, Kc, ctx);
                    else
                    {
                        status |= gr_mul(t, K, Kc, ctx);
                        status |= gr_mul(u, E, Kc, ctx);
                        status |= gr_sub(t, t, u, ctx);
                        status |= gr_pi(u, ctx);
                        status |= gr_div_ui(u, u, 2, ctx);
                        status |= gr_add(t, t, u, ctx);
                        status |= gr_div(res, t, K, ctx);
                    }
                }
                *done = 1;
            }
            GR_TMP_CLEAR5(t, K, E, Kc, u, ctx);
        }
    }

    return status;
}

/* -------------------------------------------------------------------- */
/* anchors: values linked to the generators of the context               */
/* -------------------------------------------------------------------- */

/*
    Landen's transformation: with k' = sqrt(1 - m) and k1 = (1 - k') /
    (1 + k'),

        K(m) = (1 + k1) K(k1^2),  E(m) = (1 + k') E(k1^2) - k' K(m)

    for m off the cut [1, inf). A new K(m) or E(m) is linked to a
    generator K, E at a point of the Landen chain of m, down (m -> k1^2)
    or up (m -> 4 s / (1 + s)^2 with s = +-sqrt(m), |s| < 1, the inverse),
    within LANDEN_STEPS steps: one step is taken exactly, and the value
    at the neighbour comes through the hypergeometric evaluation, which
    takes the next one. (Landen's transformation is the modular equation
    of level 2: these values are those of the duplication formulas at
    commensurable points of the modular layer.)
*/
#define LANDEN_STEPS 3
#define LANDEN_DEPTH 6

/* the direction of a chain from am reaching an anchor: -1 (down), +1 (up,
   with the sign *sg of the first root), 0 (none) */
static int
_landen_search(int * sg, const acb_t am, gr_ctx_t ctx)
{
    acb_t c, kp, s;
    slong step, prec = 128;
    int dir = 0;

    acb_init(c);
    acb_init(kp);
    acb_init(s);

    /* down */
    acb_set(c, am);
    for (step = 0; step < LANDEN_STEPS && dir == 0; step++)
    {
        acb_sub_ui(kp, c, 1, prec);
        acb_neg(kp, kp);
        acb_sqrt(kp, kp, prec);
        acb_sub_ui(s, kp, 1, prec);
        acb_neg(s, s);
        acb_add_ui(kp, kp, 1, prec);
        acb_div(c, s, kp, prec);
        acb_sqr(c, c, prec);
        if (_gr_tower_lazy_hyp_anchor_present(c, 1, ctx))
            dir = -1;
    }

    /* up: the signs of the roots as the bits of b */
    {
        ulong b;
        for (b = 0; b < (UWORD(1) << LANDEN_STEPS) && dir == 0; b++)
        {
            acb_set(c, am);
            for (step = 0; step < LANDEN_STEPS && dir == 0; step++)
            {
                mag_t r;
                acb_sqrt(s, c, prec);
                if ((b >> step) & 1)
                    acb_neg(s, s);
                mag_init(r);
                acb_get_mag(r, s);
                if (mag_cmp_2exp_si(r, 0) >= 0)
                {
                    mag_clear(r);
                    break;
                }
                mag_clear(r);
                acb_add_ui(kp, s, 1, prec);
                acb_sqr(kp, kp, prec);
                acb_mul_2exp_si(c, s, 2);
                acb_div(c, c, kp, prec);
                if (_gr_tower_lazy_hyp_anchor_present(c, 1, ctx))
                {
                    dir = 1;
                    *sg = (b & 1) ? -1 : 1;
                }
            }
        }
    }

    acb_clear(c);
    acb_clear(kp);
    acb_clear(s);
    return dir;
}

static int _elliptic(gr_ptr res, gr_srcptr m, int kind, gr_ctx_t ctx);

static int
_elliptic_landen(gr_ptr res, int * done, gr_srcptr m, int kind, gr_ctx_t ctx)
{
    acb_t am;
    int dir, sg = 1, status = GR_SUCCESS;

    *done = 0;
    if (_gr_tower_lazy_hyp_anchored(ctx, 0) >= LANDEN_DEPTH)
        return GR_SUCCESS;

    acb_init(am);
    if (gr_tower_lazy_get_acb(am, m, 128, ctx) != GR_SUCCESS || !acb_is_finite(am))
    {
        acb_clear(am);
        return GR_SUCCESS;
    }
    dir = _landen_search(&sg, am, ctx);

    if (dir < 0)
    {
        /* off the cut [1, inf) */
        arb_t t;
        arb_init(t);
        arb_sub_ui(t, acb_realref(am), 1, 128);
        if (arb_contains_zero(acb_imagref(am)) && !arb_is_negative(t))
            dir = 0;
        arb_clear(t);
    }

    if (dir != 0)
    {
        gr_ptr kp, k1, m1, K1, E1, t;
        GR_TMP_INIT3(kp, k1, m1, ctx);
        GR_TMP_INIT3(K1, E1, t, ctx);

        if (dir < 0)
        {
            /* k' = sqrt(1 - m), k1 = (1 - k') / (1 + k'), m1 = k1^2 */
            status = gr_sub_ui(kp, m, 1, ctx);
            status |= gr_neg(kp, kp, ctx);
            status |= gr_sqrt(kp, kp, ctx);
            status |= gr_sub_ui(k1, kp, 1, ctx);
            status |= gr_neg(k1, k1, ctx);
            status |= gr_add_ui(t, kp, 1, ctx);
            status |= gr_div(k1, k1, t, ctx);
            status |= gr_sqr(m1, k1, ctx);
        }
        else
        {
            /* s = +-sqrt(m) (k1 of m1), m1 = 4 s / (1 + s)^2, k' of m1 = (1 - s) / (1 + s) */
            status = gr_sqrt(k1, m, ctx);
            if (sg < 0)
                status |= gr_neg(k1, k1, ctx);
            status |= gr_add_ui(t, k1, 1, ctx);
            status |= gr_sqr(m1, t, ctx);
            status |= gr_div(m1, k1, m1, ctx);
            status |= gr_mul_2exp_si(m1, m1, 2, ctx);
            status |= gr_sub_ui(kp, k1, 1, ctx);
            status |= gr_neg(kp, kp, ctx);
            status |= gr_div(kp, kp, t, ctx);
        }

        if (status == GR_SUCCESS)
        {
            _gr_tower_lazy_hyp_anchored(ctx, 1);
            status = _elliptic(K1, m1, GR_TOWER_ELLIPTIC_K, ctx);
            if (status == GR_SUCCESS && kind == GR_TOWER_ELLIPTIC_E)
                status = _elliptic(E1, m1, GR_TOWER_ELLIPTIC_E, ctx);
            _gr_tower_lazy_hyp_anchored(ctx, -1);
        }

        if (status == GR_SUCCESS)
        {
            if (dir < 0)
            {
                /* K(m) = (1 + k1) K(m1), E(m) = (1 + k') E(m1) - k' K(m) */
                status = gr_add_ui(t, k1, 1, ctx);
                status |= gr_mul(K1, K1, t, ctx);
                if (kind == GR_TOWER_ELLIPTIC_K)
                    status |= gr_set(res, K1, ctx);
                else
                {
                    status |= gr_add_ui(t, kp, 1, ctx);
                    status |= gr_mul(E1, E1, t, ctx);
                    status |= gr_mul(t, kp, K1, ctx);
                    status |= gr_sub(res, E1, t, ctx);
                }
            }
            else
            {
                /* K(m1) = (1 + s) K(m), E(m1) = (1 + k'_1) E(m) - k'_1 K(m1) */
                if (kind == GR_TOWER_ELLIPTIC_K)
                {
                    status = gr_add_ui(t, k1, 1, ctx);
                    status |= gr_div(res, K1, t, ctx);
                }
                else
                {
                    status = gr_mul(t, kp, K1, ctx);
                    status |= gr_add(E1, E1, t, ctx);
                    status |= gr_add_ui(t, kp, 1, ctx);
                    status |= gr_div(res, E1, t, ctx);
                }
            }
            if (status == GR_SUCCESS)
                *done = 1;
        }

        /* (a failure leaves the generator to be created) */
        if (status != GR_SUCCESS)
            status = GR_SUCCESS;

        GR_TMP_CLEAR3(kp, k1, m1, ctx);
        GR_TMP_CLEAR3(K1, E1, t, ctx);
    }

    acb_clear(am);
    return status;
}

/*
    K(m), E(m) are evaluated as hypergeometric functions (lazy_hypgeom.c):
    the special values K(0), E(0), E(1), K(1/2), E(1/2) are summation
    theorems, the imaginary-modulus transformation (to a canonical
    argument with |m - 1| <= 1, Im(m) >= 0 on its boundary) is Pfaff's
    transformation, and 2F1 with parameters in (1/2, 1/2, 1) + Z^3 reduce
    to K and E. The generators themselves are created here, at the
    canonical arguments, after the singular values.
*/
static int
_elliptic(gr_ptr res, gr_srcptr m, int kind, gr_ctx_t ctx)
{
    /* at algebraic m, the singular values (Chowla-Selberg) first: the
       hypergeometric normal form would first compare the transforms of m
       with the arguments of the context, which costs much at algebraic
       numbers of high degree (lambda at a CM point, say) */
    if (_gr_tower_lazy_is_algebraic_repr_locked(m, ctx) == T_TRUE)
    {
        int done = 0, status;
        _gr_tower_lazy_lock(ctx);
        status = _elliptic_singular(res, &done, m, kind, ctx);
        if (status == GR_SUCCESS && !done)
            status = _gr_tower_lazy_elliptic_cm(res, &done, m, kind, ctx);
        _gr_tower_lazy_unlock(ctx);
        if (status != GR_SUCCESS || done)
            return status;
    }
    return _gr_tower_lazy_elliptic_hypgeom(res, m, kind, ctx);
}

int
_gr_tower_lazy_elliptic_gen(gr_ptr res, gr_srcptr m, int kind, gr_ctx_t ctx)
{
    int status, done = 0;

    _gr_tower_lazy_lock(ctx);
    status = _elliptic_singular(res, &done, m, kind, ctx);
    /* K at the other CM moduli of class number one (Chowla-Selberg) */
    if (status == GR_SUCCESS && !done)
        status = _gr_tower_lazy_elliptic_cm(res, &done, m, kind, ctx);
    /* through Landen's transformation from a generator of the context */
    if (status == GR_SUCCESS && !done)
        status = _elliptic_landen(res, &done, m, kind, ctx);
    if (status == GR_SUCCESS && !done)
        status = _gr_tower_lazy_special_gen_locked(res, m, kind, 0, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_elliptic_k(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _gr_tower_lazy_view_finish(_elliptic(res, x, GR_TOWER_ELLIPTIC_K, ctx), res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_elliptic_e(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status, real, alg;
    _gr_tower_lazy_lock(ctx);
    real = REAL(ctx);
    alg = ALG(ctx);
    status = _gr_tower_lazy_view_finish(_elliptic(res, x, GR_TOWER_ELLIPTIC_E, ctx), res, real, alg, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* constants                                                             */
/* -------------------------------------------------------------------- */

int
gr_tower_lazy_euler(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
{
    int status;
    _gr_tower_lazy_lock(ctx);
    status = ALG(ctx) ? GR_UNABLE : _gr_tower_lazy_special_gen_locked(res, NULL, GR_TOWER_CONSTANT, GR_TOWER_CONST_EULER, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_catalan(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
{
    int status;
    _gr_tower_lazy_lock(ctx);
    status = ALG(ctx) ? GR_UNABLE : _gr_tower_lazy_special_gen_locked(res, NULL, GR_TOWER_CONSTANT, GR_TOWER_CONST_CATALAN, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* -------------------------------------------------------------------- */
/* dispatch (conjugation of generators)                                  */
/* -------------------------------------------------------------------- */

int
_gr_tower_lazy_special(gr_tower_lazy_elem_t res, int kind, slong param, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status;

    _gr_tower_lazy_lock(ctx);
    switch (kind)
    {
        case GR_TOWER_GAMMA: status = _gamma(res, x, ctx); break;
        case GR_TOWER_ERF:
            if (param == 0)
                status = _erf(res, x, ctx);
            else
            {
                gr_ptr i, t;
                GR_TMP_INIT2(i, t, ctx);
                status = gr_i(i, ctx);
                status |= gr_mul(t, x, i, ctx);
                if (status == GR_SUCCESS)
                    status = _erf(t, t, ctx);
                status |= gr_mul(res, t, i, ctx);
                status |= gr_neg(res, res, ctx);
                GR_TMP_CLEAR2(i, t, ctx);
            }
            break;
        case GR_TOWER_LAMBERTW: status = _lambertw(res, x, param, ctx); break;
        case GR_TOWER_POLYGAMMA: status = _polygamma(res, param, x, ctx); break;
        case GR_TOWER_POLYLOG: status = _polylog(res, param, x, ctx); break;
        case GR_TOWER_ZETA: status = _zeta(res, x, ctx); break;
        case GR_TOWER_DIRICHLET_L:
            status = _gr_tower_lazy_dirichlet_l_prim(res, x, GR_TOWER_DIRICHLET_Q(param), GR_TOWER_DIRICHLET_K(param), ctx);
            break;
        case GR_TOWER_ELLIPTIC_K:
        case GR_TOWER_ELLIPTIC_E: status = _elliptic(res, x, kind, ctx); break;
        case GR_TOWER_MODULAR_LAMBDA: status = _gr_tower_lazy_modular_lambda(res, x, ctx); break;
        case GR_TOWER_CONSTANT: status = _gr_tower_lazy_special_gen_locked(res, NULL, kind, param, ctx); break;
        default: status = GR_UNABLE;
    }
    _gr_tower_lazy_unlock(ctx);
    return status;
}

POP_OPTIONS
