/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Numerical evaluation of special function generators. The canonical
    choice of arguments and the relations between values (functional
    equations, special values) are the responsibility of the callers
    (see lazy_special.c); this file only provides enclosures.
*/

#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_elliptic.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "ulong_extras.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"
#include "bernoulli.h"
#include "fmpz_poly.h"

/* P_j with (d/dx)^j cot(x) = P_j(cot(x)): P_0 = c, P_{j+1} = -(1 + c^2) P_j' */
void
_gr_tower_cot_derivative_poly(fmpz_poly_t P, ulong j)
{
    fmpz_poly_t R, D;
    ulong i;
    fmpz_poly_init(R);
    fmpz_poly_init(D);
    fmpz_poly_set_coeff_si(R, 0, -1);
    fmpz_poly_set_coeff_si(R, 2, -1);
    fmpz_poly_zero(P);
    fmpz_poly_set_coeff_si(P, 1, 1);
    for (i = 0; i < j; i++)
    {
        fmpz_poly_derivative(D, P);
        fmpz_poly_mul(P, D, R);
    }
    fmpz_poly_clear(R);
    fmpz_poly_clear(D);
}

/* zeta(s) / pi^s for even s >= 2: (-1)^(s/2+1) B_s 2^(s-1) / s! */
void
_gr_tower_zeta_even_over_pi(fmpq_t res, ulong s)
{
    fmpz_t f;
    fmpz_init(f);
    bernoulli_fmpq_ui(res, s);
    fmpz_fac_ui(f, s);
    fmpq_div_fmpz(res, res, f);
    fmpq_mul_2exp(res, res, s - 1);
    if ((s / 2) % 2 == 0)
        fmpq_neg(res, res);
    fmpz_clear(f);
}

int
_gr_tower_special_eval(acb_t res, int kind, slong param, const acb_t u, slong prec)
{
    switch (kind)
    {
        case GR_TOWER_GAMMA:
            acb_gamma(res, u, prec);
            break;
        case GR_TOWER_ERF:
            if (param == 0)
                acb_hypgeom_erf(res, u, prec);
            else
                acb_hypgeom_erfi(res, u, prec);
            break;
        case GR_TOWER_LAMBERTW:
            {
                fmpz_t k;
                fmpz_init_set_si(k, param);
                acb_lambertw(res, u, k, 0, prec);
                fmpz_clear(k);
            }
            break;
        case GR_TOWER_POLYGAMMA:
            {
                acb_t s;
                acb_init(s);
                acb_set_si(s, param);
                acb_polygamma(res, s, u, prec);
                acb_clear(s);
            }
            break;
        case GR_TOWER_POLYLOG:
            acb_polylog_si(res, param, u, prec);
            break;
        case GR_TOWER_ZETA:
            acb_zeta(res, u, prec);
            break;
        case GR_TOWER_ELLIPTIC_K:
            acb_elliptic_k(res, u, prec);
            break;
        case GR_TOWER_ELLIPTIC_E:
            acb_elliptic_e(res, u, prec);
            break;
        case GR_TOWER_CONSTANT:
            if (param == GR_TOWER_CONST_EULER)
                arb_const_euler(acb_realref(res), prec);
            else if (param == GR_TOWER_CONST_CATALAN)
                arb_const_catalan(acb_realref(res), prec);
            else
                return GR_UNABLE;
            arb_zero(acb_imagref(res));
            break;
        default:
            return GR_UNABLE;
    }

    return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
}

int
_gr_tower_special_real_at(int kind, slong param, const acb_t u, slong prec)
{
    arb_t t;
    int res;

    if (kind == GR_TOWER_CONSTANT)
        return 1;

    if (!arb_is_zero(acb_imagref(u)))
        return 0;

    arb_init(t);

    switch (kind)
    {
        case GR_TOWER_GAMMA:
        case GR_TOWER_ERF:
        case GR_TOWER_POLYGAMMA:
        case GR_TOWER_ZETA:
            res = 1;
            break;
        case GR_TOWER_POLYLOG:
        case GR_TOWER_ELLIPTIC_K:
        case GR_TOWER_ELLIPTIC_E:
            /* real on (-inf, 1) */
            arb_sub_ui(t, acb_realref(u), 1, prec);
            res = arb_is_negative(t);
            break;
        case GR_TOWER_LAMBERTW:
            /* W_0 is real on (-1/e, inf), W_{-1} on (-1/e, 0) */
            arb_const_e(t, prec);
            arb_inv(t, t, prec);
            arb_add(t, t, acb_realref(u), prec);    /* u + 1/e */
            if (param == 0)
                res = arb_is_positive(t);
            else if (param == -1)
                res = arb_is_positive(t) && arb_is_negative(acb_realref(u));
            else
                res = 0;
            break;
        default:
            res = 0;
    }

    arb_clear(t);
    return res;
}

/* Writes the definition of a special function value, with the argument
   given as a string (ignored for constants). */
int
_gr_tower_special_write(gr_stream_t out, int kind, slong param, const char * arg)
{
    int status = GR_SUCCESS;

    switch (kind)
    {
        case GR_TOWER_CONSTANT:
            return gr_stream_write(out, (param == GR_TOWER_CONST_EULER) ? "euler" :
                                        (param == GR_TOWER_CONST_CATALAN) ? "catalan" : "?");
        case GR_TOWER_GAMMA: status |= gr_stream_write(out, "gamma("); break;
        case GR_TOWER_ERF: status |= gr_stream_write(out, (param == 0) ? "erf(" : "erfi("); break;
        case GR_TOWER_ZETA: status |= gr_stream_write(out, "zeta("); break;
        case GR_TOWER_ELLIPTIC_K: status |= gr_stream_write(out, "elliptic_k("); break;
        case GR_TOWER_ELLIPTIC_E: status |= gr_stream_write(out, "elliptic_e("); break;
        case GR_TOWER_LAMBERTW: status |= gr_stream_write(out, "lambertw("); break;
        case GR_TOWER_POLYGAMMA:
            if (param == 0)
                status |= gr_stream_write(out, "digamma(");
            else
            {
                status |= gr_stream_write(out, "polygamma(");
                status |= gr_stream_write_si(out, param);
                status |= gr_stream_write(out, ", ");
            }
            break;
        case GR_TOWER_POLYLOG:
            status |= gr_stream_write(out, "polylog(");
            status |= gr_stream_write_si(out, param);
            status |= gr_stream_write(out, ", ");
            break;
        default:
            return gr_stream_write(out, "?");
    }

    status |= gr_stream_write(out, arg);

    if (kind == GR_TOWER_LAMBERTW && param != 0)
    {
        status |= gr_stream_write(out, ", ");
        status |= gr_stream_write_si(out, param);
    }

    status |= gr_stream_write(out, ")");
    return status;
}

/*
    The distribution relations of the gamma function at level q (see
    lazy_special.c): with x_k = log Gamma(k/q),

        x_k + x_{q-k} + log sin(pi k/q) - log pi = 0,
        sum_{j<n} x_{k + j q/n} - x_{nk} - (n-1)/2 (log 2 + log pi)
            - (1/2 - nk/q) log n = 0      (n | q, 0 < nk < q).

    Column layout (impl.h): x_k in column col_of_k[k] (0 < k < q), then
    log sin(pi k/q) for 1 <= k < q/2, log n for 2 <= n <= q (only the
    columns of primes are used: log n is written as a combination of
    logarithms of primes), and log pi. A reduced row echelon form thus
    prefers expressing constants through pi, then integers, then sines
    (the sines are pivots of the relations between constants, such as
    prod_k sin(pi k/q) = q / 2^(q-1), whenever possible), which keeps the
    algebraic factors simple. The matrix must have
    _gr_tower_gamma_relations_rows(q) rows (zero on input) and
    GR_TOWER_GAMMA_NCOLS(q) columns; returns the number of rows used.
*/
slong
_gr_tower_gamma_relations_rows(slong q)
{
    slong n, rows = q / 2 + 1;
    for (n = 2; n <= q; n++)
        if (q % n == 0)
            rows += q / n;
    return rows;
}

slong
_gr_tower_gamma_relations(fmpq_mat_t A, const slong * col_of_k, slong q)
{
    slong cpi = GR_TOWER_GAMMA_CPI(q), cint = GR_TOWER_GAMMA_CINT(q), csin = GR_TOWER_GAMMA_CSIN(q);
    slong nrows = 0, k, n, j;
    fmpq_t t;

    fmpq_init(t);

    /* reflection */
    for (k = 1; 2 * k <= q; k++)
    {
        fmpq_add_si(fmpq_mat_entry(A, nrows, col_of_k[k]), fmpq_mat_entry(A, nrows, col_of_k[k]), 1);
        fmpq_add_si(fmpq_mat_entry(A, nrows, col_of_k[q - k]), fmpq_mat_entry(A, nrows, col_of_k[q - k]), 1);
        if (2 * k < q)
            fmpq_set_si(fmpq_mat_entry(A, nrows, csin + k - 1), 1, 1);
        fmpq_set_si(fmpq_mat_entry(A, nrows, cpi), -1, 1);
        nrows++;
    }

    /* multiplication */
    for (n = 2; n <= q; n++)
    {
        if (q % n != 0)
            continue;
        for (k = 1; k * n < q; k++)
        {
            for (j = 0; j < n; j++)
            {
                slong kk = k + j * (q / n);
                fmpq_add_si(fmpq_mat_entry(A, nrows, col_of_k[kk]), fmpq_mat_entry(A, nrows, col_of_k[kk]), 1);
            }
            fmpq_sub_si(fmpq_mat_entry(A, nrows, col_of_k[n * k]), fmpq_mat_entry(A, nrows, col_of_k[n * k]), 1);
            fmpq_set_si(t, -(n - 1), 2);
            fmpq_add(fmpq_mat_entry(A, nrows, cpi), fmpq_mat_entry(A, nrows, cpi), t);
            fmpq_add(fmpq_mat_entry(A, nrows, cint), fmpq_mat_entry(A, nrows, cint), t);
            /* -(1/2 - nk/q) log n = (2nk - q) / (2q) sum_p v_p(n) log p */
            {
                slong m = n, p;
                for (p = 2; m > 1; p++)
                {
                    slong v = 0;
                    while (m % p == 0)
                    {
                        m /= p;
                        v++;
                    }
                    if (v != 0)
                    {
                        fmpq_set_si(t, (2 * n * k - q) * v, 2 * q);
                        fmpq_add(fmpq_mat_entry(A, nrows, cint + p - 2), fmpq_mat_entry(A, nrows, cint + p - 2), t);
                    }
                }
            }
            nrows++;
        }
    }

    fmpq_clear(t);
    return nrows;
}

/*
    The relations of the Hurwitz zeta function zeta(s, x) (s >= 2) at
    level q (see lazy_special.c): with x_k = zeta(s, k/q), 0 < k <= q
    (x_q = zeta(s)),

        x_k + (-1)^s x_{q-k} - E_k = 0      (1 <= k <= q/2),
        sum_{j<m} x_{b + j q/m} - m^s x_{bm} = 0      (m | q prime, 1 <= b <= q/m),
        x_q - Z = 0                          (s even),

    where E_k = zeta(s, k/q) + (-1)^s zeta(s, 1 - k/q) =
    (-1)^(s-1) pi^s / (s-1)! P_{s-1}(cot(pi k/q)) (reflection, elementary)
    and Z = zeta(s) (a rational multiple of pi^s for even s). For odd s
    and k = q/2 the reflection says E_k = 0 and is omitted.

    Column layout (impl.h): x_k in column col_of_k[k] (0 < k <= q), then
    E_k for 1 <= k <= q/2, then Z. The matrix must have
    _gr_tower_hurwitz_relations_rows(q) rows (zero on input) and
    GR_TOWER_HURWITZ_NCOLS(q) columns; returns the number of rows used.
*/
slong
_gr_tower_hurwitz_relations_rows(slong q)
{
    slong m, rows = q / 2 + 1;
    for (m = 2; m <= q; m++)
        if (q % m == 0 && n_is_prime(m))
            rows += q / m;
    return rows;
}

slong
_gr_tower_hurwitz_relations(fmpq_mat_t A, const slong * col_of_k, slong s, slong q)
{
    slong ce = GR_TOWER_HURWITZ_CE(q), cz = GR_TOWER_HURWITZ_CZ(q);
    slong nrows = 0, k, m, b, j;
    fmpz_t ms;

    fmpz_init(ms);

    /* reflection */
    for (k = 1; 2 * k <= q; k++)
    {
        if (2 * k == q && s % 2 == 1)
            continue;
        fmpq_add_si(fmpq_mat_entry(A, nrows, col_of_k[k]), fmpq_mat_entry(A, nrows, col_of_k[k]), 1);
        fmpq_add_si(fmpq_mat_entry(A, nrows, col_of_k[q - k]), fmpq_mat_entry(A, nrows, col_of_k[q - k]), (s % 2 == 0) ? 1 : -1);
        fmpq_set_si(fmpq_mat_entry(A, nrows, ce + k - 1), -1, 1);
        nrows++;
    }

    /* zeta(s) for even s */
    if (s % 2 == 0)
    {
        fmpq_set_si(fmpq_mat_entry(A, nrows, col_of_k[q]), 1, 1);
        fmpq_set_si(fmpq_mat_entry(A, nrows, cz), -1, 1);
        nrows++;
    }

    /* distribution */
    for (m = 2; m <= q; m++)
    {
        if (q % m != 0 || !n_is_prime(m))
            continue;
        fmpz_set_ui(ms, m);
        fmpz_pow_ui(ms, ms, s);
        for (b = 1; b <= q / m; b++)
        {
            for (j = 0; j < m; j++)
            {
                slong kk = b + j * (q / m);
                fmpq_add_si(fmpq_mat_entry(A, nrows, col_of_k[kk]), fmpq_mat_entry(A, nrows, col_of_k[kk]), 1);
            }
            fmpz_sub(fmpq_numref(fmpq_mat_entry(A, nrows, col_of_k[b * m])),
                     fmpq_numref(fmpq_mat_entry(A, nrows, col_of_k[b * m])), ms);
            nrows++;
        }
    }

    fmpz_clear(ms);
    return nrows;
}
