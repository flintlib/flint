/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz_mat.h"
#include "impl.h"

/* Shell handler: Fincke-Pohst enumeration of the lattice points (any
   positive definite F), into a table of all counts. */

typedef struct
{
    ulong * C;
    ulong len;
}
enum_struct;

/* F(x) = sum_i q_ii (x_i + sum_{j>i} q_ij x_j)^2 in double (for the
   enumeration bounds only); returns 0 if not positive */
int
_acb_mzeta_ldl_d(double * q, const fmpz_mat_t F)
{
    slong d = fmpz_mat_nrows(F), i, j, k;
    for (i = 0; i < d; i++)
    {
        double t = fmpz_get_d(fmpz_mat_entry(F, i, i));
        for (k = 0; k < i; k++)
            t -= q[k * d + k] * q[k * d + i] * q[k * d + i];
        if (!(t > 0))
            return 0;
        q[i * d + i] = t;
        for (j = i + 1; j < d; j++)
        {
            double u = fmpz_get_d(fmpz_mat_entry(F, i, j));
            for (k = 0; k < i; k++)
                u -= q[k * d + k] * q[k * d + i] * q[k * d + j];
            q[i * d + j] = u / t;
        }
    }
    return 1;
}

static double
enum_cost(const void * data, const _acb_mzeta_shells_struct * T, ulong len)
{
    slong d = T->d, i;
    double q[64], det = 1.0, Vmax = (double) T->g * len + T->r0, vol;
    /* memory for the counts */
    if ((double) len * T->N > 2e8 || !_acb_mzeta_ldl_d(q, T->F))
        return INFINITY;
    for (i = 0; i < d; i++)
        det *= q[i * d + i];
    /* omega_d Vmax^(d/2) / (sqrt(det F) S^d) */
    vol = pow(3.141592653589793, d / 2.0) / tgamma(d / 2.0 + 1) * pow(Vmax, d / 2.0)
        / sqrt(det) / pow((double) T->S, d);
    return 6.0 * (2.0 * vol + 100.0) + 50.0 * pow(2.0 * vol + 100.0, (d - 1.0) / d);
}

/* counts by Fincke-Pohst enumeration of x = S n - s with exact (slong)
   shell values; the double bounds only need to contain the points
   (generous margins, exact filtering).  Returns 0 if the values could
   overflow. */
static int
counts_enum(ulong * C, const _acb_mzeta_shells_struct * T, ulong len)
{
    slong d = T->d, i, j;
    double q[64], c[8], P[9], xmax[8];
    slong n[8], lo[8], hi[8], Fi[64], si[8], li[8];
    ulong N = T->N, g = T->g, r0 = T->r0, S = T->S;
    double Vmax = (double) g * (len - 1) + r0, bnd;
    slong Vm = (slong) (g * (len - 1) + r0);
    fmpz_mat_t Finv;
    fmpz_t den;

    if (!_acb_mzeta_ldl_d(q, T->F))
        return 0;

    /* |x_i| <= sqrt(Vmax (F^-1)_ii); require sum |F_ij| xmax_i xmax_j
       < 2^(FLINT_BITS - 3) for the exact (slong) values */
    fmpz_mat_init(Finv, d, d);
    fmpz_init(den);
    fmpz_mat_inv(Finv, den, T->F);
    bnd = 0.0;
    for (i = 0; i < d; i++)
        xmax[i] = sqrt(Vmax * fabs(fmpz_get_d(fmpz_mat_entry(Finv, i, i)) / fmpz_get_d(den))) * 1.01 + S + 2;
    for (i = 0; i < d; i++)
        for (j = 0; j < d; j++)
        {
            if (fmpz_bits(fmpz_mat_entry(T->F, i, j)) > FLINT_BITS - 24)
                bnd = INFINITY;
            bnd += fabs(fmpz_get_d(fmpz_mat_entry(T->F, i, j))) * xmax[i] * xmax[j];
        }
    fmpz_mat_clear(Finv);
    fmpz_clear(den);
    if (!(bnd < ldexp(1.0, FLINT_BITS - 3)))
        return 0;

    for (i = 0; i < d; i++)
    {
        for (j = 0; j < d; j++)
            Fi[i * d + j] = fmpz_get_si(fmpz_mat_entry(T->F, i, j));
        si[i] = fmpz_get_si(T->s + i);
        li[i] = fmpz_get_si(T->ell + i);
    }

    /* levels i = d-1, ..., 0; P[i] = partial sum of the levels > i - 1 */
    P[d] = 0.0;
    i = d - 1;
    goto enter;

    while (1)
    {
        /* next value at level i, or go up */
        n[i]++;
        if (n[i] > hi[i])
        {
            i++;
            if (i >= d)
                break;
            continue;
        }

        if (i > 0)
        {
            double xi = (double) S * n[i] - si[i], y = xi - c[i];
            P[i] = P[i + 1] + q[i * d + i] * y * y;
            i--;
enter:
            /* the range of x_i: |x_i - c_i| <= sqrt((Vmax - P[i+1]) / q_ii) */
            {
                double r2, r, xl, xh;
                c[i] = 0.0;
                for (j = i + 1; j < d; j++)
                    c[i] -= q[i * d + j] * ((double) S * n[j] - si[j]);
                r2 = (Vmax - P[i + 1]) / q[i * d + i];
                r = sqrt(FLINT_MAX(r2, 0.0)) * (1 + 1e-9) + 1e-6 * (1 + fabs(c[i]));
                xl = c[i] - r;
                xh = c[i] + r;
                lo[i] = (slong) ceil((xl + si[i]) / S);
                hi[i] = (slong) floor((xh + si[i]) / S);
                n[i] = lo[i] - 1;
            }
            if (i > 0)
                continue;
        }

        /* innermost: the whole line n_0 in [lo_0, hi_0] with exact values
           V = F00 x0^2 + 2 x0 L + R */
        {
            slong L = 0, R = 0, cl = 0, x0, V, k;
            for (j = 1; j < d; j++)
            {
                slong xj = (slong) S * n[j] - si[j];
                L += Fi[j] * xj;
                cl += li[j] * n[j];
                for (k = 1; k < d; k++)
                    R += Fi[j * d + k] * xj * ((slong) S * n[k] - si[k]);
            }
            for (n[0] = lo[0]; n[0] <= hi[0]; n[0]++)
            {
                x0 = (slong) S * n[0] - si[0];
                V = Fi[0] * x0 * x0 + 2 * x0 * L + R;
                if (V <= Vm)
                {
                    ulong M = ((ulong) V - r0) / g;
                    ulong r = (ulong) ((((cl + li[0] * n[0]) % (slong) N) + (slong) N) % (slong) N);
                    C[M * N + r]++;
                }
            }
            n[0] = hi[0];
            i = 0;
            /* go up */
            i++;
            if (i >= d)
                break;
        }
    }

    return 1;
}

static void *
enum_init(const _acb_mzeta_shells_struct * T)
{
    return flint_calloc(1, sizeof(enum_struct));
}

static void
enum_clear(void * data)
{
    enum_struct * E = data;
    flint_free(E->C);
    flint_free(E);
}

static int
enum_compute(void * data, const _acb_mzeta_shells_struct * T, ulong len)
{
    enum_struct * E = data;
    flint_free(E->C);
    E->C = flint_calloc(len * T->N, sizeof(ulong));
    E->len = len;
    return counts_enum(E->C, T, len);
}

static void
enum_block(ulong * out, const void * data, const _acb_mzeta_shells_struct * T, ulong M0, ulong n)
{
    _acb_mzeta_shells_table_block(out, ((const enum_struct *) data)->C, T->N, M0, n);
}

static ulong
enum_max_count(const void * data, const _acb_mzeta_shells_struct * T)
{
    const enum_struct * E = data;
    return _acb_mzeta_shells_table_max(E->C, E->len, T->N);
}

static void
enum_residues(char * mask, const void * data, const _acb_mzeta_shells_struct * T)
{
    const enum_struct * E = data;
    _acb_mzeta_shells_table_residues(mask, E->C, E->len, T->N);
}

const _acb_mzeta_shells_handler_struct _acb_mzeta_shells_enum =
{
    "enumeration", enum_init, enum_clear, enum_cost, enum_compute, enum_block,
    enum_max_count, enum_residues
};
