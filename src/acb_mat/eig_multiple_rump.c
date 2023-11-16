/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

static void
acb_approx_mag(mag_t res, const acb_t x)
{
    mag_t t;
    mag_init(t);
    arf_get_mag(res, arb_midref(acb_realref(x)));
    arf_get_mag(t, arb_midref(acb_imagref(x)));
    mag_hypot(res, res, t);
    mag_clear(t);
}

static int
close(const acb_t x, const acb_t y, const mag_t eps)
{
    arf_t t;
    mag_t a, b;
    int result;

    mag_init(a);
    mag_init(b);
    arf_init(t);
    arf_sub(t, arb_midref(acb_realref(x)), arb_midref(acb_realref(y)), MAG_BITS, ARF_RND_UP);
    arf_get_mag(a, t);
    arf_sub(t, arb_midref(acb_imagref(x)), arb_midref(acb_imagref(y)), MAG_BITS, ARF_RND_UP);
    arf_get_mag(b, t);

    mag_hypot(a, a, b);
    result = (mag_cmp(a, eps) <= 0);

    mag_clear(a);
    mag_clear(b);
    arf_clear(t);

    return result;
}

int
acb_mat_eig_multiple_rump(acb_ptr E, const acb_mat_t A, acb_srcptr E_approx, const acb_mat_t R_approx, slong prec)
{
    slong c, i, j, k, n;
    acb_mat_t X;
    acb_ptr F;
    int result;
    slong iter;
    mag_t escale, eps, tm, um;
    slong ** cluster;
    slong * cluster_size;
    slong num_clusters;

    n = acb_mat_nrows(A);

    if (n == 0)
        return 1;

    cluster = flint_malloc(sizeof(slong *) * n);
    for (i = 0; i < n; i++)
        cluster[i] = flint_malloc(sizeof(slong) * n);
    cluster_size = flint_malloc(sizeof(slong) * n);

    mag_init(eps);
    mag_init(escale);
    mag_init(tm);
    mag_init(um);

    mag_zero(escale);
    for (i = 0; i < n; i++)
    {
        acb_approx_mag(tm, E_approx + i);
        mag_max(escale, escale, tm);
    }

    /* todo: when num_clusters = 1, could use fallback global enclosure */

    /* todo: initial clustering could be allowed to be zero */
    /* M 2^(-0.75p), M 2^(-0.5p), M 2^(-0.25p), ... */
    for (iter = 0; iter < 2; iter++)
    {
        mag_mul_2exp_si(eps, escale, -prec + (iter + 1) * prec/4);

        /* Group the eigenvalue approximations. */
        num_clusters = 0;
        for (i = 0; i < n; i++)
        {
            int new_cluster = 1;

            for (j = 0; j < num_clusters && new_cluster; j++)
            {
                if (close(E_approx + i, E_approx + cluster[j][0], eps))
                {
                    cluster[j][cluster_size[j]] = i;
                    cluster_size[j]++;
                    new_cluster = 0;
                }
            }

            if (new_cluster)
            {
                cluster[num_clusters][0] = i;
                cluster_size[num_clusters] = 1;
                num_clusters++;
            }
        }

        result = 1;

        F = _acb_vec_init(num_clusters);

        for (c = 0; c < num_clusters && result; c++)
        {
            k = cluster_size[c];

            acb_mat_init(X, n, k);

            for (i = 0; i < n; i++)
                for (j = 0; j < k; j++)
                    acb_set(acb_mat_entry(X, i, j), acb_mat_entry(R_approx, i, cluster[c][j]));

            acb_mat_eig_enclosure_rump(F + c, NULL, X, A, E_approx + cluster[c][0], X, prec);

            if (!acb_is_finite(F + c))
                result = 0;

            acb_mat_clear(X);
        }

        for (i = 0; i < num_clusters; i++)
        {
            for (j = i + 1; j < num_clusters; j++)
            {
                if (acb_overlaps(F + i, F + j))
                    result = 0;
            }
        }

        if (result)
        {
            i = 0;
            for (c = 0; c < num_clusters; c++)
            {
                for (j = 0; j < cluster_size[c]; j++)
                {
                    acb_set(E + i, F + c);
                    i++;
                }
            }
        }

        _acb_vec_clear(F, num_clusters);

        if (result)
            break;
    }

    if (!result)
        _acb_vec_indeterminate(E, n);

    for (i = 0; i < n; i++)
        flint_free(cluster[i]);
    flint_free(cluster);
    flint_free(cluster_size);

    mag_clear(eps);
    mag_clear(escale);
    mag_clear(tm);
    mag_clear(um);

    return result;
}
