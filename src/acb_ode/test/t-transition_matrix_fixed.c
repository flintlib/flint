/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_types.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_ode.h"
#include "gr.h"
#include "gr_ore_poly.h"
#include "gr_poly.h"
#include "gr_vec.h"

/* generated using Claude Code (claude-opus-5-5) */

/* The leading coefficient of the operator is a multiple of
   (x - 8)*(x^2 + 64), so that its singular points are 8, 8i, -8i, and, for
   Euler operators, 0. Paths are walks on the vertices
   c + oct_re[k] + oct_im[k]*i of an octagon around a center c (encoded as -1)
   chosen among these four points. */
static const slong lc_coeffs[4] = { -512, 64, -8, 1 };
static const slong ctr_re[4] = { 0, 8, 0,  0 };
static const slong ctr_im[4] = { 0, 0, 8, -8 };
static const slong oct_re[8] = { 3, 2, 0, -2, -3, -2,  0,  2 };
static const slong oct_im[8] = { 0, 2, 3,  2,  0, -2, -3, -2 };
#define CENTER -1
#define MAX_PATH_LEN 16

/* a vertex, restricted to the real axis if real is set */
static slong
random_vertex(flint_rand_t state, int real)
{
    return real ? 4 * (slong) n_randint(state, 2) : (slong) n_randint(state, 8);
}

/* A few random moves (between the center and any vertex, or between adjacent
   vertices) followed by the shortest route to end, possibly with a detour
   through the center or around the octagon. Returns the length of the walk. */
static slong
random_walk(slong * vtx, slong start, slong end, int real, flint_rand_t state)
{
    slong n = 0, cur = start;
    slong nsteps = n_randint(state, 5);

    vtx[n++] = cur;
    for (slong i = 0; i < nsteps; i++)
    {
        if (cur == CENTER)
            cur = random_vertex(state, real);
        else if (real || !n_randint(state, 3))
            cur = CENTER;
        else
            cur = (cur + (n_randint(state, 2) ? 1 : 7)) % 8;
        vtx[n++] = cur;
    }

    if (cur != CENTER && (end == CENTER || real || n_randint(state, 2)))
        vtx[n++] = cur = CENTER;

    if (cur == CENTER && end != CENTER)
        vtx[n++] = cur = end;

    slong dir = n_randint(state, 2) ? 1 : 7;
    while (cur != end)
        vtx[n++] = cur = (cur + dir) % 8;

    FLINT_ASSERT(n <= MAX_PATH_LEN);
    return n;
}

/* With the principal branch of the logarithm at the center (used for the
   local basis there on both sides of a visit to the center), the homotopy
   class of a walk around a singular center is determined by the number of
   times it crosses the branch cut, i.e., the edge between vertices 4 and 5
   (vertex 4, with argument pi, lies on the cut). */
static slong
winding_number(const slong * vtx, slong len)
{
    slong w = 0;
    for (slong i = 0; i + 1 < len; i++)
        w += (vtx[i] == 4 && vtx[i + 1] == 5) - (vtx[i] == 5 && vtx[i + 1] == 4);
    return w;
}

static int
gaussian_integer(gr_ptr res, slong re, slong im, gr_ctx_t ctx, gr_ctx_t QQbar)
{
    gr_ptr t;
    int status = GR_SUCCESS;
    GR_TMP_INIT(t, QQbar);
    status |= gr_i(t, QQbar);
    status |= gr_mul_si(t, t, im, QQbar);
    status |= gr_add_si(t, t, re, QQbar);
    status |= gr_set_other(res, t, QQbar, ctx);
    GR_TMP_CLEAR(t, QQbar);
    return status;
}

TEST_GR_FUNCTION_START(acb_ode_transition_matrix_fixed, state, count_success, count_domain, count_unable)
{
    gr_ctx_t QQbar;
    gr_ctx_init_complex_qqbar(QQbar);

    for (slong iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        ulong seed1, seed2;
        flint_rand_get_seed(&seed1, &seed2, state);

        int status = GR_SUCCESS;
        slong prec = 2 + 50 * n_randint(state, 3) + n_randint(state, 30);

        /* operator over ZZ (usually), CC, or a random ring */

        gr_ctx_t Cst, Pol, Dop, _Pts;
        gr_ctx_struct * Pts;
        slong rnd = n_randint(state, 8);
        if (rnd == 0)
            gr_ctx_init_random(Cst, state);
        else if (rnd == 1)
            gr_ctx_init_complex_acb(Cst, prec + n_randint(state, 20));
        else
            gr_ctx_init_fmpz(Cst);
        gr_ctx_init_gr_poly(Pol, Cst);
        int euler = n_randint(state, 2);
        gr_ore_poly_ctx_init(Dop, Pol, 0, euler ? ORE_ALGEBRA_EULER_DERIVATIVE
                                                : ORE_ALGEBRA_DERIVATIVE);

        gr_ore_poly_t dop;
        gr_ore_poly_init(dop, Dop);
        slong dop_len = 2 + n_randint(state, 3);
        slong order = dop_len - 1;
        gr_ore_poly_fit_length(dop, dop_len, Dop);
        _gr_ore_poly_set_length(dop, dop_len, Dop);

        /* small coefficients to avoid huge exponents (slow) */
        for (slong i = 0; i < order; i++)
            for (slong j = 0; j < 3; j++)
                status |= gr_poly_set_coeff_si(
                        gr_ore_poly_coeff_ptr(dop, i, Dop), j,
                        (slong) n_randint(state, 17) - 8, Cst);

        slong mul = 1 + n_randint(state, 16);
        for (slong j = 0; j < 4; j++)
            status |= gr_poly_set_coeff_si(
                    gr_ore_poly_coeff_ptr(dop, order, Dop), j,
                    mul * lc_coeffs[j], Cst);
        _gr_ore_poly_normalise(dop, Dop);

        slong ctr = n_randint(state, 4);
        slong cre = ctr_re[ctr], cim = ctr_im[ctr];
        int ordinary = (ctr == 0 && !euler);

        /* when the center is the origin, fix the order of the local basis
           there to allow comparing with acb_ode_fundamental_matrix */
        gr_ctx_t CC;
        acb_ode_exponents_t expos;
        gr_ctx_init_complex_acb(CC, prec);
        acb_ode_exponents_init(expos);
        int have_expos = (cre == 0 && cim == 0 && status == GR_SUCCESS
                          && acb_ode_exponents(expos, dop, Dop, CC) == GR_SUCCESS);

        /* path points in the constant ring, CC, QQbar, or a random ring */

        rnd = n_randint(state, 16);
        Pts = _Pts;
        if (rnd == 0)
            gr_ctx_init_random(Pts, state);
        else if (rnd <= 2)
            gr_ctx_init_complex_qqbar(Pts);
        else if (rnd <= 6 || cim != 0)
            gr_ctx_init_complex_acb(Pts, 2 + n_randint(state, 200));
        else
            Pts = Cst;

        int real = (cim == 0)
                   && (Pts->which_ring == GR_CTX_FMPZ || !n_randint(state, 8));

        /* walks */

        slong start = n_randint(state, 3) ? random_vertex(state, real) : CENTER;
        slong end = n_randint(state, 3) ? random_vertex(state, real) : CENTER;

        slong vtx[2][MAX_PATH_LEN], len[2], wind[2];
        gr_vec_t path[2];
        acb_mat_t mat[2], ref;

        for (slong p = 0; p < 2; p++)
        {
            slong n = len[p] = random_walk(vtx[p], start, end, real, state);
            wind[p] = winding_number(vtx[p], n);

            gr_vec_init(path[p], n, Pts);
            for (slong i = 0; i < n; i++)
            {
                slong k = vtx[p][i];
                slong re = cre + ((k == CENTER) ? 0 : oct_re[k]);
                slong im = cim + ((k == CENTER) ? 0 : oct_im[k]);
                status |= gaussian_integer(gr_vec_entry_ptr(path[p], i, Pts),
                                           re, im, Pts, QQbar);
            }

            acb_mat_init(mat[p], order, order);
            acb_mat_randtest(mat[p], state, prec, 3);
        }
        acb_mat_init(ref, order, order);

        /* test */

        if (status != GR_SUCCESS)
            goto end;

        for (slong p = 0; p < 2; p++)
            status |= acb_ode_transition_matrix_fixed(
                    mat[p], dop, Dop, path[p], Pts,
                    (have_expos && start == CENTER) ? expos : NULL,
                    (have_expos && end == CENTER) ? expos : NULL,
                    0, prec);

        if (status != GR_SUCCESS)
            goto end;

        if ((ordinary || wind[0] == wind[1])
            && !acb_mat_overlaps(mat[0], mat[1]))
        {
            flint_printf("FAIL: homotopic paths\n\n");
            goto fail;
        }

        /* single steps starting at the origin */
        for (slong p = 0; p < 2; p++)
        {
            acb_t pt;
            slong k = vtx[p][1];

            if (!have_expos || len[p] != 2 || vtx[p][0] != CENTER)
                continue;

            acb_init(pt);
            acb_set_si_si(pt, oct_re[k], oct_im[k]);
            status |= acb_ode_fundamental_matrix(ref, dop, Dop, expos, NULL, pt,
                                                 ACB_ODE_BASIS_ECHELON, 0, prec);
            acb_clear(pt);

            if (status == GR_SUCCESS && !acb_mat_overlaps(mat[p], ref))
            {
                flint_printf("FAIL: single step\n\n");
                flint_printf("ref = %{acb_mat}\n", ref);
                goto fail;
            }
        }

        goto end;

fail:
        flint_printf("flint_rand_set_seed(state, 0x%wx, 0x%wx);\n", seed1, seed2);
        flint_printf("Dop = %{gr_ctx}\n", Dop);
        flint_printf("dop = %{gr}\n", dop, Dop);
        flint_printf("Pts = %{gr_ctx}\n", Pts);
        flint_printf("ordinary = %d, prec = %wd\n", ordinary, prec);
        for (slong p = 0; p < 2; p++)
        {
            flint_printf("path%wd = ", p);
            GR_MUST_SUCCEED(gr_vec_print(path[p], Pts));
            flint_printf("\nwinding%wd = %wd\n", p, wind[p]);
            flint_printf("mat%wd = %{acb_mat}\n", p, mat[p]);
        }
        flint_abort();

end:
        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        acb_mat_clear(ref);
        for (slong p = 0; p < 2; p++)
        {
            acb_mat_clear(mat[p]);
            gr_vec_clear(path[p], Pts);
        }
        gr_ore_poly_clear(dop, Dop);
        acb_ode_exponents_clear(expos);
        gr_ctx_clear(CC);
        if (Pts != Cst)
            gr_ctx_clear(Pts);
        gr_ctx_clear(Dop);
        gr_ctx_clear(Pol);
        gr_ctx_clear(Cst);
    }

    gr_ctx_clear(QQbar);

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
