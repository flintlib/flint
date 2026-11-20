/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Code (claude-opus-5-5) */

#include <stdlib.h>
#include <string.h>
#include <flint/acb_types.h>
#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/acb_ode.h>
#include <flint/gr.h>
#include <flint/gr_ore_poly.h>
#include <flint/gr_vec.h>

static void
usage(void)
{
    flint_printf(
        "usage: ode_transition_matrix [options] dop pt0 pt1 ...\n"
        "\n"
        "Compute the transition matrix of the differential operator dop\n"
        "(a polynomial in x and Dx = d/dx) along the piecewise linear path\n"
        "pt0 -> pt1 -> ...\n"
        "\n"
        "options (may also be written with two dashes):\n"
        "  -prec n         target precision in bits (default 64)\n"
        "  -digits n       number of digits to print (default: from prec)\n"
        "  -coeffs ring    constants of dop: zz (default), qq, qqbar, ca, acb\n"
        "  -points ring    path vertices: zz, qq, qqbar, ca, acb (default)\n"
        "\n"
        "Products are parsed commutatively: write polynomial factors to the\n"
        "left of powers of Dx.\n"
        "\n"
        "example: ode_transition_matrix \"(x^2 + 1)*Dx^2 + 2*x*Dx\" 0 \"1/2 + 1/4*I\" 1 2\n");
}

static int
init_ring(gr_ctx_t ctx, const char * name, slong prec)
{
    if (!strcmp(name, "zz"))
        gr_ctx_init_fmpz(ctx);
    else if (!strcmp(name, "qq"))
        gr_ctx_init_fmpq(ctx);
    else if (!strcmp(name, "qqbar"))
        gr_ctx_init_complex_qqbar(ctx);
    else if (!strcmp(name, "ca"))
        gr_ctx_init_complex_ca(ctx);
    else if (!strcmp(name, "acb"))
        gr_ctx_init_complex_acb(ctx, prec);
    else
        return 0;
    return 1;
}

static void
print_status(int status)
{
    if (status == GR_SUCCESS)
        flint_printf("GR_SUCCESS");
    else if (status == GR_DOMAIN)
        flint_printf("GR_DOMAIN");
    else if (status == GR_UNABLE)
        flint_printf("GR_UNABLE");
    else if (status == (GR_DOMAIN | GR_UNABLE))
        flint_printf("GR_DOMAIN | GR_UNABLE");
    else
        flint_printf("%d", status);
}

static void
print_mat(const acb_mat_t mat, slong digits)
{
    for (slong i = 0; i < acb_mat_nrows(mat); i++)
    {
        flint_printf("[");
        for (slong j = 0; j < acb_mat_ncols(mat); j++)
        {
            if (j > 0)
                flint_printf(", ");
            acb_printn(acb_mat_entry(mat, i, j), digits, 0);
        }
        flint_printf("]\n");
    }
}

int
main(int argc, char * argv[])
{
    gr_ctx_t Cst, Pol, Dop, Pts;
    gr_ptr dop;
    gr_vec_t path;
    acb_mat_t mat;
    const char * coeffs = "zz";
    const char * points = "acb";
    slong prec = 64, digits = -1;
    int i, status = GR_SUCCESS;

    for (i = 1; i + 1 < argc; i += 2)
    {
        const char * opt = argv[i];

        if (opt[0] == '-' && opt[1] == '-')
            opt++;

        if (!strcmp(opt, "-prec"))
            prec = atol(argv[i + 1]);
        else if (!strcmp(opt, "-digits"))
            digits = atol(argv[i + 1]);
        else if (!strcmp(opt, "-coeffs"))
            coeffs = argv[i + 1];
        else if (!strcmp(opt, "-points"))
            points = argv[i + 1];
        else
            break;
    }

    if (argc - i < 2 || prec < 2)
    {
        usage();
        return 1;
    }

    if (digits < 0)
        digits = prec * 0.30103 + 1;

    if (!init_ring(Cst, coeffs, prec + 10))
    {
        flint_printf("unsupported ring for -coeffs: %s\n", coeffs);
        return 1;
    }
    if (!init_ring(Pts, points, prec + 10))
    {
        flint_printf("unsupported ring for -points: %s\n", points);
        gr_ctx_clear(Cst);
        return 1;
    }

    gr_ctx_init_gr_poly(Pol, Cst);
    gr_ore_poly_ctx_init(Dop, Pol, 0, ORE_ALGEBRA_DERIVATIVE);
    status |= gr_ctx_set_gen_name(Pol, "x");
    status |= gr_ctx_set_gen_name(Dop, "Dx");

    GR_TMP_INIT(dop, Dop);
    gr_vec_init(path, argc - i - 1, Pts);
    acb_mat_init(mat, 0, 0);

    if (gr_set_str(dop, argv[i], Dop) != GR_SUCCESS)
    {
        flint_printf("unable to parse the operator: %s\n", argv[i]);
        status = GR_UNABLE;
        goto cleanup;
    }

    for (slong j = 0; j < path->length; j++)
    {
        if (gr_set_str(gr_vec_entry_ptr(path, j, Pts), argv[i + 1 + j], Pts)
                != GR_SUCCESS)
        {
            flint_printf("unable to parse the point: %s\n", argv[i + 1 + j]);
            status = GR_UNABLE;
            goto cleanup;
        }
    }

    slong order = gr_ore_poly_length(dop, Dop) - 1;
    if (order < 0)
    {
        flint_printf("the operator must be nonzero\n");
        status = GR_DOMAIN;
        goto cleanup;
    }

    flint_printf("dop = %{gr}\n", dop, Dop);
    flint_printf("path = ");
    status |= gr_vec_print(path, Pts);
    flint_printf("\n\n");

    acb_mat_clear(mat);
    acb_mat_init(mat, order, order);

    status |= acb_ode_transition_matrix_fixed(mat, dop, Dop, path, Pts,
                                              NULL, NULL, 0, prec);

    if (status == GR_SUCCESS)
        print_mat(mat, digits);
    else
    {
        flint_printf("failed (status = ");
        print_status(status);
        flint_printf(")\n");
    }

cleanup:
    acb_mat_clear(mat);
    gr_vec_clear(path, Pts);
    GR_TMP_CLEAR(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(Pts);
    gr_ctx_clear(Cst);

    flint_cleanup_master();

    return (status == GR_SUCCESS) ? 0 : 1;
}
