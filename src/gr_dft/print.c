/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "gr.h"
#include "gr_dft.h"

/* Diagnostic printer: writes a description of the plan (algorithm,
   decomposition, table state, threading configuration) with one
   indented block per sub-plan. */

static const char *
_alg_name(int alg)
{
    switch (alg)
    {
        case GR_DFT_ALG_AUTO:       return "auto";
        case GR_DFT_ALG_NAIVE:      return "naive";
        case GR_DFT_ALG_CT:         return "ct";
        case GR_DFT_ALG_BAILEY:     return "bailey";
        case GR_DFT_ALG_SPLIT:      return "split";
        case GR_DFT_ALG_PFA:        return "pfa";
        case GR_DFT_ALG_MIXED:      return "mixed";
        case GR_DFT_ALG_BLUESTEIN:  return "bluestein";
        default:                    return "?";
    }
}

static void
_print_rec(FILE * out, const gr_dft_pre_t P, int indent, const char * label)
{
    slong i;

    for (i = 0; i < indent; i++)
        flint_fprintf(out, "  ");
    flint_fprintf(out, "%s: n = %wu, alg = %s", label, P->n, _alg_name(P->alg));

    if (P->depth >= 0)
        flint_fprintf(out, ", depth = %d", P->depth);
    if (P->flags & GR_DFT_SCRAMBLED)
        flint_fprintf(out, ", scrambled");
    if (P->ctx == NULL)
        flint_fprintf(out, ", layout only");
    flint_fprintf(out, "\n");

    for (i = 0; i < indent; i++)
        flint_fprintf(out, "  ");
    flint_fprintf(out, "  tables: roots %s", (P->roots != NULL) ? "yes" :
            (P->stage_tab != NULL) ? "packed stages" : "no");
    if (P->wclass != NULL)
        flint_fprintf(out, ", complex mode%s",
                (P->wtab != NULL) ? " (karatsuba)" : "");
    if (P->nfixed_root_err > 0.0)
        flint_fprintf(out, ", root err %.3g ulp%s", P->nfixed_root_err,
                (P->roots == NULL) ? " (est.)" : "");
    if (P->bl_kern != NULL)
        flint_fprintf(out, ", kernel (conv_len %wu%s)", P->conv_len,
                P->bl_shifted ? ", shifted" : "");
    flint_fprintf(out, "\n");

    if (P->alg == GR_DFT_ALG_BAILEY || P->alg == GR_DFT_ALG_PFA)
    {
        for (i = 0; i < indent; i++)
            flint_fprintf(out, "  ");
        flint_fprintf(out, "  n1 = %wu, n2 = %wu", P->n1, P->n2);
        if (P->alg == GR_DFT_ALG_PFA)
            flint_fprintf(out, ", crt coefficients %wu, %wu", P->pfa_a, P->pfa_b);
        flint_fprintf(out, "\n");
    }

    if (P->radices != NULL)
    {
        for (i = 0; i < indent; i++)
            flint_fprintf(out, "  ");
        flint_fprintf(out, "  radices:");
        for (i = 0; i < P->num_radices; i++)
            flint_fprintf(out, " %wu", P->radices[i]);
        flint_fprintf(out, "\n");
    }

    if (P->num_threads > 0 || P->serial_block > 0)
    {
        for (i = 0; i < indent; i++)
            flint_fprintf(out, "  ");
        flint_fprintf(out, "  threading: %wd attached handles, serial block %wd\n",
                P->num_threads, (P->serial_block > 0) ? P->serial_block :
                (slong) GR_DFT_SERIAL_BLOCK_DEFAULT);
    }

    if (P->P1 != NULL)
        _print_rec(out, P->P1,
                indent + 1,
                (P->alg == GR_DFT_ALG_MIXED ||
                 P->alg == GR_DFT_ALG_BLUESTEIN) ? "child" : "sub-plan 1");
    if (P->P2 != NULL)
        _print_rec(out, P->P2, indent + 1, "sub-plan 2");
}

void
gr_dft_precomp_fprint(FILE * out, const gr_dft_pre_t P)
{
    _print_rec(out, P, 0, "plan");
}

void
gr_dft_precomp_print(const gr_dft_pre_t P)
{
    gr_dft_precomp_fprint(stdout, P);
}
