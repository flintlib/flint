/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Lazy fields: printing, symbolic expressions and parsing. */

/* (for recursive mutexes: pthread_mutexattr_settype) */
#define _GNU_SOURCE

#include "gr_tower/lazy_impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

/* the number of decimal digits that separate the generator g from the
   other roots of its origin polynomial (0 if there is none, or if the
   roots cannot be isolated at a moderate precision) */
static slong
_root_separation_digits(const gr_tower_gen_struct * g, gr_tower_t T)
{
    acb_ptr roots;
    acb_t z;
    arb_t d, dmin;
    slong n, i, prec, digits = 0;

    if (g->origin == NULL || g->origin->length < 3)
        return 0;
    n = g->origin->length - 1;
    roots = _acb_vec_init(n);
    acb_init(z);
    arb_init(d);
    arb_init(dmin);
    for (prec = 64; prec <= 1024 && digits == 0; prec *= 2)
    {
        slong found = 0;
        arb_fmpz_poly_complex_roots(roots, g->origin, 0, prec);
        if (gr_tower_step_get_acb(z, T, g->index, prec) != GR_SUCCESS)
            break;
        arb_indeterminate(dmin);
        for (i = 0; i < n; i++)
        {
            if (acb_overlaps(roots + i, z))
            {
                found++;
                continue;
            }
            acb_sub(roots + i, roots + i, z, prec);
            acb_abs(d, roots + i, prec);
            if (!arb_is_finite(dmin) || arb_lt(d, dmin))
                arb_set(dmin, d);
        }
        if (found != 1 || !arb_is_finite(dmin) || !arb_is_positive(dmin))
            continue;
        /* digits with 10^-digits |z| < dmin / 4 */
        {
            double lz = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR), li;
            double lmin = arf_get_d(arb_midref(dmin), ARF_RND_NEAR);
            li = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
            lz = sqrt(lz * lz + li * li);
            if (lz == 0 || lmin == 0)
                break;
            digits = (slong) ceil(log10(4.0 * lz / lmin)) + 1;
            if (digits < 1)
                digits = 1;
        }
    }
    _acb_vec_clear(roots, n);
    acb_clear(z);
    arb_clear(d);
    arb_clear(dmin);
    return digits;
}

/* the definition of a generator, as an expression in the other generators */
static int
_write_gen_def(gr_stream_t out, const gr_tower_gen_struct * g, gr_tower_t T, slong digits)
{
    int status = GR_SUCCESS;

    if (g->def_kind == GR_TOWER_PI)
        return gr_stream_write(out, "pi");

    if (g->def_kind == GR_TOWER_ROOT_OF_UNITY)
    {
        status |= gr_stream_write(out, "exp(2*pi*i/");
        status |= gr_stream_write_si(out, g->def_param);
        status |= gr_stream_write(out, ")");
        return status;
    }

    if (g->def_kind == GR_TOWER_TAN_PI)
    {
        status |= gr_stream_write(out, "tan(pi/");
        status |= gr_stream_write_si(out, g->def_param);
        status |= gr_stream_write(out, ")");
        return status;
    }

    if (GR_TOWER_KIND_IS_SPECIAL(g->def_kind))
    {
        char * s = _gr_tower_gen_args_str(g, T);
        status |= _gr_tower_special_write(out, g->def_kind, g->def_param, (s != NULL) ? s : "");
        flint_free(s);
        return status;
    }

    if ((g->def_kind == GR_TOWER_EXP || g->def_kind == GR_TOWER_LOG || g->def_kind == GR_TOWER_ROOT ||
         g->def_kind == GR_TOWER_TAN || g->def_kind == GR_TOWER_ATAN) && g->arg.mctx != NULL)
    {
        if (g->def_kind == GR_TOWER_TAN)
            status |= gr_stream_write(out, "tan(");
        else if (g->def_kind == GR_TOWER_ATAN)
            status |= gr_stream_write(out, "atan(");
        else if (g->def_kind == GR_TOWER_EXP)
            status |= gr_stream_write(out, "exp(");
        else if (g->def_kind == GR_TOWER_LOG)
            status |= gr_stream_write(out, "log(");
        else if (g->def_param == 2)
            status |= gr_stream_write(out, "sqrt(");
        else
            status |= gr_stream_write(out, "root(");
        status |= gr_stream_write_free(out, _gr_tower_flat_get_str(&g->arg.data, g->arg.mctx, T));
        if (g->def_kind == GR_TOWER_ROOT && g->def_param != 2)
        {
            status |= gr_stream_write(out, ", ");
            status |= gr_stream_write_si(out, g->def_param);
        }
        status |= gr_stream_write(out, ")");
        return status;
    }

    if (g->kind == GR_TOWER_ALGEBRAIC)
    {
        /* a root of its minimal polynomial, identified numerically: with
           enough digits to tell it from the other roots of its origin
           polynomial (two roots of x^20 - 2 (101 x - 1)^2 agree to 22
           digits), when that polynomial is known */
        slong k = g->index;
        char * ps;
        acb_t z;

        status |= gr_stream_write(out, "root(");
        if (gr_poly_get_str(&ps, gr_tower_step_minpoly(T, k), g->name, gr_tower_field_at(T, k - 1)) == GR_SUCCESS)
            status |= gr_stream_write_free(out, ps);
        else
            status |= gr_stream_write(out, "?");
        status |= gr_stream_write(out, ", ");
        acb_init(z);
        digits = FLINT_MAX(digits, _root_separation_digits(g, T));
        if (gr_tower_step_get_acb(z, T, k, (slong) (digits * 3.33) + 30) != GR_SUCCESS)
            acb_set(z, &g->enclosure);
        status |= gr_stream_write_free(out, arb_get_str(acb_realref(z), digits, ARB_STR_NO_RADIUS));
        /* (the enclosure of a real root may have an imaginary part
           containing zero; a nonreal root is separated from its
           conjugate at this precision) */
        if (!arb_contains_zero(acb_imagref(z)))
        {
            status |= gr_stream_write(out, (arf_sgn(arb_midref(acb_imagref(z))) < 0) ? " - " : " + ");
            arb_abs(acb_imagref(z), acb_imagref(z));
            status |= gr_stream_write_free(out, arb_get_str(acb_imagref(z), digits, ARB_STR_NO_RADIUS));
            status |= gr_stream_write(out, "*i");
        }
        acb_clear(z);
        status |= gr_stream_write(out, ")");
        return status;
    }

    return gr_stream_write(out, "?");
}

/* a numerical value with the given number of digits */
static int
_write_numeric(gr_stream_t out, const gr_tower_lazy_elem_t x, slong digits, gr_ctx_t ctx)
{
    acb_t z;
    slong prec = digits * 3.33 + 30;
    int status = GR_SUCCESS;

    acb_init(z);
    if (_gr_tower_lazy_get_acb_impl(z, x, prec, ctx) != GR_SUCCESS)
    {
        status = gr_stream_write(out, "?");
    }
    else
    {
        int show_imag = !arb_is_zero(acb_imagref(z)) &&
            !(arb_contains_zero(acb_imagref(z)) && mag_cmp_2exp_si(arb_radref(acb_imagref(z)), -prec / 2) < 0);
        /* (a real part which is zero to within the working precision
           is not shown either, like the imaginary part above) */
        int show_real = !show_imag || !(arb_is_zero(acb_realref(z)) ||
            (arb_contains_zero(acb_realref(z)) && mag_cmp_2exp_si(arb_radref(acb_realref(z)), -prec / 2) < 0));

        if (show_real)
            status |= gr_stream_write_free(out, arb_get_str(acb_realref(z), digits, ARB_STR_NO_RADIUS));
        if (show_imag)
        {
            if (show_real)
                status |= gr_stream_write(out, (arf_sgn(arb_midref(acb_imagref(z))) < 0) ? " - " : " + ");
            else if (arf_sgn(arb_midref(acb_imagref(z))) < 0)
                status |= gr_stream_write(out, "-");
            arb_abs(acb_imagref(z), acb_imagref(z));
            status |= gr_stream_write_free(out, arb_get_str(acb_imagref(z), digits, ARB_STR_NO_RADIUS));
            status |= gr_stream_write(out, "*i");
        }
    }
    acb_clear(z);
    return status;
}

int
_gr_tower_lazy_write(gr_stream_t out, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_struct * T = x->F->T;
    int flags = L->options[GR_TOWER_OPT_PRINT_FLAGS];
    int rational, status = GR_SUCCESS;
    int * mark = NULL;
    slong d, ndefs = 0;

    x = _gr_tower_lazy_flat_view(x);

    /* generators created inside the tower (by the relation search) are
       named here at the latest: their default names are those of the
       tower, which may clash with the names in the context */
    _gr_tower_lazy_assign_def_ids(x->F, ctx);

    rational = fmpz_mpoly_q_is_fmpq(&x->elem.flat.data, x->elem.flat.mctx);

    if (rational && (flags & GR_TOWER_PRINT_SYMBOLIC))
        flags = GR_TOWER_PRINT_SYMBOLIC;
    if (!(flags & (GR_TOWER_PRINT_NUMERIC | GR_TOWER_PRINT_SYMBOLIC)))
        flags |= GR_TOWER_PRINT_SYMBOLIC;

    if ((flags & GR_TOWER_PRINT_DEFS) && (flags & GR_TOWER_PRINT_SYMBOLIC))
    {
        mark = flint_calloc(FLINT_MAX(T->num_gens, 1), sizeof(int));
        _gr_tower_lazy_involved_gens(mark, x, T);
        for (d = 0; d < T->num_gens; d++)
            if (mark[d] && !_gr_tower_lazy_gen_is_constant_symbol(GR_TOWER_GEN(T, d)))
                ndefs++;
    }

    if (flags & GR_TOWER_PRINT_NUMERIC)
    {
        status |= _write_numeric(out, x, L->options[GR_TOWER_OPT_PRINT_DIGITS], ctx);
        if (flags & GR_TOWER_PRINT_SYMBOLIC)
            status |= gr_stream_write(out, " {");
    }

    if (flags & GR_TOWER_PRINT_SYMBOLIC)
    {
        status |= gr_stream_write_free(out, _gr_tower_flat_get_str(&x->elem.flat.data, x->elem.flat.mctx, T));

        if (ndefs > 0)
        {
            int first = 1;
            slong n, i;
            slong * order = _gr_tower_lazy_marked_by_creation(&n, mark, T);
            status |= gr_stream_write(out, (flags & GR_TOWER_PRINT_NUMERIC) ? " where " : " {");
            for (i = 0; i < n; i++)
            {
                const gr_tower_gen_struct * g = GR_TOWER_GEN(T, order[i]);
                if (_gr_tower_lazy_gen_is_constant_symbol(g))
                    continue;
                if (!first)
                    status |= gr_stream_write(out, "; ");
                first = 0;
                status |= gr_stream_write(out, g->name);
                status |= gr_stream_write(out, " = ");
                status |= _write_gen_def(out, g, T, L->options[GR_TOWER_OPT_PRINT_DIGITS]);
            }
            if (!(flags & GR_TOWER_PRINT_NUMERIC))
                status |= gr_stream_write(out, "}");
            flint_free(order);
        }

        if (flags & GR_TOWER_PRINT_NUMERIC)
            status |= gr_stream_write(out, "}");
    }

    flint_free(mark);
    return status;
}

/* -------------------------------------------------------------------- */
/* symbolic expressions                                                  */
/* -------------------------------------------------------------------- */

/* the symbol of a generator */
static void
_gen_symbol(fexpr_t res, const gr_tower_gen_struct * g)
{
    if (g->def_kind == GR_TOWER_PI)
        fexpr_set_symbol_builtin(res, FEXPR_Pi);
    else if (g->def_kind == GR_TOWER_ROOT_OF_UNITY && g->def_param == 4)
        fexpr_set_symbol_builtin(res, FEXPR_NumberI);
    else
    {
        /* a1 -> a_1 (the LaTeX form of the symbol) */
        const char * name = g->name;
        slong i = 0, n = strlen(name);
        while (i < n && !(name[i] >= '0' && name[i] <= '9'))
            i++;
        if (i > 0 && i < n)
        {
            char * t = flint_malloc(n + 2);
            memcpy(t, name, i);
            t[i] = '_';
            memcpy(t + i + 1, name + i, n - i + 1);
            fexpr_set_symbol_str(res, t);
            flint_free(t);
        }
        else
            fexpr_set_symbol_str(res, name);
    }
}

/* a flat element in the context mctx of the tower T, in the symbols of
   the generators */
static void
_flat_get_fexpr(fexpr_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_struct * mctx, gr_tower_t T)
{
    slong cap, v;
    slong * var_gid;
    fexpr_vec_t vars;

    if (!_gr_tower_flat_find_layout(&cap, &var_gid, mctx, &T->flat))
    {
        cap = mctx->minfo->nvars;
        var_gid = NULL;
    }

    fexpr_vec_init(vars, cap);
    for (v = 0; v < cap; v++)
    {
        slong d = (var_gid == NULL || var_gid[v] < 0) ? -1 : gr_tower_gid_order(T, var_gid[v]);
        if (d >= 0)
            _gen_symbol(vars->entries + v, GR_TOWER_GEN(T, d));
        else
            fexpr_set_symbol_builtin(vars->entries + v, FEXPR_Unknown);
    }

    fexpr_set_fmpz_mpoly_q(res, x, vars, mctx);
    fexpr_vec_clear(vars);
}

/* the definition of a generator as an expression */
static void
_gen_def_get_fexpr(fexpr_t res, const gr_tower_gen_struct * g, gr_tower_t T, slong digits)
{
    fexpr_t t, u;

    fexpr_init(t);
    fexpr_init(u);

    if (g->def_kind == GR_TOWER_ROOT_OF_UNITY)
    {
        /* Exp(Div(Mul(2, Pi, I), n)) */
        fexpr_t two, pi, I;
        fexpr_init(two);
        fexpr_init(pi);
        fexpr_init(I);
        fexpr_set_ui(two, 2);
        fexpr_set_symbol_builtin(pi, FEXPR_Pi);
        fexpr_set_symbol_builtin(I, FEXPR_NumberI);
        fexpr_mul(t, two, pi);
        fexpr_mul(t, t, I);
        fexpr_set_si(u, g->def_param);
        fexpr_div(t, t, u);
        fexpr_call_builtin1(res, FEXPR_Exp, t);
        fexpr_clear(two);
        fexpr_clear(pi);
        fexpr_clear(I);
    }
    else if (g->def_kind == GR_TOWER_TAN_PI)
    {
        /* Tan(Div(Pi, n)) */
        fexpr_set_symbol_builtin(t, FEXPR_Pi);
        fexpr_set_si(u, g->def_param);
        fexpr_div(t, t, u);
        fexpr_call_builtin1(res, FEXPR_Tan, t);
    }
    else if (g->def_kind == GR_TOWER_CONSTANT)
    {
        fexpr_set_symbol_builtin(res, (g->def_param == GR_TOWER_CONST_CATALAN) ? FEXPR_CatalanConstant : FEXPR_Euler);
    }
    else if (g->def_kind == GR_TOWER_HYPGEOM && g->arg.mctx != NULL)
    {
        /* Hypergeometric0F1(b, z), ..., Hypergeometric3F2(a1, a2, a3, b1, b2, z) */
        slong p = GR_TOWER_HYPGEOM_P(g->def_param), q = GR_TOWER_HYPGEOM_Q(g->def_param), i, n = 1 + p + q;
        fexpr_ptr args;
        fexpr_vec_t v;
        ulong head;

        fexpr_vec_init(v, n);
        args = v->entries;
        for (i = 0; i < p + q; i++)
            _flat_get_fexpr(args + i, &g->xargs[i].data, g->xargs[i].mctx, T);
        _flat_get_fexpr(args + p + q, &g->arg.data, g->arg.mctx, T);

        head = (p == 0 && q == 1) ? FEXPR_Hypergeometric0F1 : (p == 1 && q == 1) ? FEXPR_Hypergeometric1F1 :
               (p == 1 && q == 2) ? FEXPR_Hypergeometric1F2 : (p == 2 && q == 0) ? FEXPR_Hypergeometric2F0 :
               (p == 2 && q == 1) ? FEXPR_Hypergeometric2F1 : (p == 2 && q == 2) ? FEXPR_Hypergeometric2F2 :
               (p == 3 && q == 2) ? FEXPR_Hypergeometric3F2 : FEXPR_Unknown;

        if (head != FEXPR_Unknown)
        {
            fexpr_t f;
            fexpr_init(f);
            fexpr_set_symbol_builtin(f, head);
            fexpr_call_vec(res, f, args, n);
            fexpr_clear(f);
        }
        else
            fexpr_set_symbol_builtin(res, FEXPR_Unknown);

        fexpr_vec_clear(v);
    }
    else if (g->def_kind == GR_TOWER_HURWITZ_ZETA && g->arg.mctx != NULL && g->num_xargs == 1)
    {
        _flat_get_fexpr(t, &g->arg.data, g->arg.mctx, T);
        _flat_get_fexpr(u, &g->xargs[0].data, g->xargs[0].mctx, T);
        fexpr_call_builtin2(res, FEXPR_HurwitzZeta, t, u);
    }
    else if (g->def_kind == GR_TOWER_JACOBI_THETA && g->arg.mctx != NULL && g->num_xargs == 1)
    {
        /* JacobiTheta(j, z, tau) */
        fexpr_t f, jj;
        fexpr_init(f);
        fexpr_init(jj);
        _flat_get_fexpr(t, &g->arg.data, g->arg.mctx, T);
        _flat_get_fexpr(u, &g->xargs[0].data, g->xargs[0].mctx, T);
        fexpr_set_si(jj, g->def_param);
        fexpr_set_symbol_builtin(f, FEXPR_JacobiTheta);
        fexpr_call3(res, f, jj, t, u);
        fexpr_clear(f);
        fexpr_clear(jj);
    }
    else if (GR_TOWER_KIND_IS_SPECIAL(g->def_kind) && g->arg.mctx != NULL)
    {
        _flat_get_fexpr(t, &g->arg.data, g->arg.mctx, T);
        fexpr_set_si(u, g->def_param);
        switch (g->def_kind)
        {
            case GR_TOWER_GAMMA: fexpr_call_builtin1(res, FEXPR_Gamma, t); break;
            case GR_TOWER_ERF: fexpr_call_builtin1(res, (g->def_param == 0) ? FEXPR_Erf : FEXPR_Erfi, t); break;
            case GR_TOWER_ZETA: fexpr_call_builtin1(res, FEXPR_RiemannZeta, t); break;
            case GR_TOWER_ELLIPTIC_K: fexpr_call_builtin1(res, FEXPR_EllipticK, t); break;
            case GR_TOWER_ELLIPTIC_E: fexpr_call_builtin1(res, FEXPR_EllipticE, t); break;
            case GR_TOWER_MODULAR_LAMBDA: fexpr_call_builtin1(res, FEXPR_ModularLambda, t); break;
            case GR_TOWER_LAMBERTW:
                if (g->def_param == 0)
                    fexpr_call_builtin1(res, FEXPR_LambertW, t);
                else
                    fexpr_call_builtin2(res, FEXPR_LambertW, t, u);
                break;
            case GR_TOWER_POLYGAMMA:
                if (g->def_param == 0)
                    fexpr_call_builtin1(res, FEXPR_DigammaFunction, t);
                else
                    fexpr_call_builtin2(res, FEXPR_DigammaFunction, t, u);
                break;
            case GR_TOWER_POLYLOG: fexpr_call_builtin2(res, FEXPR_PolyLog, u, t); break;
            case GR_TOWER_DIRICHLET_L:
                {
                    fexpr_t qq, kk, ch;
                    fexpr_init(qq);
                    fexpr_init(kk);
                    fexpr_init(ch);
                    fexpr_set_ui(qq, GR_TOWER_DIRICHLET_Q(g->def_param));
                    fexpr_set_ui(kk, GR_TOWER_DIRICHLET_K(g->def_param));
                    fexpr_call_builtin2(ch, FEXPR_DirichletCharacter, qq, kk);
                    fexpr_call_builtin2(res, FEXPR_DirichletL, t, ch);
                    fexpr_clear(qq);
                    fexpr_clear(kk);
                    fexpr_clear(ch);
                }
                break;
            default: fexpr_set_symbol_builtin(res, FEXPR_Unknown);
        }
    }
    else if ((g->def_kind == GR_TOWER_EXP || g->def_kind == GR_TOWER_LOG || g->def_kind == GR_TOWER_ROOT ||
              g->def_kind == GR_TOWER_TAN || g->def_kind == GR_TOWER_ATAN) && g->arg.mctx != NULL)
    {
        _flat_get_fexpr(t, &g->arg.data, g->arg.mctx, T);
        if (g->def_kind == GR_TOWER_TAN)
            fexpr_call_builtin1(res, FEXPR_Tan, t);
        else if (g->def_kind == GR_TOWER_ATAN)
            fexpr_call_builtin1(res, FEXPR_Atan, t);
        else if (g->def_kind == GR_TOWER_EXP)
            fexpr_call_builtin1(res, FEXPR_Exp, t);
        else if (g->def_kind == GR_TOWER_LOG)
            fexpr_call_builtin1(res, FEXPR_Log, t);
        else if (g->def_param == 2)
            fexpr_call_builtin1(res, FEXPR_Sqrt, t);
        else
        {
            fexpr_t one, n;
            fexpr_init(one);
            fexpr_init(n);
            fexpr_set_ui(one, 1);
            fexpr_set_si(n, g->def_param);
            fexpr_div(u, one, n);
            fexpr_pow(res, t, u);
            fexpr_clear(one);
            fexpr_clear(n);
        }
    }
    else if (g->kind == GR_TOWER_ALGEBRAIC)
    {
        /* PolynomialRootNearest(List(coefficients), approximation) */
        slong k = g->index, j, len;
        fexpr_struct * coeffs;
        fexpr_t L;

        if (g->origin != NULL)
        {
            len = g->origin->length;
            coeffs = _fexpr_vec_init(len);
            for (j = 0; j < len; j++)
                fexpr_set_fmpz(coeffs + j, g->origin->coeffs + j);
        }
        else
        {
            const gr_poly_struct * m = gr_tower_step_minpoly(T, k);
            gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
            gr_tower_flat_struct * F = &T->flat;
            fmpz_mpoly_q_t c;

            gr_tower_flat_ensure(F);
            fmpz_mpoly_q_init(c, F->mctx);
            len = m->length;
            coeffs = _fexpr_vec_init(len);
            for (j = 0; j < len; j++)
            {
                if (gr_tower_flat_set_nested_at(c, gr_poly_coeff_srcptr(m, j, below), k - 1, F) == GR_SUCCESS)
                    _flat_get_fexpr(coeffs + j, c, F->mctx, T);
                else
                    fexpr_set_symbol_builtin(coeffs + j, FEXPR_Unknown);
            }
            fmpz_mpoly_q_clear(c, F->mctx);
        }

        fexpr_init(L);
        fexpr_set_symbol_builtin(t, FEXPR_List);
        fexpr_call_vec(L, t, coeffs, len);
        fexpr_set_acb_decimal(u, &g->enclosure, digits);
        fexpr_call_builtin2(res, FEXPR_PolynomialRootNearest, L, u);
        fexpr_clear(L);
        _fexpr_vec_clear(coeffs, len);
    }
    else
        fexpr_set_symbol_builtin(res, FEXPR_Unknown);

    fexpr_clear(t);
    fexpr_clear(u);
}

/*
    Where(expression, Def(a1, ...), Def(t1, ...), ...): the element as an
    expression in the generators it involves, with their definitions in
    order (pi and i are the symbols Pi and NumberI).
*/
static int
_gr_tower_lazy_get_fexpr(fexpr_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_struct * T = x->F->T;
    int * mark;
    slong d, ndefs = 0, i;
    fexpr_struct * args;
    fexpr_t t;

    x = _gr_tower_lazy_flat_view(x);

    if (fmpz_mpoly_q_is_fmpq(&x->elem.flat.data, x->elem.flat.mctx))
    {
        fmpq_t c;
        fmpq_init(c);
        (void) fmpz_mpoly_q_get_fmpq(c, &x->elem.flat.data, x->elem.flat.mctx);
        fexpr_set_fmpq(res, c);
        fmpq_clear(c);
        return GR_SUCCESS;
    }

    mark = flint_calloc(FLINT_MAX(T->num_gens, 1), sizeof(int));
    _gr_tower_lazy_involved_gens(mark, x, T);
    for (d = 0; d < T->num_gens; d++)
        if (mark[d] && !_gr_tower_lazy_gen_is_constant_symbol(GR_TOWER_GEN(T, d)))
            ndefs++;

    args = _fexpr_vec_init(ndefs + 1);
    _flat_get_fexpr(args + 0, &x->elem.flat.data, x->elem.flat.mctx, T);

    fexpr_init(t);
    i = 1;
    {
        slong n, j;
        slong * order = _gr_tower_lazy_marked_by_creation(&n, mark, T);
        for (j = 0; j < n; j++)
        {
            const gr_tower_gen_struct * g = GR_TOWER_GEN(T, order[j]);
            fexpr_t sym;
            if (_gr_tower_lazy_gen_is_constant_symbol(g))
                continue;
            fexpr_init(sym);
            _gen_symbol(sym, g);
            _gen_def_get_fexpr(t, g, T, L->options[GR_TOWER_OPT_PRINT_DIGITS]);
            fexpr_call_builtin2(args + i, FEXPR_Def, sym, t);
            fexpr_clear(sym);
            i++;
        }
        flint_free(order);
    }

    if (ndefs == 0)
        fexpr_set(res, args + 0);
    else
    {
        fexpr_set_symbol_builtin(t, FEXPR_Where);
        fexpr_call_vec(res, t, args, ndefs + 1);
    }

    fexpr_clear(t);
    _fexpr_vec_clear(args, ndefs + 1);
    flint_free(mark);
    return GR_SUCCESS;
}

/* -------------------------------------------------------------------- */
/* generators of the context                                             */
/* -------------------------------------------------------------------- */

/*
    The generators of the context: one element per definition (the
    generator in the latest tower containing it, or the alias). Distinct
    definitions with the same name (pi and i, created independently in
    towers which were never merged) are listed once, by the generator in
    the smallest tower, so that the parser, which must tell equal names
    apart by their values, does not need to compare elements of large
    towers.
*/
int
_gr_tower_lazy_gens(gr_vec_t vec, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    ulong id;
    slong n = 0, i, j;
    const char ** names;
    slong * sizes;

    /* (generators created inside the towers, unnamed so far) */
    for (i = 0; i < L->num_towers; i++)
        _gr_tower_lazy_assign_def_ids(L->towers[i], ctx);

    gr_vec_set_length(vec, L->next_def_id, ctx);
    names = flint_malloc(sizeof(const char *) * (L->next_def_id + 1));
    sizes = flint_malloc(sizeof(slong) * (L->next_def_id + 1));

    for (id = 1; id <= L->next_def_id; id++)
    {
        gr_tower_lazy_elem_struct * res = (gr_tower_lazy_elem_struct *) gr_vec_entry_ptr(vec, n, ctx);
        int found = 0;

        names[n] = NULL;
        sizes[n] = 0;
        for (i = L->num_towers - 1; i >= 0 && !found; i--)
        {
            gr_tower_flat_struct * F = L->towers[i];
            slong d = gr_tower_find_def_order(F->T, id);
            if (d >= 0)
            {
                names[n] = GR_TOWER_GEN(F->T, d)->name;
                sizes[n] = F->T->num_gens;
                /* a definition of the same name already listed: the one
                   in the smaller tower */
                for (j = 0; j < n && names[n] != NULL; j++)
                {
                    if (names[j] != NULL && strcmp(names[j], names[n]) == 0)
                    {
                        if (sizes[n] < sizes[j])
                        {
                            gr_tower_lazy_elem_struct * prev = (gr_tower_lazy_elem_struct *) gr_vec_entry_ptr(vec, j, ctx);
                            _gr_tower_lazy_set_gen_d(prev, F, d, ctx);
                            names[j] = names[n];
                            sizes[j] = sizes[n];
                        }
                        break;
                    }
                }
                if (j < n && names[n] != NULL)
                {
                    /* (listed already) */
                    found = 0;
                    break;
                }
                _gr_tower_lazy_set_gen_d(res, F, d, ctx);
                found = 1;
            }
        }
        if (!found && names[n] != NULL)
            continue;

        for (i = 0; i < L->num_aliases && !found; i++)
        {
            if (L->aliases[i].def_id == id)
            {
                fmpz_mpoly_q_t t;
                gr_tower_flat_struct * F = L->aliases[i].F;
                gr_tower_flat_ensure(F);
                fmpz_mpoly_q_init(t, F->mctx);
                gr_tower_flat_convert(t, &L->aliases[i].data, L->aliases[i].mctx, F);
                _gr_tower_lazy_install(res, F, t, ctx);
                fmpz_mpoly_q_clear(t, F->mctx);
                found = 1;
            }
        }

        /* (in a real or algebraic view, the definitions in the field) */
        if (found && VIEW(ctx)->field_flags != 0 && _gr_tower_lazy_check_member(res, ctx) != GR_SUCCESS)
            found = 0;

        if (found)
            n++;
    }

    flint_free(names);
    flint_free(sizes);
    gr_vec_set_length(vec, n, ctx);
    return GR_SUCCESS;
}

static double
_decimal_ulp(const char * s)
{
    double best = 0;
    const char * p = s;
    while (*p)
    {
        slong frac = 0, ex = 0, exsign = 1, seen = 0;
        int in_frac = 0;
        while (*p && !(*p >= '0' && *p <= '9') && *p != '.')
            p++;
        if (!*p)
            break;
        while ((*p >= '0' && *p <= '9') || *p == '.')
        {
            if (*p == '.')
                in_frac = 1;
            else
            {
                seen = 1;
                if (in_frac)
                    frac++;
            }
            p++;
        }
        if (*p == 'e' || *p == 'E')
        {
            p++;
            if (*p == '-') { exsign = -1; p++; }
            else if (*p == '+') p++;
            while (*p >= '0' && *p <= '9')
                ex = 10 * ex + (*p++ - '0');
        }
        if (seen)
        {
            double u = pow(10.0, (double) (exsign * ex - frac));
            if (best == 0 || u > best)
                best = u;
        }
    }
    return (best == 0) ? 1e-6 : best;
}

static int
_gr_tower_lazy_parse_root(gr_tower_lazy_elem_t res, const char * s, const char * name,
    const char ** names, const gr_tower_lazy_elem_t values, slong num, gr_ctx_t ctx)
{
    slong len = strlen(s), i, depth = 0, split = -1, sprec;
    char * polystr, * approx, * approx_copy;
    gr_ctx_t P;
    gr_poly_t poly;
    gr_vec_t roots;
    fmpz_vec_t mult;
    acb_t z, w;
    arb_t t;
    int status = GR_SUCCESS;

    /* root( ... , ... ) : the last top-level comma separates the parts */
    if (len < 8 || s[len - 1] != ')')
        return GR_UNABLE;
    for (i = 5; i < len - 1; i++)
    {
        if (s[i] == '(') depth++;
        else if (s[i] == ')') depth--;
        else if (s[i] == ',' && depth == 0) split = i;
    }
    if (split < 0)
        return GR_UNABLE;

    polystr = flint_malloc(len);
    approx = flint_malloc(len);
    memcpy(polystr, s + 5, split - 5);
    polystr[split - 5] = 0;
    i = split + 1;
    while (s[i] == ' ') i++;
    memcpy(approx, s + i, len - 1 - i);
    approx[len - 1 - i] = 0;

    /* the approximation */
    acb_init(z);
    acb_init(w);
    arb_init(t);
    approx_copy = flint_malloc(len);
    strcpy(approx_copy, approx);
    sprec = FLINT_MAX(64, 4 * (slong) strlen(approx) + 32);
    {
        char * star = strstr(approx, "*i");
        if (star == NULL)
        {
            if (arb_set_str(acb_realref(z), approx, sprec) != 0)
                status = GR_UNABLE;
        }
        else
        {
            /* [re +|- ] im*i */
            char * sep = NULL, * q;
            int neg = 0;
            *star = 0;
            for (q = approx + 1; *q; q++)
                if ((q[0] == ' ' && (q[1] == '+' || q[1] == '-') && q[2] == ' '))
                    sep = q;
            if (sep != NULL)
            {
                *sep = 0;
                neg = (sep[1] == '-');
                if (arb_set_str(acb_realref(z), approx, sprec) != 0 || arb_set_str(acb_imagref(z), sep + 3, sprec) != 0)
                    status = GR_UNABLE;
            }
            else if (arb_set_str(acb_imagref(z), approx, sprec) != 0)
                status = GR_UNABLE;
            if (neg)
                arb_neg(acb_imagref(z), acb_imagref(z));
        }
    }

    /* the polynomial, over the field with the name as the variable */
    gr_ctx_init_gr_poly(P, ctx);
    GR_MUST_SUCCEED(gr_ctx_set_gen_name(P, name));
    gr_poly_init(poly, ctx);
    gr_vec_init(roots, 0, ctx);
    fmpz_vec_init(mult, 0);

    if (status == GR_SUCCESS)
    {
        /* (the terminals as constant polynomials, and the variable; the
           generators of the context are not used: the names are local) */
        gr_vec_t cvalues;
        const char ** pnames = flint_malloc((num + 1) * sizeof(char *));
        gr_vec_init(cvalues, num + 1, P);
        for (i = 0; i < num && status == GR_SUCCESS; i++)
        {
            pnames[i] = names[i];
            status = gr_poly_set_scalar(gr_vec_entry_ptr(cvalues, i, P), GR_ENTRY(values, i, ctx->sizeof_elem), ctx);
        }
        pnames[num] = name;
        status |= gr_gen(gr_vec_entry_ptr(cvalues, num, P), P);
        if (status == GR_SUCCESS)
            status = gr_generic_set_str_expr_terminals(poly, polystr, 0, pnames, cvalues->entries, num + 1, 0, P);
        gr_vec_clear(cvalues, P);
        flint_free(pnames);
    }

    if (status == GR_SUCCESS)
        status = gr_poly_roots(roots, mult, poly, 0, ctx);

    if (status == GR_SUCCESS)
    {
        /* the unique root within the rounding error of the decimal
           approximation (z, its parts enlarged by the ulp of the printed
           digits); two roots within that error are not told apart
           (GR_UNABLE), rather than the nearest one taken */
        slong best = -1, prec;
        mag_t ulp;
        mag_init(ulp);
        mag_set_d(ulp, _decimal_ulp(approx_copy));
        mag_mul_2exp_si(ulp, ulp, 1);
        arb_add_error_mag(acb_realref(z), ulp);
        arb_add_error_mag(acb_imagref(z), ulp);
        for (prec = 64; prec <= 4096 && best < 0; prec *= 2)
        {
            slong hits = 0, last = -1, tight = 1;
            for (i = 0; i < roots->length; i++)
            {
                if (gr_tower_lazy_get_acb(w, gr_vec_entry_ptr(roots, i, ctx), prec, ctx) != GR_SUCCESS)
                {
                    hits = -1;
                    break;
                }
                if (acb_overlaps(w, z))
                {
                    hits++;
                    last = i;
                    if (mag_cmp(arb_radref(acb_realref(w)), ulp) > 0 || mag_cmp(arb_radref(acb_imagref(w)), ulp) > 0)
                        tight = 0;
                }
            }
            if (hits == 1)
                best = last;
            else if (hits <= 0 || tight)
                break;
        }
        mag_clear(ulp);
        if (best < 0)
            status = GR_UNABLE;
        else
            status = _gr_tower_lazy_set(res, gr_vec_entry_ptr(roots, best, ctx), ctx);
    }

    gr_vec_clear(roots, ctx);
    fmpz_vec_clear(mult);
    gr_poly_clear(poly, ctx);
    gr_ctx_clear(P);
    acb_clear(z);
    acb_clear(w);
    arb_clear(t);
    flint_free(polystr);
    flint_free(approx);
    flint_free(approx_copy);
    return status;
}

/*
    Parsing. A string "expr {name = def; name = def; ...}" as printed
    with the definitions is read by evaluating the definitions in order,
    each in the names defined before it (the names are local to the
    string: they need not be those of this context, so that the printed
    form of an element is persistent across contexts and sessions), then
    the expression in all of them. A string without definitions is read
    in the generators of the context (the generic parser knows them by
    their printed form, which must be the bare names here), and in pi,
    i and the elementary functions.
*/
int
gr_tower_lazy_set_str(gr_tower_lazy_elem_t res, const char * s, gr_ctx_t ctx)
{
    int status = GR_SUCCESS, flags;
    slong len = strlen(s);
    const char * brace;

    _gr_tower_lazy_lock(ctx);
    flags = LAZY(ctx)->options[GR_TOWER_OPT_PRINT_FLAGS];
    LAZY(ctx)->options[GR_TOWER_OPT_PRINT_FLAGS] = GR_TOWER_PRINT_SYMBOLIC;

    brace = (len > 0 && s[len - 1] == '}') ? strstr(s, " {") : NULL;

    if (brace == NULL)
    {
        status = gr_generic_set_str_ring_exponents(res, s, ctx);
    }
    else
    {
        /* the definitions */
        char * expr = flint_malloc(len + 1);
        char * defs = flint_malloc(len + 1);
        char ** names = NULL;
        gr_vec_t values;
        slong num = 0, i;
        char * p;

        memcpy(expr, s, brace - s);
        expr[brace - s] = 0;
        memcpy(defs, brace + 2, len - (brace - s) - 3);
        defs[len - (brace - s) - 3] = 0;

        gr_vec_init(values, 0, ctx);

        p = defs;
        while (*p && status == GR_SUCCESS)
        {
            char * end = strstr(p, "; ");
            char * eq;
            slong nlen;

            if (end != NULL)
                *end = 0;
            eq = strstr(p, " = ");
            if (eq == NULL)
            {
                status = GR_UNABLE;
                break;
            }
            nlen = eq - p;
            names = flint_realloc(names, (num + 1) * sizeof(char *));
            names[num] = flint_malloc(nlen + 1);
            memcpy(names[num], p, nlen);
            names[num][nlen] = 0;
            gr_vec_set_length(values, num + 1, ctx);
            if (strncmp(eq + 3, "root(", 5) == 0 && strrchr(eq + 3, ',') != NULL && strchr(strrchr(eq + 3, ','), '.') != NULL)
                /* root(m, approximation) */
                status = _gr_tower_lazy_parse_root(gr_vec_entry_ptr(values, num, ctx), eq + 3, names[num],
                    (const char **) names, values->entries, num, ctx);
            else if (strncmp(eq + 3, "root(", 5) == 0 && strrchr(eq + 3, ',') != NULL)
            {
                /* root(x, n): the principal n-th root */
                char * comma = strrchr(eq + 3, ',');
                slong n = atol(comma + 1), elen = strlen(eq + 3);
                gr_tower_lazy_elem_struct t;
                if (n < 1 || eq[3 + elen - 1] != ')')
                    status = GR_UNABLE;
                else
                {
                    *comma = 0;
                    _gr_tower_lazy_init(&t, ctx);
                    status = gr_generic_set_str_expr_terminals(&t, eq + 8, GR_PARSE_RING_EXPONENTS,
                        (const char **) names, values->entries, num, 0, ctx);
                    if (status == GR_SUCCESS)
                        status = _gr_tower_lazy_root_ui(gr_vec_entry_ptr(values, num, ctx), &t, n, ctx);
                    _gr_tower_lazy_clear(&t, ctx);
                }
            }
            else if (strncmp(eq + 3, "polylog(", 8) == 0 || strncmp(eq + 3, "polygamma(", 10) == 0 ||
                     (strncmp(eq + 3, "lambertw(", 9) == 0 && strchr(eq + 3, ',') != NULL))
            {
                /* polylog(s, x), polygamma(m, x), lambertw(x, k) */
                char * open = strchr(eq + 3, '(');
                char * comma = strchr(open, ',');
                slong elen = strlen(eq + 3);
                int is_w = (eq[3] == 'l');
                gr_tower_lazy_elem_struct t, u;
                if (comma == NULL || eq[3 + elen - 1] != ')')
                    status = GR_UNABLE;
                else
                {
                    char * arg;
                    slong param;
                    _gr_tower_lazy_init(&t, ctx);
                    _gr_tower_lazy_init(&u, ctx);
                    eq[3 + elen - 1] = 0;
                    *comma = 0;
                    param = is_w ? atol(comma + 1) : atol(open + 1);
                    arg = is_w ? open + 1 : comma + 1;
                    while (*arg == ' ')
                        arg++;
                    status = gr_generic_set_str_expr_terminals(&t, arg, GR_PARSE_RING_EXPONENTS,
                        (const char **) names, values->entries, num, 0, ctx);
                    if (status == GR_SUCCESS)
                    {
                        if (is_w)
                        {
                            fmpz_t k;
                            fmpz_init_set_si(k, param);
                            status = gr_lambertw_fmpz(gr_vec_entry_ptr(values, num, ctx), &t, k, ctx);
                            fmpz_clear(k);
                        }
                        else
                        {
                            status = gr_set_si(&u, param, ctx);
                            if (strncmp(eq + 3, "polylog(", 8) == 0)
                                status |= gr_polylog(gr_vec_entry_ptr(values, num, ctx), &u, &t, ctx);
                            else
                                status |= gr_polygamma(gr_vec_entry_ptr(values, num, ctx), &u, &t, ctx);
                        }
                    }
                    _gr_tower_lazy_clear(&t, ctx);
                    _gr_tower_lazy_clear(&u, ctx);
                }
            }
            else
                status = gr_generic_set_str_expr_terminals(gr_vec_entry_ptr(values, num, ctx), eq + 3,
                    GR_PARSE_RING_EXPONENTS, (const char **) names, values->entries, num, 0, ctx);
            num++;
            p = (end == NULL) ? p + strlen(p) : end + 2;
        }

        if (status == GR_SUCCESS)
            status = gr_generic_set_str_expr_terminals(res, expr, GR_PARSE_RING_EXPONENTS,
                (const char **) names, values->entries, num, 0, ctx);

        for (i = 0; i < num; i++)
            flint_free(names[i]);
        flint_free(names);
        gr_vec_clear(values, ctx);
        flint_free(expr);
        flint_free(defs);
    }

    /* (auxiliary generators in the definitions may lie outside the
       field; the value may not) */
    if (status == GR_SUCCESS && _gr_tower_lazy_outermost(ctx))
        status = _gr_tower_lazy_check_member(res, ctx);
    if (status != GR_SUCCESS)
        GR_MUST_SUCCEED(_gr_tower_lazy_zero(res, ctx));

    LAZY(ctx)->options[GR_TOWER_OPT_PRINT_FLAGS] = flags;
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int gr_tower_lazy_get_fexpr(fexpr_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_get_fexpr(res, x, ctx)) }

int gr_tower_lazy_get_qqbar(qqbar_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_get_qqbar_impl(res, x, ctx)) }
int gr_tower_lazy_get_acb(acb_t res, const gr_tower_lazy_elem_t x, slong prec, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_get_acb_impl(res, x, prec, ctx)) }
gr_tower_struct * gr_tower_lazy_get_tower(slong * level, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED_T(gr_tower_struct *, _gr_tower_lazy_get_tower_impl(level, x, ctx)) }
const fmpz_mpoly_q_struct * gr_tower_lazy_get_data(const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED_T(const fmpz_mpoly_q_struct *, _gr_tower_lazy_get_data_impl(x, ctx)) }
int gr_tower_lazy_get_fmpq_poly(fmpq_poly_t res, fmpz_poly_t modulus, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_get_fmpq_poly_impl(res, modulus, (gr_tower_lazy_elem_struct *) x, ctx)) }
/* (the representation tag is published with a release store: no lock) */
int gr_tower_lazy_repr(const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { return (int) LAZY_REPR((const gr_tower_lazy_elem_struct *) x); }
void gr_tower_lazy_ctx_stats(gr_ctx_t ctx) { LOCKED_V(_gr_tower_lazy_ctx_stats_impl(ctx)) }

POP_OPTIONS
