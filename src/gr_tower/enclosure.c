/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_poly.h"
#include "acb_poly.h"
#include "arb_fmpz_poly.h"
#include "fmpz_mpoly_q.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

/* Numerical evaluation of an element x of F_k. */
int
gr_tower_get_acb_at(acb_t res, gr_srcptr x, slong k, slong prec, gr_tower_t T)
{
    int status = GR_SUCCESS;

    if (k == 0)
    {
        if (T->base->which_ring == GR_CTX_FMPZ_MPOLY_Q)
        {
            /* rational function in the transcendental generators */
            const fmpz_mpoly_ctx_struct * bctx = gr_ctx_fmpz_mpoly_q_mctx(T->base);
            slong r = T->num_trans, j;
            acb_ptr vals = _acb_vec_init(r);
            int * used = flint_malloc(sizeof(int) * r);

            fmpz_mpoly_q_used_vars(used, x, bctx);
            for (j = 0; j < r && status == GR_SUCCESS; j++)
                if (used[j])
                    status |= gr_tower_trans_get_acb(vals + j, T, j + 1, prec);

            if (status == GR_SUCCESS)
                fmpz_mpoly_q_evaluate_acb(res, x, vals, prec, bctx);

            _acb_vec_clear(vals, r);
            flint_free(used);
        }
        else
        {
            gr_ctx_t C;
            gr_ctx_init_complex_acb(C, prec);
            status = gr_set_other(res, x, T->base, C);
            gr_ctx_clear(C);
        }
    }
    else
    {
        const gr_poly_struct * poly = x;
        gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
        slong i, len = poly->length;
        acb_t z, t;

        if (len == 0)
        {
            acb_zero(res);
            return GR_SUCCESS;
        }

        /* a constant: no need to evaluate the generator */
        if (len == 1)
            return gr_tower_get_acb_at(res, gr_poly_coeff_srcptr(poly, 0, below), k - 1, prec, T);

        acb_init(z);
        acb_init(t);

        status |= gr_tower_step_get_acb(z, T, k, prec);

        if (status == GR_SUCCESS)
        {
            status |= gr_tower_get_acb_at(res, gr_poly_coeff_srcptr(poly, len - 1, below), k - 1, prec, T);

            for (i = len - 2; i >= 0 && status == GR_SUCCESS; i--)
            {
                acb_mul(res, res, z, prec);
                status |= gr_tower_get_acb_at(t, gr_poly_coeff_srcptr(poly, i, below), k - 1, prec, T);
                acb_add(res, res, t, prec);
            }
        }

        acb_clear(z);
        acb_clear(t);
    }

    return status;
}

/* Numerical evaluation of an element x of the top field with relative
   accuracy (about) prec bits, increasing the working precision as
   needed (an element with large coefficients, of which the value is
   the result of cancellation, needs more than prec bits); GR_UNABLE if
   the accuracy is not reached. */
int
_gr_tower_get_acb_accurate(acb_t res, gr_srcptr x, slong prec, gr_tower_t T)
{
    slong wp;
    int status = GR_UNABLE;

    for (wp = prec; wp <= 64 * prec; wp *= 2)
    {
        status = gr_tower_get_acb(res, x, wp, T);
        if (status != GR_SUCCESS)
            return status;
        if (acb_is_finite(res) && (acb_is_exact(res) || acb_rel_accuracy_bits(res) >= prec))
            return GR_SUCCESS;
        status = GR_UNABLE;
    }

    return status;
}

/* An enclosure of the principal n-th root of the top-field element x
   which isolates it from the other n-th roots: the enclosure of x is
   refined until it does not meet the branch cut (the negative real axis;
   an exactly real negative x -- which the caller may assert through an
   exactly real *given* enclosure -- gives the root with argument pi/n),
   and the radius of the root is below |root| sin(pi/n) / 2. Returns
   GR_UNABLE when this is not reached within the precision limit (an x
   which is not known to be real and whose enclosure straddles the cut). */
int
_gr_tower_principal_root_enclosure(acb_t res, gr_srcptr x, ulong n, const acb_t given, gr_tower_t T)
{
    acb_t w;
    mag_t rad, lo;
    slong prec, limit;
    int exact_real = (given != NULL) && arb_is_zero(acb_imagref(given));
    int status = GR_UNABLE;

    if (n == 0)
        return GR_DOMAIN;
    if (n == 1)
    {
        if (given != NULL)
        {
            acb_set(res, given);
            return GR_SUCCESS;
        }
        return _gr_tower_get_acb_accurate(res, x, GR_TOWER_DEFAULT_PREC, T);
    }

    acb_init(w);
    mag_init(rad);
    mag_init(lo);
    limit = GR_TOWER_OPTION(T, GR_TOWER_OPT_CERTIFY_PREC_LIMIT);

    for (prec = GR_TOWER_DEFAULT_PREC; prec <= limit; prec *= 2)
    {
        if (given != NULL && prec == GR_TOWER_DEFAULT_PREC)
            acb_set(w, given);
        else
        {
            status = _gr_tower_get_acb_accurate(w, x, prec, T);
            if (status != GR_SUCCESS)
                break;
            status = GR_UNABLE;
            if (exact_real)
                arb_zero(acb_imagref(w));
        }

        if (acb_contains_zero(w))
            continue;

        /* the ball must not straddle the cut */
        if (!arb_is_zero(acb_imagref(w)) && arb_contains_zero(acb_imagref(w)) && !arb_is_positive(acb_realref(w)))
            continue;

        acb_root_ui(res, w, n, prec);

        /* isolation: rad * n < |res| implies rad < |res| sin(pi/n) / 2 */
        acb_get_mag_lower(lo, res);
        mag_hypot(rad, arb_radref(acb_realref(res)), arb_radref(acb_imagref(res)));
        mag_mul_ui(rad, rad, n);
        if (mag_cmp(rad, lo) < 0)
        {
            status = GR_SUCCESS;
            break;
        }
    }

    acb_clear(w);
    mag_clear(rad);
    mag_clear(lo);
    return status;
}

/*
    Selection of the factor vanishing at a root: among the factors of a
    polynomial (integer polynomials, or polynomials over F_k of the
    tower), the index of the one vanishing at the value z, computed by
    get_z at the precisions DEFAULT, 2 DEFAULT, ... up to the certification
    precision limit until exactly one factor vanishes; -1 when none
    does (z is then not a root of the product) or when the precision
    limit is reached with several candidates.
*/
slong
_gr_tower_select_fmpz_factor(const fmpz_poly_factor_t fac, int (*get_z)(acb_t, slong, void *), void * arg, gr_tower_t T)
{
    acb_t z, v;
    slong prec, i, res = -1, limit = GR_TOWER_OPTION(T, GR_TOWER_OPT_CERTIFY_PREC_LIMIT);

    if (fac->num == 1)
        return 0;

    acb_init(z);
    acb_init(v);
    for (prec = GR_TOWER_DEFAULT_PREC; prec <= limit && res == -1; prec *= 2)
    {
        slong count = 0, cand = -1;
        if (get_z(z, prec, arg) != GR_SUCCESS)
            break;
        for (i = 0; i < fac->num; i++)
        {
            arb_fmpz_poly_evaluate_acb(v, fac->p + i, z, prec);
            if (acb_contains_zero(v))
                count++, cand = i;
        }
        if (count == 1)
            res = cand;
        else if (count == 0)
            break;
    }
    acb_clear(z);
    acb_clear(v);
    return res;
}

slong
_gr_tower_select_poly_factor(gr_vec_t fac, slong k, int (*get_z)(acb_t, slong, void *), void * arg, gr_tower_t T)
{
    gr_ctx_t pctx;
    acb_t z, v;
    acb_poly_t hz;
    slong prec, i, res = -1, limit = GR_TOWER_OPTION(T, GR_TOWER_OPT_CERTIFY_PREC_LIMIT);

    if (fac->length == 1)
        return 0;

    gr_ctx_init_gr_poly(pctx, gr_tower_field_at(T, k));
    acb_init(z);
    acb_init(v);
    acb_poly_init(hz);
    for (prec = GR_TOWER_DEFAULT_PREC; prec <= limit && res == -1; prec *= 2)
    {
        slong count = 0, cand = -1;
        if (get_z(z, prec, arg) != GR_SUCCESS)
            break;
        for (i = 0; i < fac->length; i++)
        {
            if (_gr_tower_poly_get_acb_poly(hz, gr_vec_entry_ptr(fac, i, pctx), k, prec, T) != GR_SUCCESS)
            {
                count = -1;
                break;
            }
            acb_poly_evaluate(v, hz, z, prec);
            if (acb_contains_zero(v))
                count++, cand = i;
        }
        if (count == 1)
            res = cand;
        else if (count <= 0)
            break;
    }
    acb_poly_clear(hz);
    acb_clear(z);
    acb_clear(v);
    gr_ctx_clear(pctx);
    return res;
}

/*
    Selection between two candidate values by their enclosures: which of
    the candidates 0 and 1 (enclosures from get(z, which, prec, arg))
    overlaps the target (which = -1) when the other does not, with the
    precision doubled from DEFAULT up to the certification limit; -1 when
    neither overlaps (inconsistent data) or at the precision limit.
*/
int
_gr_tower_select_candidate(int (*get)(acb_t, int, slong, void *), void * arg, gr_tower_t T)
{
    acb_t z, z0, z1;
    slong prec, limit = GR_TOWER_OPTION(T, GR_TOWER_OPT_CERTIFY_PREC_LIMIT);
    int res = -1;

    acb_init(z);
    acb_init(z0);
    acb_init(z1);
    for (prec = GR_TOWER_DEFAULT_PREC; prec <= limit; prec *= 2)
    {
        int o0, o1;
        if (get(z, -1, prec, arg) != GR_SUCCESS || get(z0, 0, prec, arg) != GR_SUCCESS || get(z1, 1, prec, arg) != GR_SUCCESS)
            break;
        o0 = acb_overlaps(z, z0);
        o1 = acb_overlaps(z, z1);
        if (o0 != o1)
        {
            res = o1;
            break;
        }
        if (!o0)
            break;
    }
    acb_clear(z);
    acb_clear(z0);
    acb_clear(z1);
    return res;
}

/* get_z callback: a fixed enclosure */
int
_gr_tower_get_z_fixed(acb_t z, slong prec, void * arg)
{
    acb_set(z, (const acb_struct *) arg);
    return GR_SUCCESS;
}

/* Converts a polynomial with coefficients in F_k to an acb_poly. */
int
_gr_tower_poly_get_acb_poly(acb_poly_t res, const gr_poly_t f, slong k, slong prec, gr_tower_t T)
{
    int status = GR_SUCCESS;
    slong i, len = f->length;
    gr_ctx_struct * ctx = gr_tower_field_at(T, k);

    acb_poly_fit_length(res, len);

    for (i = 0; i < len && status == GR_SUCCESS; i++)
        status |= gr_tower_get_acb_at(res->coeffs + i, gr_poly_coeff_srcptr(f, i, ctx), k, prec, T);

    _acb_poly_set_length(res, len);
    return status;
}

/*
    Interval Newton refinement of a root of f (a polynomial with ball
    coefficients, assumed accurate to roughly wp bits). Given an enclosure
    z, computes t = mid(z) - f(mid(z)) / f'(z). If t is contained in z,
    then z contains exactly one root, which is also contained in t.

    Returns 1 and sets res = t on success, 0 if the inclusion fails.

    For a real polynomial and a z straddling the real axis, the unique
    root (certified in the complex rectangle z) is real when a real
    Newton step certifies a real root in Re(z): res is then given an
    exactly zero imaginary part (which matters for branch cuts of
    logarithms of the root).
*/
int
_gr_tower_newton_step(acb_t res, const acb_poly_t f, const acb_t z_in, slong wp)
{
    acb_t zmid, t, u, z;
    acb_poly_t df;
    int success;

    acb_init(zmid);
    acb_init(t);
    acb_init(u);
    acb_poly_init(df);
    acb_init(z);
    acb_set(z, z_in);   /* (res may alias z_in) */

    acb_poly_derivative(df, f, wp);
    acb_get_mid(zmid, z);
    acb_poly_evaluate(t, f, zmid, wp);
    acb_poly_evaluate(u, df, z, wp);

    if (acb_contains_zero(u))
    {
        success = 0;
    }
    else
    {
        acb_div(t, t, u, wp);
        acb_sub(t, zmid, t, wp);
        success = acb_contains(z, t);
        if (success)
            acb_set(res, t);
    }

    if (success && arb_contains_zero(acb_imagref(z)) && !arb_is_zero(acb_imagref(z)))
    {
        slong i;
        int real = 1;

        for (i = 0; i < f->length; i++)
            if (!arb_is_zero(acb_imagref(f->coeffs + i)))
                real = 0;

        if (real)
        {
            arb_poly_t fr, dfr;
            arb_t xmid, tr, ur;

            arb_poly_init(fr);
            arb_poly_init(dfr);
            arb_init(xmid);
            arb_init(tr);
            arb_init(ur);

            arb_poly_fit_length(fr, f->length);
            for (i = 0; i < f->length; i++)
                arb_set(fr->coeffs + i, acb_realref(f->coeffs + i));
            _arb_poly_set_length(fr, f->length);
            arb_poly_derivative(dfr, fr, wp);

            arb_get_mid_arb(xmid, acb_realref(z));
            arb_poly_evaluate(tr, fr, xmid, wp);
            arb_poly_evaluate(ur, dfr, acb_realref(z), wp);

            if (!arb_contains_zero(ur))
            {
                arb_div(tr, tr, ur, wp);
                arb_sub(tr, xmid, tr, wp);
                if (arb_contains(acb_realref(z), tr))
                {
                    /* the real root in Re(z) is the unique root in z */
                    arb_set(acb_realref(res), tr);
                    arb_zero(acb_imagref(res));
                }
            }

            arb_poly_clear(fr);
            arb_poly_clear(dfr);
            arb_clear(xmid);
            arb_clear(tr);
            arb_clear(ur);
        }
    }

    acb_clear(zmid);
    acb_clear(t);
    acb_clear(u);
    acb_clear(z);
    acb_poly_clear(df);

    return success;
}

/*
    Certifies that a disk around the midpoint of z contains exactly one
    root of m (a monic polynomial over F_k) and returns a refined
    enclosure of that root. The disk is z itself, or, when z is too
    tight (an exact approximation, say), z inflated by |z| 2^-j for
    decreasing j down to 2; the disks being nested, the certified root
    is the root of m nearest to the midpoint of z (a root nearer to it
    would lie in the certified disk). If no disk can be certified
    (several roots within |z|/4 of the midpoint, or insufficient
    precision), returns GR_UNABLE.
*/
int
_gr_tower_certify_root(acb_t res, const gr_poly_t m, slong k, const acb_t z, slong prec, gr_tower_t T)
{
    acb_poly_t f;
    acb_t zz;
    mag_t r, s;
    slong wp, j, acc;
    int status = GR_UNABLE;

    if (m->length == 2)
    {
        /* linear: the root is -m[0], exactly */
        gr_ctx_struct * ctx = gr_tower_field_at(T, k);
        int st = gr_tower_get_acb_at(res, gr_poly_coeff_srcptr(m, 0, ctx), k, prec, T);
        acb_neg(res, res);
        return st;
    }

    acb_poly_init(f);
    acb_init(zz);
    mag_init(r);
    mag_init(s);

    acc = acb_rel_accuracy_bits(z);
    if (acc < 0) acc = 0;
    if (acc > prec) acc = prec;

    for (wp = FLINT_MAX(64, prec); wp <= 16 * FLINT_MAX(64, prec); wp *= 2)
    {
        if (_gr_tower_poly_get_acb_poly(f, m, k, wp, T) != GR_SUCCESS)
            break;

        /* first try the disk as given */
        if (_gr_tower_newton_step(zz, f, z, wp))
        {
            status = GR_SUCCESS;
            break;
        }

        /* then inflate the radius, starting from roughly the accuracy of z */
        acb_get_mag(s, z);
        if (mag_is_zero(s))
            mag_one(s);

        for (j = FLINT_MIN(acc, wp); j >= 2; j /= 2)
        {
            acb_set(zz, z);
            mag_mul_2exp_si(r, s, -j);
            mag_add(arb_radref(acb_realref(zz)), arb_radref(acb_realref(zz)), r);
            mag_add(arb_radref(acb_imagref(zz)), arb_radref(acb_imagref(zz)), r);

            if (_gr_tower_newton_step(zz, f, zz, wp))
            {
                status = GR_SUCCESS;
                break;
            }
        }

        if (status == GR_SUCCESS)
            break;
    }

    if (status == GR_SUCCESS)
        acb_set(res, zz);

    acb_poly_clear(f);
    acb_clear(zz);
    mag_clear(r);
    mag_clear(s);

    return status;
}

/* Enclosure of the generator a_k, refined to (relative) precision prec. */
int
gr_tower_step_get_acb(acb_t res, gr_tower_t T, slong k, slong prec)
{
    gr_tower_gen_struct * step = GR_TOWER_STEP(T, k - 1);
    acb_t z, t;
    acb_poly_t f;
    slong wp, iter;
    int status = GR_SUCCESS;

    if (acb_is_exact(&step->enclosure) || acb_rel_accuracy_bits(&step->enclosure) >= prec)
    {
        /* the cached enclosure may be far more precise than requested;
           a rounded copy keeps subsequent arithmetic at the requested cost */
        acb_set_round(res, &step->enclosure, prec + 10);
        return GR_SUCCESS;
    }

    acb_init(z);
    acb_init(t);
    acb_poly_init(f);

    acb_set(z, &step->enclosure);
    wp = FLINT_MAX(acb_rel_accuracy_bits(z), 32) + 20;

    /* a generator with a linear modulus X - c is the element c: no
       Newton iteration (which requires the enclosure of c to shrink
       inside the cached enclosure) */
    if (gr_poly_quotient_ctx_modulus(step->ctx)->length == 2)
    {
        gr_ctx_struct * below = gr_tower_field_at(T, k - 1);
        gr_srcptr c0 = gr_poly_coeff_srcptr(gr_poly_quotient_ctx_modulus(step->ctx), 0, below);
        slong last = -1;

        wp = prec + 20;
        for (iter = 0; iter < 8; iter++)
        {
            status = gr_tower_get_acb_at(t, c0, k - 1, wp, T);
            if (status != GR_SUCCESS)
                break;
            acb_neg(t, t);
            if (acb_rel_accuracy_bits(t) >= prec || acb_rel_accuracy_bits(t) <= last)
                break;
            last = acb_rel_accuracy_bits(t);
            wp *= 2;
        }

        if (status == GR_SUCCESS)
        {
            if (acb_rel_accuracy_bits(t) > acb_rel_accuracy_bits(&step->enclosure))
            {
                acb_set(&step->enclosure, t);
                step->enclosure_prec = prec;
            }
            acb_set_round(res, t, prec + 10);
        }

        acb_clear(z);
        acb_clear(t);
        acb_poly_clear(f);
        return status;
    }

    for (iter = 0; ; iter++)
    {
        /* the working precision needed is a small multiple of prec
           (plus the loss in the coefficients); far beyond that, the
           enclosure is not going to improve */
        if (iter > 60 || wp > 64 * prec + 65536)
        {
            status = GR_UNABLE;
            break;
        }

        wp *= 2;

        status = _gr_tower_poly_get_acb_poly(f, gr_poly_quotient_ctx_modulus(step->ctx), k - 1, wp, T);
        if (status != GR_SUCCESS)
            break;

        if (_gr_tower_newton_step(t, f, z, wp))
        {
            /* Newton converges quadratically once it converges at all;
               accept t if it is better than z */
            if (acb_rel_accuracy_bits(t) > acb_rel_accuracy_bits(z))
                acb_set(z, t);

            if (acb_rel_accuracy_bits(z) >= prec)
                break;
        }
        else if (arb_is_zero(acb_imagref(z)) || arb_is_zero(acb_realref(z)))
        {
            /* an enclosure with an exactly zero part (a real root found
               in a real tower, say) cannot contain the Newton image once
               the coefficients are evaluated with rounding errors in
               that part (the tower having grown by complex generators):
               widen the part slightly */
            mag_t r;
            mag_init(r);
            acb_get_mag(r, z);
            mag_mul_2exp_si(r, r, -FLINT_MAX(acb_rel_accuracy_bits(z), 32));
            if (arb_is_zero(acb_imagref(z)))
                mag_set(arb_radref(acb_imagref(z)), r);
            if (arb_is_zero(acb_realref(z)))
                mag_set(arb_radref(acb_realref(z)), r);
            mag_clear(r);
        }
        /* otherwise: increase the working precision and retry */
    }

    if (status == GR_SUCCESS)
    {
        acb_set(&step->enclosure, z);
        step->enclosure_prec = prec;
        acb_set_round(res, z, prec + 10);
    }

    acb_clear(z);
    acb_clear(t);
    acb_poly_clear(f);

    return status;
}

/* Enclosure of the transcendental generator t_j, computed from its definition. */
int
gr_tower_trans_get_acb(acb_t res, gr_tower_t T, slong j, slong prec)
{
    gr_tower_gen_struct * t = GR_TOWER_TRANS(T, j - 1);
    int status = GR_SUCCESS;

    if (t->enclosure_prec >= prec && (acb_is_exact(&t->enclosure) || acb_rel_accuracy_bits(&t->enclosure) >= prec))
    {
        acb_set_round(res, &t->enclosure, prec + 10);
        return GR_SUCCESS;
    }

    if (t->kind == GR_TOWER_FREE)
    {
        /* (a formal variable has no value) */
        acb_indeterminate(res);
        return GR_UNABLE;
    }
    else if (t->kind == GR_TOWER_PI)
    {
        acb_const_pi(res, prec);
    }
    else if (t->kind == GR_TOWER_CONSTANT)
    {
        slong wp;
        for (wp = prec + 20; ; wp *= 2)
        {
            status = _gr_tower_special_eval(res, t->kind, t->def_param, NULL, wp);
            if (status != GR_SUCCESS || acb_rel_accuracy_bits(res) >= prec)
                break;
        }
    }
    else if (t->num_xargs > 0)
    {
        /* a function of several arguments */
        slong wp, n = _gr_tower_gen_num_args(t);
        acb_ptr u = _acb_vec_init(n);
        int flags = (t->kind == GR_TOWER_HYPGEOM) ? _gr_tower_gen_hypgeom_flags(t, &T->flat) : 0;

        for (wp = prec + 20; ; wp *= 2)
        {
            status = _gr_tower_gen_args_get_acb(u, t, wp, &T->flat);
            if (status != GR_SUCCESS)
                break;
            status = _gr_tower_special_eval_multi_flags(res, t->kind, t->def_param, u, n, flags, wp);
            if (status != GR_SUCCESS)
            {
                status = GR_SUCCESS;
                acb_indeterminate(res);
            }
            else if (_gr_tower_special_real_at_multi(t->kind, t->def_param, u, n, wp))
                arb_zero(acb_imagref(res));
            if (acb_is_finite(res) && (acb_is_exact(res) || acb_rel_accuracy_bits(res) >= prec))
                break;
            if (wp > 16 * prec + 4096)
            {
                status = GR_UNABLE;
                break;
            }
        }

        _acb_vec_clear(u, n);
    }
    else if (GR_TOWER_KIND_HAS_ARG(t->kind))
    {
        acb_t u;
        fmpz_mpoly_q_t arg;
        gr_tower_flat_struct * F = &T->flat;

        acb_init(u);
        gr_tower_flat_ensure(F);
        fmpz_mpoly_q_init(arg, F->mctx);
        gr_tower_flat_convert(arg, &t->arg.data, t->arg.mctx, F);

        /* the argument may cancel catastrophically, or be huge (exp needs
           it to absolute precision): evaluate with increasing working
           precision until the result is accurate */
        {
            slong wp;
            for (wp = prec + 20; ; wp *= 2)
            {
                status = gr_tower_flat_get_acb(u, arg, wp, F);
                if (status != GR_SUCCESS)
                    break;
                if (t->kind == GR_TOWER_EXP)
                    acb_exp(res, u, wp);
                else if (t->kind == GR_TOWER_LOG)
                    acb_log(res, u, wp);
                else if (GR_TOWER_KIND_IS_SPECIAL(t->kind))
                {
                    /* an argument whose imaginary part is exactly zero
                       is real (real generators have such enclosures);
                       the value is then real where the function is */
                    status = _gr_tower_special_eval(res, t->kind, t->def_param, u, wp);
                    if (status != GR_SUCCESS)
                    {
                        /* (e.g. the argument is too wide near a pole) */
                        status = GR_SUCCESS;
                        acb_indeterminate(res);
                    }
                    else if (_gr_tower_special_real_at(t->kind, t->def_param, u, wp))
                        arb_zero(acb_imagref(res));
                }
                else
                {
                    /* a generator asserted real (tan(u), atan(u) of the
                       real u the caller promised): a real value, with an
                       exactly zero imaginary part as for the other real
                       generators */
                    if (t->real && arb_contains_zero(acb_imagref(u)))
                        arb_zero(acb_imagref(u));
                    if (t->kind == GR_TOWER_TAN)
                        acb_tan(res, u, wp);
                    else
                        acb_atan(res, u, wp);
                    if (t->real)
                        arb_zero(acb_imagref(res));
                }
                if (acb_is_finite(res) && (acb_is_exact(res) || acb_rel_accuracy_bits(res) >= prec))
                    break;
                if (wp > 16 * prec + 4096)
                {
                    status = GR_UNABLE;
                    break;
                }
            }
        }

        fmpz_mpoly_q_clear(arg, F->mctx);
        acb_clear(u);
    }
    else
    {
        return GR_UNABLE;
    }

    if (status == GR_SUCCESS)
    {
        acb_set(&t->enclosure, res);
        t->enclosure_prec = prec;
        acb_set_round(res, res, prec + 10);
    }

    return status;
}

POP_OPTIONS
