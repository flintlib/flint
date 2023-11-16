/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_calc.h"
#include "acb_calc.h"

static void
quad_simple(acb_t res, acb_calc_func_t f, void * param,
        const acb_t a, const acb_t b, slong prec)
{
    acb_t mid, delta, wide;
    mag_t tmpm;

    acb_init(mid);
    acb_init(delta);
    acb_init(wide);
    mag_init(tmpm);

    /* delta = (b-a)/2 */
    acb_sub(delta, b, a, prec);
    acb_mul_2exp_si(delta, delta, -1);

    /* mid = (a+b)/2 */
    acb_add(mid, a, b, prec);
    acb_mul_2exp_si(mid, mid, -1);

    /* wide = mid +- [delta] */
    acb_set(wide, mid);
    arb_get_mag(tmpm, acb_realref(delta));
    arb_add_error_mag(acb_realref(wide), tmpm);
    arb_get_mag(tmpm, acb_imagref(delta));
    arb_add_error_mag(acb_imagref(wide), tmpm);

    /* Direct evaluation: integral = (b-a) * f([a,b]). */
    f(res, wide, param, 0, prec);
    acb_mul(res, res, delta, prec);
    acb_mul_2exp_si(res, res, 1);

    acb_clear(mid);
    acb_clear(delta);
    acb_clear(wide);
    mag_clear(tmpm);
}

static void
heap_up(acb_ptr as, acb_ptr bs, acb_ptr vs, mag_ptr ms, slong n)
{
    slong i, max, l, r;
    i = 0;
    for (;;)
    {
        max = i;
        l = 2 * i + 1;
        r = 2 * i + 2;
        if (l < n && mag_cmp(ms + l, ms + max) > 0)
            max = l;
        if (r < n && mag_cmp(ms + r, ms + max) > 0)
            max = r;
        if (max != i)
        {
            acb_swap(as + i, as + max);
            acb_swap(bs + i, bs + max);
            acb_swap(vs + i, vs + max);
            mag_swap(ms + i, ms + max);
            i = max;
        }
        else
        {
            break;
        }
    }
}

static void
heap_down(acb_ptr as, acb_ptr bs, acb_ptr vs, mag_ptr ms, slong n)
{
    slong j, k;

    k = n - 1;
    j = (k - 1) / 2;

    while (k > 0 && mag_cmp(ms + j, ms + k) < 0)
    {
        acb_swap(as + j, as + k);
        acb_swap(bs + j, bs + k);
        acb_swap(vs + j, vs + k);
        mag_swap(ms + j, ms + k);
        k = j;
        j = (j - 1) / 2;
    }
}

static int
_acb_overlaps(acb_t tmp, const acb_t a, const acb_t b, slong prec)
{
    acb_sub(tmp, a, b, prec);
    return acb_contains_zero(tmp);
}

int
acb_calc_integrate(acb_t res, acb_calc_func_t f, void * param,
    const acb_t a, const acb_t b,
    slong goal, const mag_t tol,
    const acb_calc_integrate_opt_t options,
    slong prec)
{
    acb_ptr as, bs, vs;
    mag_ptr ms;
    acb_t s, t, u;
    mag_t tmpm, tmpn, new_tol;
    slong depth_limit, eval_limit, deg_limit;
    slong depth, depth_max, eval, feval, top;
    slong leaf_interval_count;
    slong alloc;
    int stopping, real_error, use_heap, status, gl_status, verbose;

    if (options == NULL)
    {
        acb_calc_integrate_opt_t opt;
        acb_calc_integrate_opt_init(opt);
        return acb_calc_integrate(res, f, param, a, b, goal, tol, opt, prec);
    }

    status = ARB_CALC_SUCCESS;

    acb_init(s);
    acb_init(t);
    acb_init(u);
    mag_init(tmpm);
    mag_init(tmpn);
    mag_init(new_tol);

    depth_limit = options->depth_limit;
    if (depth_limit <= 0)
        depth_limit = 2 * prec;
    depth_limit = FLINT_MAX(depth_limit, 1);

    eval_limit = options->eval_limit;
    if (eval_limit <= 0)
        eval_limit = 1000 * prec + prec * prec;
    eval_limit = FLINT_MAX(eval_limit, 1);

    goal = FLINT_MAX(goal, 0);
    deg_limit = options->deg_limit;
    if (deg_limit <= 0)
        deg_limit = 0.5 * FLINT_MIN(goal, prec) + 60;

    verbose = options->verbose;
    use_heap = options->use_heap;

    alloc = 4;
    as = _acb_vec_init(alloc);
    bs = _acb_vec_init(alloc);
    vs = _acb_vec_init(alloc);
    ms = _mag_vec_init(alloc);

    /* Compute initial crude estimate for the whole interval. */
    acb_set(as, a);
    acb_set(bs, b);
    quad_simple(vs, f, param, as, bs, prec);
    mag_hypot(ms, arb_radref(acb_realref(vs)), arb_radref(acb_imagref(vs)));

    depth = depth_max = 1;
    eval = 1;
    stopping = 0;
    leaf_interval_count = 0;

    /* Adjust absolute tolerance based on new information. */
    acb_get_mag_lower(tmpm, vs);
    mag_mul_2exp_si(tmpm, tmpm, -goal);
    mag_max(new_tol, tol, tmpm);

    acb_zero(s);

    while (depth >= 1)
    {
        if (stopping == 0 && eval >= eval_limit - 1)
        {
            if (verbose > 0)
                flint_printf("stopping at eval_limit %wd\n", eval_limit);
            status = ARB_CALC_NO_CONVERGENCE;
            stopping = 1;
            continue;
        }

        if (use_heap)
            top = 0;
        else
            top = depth - 1;

        /* We are done with this subinterval. */
        if (mag_cmp(ms + top, new_tol) < 0 ||
            _acb_overlaps(u, as + top, bs + top, prec) || stopping)
        {
            acb_add(s, s, vs + top, prec);
            leaf_interval_count++;

            depth--;
            if (use_heap && depth > 0)
            {
                acb_swap(as, as + depth);
                acb_swap(bs, bs + depth);
                acb_swap(vs, vs + depth);
                mag_swap(ms, ms + depth);
                heap_up(as, bs, vs, ms, depth);
            }
            continue;
        }

        /* Attempt using Gauss-Legendre rule. */
        if (acb_is_finite(vs + top))
        {
            gl_status = acb_calc_integrate_gl_auto_deg(u, &feval, f, param,
                as + top, bs + top, new_tol, deg_limit, verbose > 1, prec);
            eval += feval;

            /* We are done with this subinterval. */
            if (gl_status == ARB_CALC_SUCCESS)
            {
                /* We know that the result is real. */
                real_error = acb_is_finite(vs + top) && acb_is_real(vs + top);

                if (real_error)
                    arb_zero(acb_imagref(u));

                acb_add(s, s, u, prec);
                leaf_interval_count++;

                /* Adjust absolute tolerance based on new information. */
                acb_get_mag_lower(tmpm, u);
                mag_mul_2exp_si(tmpm, tmpm, -goal);
                mag_max(new_tol, new_tol, tmpm);

                depth--;
                if (use_heap && depth > 0)
                {
                    acb_swap(as, as + depth);
                    acb_swap(bs, bs + depth);
                    acb_swap(vs, vs + depth);
                    mag_swap(ms, ms + depth);
                    heap_up(as, bs, vs, ms, depth);
                }
                continue;
            }
        }

        if (depth >= depth_limit - 1)
        {
            if (verbose > 0)
                flint_printf("stopping at depth_limit %wd\n", depth_limit);
            status = ARB_CALC_NO_CONVERGENCE;
            stopping = 1;
            continue;
        }

        if (depth >= alloc - 1)
        {
            slong k;
            as = flint_realloc(as, 2 * alloc * sizeof(acb_struct));
            bs = flint_realloc(bs, 2 * alloc * sizeof(acb_struct));
            vs = flint_realloc(vs, 2 * alloc * sizeof(acb_struct));
            ms = flint_realloc(ms, 2 * alloc * sizeof(mag_struct));
            for (k = alloc; k < 2 * alloc; k++)
            {
                acb_init(as + k);
                acb_init(bs + k);
                acb_init(vs + k);
                mag_init(ms + k);
            }
            alloc *= 2;
        }

        /* Bisection. */
        /* Interval [depth] becomes [mid, b]. */
        acb_set(bs + depth, bs + top);
        acb_add(as + depth, as + top, bs + top, prec);
        acb_mul_2exp_si(as + depth, as + depth, -1);

        /* Interval [top] becomes [a, mid]. */
        acb_set(bs + top, as + depth);

        /* Evaluate on [a, mid] */
        quad_simple(vs + top, f, param, as + top, bs + top, prec);
        mag_hypot(ms + top, arb_radref(acb_realref(vs + top)), arb_radref(acb_imagref(vs + top)));
        eval++;
        /* Adjust absolute tolerance based on new information. */
        acb_get_mag_lower(tmpm, vs + top);
        mag_mul_2exp_si(tmpm, tmpm, -goal);
        mag_max(new_tol, new_tol, tmpm);

        /* Evaluate on [mid, b] */
        quad_simple(vs + depth, f, param, as + depth, bs + depth, prec);
        mag_hypot(ms + depth, arb_radref(acb_realref(vs + depth)), arb_radref(acb_imagref(vs + depth)));
        eval++;
        /* Adjust absolute tolerance based on new information. */
        acb_get_mag_lower(tmpm, vs + depth);
        mag_mul_2exp_si(tmpm, tmpm, -goal);
        mag_max(new_tol, new_tol, tmpm);

        /* Make the interval with the larger error the priority. */
        if (mag_cmp(ms + top, ms + depth) < 0)
        {
            acb_swap(as + top, as + depth);
            acb_swap(bs + top, bs + depth);
            acb_swap(vs + top, vs + depth);
            mag_swap(ms + top, ms + depth);
        }

        if (use_heap)
        {
            heap_up(as, bs, vs, ms, depth);
            heap_down(as, bs, vs, ms, depth + 1);
        }

        depth++;
        depth_max = FLINT_MAX(depth, depth_max);
    }

    if (verbose > 0)
    {
        flint_printf("depth %wd/%wd, eval %wd/%wd, %wd leaf intervals\n",
            depth_max, depth_limit, eval, eval_limit, leaf_interval_count);
    }

    acb_set(res, s);

    _acb_vec_clear(as, alloc);
    _acb_vec_clear(bs, alloc);
    _acb_vec_clear(vs, alloc);
    _mag_vec_clear(ms, alloc);
    acb_clear(s);
    acb_clear(t);
    acb_clear(u);
    mag_clear(tmpm);
    mag_clear(tmpn);
    mag_clear(new_tol);

    return status;
}

