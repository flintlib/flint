/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <float.h>
#include "ulong_extras.h"
#include "arf.h"
#include "arb.h"
#include "acb.h"
#include "gr.h"
#include "gr_vec.h"
#include "thread_pool.h"
#include "thread_support.h"
#include "gr_dft.h"

/* Drop-in replacements for acb_dft and acb_dft_inverse, computing DFTs
   of acb vectors through the generic plans of this module, using
   either acb/arb ball arithmetic or the fixed-point contexts of
   nfixed.c internally.

   In the fixed-point path, only the input midpoints enter the
   transform: they are scaled by an exact power of two chosen from the
   plan's magnitude bound so that no intermediate value can overflow
   the |t| < 1 representation, and truncated to fixed point (error
   below 1 ulp per component). The output radii combine two
   independent contributions:

   1. the fixed-point roundoff and root-table errors, from
      gr_dft_precomp_nfixed_bound applied to the actual input scale
      and the 1 ulp conversion error;

   2. the effect of the input radii, which is not propagated through
      the transform at all but bounded directly by the standard DFT
      perturbation result: since the DFT matrix (and the unscaled
      inverse) has unit-modulus entries,

          |Delta X_k| <= sum_j |Delta x_j|,

      so a single magnitude R = sum_j hypot(rad re_j, rad im_j),
      computed with upward rounding, bounds the output perturbation of
      every component. (This is sharper than pushing intervals through
      the computation, and independent of the fixed-point precision.)

   The number of limbs is chosen from the requested precision plus
   guard bits determined by the error bound of the plan, so that the
   absolute output error from roundoff is below 2^-prec times the
   largest input magnitude; outputs of that scale therefore carry
   roughly prec relative bits. */

/* perturbation bound R = sum_j hypot(rad re_j, rad im_j) */
static int
_gr_dft_acb_ball(acb_ptr w, acb_srcptr v, slong n, int inverse, slong prec)
{
    int status = GR_SUCCESS;
    gr_ctx_t rctx, cctx;
    gr_dft_pre_t P;

    gr_ctx_init_real_arb(rctx, prec);
    gr_ctx_init_complex_acb(cctx, prec);

    status |= gr_dft_precomp_init_karatsuba(P, n, GR_DFT_ALG_AUTO, 0,
            rctx, cctx);

    if (status == GR_SUCCESS)
    {
        if (inverse)
            status |= gr_dft_inverse_precomp((gr_ptr) w, (gr_srcptr) v,
                    P, cctx);
        else
            status |= gr_dft_precomp((gr_ptr) w, (gr_srcptr) v, P, cctx);
    }

    gr_dft_precomp_clear(P);
    gr_ctx_clear(rctx);
    gr_ctx_clear(cctx);

    return status;
}

/* Scaling and limb-count selection shared by the cyclic and product
   fixed-point wrappers, parameterized by the plan's error-bound
   evaluator. On failure the caller owns cleaning up its plan. */
typedef void (*_gr_dft_nfixed_bound_fn)(double * peak, double * err_ulps,
        double in_mag, double in_err, const void * plan);

static int
_gr_dft_acb_select_core(slong * nl_out, slong * p1e_out, double * in_mag_out,
        double * errulps_out, _gr_dft_nfixed_bound_fn bound,
        const void * plan, slong prec)
{
    double peak1, err0, peak, in_mag, errulps;
    slong nl, p1e;

    {
        int p1exp;
        bound(&peak1, &err0, 1.0, 0.0, plan);

        /* the bound of a unit input is at least 1 and, for sane
           transform lengths, far below 2^500; refuse defensively if
           the double computation escaped that range (the comparison
           is written so that NaN also fails). Besides protecting the
           frexp below, the cap guarantees in_mag >= 2^-501, so that
           every intermediate quantity inside the subsequent bound
           evaluations stays comfortably within the normal double
           range: an underflow to zero there would silently understate
           the bound. */
        if (!(peak1 >= 1.0 && peak1 <= ldexp(1.0, 500)))
            return GR_UNABLE;

        frexp(peak1, &p1exp);   /* peak1 < 2^p1exp */
        p1e = p1exp;
    }
    /* p1exp lies in [1, 1024] here, so this neither under- nor
       overflows */
    in_mag = ldexp(1.0, -p1e - 1);

    bound(&peak, &errulps, in_mag, 1.0, plan);
    if (!(peak < 1.0) || !(errulps >= 1.0 && errulps <= DBL_MAX))
        return GR_UNABLE;

    /* guard bits: absolute roundoff errulps 2^(-64 nl) should stay
       below 2^-(prec+2) in_mag */
    nl = (prec + FLINT_BITS - 1) / FLINT_BITS;
    {
        double need = prec + 2.0 + log2(errulps) - log2(in_mag);
        slong nl2 = (slong) (need / FLINT_BITS) + 1;
        nl = FLINT_MAX(nl, nl2);
    }
    if (nl > GR_DFT_NFIXED_MAX_NLIMBS)
        return GR_UNABLE;

    *nl_out = nl;
    *p1e_out = p1e;
    *in_mag_out = in_mag;
    *errulps_out = errulps;
    return GR_SUCCESS;
}

static void
_gr_dft_nfixed_bound_cyc(double * peak, double * err_ulps, double in_mag,
        double in_err, const void * plan)
{
    gr_dft_precomp_nfixed_bound(peak, err_ulps, in_mag, in_err,
            (const gr_dft_pre_struct *) plan);
}

static void
_gr_dft_nfixed_bound_prod(double * peak, double * err_ulps, double in_mag,
        double in_err, const void * plan)
{
    gr_dft_prod_precomp_nfixed_bound(peak, err_ulps, in_mag, in_err,
            (const gr_dft_prod_pre_struct *) plan);
}

static int
_gr_dft_acb_nfixed_select(gr_dft_acb_pre_struct * Q, slong n, slong prec)
{
    int status = GR_SUCCESS;
    double peak, errulps2;

    prec = FLINT_MAX(prec, 2);

    status = _gr_dft_precomp_init_layout(Q->P, n, GR_DFT_ALG_AUTO, 0, 1);
    if (status != GR_SUCCESS)
        return status;

    status = _gr_dft_acb_select_core(&Q->nl, &Q->p1e, &Q->in_mag,
            &Q->errulps, _gr_dft_nfixed_bound_cyc, Q->P, prec);
    if (status != GR_SUCCESS)
    {
        gr_dft_precomp_clear(Q->P);
        return status;
    }

    status |= gr_dft_ctx_init_nfixed(Q->rctx, Q->nl);
    status |= gr_dft_ctx_init_nfixed_complex(Q->cctx, Q->nl);
    status |= _gr_dft_precomp_realize(Q->P, Q->rctx, Q->cctx);
    if (status != GR_SUCCESS)
    {
        /* realize cleared the plan already */
        gr_ctx_clear(Q->rctx);
        gr_ctx_clear(Q->cctx);
        return status;
    }

    /* re-evaluate the error bound with the realized plan (actual
       complex-multiplication constant and root-table error), for
       tighter output radii */
    gr_dft_precomp_nfixed_bound(&peak, &errulps2, Q->in_mag, 1.0, Q->P);
    Q->errulps = errulps2;
    if (!(peak < 1.0) ||
        !(Q->errulps >= 1.0 && Q->errulps <= DBL_MAX))
    {
        gr_dft_precomp_clear(Q->P);
        gr_ctx_clear(Q->rctx);
        gr_ctx_clear(Q->cctx);
        return GR_UNABLE;
    }

    return GR_SUCCESS;
}

/* the per-call part of the fixed-point transform: scale, convert,
   transform, assemble. Returns GR_UNABLE for inputs the fixed-point
   path cannot handle (non-finite or extreme midpoints). */
/* Fast conversions between arf midpoints and nfixed components
   (sign limb + nl fraction limbs), following the pattern of
   _nfloat_get_nfixed / _nfloat_set_nfixed: the input direction
   extracts the integer part of |mid| 2^(64 nl - e) directly into the
   fraction limbs with _arf_get_integer_mpn (truncation toward zero,
   error below 1 ulp, matching gr_dft_nfixed_set_arf), and the output
   direction builds the rounded arf at the target precision in one
   step with _arf_set_mpn_fixed, so no temporary arf, no separate
   rounding pass, and no arb_set_arf copy are needed. */

static void
_acb_dft_arf_to_nfixed(nn_ptr res, arf_srcptr mid, slong e, slong nl)
{
    slong rel;
    nn_srcptr xp;
    slong xn;

    /* res is already zeroed (fresh gr vector) */
    if (ARF_IS_SPECIAL(mid))
        return;                     /* finite by the earlier scan: zero */

    if (COEFF_IS_MPZ(ARF_EXP(mid)))
        return;                     /* a bignum exponent must be hugely
                                       negative here (a positive one
                                       would contradict the caller's
                                       magnitude bound EM, which is
                                       checked to be far within slong
                                       range): the value underflows to
                                       zero, with error below 1 ulp */

    rel = ARF_EXP(mid) - e;         /* < 0 by choice of e */
    if (nl * FLINT_BITS + rel <= 0)
        return;                     /* underflow */

    res[0] = ARF_SGNBIT(mid);
    ARF_GET_MPN_READONLY(xp, xn, mid);
    _arf_get_integer_mpn(res + 1, xp, xn, nl * FLINT_BITS + rel);
}

/* part <- (component) 2^e / den, rounded to prec, with the rounding
   error added to the radius (which the caller has preset); den = 0
   means no division */
static void
_acb_dft_nfixed_to_arb(arb_ptr part, nn_srcptr comp, slong e, slong nl,
        ulong den, slong prec)
{
    arf_ptr mid = arb_midref(part);
    int inexact;

    if (den != 0)
    {
        /* convert exactly (the fixed-point value has at most 64 nl
           bits) so that the division performs the single rounding */
        _arf_set_mpn_fixed(mid, comp + 1, nl, nl,
                (int) comp[0], nl * FLINT_BITS, ARF_RND_DOWN);
        arf_mul_2exp_si(mid, mid, e);
        inexact = arf_div_ui(mid, mid, den, prec, ARF_RND_DOWN);
    }
    else
    {
        inexact = _arf_set_mpn_fixed(mid, comp + 1, nl, nl,
                (int) comp[0], prec, ARF_RND_DOWN);
        arf_mul_2exp_si(mid, mid, e);
    }

    if (inexact)
        arf_mag_add_ulp(arb_radref(part), arb_radref(part), mid, prec);
}

/* Threaded input/output conversion: at low precision and large n the
   linear passes cost as much as the transform itself, and once the
   transform is threaded a serial conversion would dominate the wall
   time; both loops are embarrassingly parallel over disjoint
   components. */

#define GR_DFT_ACB_CONV_SERIAL 8192

typedef struct
{
    acb_srcptr v;
    nn_ptr x;
    acb_ptr w;
    nn_srcptr y;
    slong lo;
    slong hi;             /* component (not coefficient) range */
    slong e;
    slong nl;
    slong rsz;
    ulong den;
    slong prec;
    const mag_struct * errmag;
}
_acb_dft_conv_work_t;

static void
_acb_dft_conv_in_worker(void * arg)
{
    _acb_dft_conv_work_t * a = arg;
    slong j;

    for (j = a->lo; j < a->hi; j++)
    {
        arf_srcptr mid = (j % 2)
            ? arb_midref(acb_imagref(a->v + j / 2))
            : arb_midref(acb_realref(a->v + j / 2));

        _acb_dft_arf_to_nfixed((nn_ptr) ((char *) a->x + j * a->rsz),
                mid, a->e, a->nl);
    }
}

static void
_acb_dft_conv_out_worker(void * arg)
{
    _acb_dft_conv_work_t * a = arg;
    slong j;

    for (j = a->lo; j < a->hi; j++)
    {
        arb_ptr part = (j % 2)
            ? acb_imagref(a->w + j / 2)
            : acb_realref(a->w + j / 2);

        mag_set(arb_radref(part), a->errmag);
        _acb_dft_nfixed_to_arb(part,
                (nn_srcptr) ((const char *) a->y + j * a->rsz),
                a->e, a->nl, a->den, a->prec);
    }
}

static void
_acb_dft_conv_run(int out, acb_srcptr v, nn_ptr x, acb_ptr w, nn_srcptr y,
        slong ncomp, slong e, slong nl, slong rsz, ulong den, slong prec,
        const mag_struct * errmag)
{
    thread_pool_handle * handles = NULL;
    slong num_workers = 0, nchunks, i;
    _acb_dft_conv_work_t args[16];

    if (ncomp >= 2 * GR_DFT_ACB_CONV_SERIAL && flint_get_num_threads() > 1)
        num_workers = flint_request_threads(&handles,
                FLINT_MIN(flint_get_num_threads(), 16));

    nchunks = num_workers + 1;

    for (i = 0; i < nchunks; i++)
    {
        args[i].v = v; args[i].x = x; args[i].w = w; args[i].y = y;
        args[i].lo = (ncomp * i) / nchunks;
        args[i].hi = (ncomp * (i + 1)) / nchunks;
        args[i].e = e; args[i].nl = nl; args[i].rsz = rsz;
        args[i].den = den; args[i].prec = prec; args[i].errmag = errmag;

        if (i < num_workers)
            thread_pool_wake(global_thread_pool, handles[i], 0,
                    out ? _acb_dft_conv_out_worker : _acb_dft_conv_in_worker,
                    &args[i]);
    }

    if (out)
        _acb_dft_conv_out_worker(&args[num_workers]);
    else
        _acb_dft_conv_in_worker(&args[num_workers]);

    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_give_back_threads(handles, num_workers);
}

/* Parallel input scan: per-worker partial reductions (finiteness
   flag, radius row sum, pointer to the largest-abs midpoint),
   combined by the caller. */
typedef struct
{
    acb_srcptr v;
    slong lo;
    slong hi;             /* coefficient range */
    arf_srcptr best;
    mag_t R;
    int finite;
}
_acb_dft_scan_work_t;

static void
_acb_dft_scan_worker(void * arg)
{
    _acb_dft_scan_work_t * a = arg;
    arf_srcptr best = arb_midref(acb_realref(a->v + a->lo));
    slong j;

    a->finite = 1;

    for (j = a->lo; j < a->hi; j++)
    {
        arf_srcptr re = arb_midref(acb_realref(a->v + j));
        arf_srcptr im = arb_midref(acb_imagref(a->v + j));

        if (!arf_is_finite(re) || !arf_is_finite(im))
        {
            a->finite = 0;
            return;
        }

        mag_add(a->R, a->R, arb_radref(acb_realref(a->v + j)));
        mag_add(a->R, a->R, arb_radref(acb_imagref(a->v + j)));

        if (arf_cmpabs(re, best) > 0)
            best = re;
        if (arf_cmpabs(im, best) > 0)
            best = im;
    }

    a->best = best;
}

/* executor: runs the chunks whose index is congruent to its id */
typedef struct
{
    _acb_dft_scan_work_t * chunks;
    slong nchunks;
    slong id;
    slong stride;
}
_acb_dft_scan_exec_t;

static void
_acb_dft_scan_exec(void * arg)
{
    _acb_dft_scan_exec_t * a = arg;
    slong i;

    for (i = a->id; i < a->nchunks; i += a->stride)
        _acb_dft_scan_worker(&a->chunks[i]);
}

/* Returns 0 if any coefficient is non-finite; otherwise R and t
   receive the radius row sum and the largest-abs midpoint. The
   reduction uses a fixed 16-chunk grid combined in canonical order,
   so the result (in particular the upward-rounded mag sum, which is
   not associative) is independent of the number of threads. */
static int
_acb_dft_scan(mag_t R, arf_t t, acb_srcptr v, slong n)
{
    thread_pool_handle * handles = NULL;
    slong num_workers = 0, nexec, nchunks, i;
    _acb_dft_scan_work_t args[16];
    _acb_dft_scan_exec_t execs[16];
    int finite = 1;

    nchunks = (n >= GR_DFT_ACB_CONV_SERIAL) ? 16 : 1;

    for (i = 0; i < nchunks; i++)
    {
        args[i].v = v;
        args[i].lo = (n * i) / nchunks;
        args[i].hi = (n * (i + 1)) / nchunks;
        args[i].best = NULL;
        mag_init(args[i].R);
    }

    if (nchunks > 1 && flint_get_num_threads() > 1)
        num_workers = flint_request_threads(&handles,
                FLINT_MIN(flint_get_num_threads(), nchunks));

    nexec = num_workers + 1;

    for (i = 0; i < nexec; i++)
    {
        execs[i].chunks = args;
        execs[i].nchunks = nchunks;
        execs[i].id = i;
        execs[i].stride = nexec;

        if (i < num_workers)
            thread_pool_wake(global_thread_pool, handles[i], 0,
                    _acb_dft_scan_exec, &execs[i]);
    }

    _acb_dft_scan_exec(&execs[num_workers]);

    for (i = 0; i < num_workers; i++)
        thread_pool_wait(global_thread_pool, handles[i]);

    flint_give_back_threads(handles, num_workers);

    arf_zero(t);
    for (i = 0; i < nchunks; i++)
    {
        if (!args[i].finite)
            finite = 0;
        else
        {
            mag_add(R, R, args[i].R);
            if (arf_cmpabs(args[i].best, t) > 0)
                arf_set(t, args[i].best);
        }
        mag_clear(args[i].R);
    }

    return finite;
}

static int
_gr_dft_acb_nfixed_run(acb_ptr w, acb_srcptr v, int inverse,
        slong n, slong nl, slong p1e, double errulps,
        gr_ctx_struct * rctx, gr_ctx_struct * cctx,
        int (*transform)(gr_ptr, gr_srcptr, int, const void *, gr_ctx_t),
        const void * plan, slong prec)
{
    int status = GR_SUCCESS;
    gr_ptr x, y;
    mag_t R, errmag;
    arf_t t;
    slong j, e, EM, rsz;

    mag_init(R);
    mag_init(errmag);
    arf_init(t);

    /* one (parallel) pass over the input: finiteness, radius row
       sum, and the largest-abs midpoint */
    if (!_acb_dft_scan(R, t, v, n))
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    if (arf_is_zero(t))
    {
        /* all midpoints zero: output midpoints are zero, radii come
           from the perturbation bound alone */
        if (inverse)
            mag_div_ui(R, R, n);
        for (j = 0; j < n; j++)
        {
            acb_zero(w + j);
            mag_set(arb_radref(acb_realref(w + j)), R);
            mag_set(arb_radref(acb_imagref(w + j)), R);
        }
        status = GR_SUCCESS;
        goto cleanup;
    }

    /* arf_abs_bound_lt_2exp_si saturates at +-ARF_PREC_EXACT for
       out-of-slong (bignum) exponents, so extreme and adversarial
       magnitudes are guaranteed to fail this test and be routed to
       ball arithmetic. The margin also ensures that e below, the
       relative exponents in the input conversion, and the shifts
       applied to the error magnitude all stay far from slong
       overflow, on 32-bit as well as 64-bit systems. It further
       makes the input conversion's treatment of individual bignum
       exponents safe: with EM in range, a coefficient whose
       exponent does not fit in an slong can only be tiny (a huge
       coefficient would contradict the bound EM), so truncating it
       to zero is correct. */
    EM = arf_abs_bound_lt_2exp_si(t);
    if (FLINT_ABS(EM) > WORD_MAX / 4)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    /* input scale 2^-e with (m 2^-e) peak(1) <= 1/2 */
    e = EM + p1e + 1;

    rsz = rctx->sizeof_elem;

    /* fixed-point elements are plain limb blocks: allocate directly
       (x zeroed, since the input conversion leaves zero and
       underflowing components untouched; y is fully written by the
       transform) */
    x = flint_calloc(2 * n, rsz);
    y = flint_malloc(2 * (size_t) n * rsz);

    /* scale (exactly) and truncate the midpoints, writing the
       fraction limbs directly */
    _acb_dft_conv_run(0, v, (nn_ptr) x, NULL, NULL, 2 * n, e, nl,
            rsz, 0, prec, NULL);

    status |= transform(y, x, inverse, plan, cctx);

    /* absolute roundoff error, back at the input scale */
    mag_set_d(errmag, errulps);
    mag_mul_2exp_si(errmag, errmag, -nl * FLINT_BITS + e);
    mag_add(errmag, errmag, R);

    /* the raw inverse omits the 1/n normalization: for power-of-two
       lengths it is folded into the output exponent for free, and
       otherwise performed by a single rounded division per component
       during the output conversion; the shared error magnitude is
       divided once */
    {
        slong eout = e;
        ulong den = 0;

        if (inverse)
        {
            if ((n & (n - 1)) == 0)
            {
                slong k = FLINT_BIT_COUNT(n) - 1;
                eout = e - k;
                mag_mul_2exp_si(errmag, errmag, -k);
            }
            else
            {
                den = n;
                mag_div_ui(errmag, errmag, n);
            }
        }

        _acb_dft_conv_run(1, NULL, NULL, w, (nn_srcptr) y, 2 * n, eout,
                nl, rsz, den, prec, errmag);
    }

    flint_free(x);
    flint_free(y);

cleanup:
    mag_clear(R);
    mag_clear(errmag);
    arf_clear(t);

    return status;
}


static int
_gr_dft_acb_transform_cyc(gr_ptr y, gr_srcptr x, int inverse,
        const void * plan, gr_ctx_t cctx)
{
    const gr_dft_pre_struct * P = plan;

    if (inverse)
        return _gr_dft_precomp_raw(y, x, 1, P, cctx);
    else
        return gr_dft_precomp(y, x, P, cctx);
}

static int
_gr_dft_acb_transform_prod(gr_ptr y, gr_srcptr x, int inverse,
        const void * plan, gr_ctx_t cctx)
{
    const gr_dft_prod_pre_struct * P = plan;

    return _gr_dft_prod_precomp_raw(y, x, inverse, P, cctx);
}

static int
_gr_dft_acb_nfixed_apply(acb_ptr w, acb_srcptr v, int inverse,
        const gr_dft_acb_pre_struct * Q, slong prec)
{
    return _gr_dft_acb_nfixed_run(w, v, inverse, Q->n, Q->nl, Q->p1e,
            Q->errulps, (gr_ctx_struct *) Q->rctx, (gr_ctx_struct *) Q->cctx,
            _gr_dft_acb_transform_cyc, Q->P, prec);
}

static int
_gr_dft_acb_nfixed(acb_ptr w, acb_srcptr v, slong n, int inverse, slong prec)
{
    int status;
    gr_dft_acb_pre_struct Q;

    Q.n = n;
    Q.prec = prec;
    Q.which = 2;

    status = _gr_dft_acb_nfixed_select(&Q, n, prec);
    if (status != GR_SUCCESS)
        return status;

    status = _gr_dft_acb_nfixed_apply(w, v, inverse, &Q, prec);

    gr_dft_precomp_clear(Q.P);
    gr_ctx_clear(Q.rctx);
    gr_ctx_clear(Q.cctx);

    return status;
}

/* Precomputation interface *****************************************/

int
gr_dft_acb_precomp_init(gr_dft_acb_pre_t Q, slong n, slong prec)
{
    int status = GR_SUCCESS;
    int use_nfixed = 1;

    if (n <= 0)
        return GR_DOMAIN;

    prec = FLINT_MAX(prec, 2);

    Q->n = n;
    Q->prec = prec;

    if (use_nfixed)
    {
        status = _gr_dft_acb_nfixed_select(Q, n, prec);
        if (status == GR_SUCCESS)
        {
            Q->which = 2;
            return GR_SUCCESS;
        }
    }

    /* ball arithmetic */
    gr_ctx_init_real_arb(Q->rctx, prec);
    gr_ctx_init_complex_acb(Q->cctx, prec);
    status = gr_dft_precomp_init_karatsuba(Q->P, n, GR_DFT_ALG_AUTO, 0,
            Q->rctx, Q->cctx);
    if (status != GR_SUCCESS)
    {
        gr_dft_precomp_clear(Q->P);
        gr_ctx_clear(Q->rctx);
        gr_ctx_clear(Q->cctx);
        Q->which = 0;
        return status;
    }
    Q->which = 1;

    return GR_SUCCESS;
}

void
gr_dft_acb_precomp_clear(gr_dft_acb_pre_t Q)
{
    if (Q->which != 0)
    {
        gr_dft_precomp_clear(Q->P);
        gr_ctx_clear(Q->rctx);
        gr_ctx_clear(Q->cctx);
        Q->which = 0;
    }
}

int
_gr_dft_acb_precomp(acb_ptr w, acb_srcptr v, int inverse,
        const gr_dft_acb_pre_t Q, slong prec)
{
    int status;
    slong n = Q->n;

    prec = FLINT_MAX(prec, 2);

    if (Q->which == 2)
    {
        status = _gr_dft_acb_nfixed_apply(w, v, inverse, Q, prec);

        /* inputs the fixed-point path cannot handle: fall back to a
           one-shot ball transform */
        if (status != GR_SUCCESS)
            status = _gr_dft_acb_ball(w, v, n, inverse, prec);
    }
    else if (Q->which == 1)
    {
        if (inverse)
            status = gr_dft_inverse_precomp((gr_ptr) w, (gr_srcptr) v,
                    Q->P, (gr_ctx_struct *) Q->cctx);
        else
            status = gr_dft_precomp((gr_ptr) w, (gr_srcptr) v,
                    Q->P, (gr_ctx_struct *) Q->cctx);
    }
    else
    {
        status = GR_UNABLE;
    }

    if (status != GR_SUCCESS)
    {
        slong j;
        for (j = 0; j < n; j++)
            acb_indeterminate(w + j);
    }

    return status;
}

void
gr_dft_acb_precomp(acb_ptr w, acb_srcptr v, const gr_dft_acb_pre_t Q,
        slong prec)
{
    GR_IGNORE(_gr_dft_acb_precomp(w, v, 0, Q, prec));
}

void
gr_dft_acb_inverse_precomp(acb_ptr w, acb_srcptr v,
        const gr_dft_acb_pre_t Q, slong prec)
{
    GR_IGNORE(_gr_dft_acb_precomp(w, v, 1, Q, prec));
}

/* which: 0 = automatic, 1 = force acb/arb ball arithmetic internally,
   2 = force fixed-point arithmetic internally */
int
_gr_dft_acb(acb_ptr w, acb_srcptr v, slong n, int inverse, int which,
        slong prec)
{
    int status;

    if (n <= 0)
        return (n == 0) ? GR_SUCCESS : GR_DOMAIN;

    prec = FLINT_MAX(prec, 2);

    if (which == 1)
        return _gr_dft_acb_ball(w, v, n, inverse, prec);

    status = _gr_dft_acb_nfixed(w, v, n, inverse, prec);

    if (status != GR_SUCCESS && which == 0)
        status = _gr_dft_acb_ball(w, v, n, inverse, prec);

    if (status != GR_SUCCESS)
    {
        slong j;
        for (j = 0; j < n; j++)
            acb_indeterminate(w + j);
    }

    return status;
}

void
gr_dft_acb(acb_ptr w, acb_srcptr v, slong n, slong prec)
{
    if (n > 0)
        GR_IGNORE(_gr_dft_acb(w, v, n, 0, 0, prec));
}

void
gr_dft_acb_inverse(acb_ptr w, acb_srcptr v, slong n, slong prec)
{
    if (n > 0)
        GR_IGNORE(_gr_dft_acb(w, v, n, 1, 0, prec));
}

/* Product-group variant *******************************************/

static int
_gr_dft_acb_prod_ball(acb_ptr w, acb_srcptr v, const ulong * cyc, slong num,
        int inverse, slong prec)
{
    int status;
    gr_ctx_t ctx;

    gr_ctx_init_complex_acb(ctx, prec);
    if (inverse)
        status = gr_dft_prod_inverse((gr_ptr) w, (gr_srcptr) v, cyc, num, ctx);
    else
        status = gr_dft_prod((gr_ptr) w, (gr_srcptr) v, cyc, num, ctx);
    gr_ctx_clear(ctx);

    return status;
}

static int
_gr_dft_acb_prod_nfixed_select(gr_dft_acb_prod_pre_struct * Q,
        const ulong * cyc, slong num, slong prec)
{
    int status = GR_SUCCESS;
    double peak, errulps2;

    prec = FLINT_MAX(prec, 2);

    status = _gr_dft_prod_precomp_init_layout(Q->P, cyc, num, 0, 1);
    if (status != GR_SUCCESS)
    {
        gr_dft_prod_precomp_clear(Q->P);
        return status;
    }

    status = _gr_dft_acb_select_core(&Q->nl, &Q->p1e, &Q->in_mag,
            &Q->errulps, _gr_dft_nfixed_bound_prod, Q->P, prec);
    if (status != GR_SUCCESS)
    {
        gr_dft_prod_precomp_clear(Q->P);
        return status;
    }

    status |= gr_dft_ctx_init_nfixed(Q->rctx, Q->nl);
    status |= gr_dft_ctx_init_nfixed_complex(Q->cctx, Q->nl);
    status |= _gr_dft_prod_precomp_realize(Q->P, Q->rctx, Q->cctx);
    if (status != GR_SUCCESS)
    {
        gr_dft_prod_precomp_clear(Q->P);
        gr_ctx_clear(Q->rctx);
        gr_ctx_clear(Q->cctx);
        return status;
    }

    /* re-evaluate with the realized component plans */
    gr_dft_prod_precomp_nfixed_bound(&peak, &errulps2, Q->in_mag, 1.0, Q->P);
    Q->errulps = errulps2;
    if (!(peak < 1.0) ||
        !(Q->errulps >= 1.0 && Q->errulps <= DBL_MAX))
    {
        gr_dft_prod_precomp_clear(Q->P);
        gr_ctx_clear(Q->rctx);
        gr_ctx_clear(Q->cctx);
        return GR_UNABLE;
    }

    return GR_SUCCESS;
}

static int
_gr_dft_acb_prod_nfixed_apply(acb_ptr w, acb_srcptr v, int inverse,
        const gr_dft_acb_prod_pre_struct * Q, slong prec)
{
    return _gr_dft_acb_nfixed_run(w, v, inverse, Q->n, Q->nl, Q->p1e,
            Q->errulps, (gr_ctx_struct *) Q->rctx, (gr_ctx_struct *) Q->cctx,
            _gr_dft_acb_transform_prod, Q->P, prec);
}

int
gr_dft_acb_prod_precomp_init(gr_dft_acb_prod_pre_t Q, const ulong * cyc,
        slong num, slong prec)
{
    int status = GR_SUCCESS;
    slong a, n = 1;

    for (a = 0; a < num; a++)
    {
        if (cyc[a] == 0 || (ulong) n > (WORD_MAX / 8) / cyc[a])
            return GR_DOMAIN;
        n *= (slong) cyc[a];
    }

    prec = FLINT_MAX(prec, 2);

    Q->n = n;
    Q->num = num;
    Q->prec = prec;
    Q->cyc = flint_malloc(FLINT_MAX(num, 1) * sizeof(ulong));
    for (a = 0; a < num; a++)
        Q->cyc[a] = cyc[a];

    status = _gr_dft_acb_prod_nfixed_select(Q, cyc, num, prec);
    if (status == GR_SUCCESS)
    {
        Q->which = 2;
        return GR_SUCCESS;
    }

    /* ball arithmetic */
    gr_ctx_init_real_arb(Q->rctx, prec);
    gr_ctx_init_complex_acb(Q->cctx, prec);
    status = gr_dft_prod_precomp_init(Q->P, cyc, num, 0,
            (gr_ctx_struct *) Q->cctx);
    if (status != GR_SUCCESS)
    {
        gr_dft_prod_precomp_clear(Q->P);
        gr_ctx_clear(Q->rctx);
        gr_ctx_clear(Q->cctx);
        flint_free(Q->cyc);
        Q->cyc = NULL;
        Q->which = 0;
        return status;
    }
    Q->which = 1;

    return GR_SUCCESS;
}

void
gr_dft_acb_prod_precomp_clear(gr_dft_acb_prod_pre_t Q)
{
    if (Q->which != 0)
    {
        gr_dft_prod_precomp_clear(Q->P);
        gr_ctx_clear(Q->rctx);
        gr_ctx_clear(Q->cctx);
        Q->which = 0;
    }
    flint_free(Q->cyc);
    Q->cyc = NULL;
}

int
_gr_dft_acb_prod_precomp(acb_ptr w, acb_srcptr v, int inverse,
        const gr_dft_acb_prod_pre_t Q, slong prec)
{
    int status;
    slong n = Q->n;

    prec = FLINT_MAX(prec, 2);

    if (Q->which == 2)
    {
        status = _gr_dft_acb_prod_nfixed_apply(w, v, inverse, Q, prec);

        /* inputs the fixed-point path cannot handle: fall back to a
           one-shot ball transform */
        if (status != GR_SUCCESS)
            status = _gr_dft_acb_prod_ball(w, v, Q->cyc, Q->num, inverse,
                    prec);
    }
    else if (Q->which == 1)
    {
        if (inverse)
            status = gr_dft_prod_inverse_precomp((gr_ptr) w, (gr_srcptr) v,
                    Q->P, (gr_ctx_struct *) Q->cctx);
        else
            status = gr_dft_prod_precomp((gr_ptr) w, (gr_srcptr) v,
                    Q->P, (gr_ctx_struct *) Q->cctx);
    }
    else
    {
        status = GR_UNABLE;
    }

    if (status != GR_SUCCESS)
    {
        slong j;
        for (j = 0; j < n; j++)
            acb_indeterminate(w + j);
    }

    return status;
}

void
gr_dft_acb_prod_precomp(acb_ptr w, acb_srcptr v,
        const gr_dft_acb_prod_pre_t Q, slong prec)
{
    GR_IGNORE(_gr_dft_acb_prod_precomp(w, v, 0, Q, prec));
}

void
gr_dft_acb_prod_inverse_precomp(acb_ptr w, acb_srcptr v,
        const gr_dft_acb_prod_pre_t Q, slong prec)
{
    GR_IGNORE(_gr_dft_acb_prod_precomp(w, v, 1, Q, prec));
}

int
_gr_dft_acb_prod(acb_ptr w, acb_srcptr v, const ulong * cyc, slong num,
        int inverse, int which, slong prec)
{
    int status = GR_SUCCESS;
    gr_dft_acb_prod_pre_struct Q;
    slong a, n = 1;

    for (a = 0; a < num; a++)
    {
        if (cyc[a] == 0 || (ulong) n > (WORD_MAX / 8) / cyc[a])
            return GR_DOMAIN;
        n *= (slong) cyc[a];
    }

    prec = FLINT_MAX(prec, 2);

    if (which == 0 || which == 2)
    {
        Q.n = n;
        Q.num = num;
        Q.prec = prec;
        Q.which = 2;

        status = _gr_dft_acb_prod_nfixed_select(&Q, cyc, num, prec);
        if (status == GR_SUCCESS)
        {
            status = _gr_dft_acb_prod_nfixed_apply(w, v, inverse, &Q, prec);

            gr_dft_prod_precomp_clear(Q.P);
            gr_ctx_clear(Q.rctx);
            gr_ctx_clear(Q.cctx);
        }

        if (status == GR_SUCCESS || which == 2)
        {
            if (status != GR_SUCCESS)
            {
                slong j;
                for (j = 0; j < n; j++)
                    acb_indeterminate(w + j);
            }
            return status;
        }
    }

    status = _gr_dft_acb_prod_ball(w, v, cyc, num, inverse, prec);
    if (status != GR_SUCCESS)
    {
        slong j;
        for (j = 0; j < n; j++)
            acb_indeterminate(w + j);
    }

    return status;
}

void
gr_dft_acb_prod(acb_ptr w, acb_srcptr v, const ulong * cyc, slong num,
        slong prec)
{
    GR_IGNORE(_gr_dft_acb_prod(w, v, cyc, num, 0, 0, prec));
}

void
gr_dft_acb_prod_inverse(acb_ptr w, acb_srcptr v, const ulong * cyc,
        slong num, slong prec)
{
    GR_IGNORE(_gr_dft_acb_prod(w, v, cyc, num, 1, 0, prec));
}
