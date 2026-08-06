/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "acb.h"
#include "acb_dirichlet.h"
#include "arb_poly.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "secondary_zeta.h"

/*
    Secondary zeta function

        Z(s) = sum_{n>=1} alpha_n^{-s},   rho_n = 1/2 + i alpha_n,

    following J. Arias de Reyna, "Computation of the secondary zeta
    function", arXiv:2006.04869. With a free parameter a > 0,

        Z(s) = A(s) - P(s) + E(s) - S(s).

    The error radii for A, P and S implement unconditional bounds derived
    by Claude: the nonrigorous tail estimates for A and S in Arias de Reyna
    have been replaced by rigorous bounds, and the P bound is significantly
    tightened. (The E term is trivial, being a hypergeometric series.)

    The zeros used by A() are obtained from acb_dirichlet_hardy_z_zeros
    (rigorously isolated); for heights below the
    Platt-Trudgian verification bound 3e12 they are simple and on the line.
*/

#define SECOND_ZETA_RH_VERIFIED_HEIGHT 3e12

/* Semi-private helpers */
void _acb_dirichlet_secondary_zeta_term_E(acb_t res, const acb_t s,
        const arb_t a, slong prec);
void _acb_dirichlet_secondary_zeta_term_S(acb_t res, const acb_t s,
        const arb_t a, slong prec);
void _acb_dirichlet_secondary_zeta_term_P(acb_t res, const acb_t s,
        const arb_t a, slong prec);
void _acb_dirichlet_secondary_zeta_term_A(acb_t res, const acb_t s,
        const arb_t a, slong prec);
void _acb_dirichlet_secondary_zeta_zeros(arb_ptr res, slong num, slong prec);


#define SECOND_ZETA_BOUND_PREC MAG_BITS

/*****************************************************************************/
/*                                Zeta zeros                                 */
/*****************************************************************************/

/* The secondary zeta function repeatedly needs gamma_1,...,gamma_N at some
   working precision; recomputing them via acb_dirichlet_hardy_z_zeros (which
   rigorously isolates each zero) dominates the cost so a cache is crucial.

   TODO: consider whether we want to expose this cache for other uses
   TODO: consider tracking precision per zero.
   TODO: implement refinement of existing zeros (much faster than isolating
         from scratch)
*/

FLINT_TLS_PREFIX slong _second_zeta_zeros_cached_num = 0;
FLINT_TLS_PREFIX slong _second_zeta_zeros_cached_prec = 0;
FLINT_TLS_PREFIX slong _second_zeta_zeros_alloc = 0;
FLINT_TLS_PREFIX arb_ptr _second_zeta_zeros_cache = NULL;

static void
_second_zeta_zeros_cleanup(void)
{
    if (_second_zeta_zeros_alloc != 0)
    {
        _arb_vec_clear(_second_zeta_zeros_cache, _second_zeta_zeros_alloc);
        _second_zeta_zeros_cache = NULL;
        _second_zeta_zeros_alloc = 0;
    }
    _second_zeta_zeros_cached_num = 0;
    _second_zeta_zeros_cached_prec = 0;
}

/* ensure the cache holds >= num ordinates, each to >= prec bits */
static void
_second_zeta_zeros_ensure(slong num, slong prec)
{
    fmpz_t n0;

    if (num <= _second_zeta_zeros_cached_num
            && prec <= _second_zeta_zeros_cached_prec)
        return;

    if (_second_zeta_zeros_alloc == 0)
        flint_register_cleanup_function(_second_zeta_zeros_cleanup);

    if (prec > _second_zeta_zeros_cached_prec)
    {
        /* Precision insufficient: recompute the whole (possibly larger)
          vector at the requested precision. To avoid repeated near-no-op
          recomputations from a sequence of tiny precision bumps, grow the
          cache precision by at least 1.25x the previous value. */
        slong newnum = FLINT_MAX(num, _second_zeta_zeros_cached_num);
        slong newprec = FLINT_MAX(prec,
                _second_zeta_zeros_cached_prec
                + _second_zeta_zeros_cached_prec / 4);

        if (newnum > _second_zeta_zeros_alloc)
        {
            if (_second_zeta_zeros_alloc != 0)
                _arb_vec_clear(_second_zeta_zeros_cache, _second_zeta_zeros_alloc);
            _second_zeta_zeros_cache = _arb_vec_init(newnum);
            _second_zeta_zeros_alloc = newnum;
        }

        fmpz_init(n0);
        fmpz_one(n0);   /* start at zero index 1 */
        acb_dirichlet_hardy_z_zeros(_second_zeta_zeros_cache, n0, newnum, newprec);
        fmpz_clear(n0);

        _second_zeta_zeros_cached_num = newnum;
        _second_zeta_zeros_cached_prec = newprec;
    }
    else
    {
        /* precision ok, just need more entries: append cached_num+1 .. num */
        slong have = _second_zeta_zeros_cached_num;
        slong want = num;

        if (want > _second_zeta_zeros_alloc)
        {
            slong i;
            _second_zeta_zeros_cache = flint_realloc(_second_zeta_zeros_cache,
                    want * sizeof(arb_struct));
            for (i = _second_zeta_zeros_alloc; i < want; i++)
                arb_init(_second_zeta_zeros_cache + i);
            _second_zeta_zeros_alloc = want;
        }

        fmpz_init(n0);
        fmpz_set_si(n0, have + 1);    /* next zero index */
        acb_dirichlet_hardy_z_zeros(_second_zeta_zeros_cache + have, n0,
                want - have, _second_zeta_zeros_cached_prec);
        fmpz_clear(n0);

        _second_zeta_zeros_cached_num = want;
    }
}

void
_acb_dirichlet_secondary_zeta_zeros(arb_ptr res, slong num, slong prec)
{
    slong i;

    if (num <= 0)
        return;

    _second_zeta_zeros_ensure(num, prec);

    for (i = 0; i < num; i++)
        arb_set_round(res + i, _second_zeta_zeros_cache + i, prec);
}

/*****************************************************************************/
/*                                  E term                                   */
/*****************************************************************************/

/*
    Exponential term  (Arias de Reyna eq. 47)

        E(s) = 1/Gamma(s/2) * sum_{n>=0} 1/(4^n n!) * a^{n+s/2}/(n+s/2) =
             = a^{s/2} * 1F1(s/2; s/2+1; a/4) / Gamma(s/2+1)
             = a^{s/2} * M(s/2, s/2+1, a/4),
*/

void
_acb_dirichlet_secondary_zeta_term_E(acb_t res, const acb_t s,
        const arb_t a, slong prec)
{
    acb_t shalf, b, z, apow;
    slong wp = prec + 8;

    acb_init(shalf); acb_init(b); acb_init(z); acb_init(apow);

    acb_mul_2exp_si(shalf, s, -1);          /* s/2 */
    acb_add_ui(b, shalf, 1, wp);            /* s/2 + 1 */

    acb_set_arb(z, a);
    acb_mul_2exp_si(z, z, -2);              /* a/4 */

    /* M(s/2, s/2+1, a/4) = 1F1(s/2; s/2+1; a/4) / Gamma(s/2+1) */
    acb_hypgeom_m(res, shalf, b, z, 1, wp);

    /* a^{s/2} */
    acb_set_arb(apow, a);
    acb_pow(apow, apow, shalf, wp);

    acb_mul(res, res, apow, wp);

    acb_set_round(res, res, prec);

    acb_clear(shalf); acb_clear(b); acb_clear(z); acb_clear(apow);
}

/*****************************************************************************/
/*                                  A term                                   */
/*****************************************************************************/

/*
    Zero term  (Arias de Reyna eq. 28)

      A(s) = sum_{n>=1} Gamma(s/2, a alpha_n^2)/Gamma(s/2) * alpha_n^{-s},

    where rho_n = 1/2 + i alpha_n are the zeros with Im rho_n > 0. For heights
    below 3e12 (Platt-Trudgian) the alpha_n = gamma_n are real (RH verified),
    so we sum with the cached real ordinates gamma_n; the term then is

      Gamma(s/2, a gamma_n^2)/Gamma(s/2) * gamma_n^{-s}.

    Truncating at N, the tail is bounded by the companion A-note (Thm 5.1):
    with Y = gamma_N, sigma = Re s, mu = max(sigma/2,1),
    cE = 0.112 + 0.278/log Y, and the Gaussian moments (eq. 4.2)
      G0 = 1/(2 sqrt a) Gamma(1/2, aY^2),
      G2 = (sqrt a /2) Gamma(-1/2, aY^2),
      G3 = (a/2)      Gamma(-1,   aY^2),
    h = e^{-aY^2}/Y^2, M(Y) the Trudgian majorant of N(Y),

      |A - A_N| <= 2 mu a^{sigma/2-1} e^{a/4} / |Gamma(s/2)| *
        [ h (M(Y)-N) + 1/(2pi) ((log(Y/2pi)-1/2) G2 + G0/(2Y^2)) + cE G3 ].

    Valid for sigma > -N; if sigma > 2 also requires a(Y^2-1/4) >= sigma/2.
*/

/* Simple upper bound on erfc(x), x >= 0, given a lower bound x_lo on x:
   erfc(x) <= exp(-x^2)/(sqrt(1+(c x)^2)+c x),  c = sqrt(pi)/2. */
static void
mag_erfc(mag_t res, const mag_t x_lo)
{
    mag_t cx, x2, denom, num, t, one, rs;
    mag_init(cx); mag_init(x2); mag_init(denom); mag_init(num);
    mag_init(t); mag_init(one); mag_init(rs);

    mag_set_d_lower(cx, 0.8862269);       /* c = sqrt(pi)/2 lower */
    mag_mul_lower(cx, cx, x_lo);          /* c x lower */

    mag_mul_lower(t, cx, cx);             /* (c x)^2 lower */
    mag_one(one);
    mag_add_lower(t, t, one);             /* 1 + (c x)^2 lower */
    mag_rsqrt(rs, t);                     /* 1/sqrt(t) upper */
    mag_div_lower(denom, one, rs);        /* sqrt(t) lower */
    mag_add_lower(denom, denom, cx);      /* + c x */

    mag_mul_lower(x2, x_lo, x_lo);        /* x^2 lower */
    mag_expinv(num, x2);                  /* e^{-x^2} upper */

    mag_div(res, num, denom);

    mag_clear(cx); mag_clear(x2); mag_clear(denom); mag_clear(num);
    mag_clear(t); mag_clear(one); mag_clear(rs);
}

/*
    Upper bounds on the incomplete gammas at z = a Y^2 (z supplied as a lower
    bound z_lo and the exact mag z_up where needed):
        Gamma(1/2, z) = sqrt(pi) erfc(sqrt z)        <= sqrt(pi) * mag_erfc(sqrt z_lo)
        Gamma(-1/2,z) <= z^{-3/2} e^{-z}             (s<=1: integrand decreasing)
        Gamma(-1,  z) <= z^{-2}  e^{-z}
    For the negative-order ones we need z^{-p} upper (=> z_lo^{-p}) and e^{-z}
    upper (=> z_lo in expinv).
*/
static void
mag_gamma_upper_half(mag_t res, const mag_t z_lo)   /* Gamma(1/2,z) */
{
    mag_t sq;
    mag_init(sq);
    /* sqrt(z) lower for erfc argument */
    {
        mag_t one, rs;
        mag_init(one); mag_init(rs);
        mag_one(one);
        mag_rsqrt(rs, z_lo);              /* 1/sqrt(z) upper */
        mag_div_lower(sq, one, rs);       /* sqrt(z) lower */
        mag_clear(one); mag_clear(rs);
    }
    mag_erfc(res, sq);
    /* * sqrt(pi) upper < 1.7724539 */
    {
        mag_t sp; mag_init(sp); mag_set_d(sp, 1.7724539);
        mag_mul(res, res, sp); mag_clear(sp);
    }
    mag_clear(sq);
}

static void
mag_gamma_upper_negp(mag_t res, const mag_t z_lo, ulong twice_p)
{
    /* Gamma(s,z) <= z^{-p} e^{-z}, p = twice_p/2 (twice_p = 3 -> z^{-3/2}, 4 -> z^{-2}) */
    mag_t zp, ez, one;
    mag_init(zp); mag_init(ez); mag_init(one);
    mag_one(one);

    if (twice_p == 3)
    {
        /* z^{-3/2} = (1/z) * (1/sqrt z) ; upper => use z_lo in denominators */
        mag_t inv, rs;
        mag_init(inv); mag_init(rs);
        mag_div(inv, one, z_lo);          /* 1/z upper */
        mag_rsqrt(rs, z_lo);              /* 1/sqrt z upper */
        mag_mul(zp, inv, rs);             /* z^{-3/2} upper */
        mag_clear(inv); mag_clear(rs);
    }
    else /* twice_p == 4 -> z^{-2} */
    {
        mag_div(zp, one, z_lo);
        mag_mul(zp, zp, zp);              /* z^{-2} upper */
    }
    mag_expinv(ez, z_lo);                 /* e^{-z} upper (z_lo lower) */
    mag_mul(res, zp, ez);

    mag_clear(zp); mag_clear(ez); mag_clear(one);
}

/*
    N-independent data for the A tail bound, precomputed once per (s,a). Only
    the prefactor 2 mu a^{sigma/2-1} e^{a/4} |1/Gamma(s/2)| and sigma are
    N-independent; everything else depends on Y = gamma_N.
*/
typedef struct
{
    mag_t pref;        /* 2 mu a^{sigma/2-1} e^{a/4} |1/Gamma(s/2)| */
    mag_t am;          /* a (upper) */
    mag_t am_lo;       /* a (lower) */
    arb_t a_arb;       /* a */
    arb_t sigma;       /* Re s */
    arb_t half_sigma;  /* sigma/2 (for the sigma>2 validity condition) */
}
A_pre_struct;

static void
A_pre_init(A_pre_struct *A, const acb_t s, const arb_t a, slong prec)
{
    mag_t rg, mu, apow, ea4;
    arb_t sigma;

    mag_init(A->pref); mag_init(A->am); mag_init(A->am_lo);
    arb_init(A->a_arb); arb_init(A->sigma); arb_init(A->half_sigma);
    mag_init(rg); mag_init(mu); mag_init(apow); mag_init(ea4);
    arb_init(sigma);

    arb_set(sigma, acb_realref(s));
    arb_set(A->sigma, sigma);
    arb_mul_2exp_si(A->half_sigma, sigma, -1);   /* sigma/2 */
    arb_set(A->a_arb, a);
    arb_get_mag(A->am, a);
    arb_get_mag_lower(A->am_lo, a);

    /* mu = max(sigma/2, 1) (upper); max-with-1 in arb to handle sign */
    {
        arb_t sh, one;
        arb_init(sh); arb_init(one);
        arb_mul_2exp_si(sh, sigma, -1);
        arb_one(one); arb_max(sh, sh, one, prec);
        arb_get_mag(mu, sh);
        arb_clear(sh); arb_clear(one);
    }

    /* a^{sigma/2-1} (upper): a<1 decreasing in exponent */
    {
        arb_t ax, ex;
        arb_init(ax); arb_init(ex);
        arb_mul_2exp_si(ex, sigma, -1); arb_sub_ui(ex, ex, 1, prec);
        arb_pow(ax, a, ex, prec);
        arb_get_mag(apow, ax);
        arb_clear(ax); arb_clear(ex);
    }

    /* e^{a/4} (upper) */
    mag_mul_2exp_si(ea4, A->am, -2);
    mag_exp(ea4, ea4);

    /* |1/Gamma(s/2)| (upper) */
    {
        acb_t shalf, rgv;
        acb_init(shalf); acb_init(rgv);
        acb_mul_2exp_si(shalf, s, -1);
        acb_rgamma(rgv, shalf, prec);
        acb_get_mag(rg, rgv);
        acb_clear(shalf); acb_clear(rgv);
    }

    /* pref = 2 mu a^{sigma/2-1} e^{a/4} |1/Gamma(s/2)| */
    mag_mul(A->pref, mu, apow);
    mag_mul_2exp_si(A->pref, A->pref, 1);
    mag_mul(A->pref, A->pref, ea4);
    mag_mul(A->pref, A->pref, rg);

    mag_clear(rg); mag_clear(mu); mag_clear(apow); mag_clear(ea4);
    arb_clear(sigma);
}

static void
A_pre_clear(A_pre_struct *A)
{
    mag_clear(A->pref); mag_clear(A->am); mag_clear(A->am_lo);
    arb_clear(A->a_arb); arb_clear(A->sigma); arb_clear(A->half_sigma);
}

/* rigorous tail bound for truncation at N zeros (Y = gamma_N), using
   precomputed A. Returns +inf if the validity condition (only active for
   sigma > 2) a(Y^2 - 1/4) >= sigma/2 is not rigorously satisfied. */
static void
A_bound(mag_t bound, const A_pre_struct *A, const arb_t Y, slong N, slong prec)
{
    mag_t aY2_lo, G0, G2, G3, h, inner, t, u, Ylo, Yup, cE;

    mag_init(aY2_lo); mag_init(G0); mag_init(G2); mag_init(G3); mag_init(h);
    mag_init(inner); mag_init(t); mag_init(u); mag_init(Ylo); mag_init(Yup);
    mag_init(cE);

    /* validity condition (only restrictive when sigma > 2):
       a (Y^2 - 1/4) >= sigma/2. Check rigorously; if it can fail, +inf. */
    {
        arb_t lhs, cond;
        arb_init(lhs); arb_init(cond);
        arb_mul(lhs, Y, Y, prec);              /* Y^2 */
        arb_sub_ui(cond, lhs, 0, prec);
        arb_set_d(cond, 0.25);
        arb_sub(lhs, lhs, cond, prec);          /* Y^2 - 1/4 */
        arb_mul(lhs, lhs, A->a_arb, prec);      /* a(Y^2-1/4) */
        arb_sub(cond, lhs, A->half_sigma, prec);/* a(Y^2-1/4) - sigma/2 */
        if (!arb_is_nonnegative(cond))
        {
            mag_inf(bound);
            arb_clear(lhs); arb_clear(cond);
            goto cleanup;
        }
        arb_clear(lhs); arb_clear(cond);
    }

    arb_get_mag(Yup, Y);
    arb_get_mag_lower(Ylo, Y);

    /* aY2 lower = a_lo * Y_lo^2 */
    mag_mul_lower(aY2_lo, Ylo, Ylo);
    mag_mul_lower(aY2_lo, aY2_lo, A->am_lo);

    /* cE = 0.112 + 0.278/log Y  (upper): log Y lower in divisor */
    {
        arb_t lY; mag_t lYlo;
        arb_init(lY); mag_init(lYlo);
        arb_log(lY, Y, prec);
        arb_get_mag_lower(lYlo, lY);
        mag_set_d(cE, 0.278);
        mag_div(cE, cE, lYlo);
        mag_set_d(t, 0.112);
        mag_add(cE, cE, t);
        mag_clear(lYlo); arb_clear(lY);
    }

    /* moments (upper), closed forms in aY2_lo */
    mag_gamma_upper_half(G0, aY2_lo);          /* Gamma(1/2,aY^2) */
    {
        mag_t rs; mag_init(rs);
        mag_rsqrt(rs, A->am_lo);               /* 1/sqrt a upper */
        mag_mul(G0, G0, rs);
        mag_mul_2exp_si(G0, G0, -1);           /* /2 */
        mag_clear(rs);
    }
    mag_gamma_upper_negp(G2, aY2_lo, 3);       /* Gamma(-1/2,aY^2) */
    {
        mag_t sa; mag_init(sa);
        mag_sqrt(sa, A->am);                   /* sqrt a upper */
        mag_mul(G2, G2, sa);
        mag_mul_2exp_si(G2, G2, -1);
        mag_clear(sa);
    }
    mag_gamma_upper_negp(G3, aY2_lo, 4);       /* Gamma(-1,aY^2) */
    mag_mul(G3, G3, A->am);
    mag_mul_2exp_si(G3, G3, -1);

    /* h = e^{-aY^2}/Y^2 (upper) */
    mag_expinv(h, aY2_lo);
    {
        mag_t y2; mag_init(y2);
        mag_mul_lower(y2, Ylo, Ylo);
        mag_div(h, h, y2);
        mag_clear(y2);
    }

    /* inner = h (M(Y)-N) + 1/(2pi)((log(Y/2pi)-1/2)G2 + G0/(2Y^2)) + cE G3 */
    {
        arb_t MmN, Marb, tt, uu, twopi;
        mag_t mn;
        arb_init(MmN); arb_init(Marb); arb_init(tt); arb_init(uu); arb_init(twopi);
        mag_init(mn);
        arb_const_pi(twopi, prec); arb_mul_2exp_si(twopi, twopi, 1);
        arb_div(tt, Y, twopi, prec); arb_div(uu, Y, twopi, prec);
        arb_log(uu, uu, prec); arb_mul(Marb, tt, uu, prec); arb_sub(Marb, Marb, tt, prec);
        arb_set_d(uu, 0.875); arb_add(Marb, Marb, uu, prec);
        arb_log(tt, Y, prec); arb_set_d(uu, 0.112); arb_addmul(Marb, tt, uu, prec);
        arb_log(uu, tt, prec); arb_set_d(tt, 0.278); arb_addmul(Marb, uu, tt, prec);
        arb_set_d(uu, 3.385); arb_add(Marb, Marb, uu, prec);
        arb_sub_ui(MmN, Marb, (ulong) N, prec);
        arb_get_mag(mn, MmN);
        mag_mul(inner, h, mn);
        arb_clear(MmN); arb_clear(Marb); arb_clear(tt); arb_clear(uu); arb_clear(twopi);
        mag_clear(mn);
    }
    {
        arb_t lYt, twopi; mag_t lYm, y2;
        arb_init(lYt); arb_init(twopi); mag_init(lYm); mag_init(y2);
        arb_const_pi(twopi, prec); arb_mul_2exp_si(twopi, twopi, 1);
        arb_div(lYt, Y, twopi, prec); arb_log(lYt, lYt, prec);
        arb_set_d(twopi, 0.5); arb_sub(lYt, lYt, twopi, prec);
        arb_get_mag(lYm, lYt);
        mag_mul(t, lYm, G2);
        mag_mul_lower(y2, Ylo, Ylo); mag_mul_2exp_si(y2, y2, 1);
        mag_div(u, G0, y2);
        mag_add(t, t, u);
        mag_set_d(u, 0.1591550);               /* 1/(2pi) upper */
        mag_mul(t, t, u);
        mag_add(inner, inner, t);
        arb_clear(lYt); arb_clear(twopi); mag_clear(lYm); mag_clear(y2);
    }
    mag_mul(t, cE, G3);
    mag_add(inner, inner, t);

    /* * precomputed prefactor */
    mag_mul(bound, inner, A->pref);

cleanup:
    mag_clear(aY2_lo); mag_clear(G0); mag_clear(G2); mag_clear(G3); mag_clear(h);
    mag_clear(inner); mag_clear(t); mag_clear(u); mag_clear(Ylo); mag_clear(Yup);
    mag_clear(cE);
}

/* choose number of zeros N so the tail bound < 2^{-prec}. Starts at a fixed N
   (ignoring s, avoiding any double cast of sigma) and grows by ~1/8; A_bound
   returns +inf for N that violate the validity condition, which the search
   steps past. */
static slong
A_choose_N(const A_pre_struct *A, slong prec, slong zeros_prec)
{
    slong N = 2;
    mag_t b, target;
    arb_t rh;

    mag_init(b); mag_init(target); arb_init(rh);
    mag_one(target);
    mag_mul_2exp_si(target, target, -prec);
    arb_set_d(rh, SECOND_ZETA_RH_VERIFIED_HEIGHT);

    /*
        Cheap infeasibility check. The bound is valid only for Y = gamma_N with
        a(Y^2 - 1/4) >= sigma/2, i.e. Y >= sqrt(sigma/(2a) + 1/4). If even the
        RH-verified height does not satisfy this (which happens for very large
        real sigma), no certifiable N exists; return a sentinel large N so the
        caller's RH-height guard fires immediately, without fetching zeros up to
        an uncertifiable height.
    */
    {
        arb_t need, rh2;
        arb_init(need); arb_init(rh2);
        /* need = sigma/(2a) + 1/4  (lower bound; compare to RH^2) */
        arb_mul_2exp_si(need, A->a_arb, 1);          /* 2a */
        arb_div(need, A->sigma, need, prec);         /* sigma/(2a) */
        arb_set_d(rh2, 0.25); arb_add(need, need, rh2, prec);
        arb_mul(rh2, rh, rh, prec);                  /* RH^2 */
        if (arb_gt(need, rh2))                       /* required Y^2 > RH^2 */
        {
            arb_clear(need); arb_clear(rh2);
            mag_clear(b); mag_clear(target); arb_clear(rh);
            return (WORD(1) << 30);   /* triggers RH-height guard in caller */
        }
        arb_clear(need); arb_clear(rh2);
    }

    while (N < (WORD(1) << 30))
    {
        arb_ptr g = _arb_vec_init(N);
        /* zeros at the (target) working precision so they populate the cache
           for the summation that follows; the bound itself is cheap */
        _acb_dirichlet_secondary_zeta_zeros(g, N, zeros_prec);
        A_bound(b, A, g + (N - 1), N, SECOND_ZETA_BOUND_PREC);
        if (mag_cmp(b, target) < 0)
        {
            _arb_vec_clear(g, N);
            break;
        }
        /* gamma_N at/over the RH-verified height: stop and let the caller
           return indeterminate rather than ask for an uncertifiable height. */
        if (!arb_lt(g + (N - 1), rh))
        {
            _arb_vec_clear(g, N);
            break;
        }
        _arb_vec_clear(g, N);
        N = N + N / 8 + 1;
    }
    mag_clear(b); mag_clear(target); arb_clear(rh);
    return N;
}

void
_acb_dirichlet_secondary_zeta_term_A(acb_t res, const acb_t s,
        const arb_t a, slong prec)
{
    acb_t sum, shalf, rg_half, term, gu, zarg, gpow, t;
    arb_t Y, g2, az2;
    arb_ptr gammas;
    mag_t bound;
    A_pre_struct A;
    slong n, N, wp;

    wp = prec;

    A_pre_init(&A, s, a, SECOND_ZETA_BOUND_PREC);
    N = A_choose_N(&A, wp, wp);

    /* A_choose_N returns the sentinel 2^30 when no certifiable cutoff exists
       (e.g. very large real s needs zeros above the RH-verified height). */
    if (N >= (WORD(1) << 30))
    {
        acb_indeterminate(res);
        A_pre_clear(&A);
        return;
    }

    gammas = _arb_vec_init(N);
    _acb_dirichlet_secondary_zeta_zeros(gammas, N, wp);

    acb_init(sum); acb_init(shalf); acb_init(rg_half);
    acb_init(term); acb_init(gu); acb_init(zarg); acb_init(gpow); acb_init(t);
    arb_init(Y); arb_init(g2); arb_init(az2);
    mag_init(bound);

    acb_mul_2exp_si(shalf, s, -1);
    acb_rgamma(rg_half, shalf, wp);          /* 1/Gamma(s/2) */

    /*
        RH-verified-height check: the rigor of treating alpha_n = gamma_n as
        real relies on RH up to 3e12 (Platt-Trudgian). gamma_N is the largest
        ordinate used; if its lower bound exceeds the verified height we cannot
        certify the result, so we set it to indeterminate (NaN) rather than
        return an unsound enclosure. In practice N is small (a few hundred at
        most) so this never triggers for sane inputs.
    */
    arb_set_d(Y, SECOND_ZETA_RH_VERIFIED_HEIGHT);
    if (!arb_lt(gammas + (N - 1), Y))
    {
        acb_indeterminate(res);
    }
    else
    {
        acb_zero(sum);

        for (n = 0; n < N; n++)
        {
            /* zarg = a gamma_n^2 ; gpow = gamma_n^{-s} */
            arb_sqr(g2, gammas + n, wp);
            arb_mul(az2, g2, a, wp);
            acb_set_arb(zarg, az2);

            /* Gamma(s/2, a gamma_n^2) */
            acb_hypgeom_gamma_upper(gu, shalf, zarg, 0, wp);

            /* gamma_n^{-s} = exp(-s log gamma_n) */
            arb_log(g2, gammas + n, wp);
            acb_set_arb(gpow, g2);
            acb_mul(gpow, gpow, s, wp);
            acb_neg(gpow, gpow);
            acb_exp(gpow, gpow, wp);

            acb_mul(term, gu, gpow, wp);
            acb_add(sum, sum, term, wp);
        }

        /* multiply by 1/Gamma(s/2) */
        acb_mul(res, sum, rg_half, wp);

        /* attach rigorous tail radius using Y = gamma_N */
        A_bound(bound, &A, gammas + (N - 1), N, SECOND_ZETA_BOUND_PREC);
        arb_add_error_mag(acb_realref(res), bound);
        if (!acb_is_real(s))
            arb_add_error_mag(acb_imagref(res), bound);

        acb_set_round(res, res, prec);
    }

    _arb_vec_clear(gammas, N);
    acb_clear(sum); acb_clear(shalf); acb_clear(rg_half);
    acb_clear(term); acb_clear(gu); acb_clear(zarg); acb_clear(gpow); acb_clear(t);
    arb_clear(Y); arb_clear(g2); arb_clear(az2);
    mag_clear(bound);
    A_pre_clear(&A);
}

/*****************************************************************************/
/*                                  P term                                   */
/*****************************************************************************/

/*
    Prime term  (Arias de Reyna eq. 38)

      P(s) = 1/(2 sqrt(pi)) sum_{n>=2} Lambda(n)/sqrt(n)
             * Gamma((1-s)/2, log^2 n /(4a)) / Gamma(s/2) * (2/log n)^{1-s},

    Lambda the von Mangoldt function (nonzero only at prime powers n=p^k,
    where Lambda = log p). We sum the prime-power terms up to a cutoff N
    chosen so that the rigorous remainder bound falls below the target, and
    attach that bound as the error radius.

    Tail bound (companion note). With sigma = Re s,
    delta = max(1,(1-sigma)/2), c = 1.04 (Rosser-Schoenfeld psi(x)<1.03883x),
    z = log N/(2 sqrt a) - sqrt(a)/2, w(u)=1/(sqrt u log^2 u) e^{-log^2u/4a},
    and P0 = delta * 2 a^{(1+sigma)/2} / (sqrt(pi) |Gamma(s/2)|),

      erfc form (4.1):  |P_N| <= delta * 2 a^{1+sigma/2} e^{a/4}
                                  / (|Gamma(s/2)| log N) * erfc(z),
      psi  form (4.2):  |P_N| <= P0 * [ (cN - psi(N)) w(N)
                                  + c sqrt(pi a)/log^2 N * e^{a/4} erfc(z) ].

    Both are rigorous; we take the smaller for the radius. The paper's eq. (46)
    is the erfc form after two avoidable concessions and is ~40-200 orders of
    magnitude weaker at small a.
*/

/* rigorous upper bound on |1/Gamma(s/2)| */
static void
rgamma_half_mag(mag_t m, const acb_t s, slong prec)
{
    acb_t shalf, rg;
    acb_init(shalf); acb_init(rg);
    acb_mul_2exp_si(shalf, s, -1);
    acb_rgamma(rg, shalf, prec);
    acb_get_mag(m, rg);
    acb_clear(shalf); acb_clear(rg);
}

/*
    N-independent data for the P tail bound, precomputed once per (s,a) and
    reused across the N search. The bound (erfc form; the psi form was removed
    as it gave no useful improvement) is

      |P_N| <= delta * 2 a^{1+sigma/2} e^{a/4} |1/Gamma(s/2)| / log N * erfc(z),
      z = log N/(2 sqrt a) - sqrt(a)/2,   delta = max(1, (1-sigma)/2),

    valid for sigma > 1 - log^2 N /(2a); P_bound returns +inf when that
    (rigorously checked) condition can fail, so the N search skips such N.
*/
typedef struct
{
    mag_t pref;        /* delta * 2 a^{1+sigma/2} e^{a/4} |1/Gamma(s/2)| */
    mag_t inv_2sqrta;  /* 1/(2 sqrt a)  (upper) for z's leading term */
    mag_t inv_2sqrta_lo; /* 1/(2 sqrt a) (lower) */
    mag_t halfsqrta;   /* sqrt(a)/2     (upper) subtracted in z */
    arb_t a_arb;       /* a (for log^2 N >= 2a(1-sigma) validity test) */
    arb_t sigma;       /* Re s */
}
P_pre_struct;

static void
P_pre_init(P_pre_struct *P, const acb_t s, const arb_t a, slong prec)
{
    mag_t rg, am, delta, ea4, apow;
    arb_t sigma;

    mag_init(P->pref); mag_init(P->inv_2sqrta); mag_init(P->inv_2sqrta_lo);
    mag_init(P->halfsqrta); arb_init(P->a_arb); arb_init(P->sigma);
    mag_init(rg); mag_init(am); mag_init(delta); mag_init(ea4); mag_init(apow);
    arb_init(sigma);

    arb_set(sigma, acb_realref(s));
    arb_set(P->sigma, sigma);
    arb_set(P->a_arb, a);

    rgamma_half_mag(rg, s, prec);          /* |1/Gamma(s/2)| upper */
    arb_get_mag(am, a);

    /* delta = max(1, (1-sigma)/2) (upper); max-with-1 in arb to handle sign */
    {
        arb_t d, one;
        arb_init(d); arb_init(one);
        arb_sub_ui(d, sigma, 1, prec); arb_neg(d, d); arb_mul_2exp_si(d, d, -1);
        arb_one(one); arb_max(d, d, one, prec);
        arb_get_mag(delta, d);
        arb_clear(d); arb_clear(one);
    }

    /* e^{a/4} (upper) */
    mag_mul_2exp_si(ea4, am, -2);
    mag_exp(ea4, ea4);

    /* a^{1+sigma/2} (upper): a<1 decreasing in exponent */
    {
        arb_t ax, ex;
        arb_init(ax); arb_init(ex);
        arb_mul_2exp_si(ex, sigma, -1); arb_add_ui(ex, ex, 1, prec);
        arb_pow(ax, a, ex, prec);
        arb_get_mag(apow, ax);
        arb_clear(ax); arb_clear(ex);
    }

    /* pref = delta * 2 * a^{1+sigma/2} * e^{a/4} * |1/Gamma(s/2)| */
    mag_mul(P->pref, delta, apow);
    mag_mul_2exp_si(P->pref, P->pref, 1);
    mag_mul(P->pref, P->pref, ea4);
    mag_mul(P->pref, P->pref, rg);

    /* 1/(2 sqrt a) upper & lower, and sqrt(a)/2 upper */
    {
        arb_t sq, h;
        mag_t slo, sup;
        arb_init(sq); arb_init(h); mag_init(slo); mag_init(sup);
        arb_sqrt(sq, a, prec);             /* sqrt a */
        arb_get_mag(sup, sq);              /* sqrt a upper */
        arb_get_mag_lower(slo, sq);        /* sqrt a lower */
        /* 1/(2 sqrt a): upper uses sqrt a lower; lower uses sqrt a upper */
        mag_one(P->inv_2sqrta);
        mag_div(P->inv_2sqrta, P->inv_2sqrta, slo);
        mag_mul_2exp_si(P->inv_2sqrta, P->inv_2sqrta, -1);
        mag_one(P->inv_2sqrta_lo);
        mag_div_lower(P->inv_2sqrta_lo, P->inv_2sqrta_lo, sup);
        mag_mul_2exp_si(P->inv_2sqrta_lo, P->inv_2sqrta_lo, -1);
        /* sqrt(a)/2 upper */
        mag_set(P->halfsqrta, sup);
        mag_mul_2exp_si(P->halfsqrta, P->halfsqrta, -1);
        arb_clear(sq); arb_clear(h); mag_clear(slo); mag_clear(sup);
    }

    mag_clear(rg); mag_clear(am); mag_clear(delta); mag_clear(ea4); mag_clear(apow);
    arb_clear(sigma);
}

static void
P_pre_clear(P_pre_struct *P)
{
    mag_clear(P->pref); mag_clear(P->inv_2sqrta); mag_clear(P->inv_2sqrta_lo);
    mag_clear(P->halfsqrta); arb_clear(P->a_arb); arb_clear(P->sigma);
}

/* erfc-form bound at cutoff N using precomputed P. Returns +inf if the
   validity condition sigma > 1 - log^2 N/(2a) is not rigorously satisfied. */
static void
P_bound(mag_t bound, const P_pre_struct *P, slong N, slong prec)
{
    mag_t logN_lo, z_lo, erfcz, t;
    arb_t L, valid;
    int ok;

    mag_init(logN_lo); mag_init(z_lo); mag_init(erfcz); mag_init(t);
    arb_init(L); arb_init(valid);

    arb_log_ui(L, (ulong) N, prec);        /* log N (enclosure) */

    /* validity: log^2 N /(2a) > 1 - sigma, i.e. sigma + log^2N/(2a) - 1 > 0.
       If not rigorously > 0, return +inf so the search advances. */
    arb_sqr(valid, L, prec);
    {
        arb_t twoa;
        arb_init(twoa);
        arb_mul_2exp_si(twoa, P->a_arb, 1);          /* 2a */
        arb_div(valid, valid, twoa, prec);           /* log^2N/(2a) */
        arb_clear(twoa);
    }
    arb_add(valid, valid, P->sigma, prec);
    arb_sub_ui(valid, valid, 1, prec);               /* sigma + log^2N/(2a) - 1 */
    ok = arb_is_positive(valid);

    if (!ok)
    {
        mag_inf(bound);
        goto cleanup;
    }

    arb_get_mag_lower(logN_lo, L);         /* log N lower (>0 for N>=2) */

    /* z_lo = logN_lo/(2 sqrt a)_lo - (sqrt a /2)_up   (lower bound on z) */
    {
        mag_t hi;
        mag_init(hi);
        mag_mul_lower(z_lo, logN_lo, P->inv_2sqrta_lo);  /* lower term */
        /* subtract upper bound of sqrt(a)/2; guard against underflow below 0 */
        if (mag_cmp(z_lo, P->halfsqrta) <= 0)
        {
            /* z could be <= 0: erfc <= 1, use erfc bound 1 (conservative) */
            mag_one(erfcz);
            mag_clear(hi);
            goto assemble;
        }
        mag_sub_lower(z_lo, z_lo, P->halfsqrta);  /* (a - b) lower */
        mag_clear(hi);
    }
    mag_erfc(erfcz, z_lo);

assemble:
    /* bound = pref / logN * erfc(z) ; logN lower divisor -> upper */
    mag_div(bound, P->pref, logN_lo);
    mag_mul(bound, bound, erfcz);

cleanup:
    mag_clear(logN_lo); mag_clear(z_lo); mag_clear(erfcz); mag_clear(t);
    arb_clear(L); arb_clear(valid);
}

/* choose cutoff N so that the erfc-form bound < 2^{-prec}. Starts at a fixed N
   (ignoring s) and grows by ~1/8; P_bound returns +inf for invalid N, which the
   search simply steps past. */
static slong
P_choose_N(const P_pre_struct *P, slong prec)
{
    slong N = 8;
    mag_t b, target;
    mag_init(b); mag_init(target);

    mag_one(target);
    mag_mul_2exp_si(target, target, -prec);

    while (N < (WORD(1) << 30))
    {
        P_bound(b, P, N, SECOND_ZETA_BOUND_PREC);
        if (mag_cmp(b, target) < 0)
            break;
        N = N + N / 8 + 1;     /* ~1.125x growth */
    }
    mag_clear(b); mag_clear(target);
    return N;
}

void
_acb_dirichlet_secondary_zeta_term_P(acb_t res, const acb_t s,
        const arb_t a, slong prec)
{
    acb_t sum, rg_half, onem_s, shalf_arg, term, gu, powfac;
    arb_t a4, logp, logn, z_re, t;
    mag_t bound;
    P_pre_struct P;
    n_factor_t fac;
    slong n, N, wp;

    wp = prec;

    /* N-independent bound data, then choose the cutoff */
    P_pre_init(&P, s, a, SECOND_ZETA_BOUND_PREC);
    N = P_choose_N(&P, wp);

    acb_init(sum); acb_init(rg_half); acb_init(onem_s);
    acb_init(shalf_arg); acb_init(term); acb_init(gu); acb_init(powfac);
    arb_init(a4); arb_init(logp); arb_init(logn); arb_init(z_re); arb_init(t);
    mag_init(bound);

    /* 1/Gamma(s/2) */
    acb_mul_2exp_si(shalf_arg, s, -1);
    acb_rgamma(rg_half, shalf_arg, wp);

    /* (1-s)/2  and  1-s */
    acb_one(onem_s);
    acb_sub(onem_s, onem_s, s, wp);          /* 1 - s */
    acb_mul_2exp_si(shalf_arg, onem_s, -1);  /* (1-s)/2 */

    /* a4 = 1/(4a) */
    arb_mul_ui(a4, a, 4, wp);
    arb_inv(a4, a4, wp);

    acb_zero(sum);

    for (n = 2; n <= N; n++)
    {
        ulong p;
        n_factor_init(&fac);
        n_factor(&fac, (ulong) n, 0);
        if (fac.num != 1)
            continue;                        /* not a prime power */
        p = fac.p[0];

        /* Lambda(n) = log p ;  log n = exp * log p */
        arb_log_ui(logp, p, wp);
        arb_mul_ui(logn, logp, (ulong) fac.exp[0], wp);

        /* z_re = log^2 n /(4a) */
        arb_sqr(z_re, logn, wp);
        arb_mul(z_re, z_re, a4, wp);

        /* gu = Gamma((1-s)/2, log^2 n/(4a)) */
        acb_set_arb(term, z_re);             /* reuse term as the z argument */
        acb_hypgeom_gamma_upper(gu, shalf_arg, term, 0, wp);

        /* powfac = (2/log n)^{1-s} = exp((1-s) log(2/log n)) */
        arb_set_ui(t, 2);
        arb_div(t, t, logn, wp);
        arb_log(t, t, wp);
        acb_set_arb(powfac, t);
        acb_mul(powfac, powfac, onem_s, wp);
        acb_exp(powfac, powfac, wp);

        /* term = (log p)/sqrt(n) * gu * powfac */
        acb_mul(term, gu, powfac, wp);
        arb_sqrt_ui(t, (ulong) n, wp);
        arb_div(t, logp, t, wp);             /* (log p)/sqrt n */
        acb_mul_arb(term, term, t, wp);

        acb_add(sum, sum, term, wp);
    }

    /* multiply by 1/(2 sqrt pi) and by 1/Gamma(s/2) */
    acb_mul(sum, sum, rg_half, wp);
    arb_const_pi(t, wp);
    arb_sqrt(t, t, wp);
    arb_mul_2exp_si(t, t, 1);                /* 2 sqrt pi */
    acb_div_arb(res, sum, t, wp);

    P_bound(bound, &P, N, SECOND_ZETA_BOUND_PREC);
    arb_add_error_mag(acb_realref(res), bound);
    if (!acb_is_real(s))
        arb_add_error_mag(acb_imagref(res), bound);

    acb_set_round(res, res, prec);

    acb_clear(sum); acb_clear(rg_half); acb_clear(onem_s);
    acb_clear(shalf_arg); acb_clear(term); acb_clear(gu); acb_clear(powfac);
    arb_clear(a4); arb_clear(logp); arb_clear(logn); arb_clear(z_re);
    arb_clear(t);
    mag_clear(bound);
    P_pre_clear(&P);
}

/*
    Fill bpv[0..len-1] with B_n(x)/n!, the Taylor coefficients of the Bernoulli
    polynomial generating function

        t e^{xt}/(e^t - 1) = sum_{n>=0} (B_n(x)/n!) t^n.

    Computed as a single power-series division:
        numerator   e^{xt}            (coeffs x^k/k!),
        denominator (e^t-1)/t         (coeffs 1/(k+1)!),
    so B_n(x)/n! = [t^n] (e^{xt} / ((e^t-1)/t)). bpv must be preallocated to
    length len; caller owns it.
*/
static void
_bernoulli_poly_over_factorial_vec(arb_ptr bpv, const arb_t x, slong len, slong prec)
{
    arb_ptr num, den, xt;
    slong k;

    num = _arb_vec_init(len);
    den = _arb_vec_init(len);
    xt = _arb_vec_init(2);

    /* numerator e^{xt}: exp of the degree-1 series (0 + x t) */
    arb_set(xt + 1, x);                       /* xt = x t */
    _arb_poly_exp_series(num, xt, 2, len, prec);

    /* denominator (e^t - 1)/t: coefficient k is 1/(k+1)!, built in place by
       den[k] = den[k-1]/(k+1) from den[0] = 1 */
    arb_one(den + 0);
    for (k = 1; k < len; k++)
        arb_div_ui(den + k, den + (k - 1), (ulong)(k + 1), prec);

    _arb_poly_div_series(bpv, num, len, den, len, len, prec);

    _arb_vec_clear(num, len);
    _arb_vec_clear(den, len);
    _arb_vec_clear(xt, 2);
}

/*****************************************************************************/
/*                                    S term                                 */
/*****************************************************************************/

/*
    Singular term  (Arias de Reyna eq. 53), truncated at order N:

      S(s) ~ a^{(s-1)/2}/(4 sqrt(pi) Gamma(s/2)) *
             { -2/(s-1)^2 + (C0 + log(16 pi^2 a))/(s-1)
               + sum_{n=1}^N B_n(3/4)/n! * (4 sqrt a)^n Gamma(n/2)/(s+n-1) }.

    The remainder S(s) - S_N(s) is bounded (companion note):
    for Re s = sigma > -N and any a > 0,

      |S - S_N| <= 1/(4 sqrt(pi) |Gamma(s/2)|) *
                   [ 2 c_N/(N+sigma) a^{(N+sigma)/2}
                     + (c_N'/kappa) a^{(N+sigma+1)/2} e^{-kappa/a} ],

    with kappa = pi^2/32,
      c_N  = (4 zeta(N+1)/pi) (2/pi)^N Gamma((N+1)/2),
      c_N' = (1/2) 32^{N/2} Gamma(N/2).

    Because the series is asymptotic (divergent), N is chosen to minimize the
    bound; the second (super-exponentially small) summand is negligible but
    included for full rigor. The bracket's spurious poles at s = 0,-2,... are
    killed by the 1/Gamma(s/2) prefactor; genuine poles sit at s = 1,-1,-3,...
*/

/* leading-bound magnitude for a given N (for choosing N), and full radius */
/*
    Upper bound (mag) on the S-term remainder |S(s) - S_N(s)|.
    With sigma = Re s, kappa = pi^2/32,

      |R_N| <= |1/Gamma(s/2)| / (4 sqrt pi) *
               [ 2 c_N/(N+sigma) a^{(N+sigma)/2}
                 + (c_N'/kappa) a^{(N+sigma+1)/2} e^{-kappa/a} ],

    c_N  = (4 zeta(N+1)/pi) (2/pi)^N Gamma((N+1)/2),
    c_N' = (1/2) 32^{N/2} Gamma(N/2).

    Implemented in mag_t with all N-independent quantities hoisted to mag_t
    constants so the N-selection loop is cheap.

    The Gamma factors are reduced to factorials:
      N even:  Gamma((N+1)/2) = Gamma(N/2 + 1/2) <= (N/2)! / sqrt(N/2),
               Gamma(N/2)     = (N/2 - 1)!
      N odd:   Gamma((N+1)/2) = ((N-1)/2)!,
               Gamma(N/2)     = Gamma((N-1)/2 + 1/2) <= ((N-1)/2)! / sqrt((N-1)/2)
    with Gamma(1/2) = sqrt(pi) for the X=0 edge case.
*/

/* upper bound on Gamma(m + 1/2) via Gamma(m+1/2) <= m!/sqrt(m), m>=1;
   Gamma(1/2)=sqrt(pi) for m=0 */
static void
mag_gamma_half_int(mag_t res, ulong m)
{
    if (m == 0)
    {
        /* sqrt(pi) < 1.7724539, round up */
        mag_set_d(res, 1.7724539);
    }
    else
    {
        mag_t r;
        mag_init(r);
        mag_fac_ui(res, m);            /* upper bound on m! */
        mag_set_ui_lower(r, m);        /* lower bound on m */
        mag_rsqrt(r, r);               /* upper bound on 1/sqrt(m) */
        mag_mul(res, res, r);
        mag_clear(r);
    }
}

/* upper bound on Gamma(n) = (n-1)! for integer n>=1 */
static void
mag_gamma_int(mag_t res, ulong n)
{
    mag_fac_ui(res, n - 1);
}

/*
    N-independent data for the S remainder bound, precomputed once per (s,a).
*/
typedef struct
{
    mag_t rgam;        /* |1/Gamma(s/2)| */
    mag_t four_sqrt_pi_inv;  /* 1/(4 sqrt pi) */
    mag_t kappa_inv;   /* 32/pi^2 */
    mag_t two_over_pi; /* 2/pi */
    mag_t expfac;      /* e^{-kappa/a} */
    arb_t a_arb;       /* a */
    arb_t sigma;       /* Re s */
}
S_pre_struct;

static void
S_pre_init(S_pre_struct *S, const acb_t s, const arb_t a, slong prec)
{
    mag_t am, koa;

    mag_init(S->rgam); mag_init(S->four_sqrt_pi_inv); mag_init(S->kappa_inv);
    mag_init(S->two_over_pi); mag_init(S->expfac);
    arb_init(S->a_arb); arb_init(S->sigma);
    mag_init(am); mag_init(koa);

    arb_set(S->sigma, acb_realref(s));
    arb_set(S->a_arb, a);
    arb_get_mag(am, a);

    mag_set_d(S->two_over_pi, 0.63661978);       /* 2/pi upper */
    mag_set_d(S->four_sqrt_pi_inv, 0.141047396); /* 1/(4 sqrt pi) upper */
    mag_set_d(S->kappa_inv, 3.2424114);          /* 32/pi^2 upper */

    /* |1/Gamma(s/2)| (upper) -- =0 at s=0,-2,... */
    {
        acb_t shalf, rg;
        acb_init(shalf); acb_init(rg);
        acb_mul_2exp_si(shalf, s, -1);
        acb_rgamma(rg, shalf, prec);
        acb_get_mag(S->rgam, rg);
        acb_clear(shalf); acb_clear(rg);
    }

    /* e^{-kappa/a}: kappa = pi^2/32 > 0.30843 lower; kappa/a lower -> expinv upper */
    mag_set_d_lower(koa, 0.3084251);
    mag_div_lower(koa, koa, am);
    mag_expinv(S->expfac, koa);

    mag_clear(am); mag_clear(koa);
}

static void
S_pre_clear(S_pre_struct *S)
{
    mag_clear(S->rgam); mag_clear(S->four_sqrt_pi_inv); mag_clear(S->kappa_inv);
    mag_clear(S->two_over_pi); mag_clear(S->expfac);
    arb_clear(S->a_arb); arb_clear(S->sigma);
}

/* remainder bound at order N using precomputed S. Returns +inf if the validity
   condition sigma > -N is not rigorously satisfied. */
static void
S_bound(mag_t bound, const S_pre_struct *S, slong N, slong prec)
{
    mag_t Nps_lo, am_pow, zeta_m, gamma_m, cN, cNp, term1, term2, t;
    arb_t nps_arb, ex;

    mag_init(Nps_lo); mag_init(am_pow); mag_init(zeta_m); mag_init(gamma_m);
    mag_init(cN); mag_init(cNp); mag_init(term1); mag_init(term2); mag_init(t);
    arb_init(nps_arb); arb_init(ex);

    /* N + sigma (need a rigorous lower bound for the 1/(N+sigma) divisor).
       Validity sigma > -N: if N+sigma is not rigorously positive, +inf. */
    arb_add_ui(nps_arb, S->sigma, (ulong) N, prec);
    if (!arb_is_positive(nps_arb))
    {
        mag_inf(bound);
        goto cleanup;
    }
    arb_get_mag_lower(Nps_lo, nps_arb);
    if (mag_is_zero(Nps_lo))
        mag_set_d(Nps_lo, 1e-300);

    /* zeta(N+1) upper */
    mag_hurwitz_zeta_uiui(zeta_m, (ulong)(N + 1), 1);

    /* c_N = (4/pi) zeta(N+1) (2/pi)^N Gamma((N+1)/2);  4/pi < 1.2732396 */
    mag_set_d(cN, 1.2732396);
    mag_mul(cN, cN, zeta_m);
    mag_pow_ui(t, S->two_over_pi, (ulong) N);
    mag_mul(cN, cN, t);
    if (N % 2 == 1)
        mag_gamma_int(gamma_m, (ulong)((N + 1) / 2));
    else
        mag_gamma_half_int(gamma_m, (ulong)(N / 2));
    mag_mul(cN, cN, gamma_m);

    /* c_N' = (1/2) 32^{N/2} Gamma(N/2);  sqrt 32 < 5.6568543 */
    mag_set_d(cNp, 5.6568543);
    mag_pow_ui(cNp, cNp, (ulong) N);
    if (N % 2 == 0)
        mag_gamma_int(gamma_m, (ulong)(N / 2));
    else
        mag_gamma_half_int(gamma_m, (ulong)((N - 1) / 2));
    mag_mul(cNp, cNp, gamma_m);
    mag_mul_2exp_si(cNp, cNp, -1);

    /* a-powers (a<1 decreasing in exponent -> mag of arb_pow is an upper bound) */
    {
        arb_t ax;
        arb_init(ax);
        arb_mul_2exp_si(ex, nps_arb, -1);       /* (N+sigma)/2 */
        arb_pow(ax, S->a_arb, ex, prec);
        arb_get_mag(am_pow, ax);
        arb_add_ui(ex, nps_arb, 1, prec);
        arb_mul_2exp_si(ex, ex, -1);            /* (N+sigma+1)/2 */
        arb_pow(ax, S->a_arb, ex, prec);
        arb_get_mag(t, ax);                     /* a^{(N+sigma+1)/2} */
        arb_clear(ax);
    }

    /* term1 = 2 c_N/(N+sigma) a^{(N+sigma)/2} */
    mag_mul(term1, cN, am_pow);
    mag_mul_2exp_si(term1, term1, 1);
    mag_div(term1, term1, Nps_lo);

    /* term2 = (c_N'/kappa) a^{(N+sigma+1)/2} e^{-kappa/a} */
    mag_mul(term2, cNp, S->kappa_inv);
    mag_mul(term2, term2, t);
    mag_mul(term2, term2, S->expfac);

    mag_add(term1, term1, term2);

    /* * 1/(4 sqrt pi) * |1/Gamma(s/2)| */
    mag_mul(term1, term1, S->four_sqrt_pi_inv);
    mag_mul(bound, term1, S->rgam);

cleanup:
    mag_clear(Nps_lo); mag_clear(am_pow); mag_clear(zeta_m); mag_clear(gamma_m);
    mag_clear(cN); mag_clear(cNp); mag_clear(term1); mag_clear(term2); mag_clear(t);
    arb_clear(nps_arb); arb_clear(ex);
}

/* choose N minimizing the remainder bound (asymptotic, divergent series).
   Starts at a fixed N=1 and scans upward; S_bound returns +inf for invalid N
   (sigma <= -N), which is skipped without being mistaken for an increasing
   run. Stops once the bound has genuinely increased for several finite steps. */
static slong
S_choose_N(const S_pre_struct *S, slong prec)
{
    slong N, bestN = 0, Nmax = 4 + prec;
    mag_t b, best;
    mag_init(b); mag_init(best);
    mag_inf(best);

    for (N = 1; N <= Nmax; N++)
    {
        S_bound(b, S, N, SECOND_ZETA_BOUND_PREC);
        if (mag_is_inf(b))
            continue;                  /* invalid N (sigma <= -N): skip */
        if (bestN == 0 || mag_cmp(b, best) < 0)
        {
            mag_set(best, b);
            bestN = N;
        }
        else if (N > bestN + 3)
        {
            break;                     /* bound increasing for several steps */
        }
    }
    if (bestN == 0)
        bestN = Nmax;                  /* should not happen for valid s */
    mag_clear(b); mag_clear(best);
    return bestN;
}

void
_acb_dirichlet_secondary_zeta_term_S(acb_t res, const acb_t s,
        const arb_t a, slong prec)
{
    acb_t shalf, sm1, pref, bracket, term, t, g;
    arb_t alog, sqa, sq4a, bp, gh, u;
    mag_t bound;
    S_pre_struct S;
    slong n, N, wp;

    wp = prec;

    S_pre_init(&S, s, a, SECOND_ZETA_BOUND_PREC);
    N = S_choose_N(&S, wp);

    acb_init(shalf); acb_init(sm1); acb_init(pref);
    acb_init(bracket); acb_init(term); acb_init(t); acb_init(g);
    arb_init(alog); arb_init(sqa); arb_init(sq4a);
    arb_init(bp); arb_init(gh); arb_init(u);
    mag_init(bound);

    /* sm1 = s - 1 */
    acb_sub_ui(sm1, s, 1, wp);

    /* bracket = -2/(s-1)^2 + (C0 + log(16 pi^2 a))/(s-1) */
    acb_sqr(t, sm1, wp);
    acb_set_si(bracket, -2);
    acb_div(bracket, bracket, t, wp);

    /* C0 + log(16 pi^2 a) */
    arb_const_pi(gh, wp);
    arb_sqr(gh, gh, wp);
    arb_mul_ui(gh, gh, 16, wp);
    arb_mul(gh, gh, a, wp);
    arb_log(gh, gh, wp);
    arb_const_euler(u, wp);
    arb_add(gh, gh, u, wp);
    acb_set_arb(term, gh);
    acb_div(term, term, sm1, wp);
    acb_add(bracket, bracket, term, wp);

    /* sum_{n=1}^N B_n(3/4)/n! (4 sqrt a)^n Gamma(n/2)/(s+n-1) */
    arb_sqrt(sqa, a, wp);
    arb_mul_2exp_si(sq4a, sqa, 2);        /* 4 sqrt a */

    {
        arb_ptr bpv;
        arb_t pown, ghalf, gh2, x34;
        slong len = N + 1;

        bpv = _arb_vec_init(len);
        arb_init(pown); arb_init(ghalf); arb_init(gh2); arb_init(x34);

        /* B_n(3/4)/n! for n=0..N in one shot */
        arb_set_d(x34, 0.75);
        _bernoulli_poly_over_factorial_vec(bpv, x34, len, wp);

        /* running (4 sqrt a)^n and Gamma(n/2), seeded for n=1 */
        arb_set(pown, sq4a);                         /* (4 sqrt a)^1 */
        arb_const_pi(ghalf, wp); arb_sqrt(ghalf, ghalf, wp);  /* Gamma(1/2)=sqrt pi */
        arb_one(gh2);                                /* Gamma(1) = 1 (for n=2) */

        for (n = 1; n <= N; n++)
        {
            /* coefficient B_n(3/4)/n! * (4 sqrt a)^n * Gamma(n/2) */
            arb_mul(bp, bpv + n, pown, wp);
            if (n % 2 == 1)
                arb_mul(bp, bp, ghalf, wp);          /* Gamma(n/2), n odd */
            else
                arb_mul(bp, bp, gh2, wp);            /* Gamma(n/2), n even */

            /* / (s + n - 1) */
            acb_add_ui(t, s, (ulong)(n - 1), wp);
            acb_set_arb(term, bp);
            acb_div(term, term, t, wp);
            acb_add(bracket, bracket, term, wp);

            /* advance recurrences for n -> n+1 */
            arb_mul(pown, pown, sq4a, wp);           /* (4 sqrt a)^{n+1} */
            /* advance the chain of parity n to its next member (used at n+2):
               Gamma((n+2)/2) = (n/2) Gamma(n/2). */
            arb_set_ui(gh, (ulong) n);
            arb_mul_2exp_si(gh, gh, -1);             /* n/2 */
            if (n % 2 == 1)
                arb_mul(ghalf, ghalf, gh, wp);       /* odd chain: Gamma((n+2)/2) */
            else
                arb_mul(gh2, gh2, gh, wp);           /* even chain */
        }

        _arb_vec_clear(bpv, len);
        arb_clear(pown); arb_clear(ghalf); arb_clear(gh2); arb_clear(x34);
    }

    /* prefactor a^{(s-1)/2} / (4 sqrt(pi) Gamma(s/2)) */
    arb_log(alog, a, wp);
    acb_set_arb(pref, alog);
    acb_mul(pref, pref, sm1, wp);
    acb_mul_2exp_si(pref, pref, -1);      /* (s-1)/2 * log a */
    acb_exp(pref, pref, wp);              /* a^{(s-1)/2} */

    arb_const_pi(gh, wp);
    arb_sqrt(gh, gh, wp);
    arb_mul_2exp_si(gh, gh, 2);           /* 4 sqrt pi */
    acb_div_arb(pref, pref, gh, wp);

    acb_mul_2exp_si(shalf, s, -1);
    acb_rgamma(g, shalf, wp);             /* 1/Gamma(s/2) */
    acb_mul(pref, pref, g, wp);

    acb_mul(res, pref, bracket, wp);

    S_bound(bound, &S, N, SECOND_ZETA_BOUND_PREC);
    arb_add_error_mag(acb_realref(res), bound);
    if (!acb_is_real(s))
        arb_add_error_mag(acb_imagref(res), bound);

    acb_set_round(res, res, prec);

    acb_clear(shalf); acb_clear(sm1); acb_clear(pref);
    acb_clear(bracket); acb_clear(term); acb_clear(t); acb_clear(g);
    arb_clear(alog); arb_clear(sqa); arb_clear(sq4a);
    arb_clear(bp); arb_clear(gh); arb_clear(u);
    mag_clear(bound);
    S_pre_clear(&S);
}


/*****************************************************************************/
/*                               Main routine                                */
/*****************************************************************************/

/* extra guard bits from the 1/Gamma(s/2) magnitude (A/P cancellation) */
static slong
cancellation_guard(const acb_t s, slong prec)
{
    acb_t shalf, rg;
    arb_t lg;
    mag_t m;
    slong bits;

    acb_init(shalf); acb_init(rg); arb_init(lg); mag_init(m);

    acb_mul_2exp_si(shalf, s, -1);
    acb_rgamma(rg, shalf, 32 + prec / 8);   /* cheap, low precision */
    acb_get_mag(m, rg);                      /* upper bound on 1/|Gamma(s/2)| */

    bits = 0;
    if (!mag_is_zero(m) && !mag_is_inf(m))
    {
        /* bits ~ max(0, log2(|1/Gamma(s/2)|)) */
        double l = mag_get_d_log2_approx(m);

        if (l > 1e9)
            bits = WORD_MAX;
        else if (l > 0)
            bits = (slong)(l + 1);
    }

    acb_clear(shalf); acb_clear(rg); arb_clear(lg); mag_clear(m);
    return bits;
}

/* Choose the free parameter a > 0. This needs less than about 0.444 / prec
   due to the exp(-kappa/a) factor in the error term of the asymptotic
   expansion. */

static void
choose_a(arb_t a, const acb_t s, slong prec)
{
    double aval = 0.4 / (double) FLINT_MAX(prec, 1);

    if (aval < 1e-7) aval = 1e-7;
    if (aval > 0.0625) aval = 0.0625;

    arb_set_d(a, aval);
}

void
acb_dirichlet_secondary_zeta(acb_t res, const acb_t s, slong prec)
{
    acb_t A, P, E, S;
    arb_t a;
    slong wp, guard;

    if (acb_contains_int(s))
    {
        /* Double pole at s = 1 */
        if (arb_contains_si(acb_realref(s), 1))
        {
            acb_indeterminate(res);
            return;
        }

        /* Simple poles at s = -1, -3, -5, ... */
        if (!arb_is_nonnegative(acb_realref(s)))
        {
            arb_t t;
            arb_init(t);
            arb_add_ui(t, acb_realref(s), 1, 2 * prec);
            arb_mul_2exp_si(t, t, -1);

            if (arb_contains_int(t))
            {
                acb_indeterminate(res);
                arb_clear(t);
                return;
            }

            arb_clear(t);
        }
    }

    /* Z(-2n) = (-1)^n (8 - E_{2n}) / 2^{2n+3} */
    if (acb_is_int(s) && arf_sgn(arb_midref(acb_realref(s))) <= 0)
    {
        /* Could be handled, if we really wanted to... */
        if (arf_cmp_si(arb_midref(acb_realref(s)), WORD_MIN) < 0)
        {
            acb_indeterminate(res);
            return;
        }

        slong n = -(arf_get_si(arb_midref(acb_realref(s)), ARF_RND_DOWN) / 2);
        arb_ptr val = acb_realref(res);

        arb_euler_number_ui(val, 2 * n, prec + 8);
        arb_sub_ui(val, val, 8, prec);
        arb_mul_2exp_si(val, val, -(2 * n + 3));
        if ((n & 1) == 0)
            arb_neg(val, val);

        arb_zero(acb_imagref(res));
        return;
    }

    guard = cancellation_guard(s, prec);
    /* flint_printf("guard %wd\n", guard); */

    /* Abort if the input looks unreasonable; the user can retry with
       higher precision if they really want this value. */
    if (guard > 100 * prec)
    {
        acb_indeterminate(res);
        return;
    }

    wp = prec + guard + 16 + 2 * FLINT_BIT_COUNT(prec);

    acb_init(A); acb_init(P); acb_init(E); acb_init(S);
    arb_init(a);

    choose_a(a, s, wp);

    /* Z = A - P + E - S */
    _acb_dirichlet_secondary_zeta_term_E(E, s, a, wp);
    _acb_dirichlet_secondary_zeta_term_A(A, s, a, wp);
    _acb_dirichlet_secondary_zeta_term_P(P, s, a, wp);
    _acb_dirichlet_secondary_zeta_term_S(S, s, a, wp);
    acb_sub(res, A, P, wp);
    acb_add(res, res, E, wp);
    acb_sub(res, res, S, wp);

    acb_set_round(res, res, prec);

    acb_clear(A); acb_clear(P); acb_clear(E); acb_clear(S);
    arb_clear(a);
}

