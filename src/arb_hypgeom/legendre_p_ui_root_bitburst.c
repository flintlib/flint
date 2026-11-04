/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "arb.h"
#include "arb_hypgeom.h"
#include "impl.h"

/*
    The root x_k of P_n by bit-burst Newton steps.

    Stage 0 evaluates P_n and P_(n-1) at a point x_0 = A / 2^B close to
    x_k by the three-term recurrence.  With u_j = j! 2^(Bj) P_j(x_0), the
    recurrence reads

        u_(j+1) = (2j+1) A u_j - j^2 2^(2B) u_(j-1),

    with no division.  Below LIMB_MAX_PREC bits, and when n is small
    enough that a numerator A of FLINT_BITS - 2 - bits(2n) bits leaves
    word-size multipliers, it runs on fixed-size limb windows
    (stage0_limbs), with a rigorous bound for the dropped limbs;
    otherwise exact products of blocks of steps by binary splitting are
    applied to a ball vector (stage0_blocks).

    Each following stage rounds the Newton correction -P_n/P_n' to a
    number delta = D / 2^E with about twice as many correct bits as the
    current point, and moves (P_n, P_n') from x to x + delta by the
    Taylor series P_n(x + t) = sum c_j t^j, whose coefficients satisfy,
    from the Legendre equation (1 - x^2) y'' - 2 x y' + n (n+1) y = 0,

        (1 - x^2) (j+1) (j+2) c_(j+2)
            = 2 x (j+1)^2 c_(j+1) + (j (j+1) - n (n+1)) c_j.

    With x = A / 2^B, a_j = c_j delta^j and sigma = 2^(2B) - A^2, this is

        gamma_j a_(j+2) = alpha_j a_(j+1) + beta_j a_j,
        alpha_j = 2 A D (j+1)^2 2^(E-B),
        beta_j  = D^2 (j (j+1) - n (n+1)),
        gamma_j = sigma (j+1) (j+2) 2^(2(E-B)),

    summed directly in fixed point (transfer_fixed) or by binary
    splitting (transfer_blocks), with the Cauchy bound for the terms
    beyond the K-th (_arb_hypgeom_legendre_p_ui_taylor_tail).  The size of
    the coefficients doubles from one stage to the next while the number
    of terms needed halves.

    The last stage gives P_n and P_n' at a point within about
    2^(-prec/2) of x_k, and an interval Newton step on the Taylor
    polynomial of degree 2 with its tail bound encloses x_k; containment
    in the isolating enclosure of Petras identifies the root.
*/

/* Default tuning parameters: the working precision below which stage 0
   uses limb windows, K * prec above which a transfer uses binary
   splitting, and the entry size (bits) from which the block products
   use fmpz_mat_mul_fft_small (measured on x86-64). */
#define LIMB_MAX_PREC 130000
#define TRANSFER_BLOCKS_MIN 1e10
#define MAT_FFT_BITS 10000

/* An approximation of x_k in double precision (for choosing x_0 only):
   Petras' approximation and Newton steps by the recurrence. */
static double
root_d(ulong n, ulong k)
{
    double phi, psi, t, s, x, p0, p1, p2, dp;
    ulong j;
    int it;

    t = 4.0 * n + 2.0;
    phi = (4.0 * k + 3.0) * 3.141592653589793 / t;
    s = sin(phi);
    psi = phi + 2.0 * cos(phi) / (t * t * s) * (1.0 - 11.0 / (t * t * s * s));
    x = cos(psi);

    for (it = 0; it < 3; it++)
    {
        p0 = 1.0;
        p1 = x;
        for (j = 1; j < n; j++)
        {
            p2 = ((2.0 * j + 1.0) * x * p1 - j * p0) / (j + 1.0);
            p0 = p1;
            p1 = p2;
        }
        dp = n * (p0 - x * p1) / (1.0 - x * x);
        x -= p1 / dp;
    }

    return x;
}

/* y = P_n(x), yp = P_n'(x) from enclosures of P_n, P_(n-1) at x = A/2^B */
static void
deriv_from_pair(arb_t y, arb_t yp, const arb_t pn, const arb_t pn1, ulong n, const fmpz_t A, slong B, slong prec)
{
    arb_t x, t, u;

    arb_init(x);
    arb_init(t);
    arb_init(u);

    arb_set_fmpz(x, A);
    arb_mul_2exp_si(x, x, -B);

    /* P_n' = n (P_(n-1) - x P_n) / (1 - x^2) */
    arb_mul(t, x, pn, prec);
    arb_sub(t, pn1, t, prec);
    arb_mul(u, x, x, prec);
    arb_sub_ui(u, u, 1, prec);
    arb_neg(u, u);
    arb_div(t, t, u, prec);
    arb_mul_ui(yp, t, n, prec);
    arb_set(y, pn);

    arb_clear(x);
    arb_clear(t);
    arb_clear(u);
}

/* An upper bound for 1/sin(theta) = 2^B / sqrt(2^(2B) - A^2) */
static void
inv_sin_bound(mag_t res, const fmpz_t A, slong B)
{
    fmpz_t s;
    mag_t t;

    fmpz_init(s);
    mag_init(t);

    fmpz_one(s);
    fmpz_mul_2exp(s, s, 2 * B);
    fmpz_submul(s, A, A);

    if (fmpz_sgn(s) <= 0)
    {
        mag_inf(res);
    }
    else
    {
        mag_set_fmpz_lower(t, s);
        mag_sqrt_lower(t, t);
        mag_one(res);
        mag_mul_2exp_si(res, res, B);
        mag_div(res, res, t);
    }

    fmpz_clear(s);
    mag_clear(t);
}

/* Stage 0 on limb windows.  With beta = 2^FLINT_BITS, u_(j-1) and u_j
   are kept as W-limb two's complement windows at a common limb exponent
   S (value = window * beta^S); each step computes
   u_(j+1) = c1 u_j - c2 beta u_(j-1) with c1 = (2j+1) A,
   c2 = j^2 2^(2B) / beta, and drops the low limbs as the values grow.

   Error bound.  When step j (computing u_(j+1)) drops limbs, the new
   u_(j+1) and the new u_j (the next u_(j-1)) each lose less than one
   unit beta^S of the new exponent.  In P_i = u_i / F_i, F_i = i! 2^(Bi),
   the first is an error of at most beta^S / F_(j+1) in P_(j+1); the
   second reaches u_(j+2) through the multiplier (j+1)^2 2^(2B), an error
   of at most (j+1)^2 2^(2B) beta^S / F_(j+2) <= (j+1) 2^B beta^S / F_(j+1)
   in P_(j+2); together at most (j+2) 2^B beta^S / F_(j+1).  The final
   u_(n-1) carries its own loss beta^S / F_(n-1).

   The computed sequence is the exact one plus, for each error e in P_i,
   the solution of the recurrence with y_i = e, y_(i-1) = 0.  With
   x = cos(theta), s = sin(theta) and
       E_m = m (y_m^2 - 2 x y_m y_(m-1) + y_(m-1)^2)
           = m (y_(m-1) - x y_m)^2 + m s^2 y_m^2,
   the recurrence (m+1) y_(m+1) = (2m+1) x y_m - m y_(m-1) gives
       E_(m+1) = m^2/(m+1) (y_(m-1) - x y_m)^2 + (m+1) s^2 y_m^2
              <= (m+1)/m E_m,
   so E_m <= (m/i) E_i = m e^2 and s^2 y_m^2 <= E_m / m <= e^2: every
   error reaches P_n and P_(n-1) multiplied by at most 1/s.

   Requires n >= 2, 2B >= FLINT_BITS (for the shift in c2) and
   B + bits(2n) <= FLINT_BITS - 2, so that c1 < 2n 2^B and
   c2 < 2^(2 bits(n) + 2B - FLINT_BITS) fit in a word. */
static void
stage0_limbs(arb_t pn, arb_t pn1, ulong n, ulong A, slong B, slong prec)
{
    slong L, W, S, j, i, red, d, off;
    ulong c1, c2, cy, bw, ext;
    nn_ptr buf[3], U0, U1, wv;
    slong cur, prv, nxt;
    mag_t G, sum, t;
    fmpz_t z, Af;
    arb_t f;

    L = (prec + FLINT_BITS - 1) / FLINT_BITS + 2;
    W = L + 2;

    mag_init(G);
    mag_init(sum);
    mag_init(t);
    fmpz_init(z);
    fmpz_init(Af);
    arb_init(f);

    /* three buffers of W + 3 limbs; u_(j-1), u_j are W-limb views into
       them (the limbs above a view are its sign extension) */
    for (i = 0; i < 3; i++)
        buf[i] = flint_calloc(W + 3, sizeof(ulong));

    /* u_0 = 1, u_1 = A, at the exponent S = -(W - 2) */
    S = -(W - 2);
    cur = 0;
    prv = 2;
    buf[prv][W - 2] = 1;
    buf[cur][W - 2] = A;
    U0 = buf[prv];
    U1 = buf[cur];
    d = 0;

    /* G >= 1 / F_1 = 2^-B */
    mag_one(G);
    mag_mul_2exp_si(G, G, -B);

    for (j = 1; j < (slong) n; j++)
    {
        c1 = (2 * (ulong) j + 1) * A;
        c2 = ((ulong) j * (ulong) j) << (2 * B - FLINT_BITS);
        nxt = 3 - cur - prv;
        wv = buf[nxt];

        /* wv = c1 u_j - c2 B u_(j-1)  (W + 3 limbs, two's complement) */
        cy = mpn_mul_1(wv, U1, W, c1);
        wv[W] = cy;
        wv[W + 1] = 0;
        wv[W + 2] = 0;
        if ((slong) U1[W - 1] < 0)
            mpn_sub_1(wv + W, wv + W, 3, c1);
        bw = mpn_submul_1(wv + 1, U0, W, c2);
        if ((slong) U0[W - 1] < 0)
            mpn_add_1(wv + W + 1, wv + W + 1, 2, c2);
        mpn_sub_1(wv + W + 1, wv + W + 1, 2, bw);

        /* redundant top limbs (sign extension) */
        red = 0;
        for (i = W + 2; i > 0; i--)
        {
            ext = ((slong) wv[i - 1] < 0) ? UWORD_MAX : 0;
            if (wv[i] != ext)
                break;
            red++;
        }
        d = FLINT_MAX(3 - red, 0);

        /* the new u_(j-1) is u_j shifted down by d limbs: its view moves
           up by d limbs within its buffer, above which the sign extension
           must be present */
        if (d > 0)
        {
            ext = ((slong) U1[W - 1] < 0) ? UWORD_MAX : 0;
            off = U1 - buf[cur];
            if (off + d + W > W + 3)
            {
                for (i = 0; i < W; i++)
                    buf[cur][i] = U1[i];
                U1 = buf[cur];
                off = 0;
            }
            for (i = off + W; i < off + W + d; i++)
                buf[cur][i] = ext;
            U0 = U1 + d;
        }
        else
        {
            U0 = U1;
        }

        prv = cur;
        cur = nxt;
        U1 = wv + d;
        S += d;

        /* G >= 1 / F_(j+1) */
        mag_div_ui(G, G, j + 1);
        mag_mul_2exp_si(G, G, -B);

        if (d > 0)
        {
            mag_mul_ui(t, G, j + 2);
            mag_mul_2exp_si(t, t, B + FLINT_BITS * S);
            mag_add(sum, sum, t);
        }
    }

    /* P_n = u_n B^S / F_n, P_(n-1) = u_(n-1) B^S / F_(n-1) */
    arb_fac_ui(f, n, prec);

    fmpz_set_signed_ui_array(z, U1, W);
    arb_set_fmpz(pn, z);
    arb_mul_2exp_si(pn, pn, FLINT_BITS * S - B * (slong) n);
    arb_div(pn, pn, f, prec);

    fmpz_set_signed_ui_array(z, U0, W);
    arb_set_fmpz(pn1, z);
    arb_mul_2exp_si(pn1, pn1, FLINT_BITS * S - B * ((slong) n - 1));
    arb_div_ui(f, f, n, prec);
    arb_div(pn1, pn1, f, prec);

    /* the errors, through 1/s */
    fmpz_set_ui(Af, A);
    inv_sin_bound(t, Af, B);
    mag_mul(sum, sum, t);
    arb_add_error_mag(pn, sum);
    arb_add_error_mag(pn1, sum);

    if (d > 0)
    {
        /* the final u_(n-1): B^S / F_(n-1) = n 2^B B^S / F_n */
        mag_mul_ui(t, G, n);
        mag_mul_2exp_si(t, t, B + FLINT_BITS * S);
        arb_add_error_mag(pn1, t);
    }

    for (i = 0; i < 3; i++)
        flint_free(buf[i]);

    mag_clear(G);
    mag_clear(sum);
    mag_clear(t);
    fmpz_clear(z);
    fmpz_clear(Af);
    arb_clear(f);
}

/* C = H L for 2 x 2 or 4 x 2 by 2 x 2 integer matrices (C not aliased) */
static void
mat_mul(fmpz_mat_t C, const fmpz_mat_t H, const fmpz_mat_t L, slong fft_bits)
{
    slong hb, lb;

    hb = FLINT_ABS(fmpz_mat_max_bits(H));
    lb = FLINT_ABS(fmpz_mat_max_bits(L));

    if (FLINT_MIN(hb, lb) >= fft_bits && fmpz_mat_mul_fft_small(C, H, L))
        return;

    fmpz_mat_mul(C, H, L);
}

/* R = M_(hi-1) ... M_lo, M_j = [[(2j+1) A, -j^2 2^(2B)], [1, 0]]
   (u_(j+1), u_j) = M_j (u_j, u_(j-1)), exactly */
static void
stage0_bsplit(fmpz_mat_t R, slong lo, slong hi, const fmpz_t A, slong B, slong fft_bits)
{
    if (hi - lo == 1)
    {
        ulong j = lo;
        fmpz_mul_ui(fmpz_mat_entry(R, 0, 0), A, 2 * j + 1);
        fmpz_set_ui(fmpz_mat_entry(R, 0, 1), j);
        fmpz_mul_ui(fmpz_mat_entry(R, 0, 1), fmpz_mat_entry(R, 0, 1), j);
        fmpz_mul_2exp(fmpz_mat_entry(R, 0, 1), fmpz_mat_entry(R, 0, 1), 2 * B);
        fmpz_neg(fmpz_mat_entry(R, 0, 1), fmpz_mat_entry(R, 0, 1));
        fmpz_one(fmpz_mat_entry(R, 1, 0));
        fmpz_zero(fmpz_mat_entry(R, 1, 1));
    }
    else
    {
        fmpz_mat_t H, L;
        slong mid = lo + (hi - lo) / 2;

        fmpz_mat_init(H, 2, 2);
        fmpz_mat_init(L, 2, 2);
        stage0_bsplit(L, lo, mid, A, B, fft_bits);
        stage0_bsplit(H, mid, hi, A, B, fft_bits);
        mat_mul(R, H, L, fft_bits);
        fmpz_mat_clear(H);
        fmpz_mat_clear(L);
    }
}

/* Stage 0 by exact blocks of steps of about prec bits, applied to the
   ball vector (u_j, u_(j-1)).  The radius grows by about the norm of a
   block's transition matrix per block, which is O(1/s) (see
   stage0_limbs); the working precision has extra bits for that. */
static void
stage0_blocks(arb_t pn, arb_t pn1, ulong n, const fmpz_t A, slong B, slong prec, slong fft_bits)
{
    slong h, lo, hi, nblocks, wp;
    fmpz_mat_t R;
    arb_t v0, v1, t, u, f;
    double xd, lg;

    h = FLINT_MAX(1, (prec + 64) / (B + 2 * FLINT_BIT_COUNT(n) + 4));
    nblocks = ((slong) n - 1 + h - 1) / h;

    xd = ldexp(fmpz_get_d(A), -B);
    lg = 0.5 * log2(1.0 / FLINT_MAX(1.0 - xd * xd, 1e-300)) + 3;
    wp = prec + nblocks * (slong) ceil(lg);

    fmpz_mat_init(R, 2, 2);
    arb_init(v0);
    arb_init(v1);
    arb_init(t);
    arb_init(u);
    arb_init(f);

    /* (u_1, u_0) = (A, 1) */
    arb_set_fmpz(v0, A);
    arb_one(v1);

    for (lo = 1; lo < (slong) n; lo = hi)
    {
        hi = FLINT_MIN((slong) n, lo + h);
        stage0_bsplit(R, lo, hi, A, B, fft_bits);
        arb_mul_fmpz(t, v0, fmpz_mat_entry(R, 0, 0), wp);
        arb_addmul_fmpz(t, v1, fmpz_mat_entry(R, 0, 1), wp);
        arb_mul_fmpz(u, v0, fmpz_mat_entry(R, 1, 0), wp);
        arb_addmul_fmpz(u, v1, fmpz_mat_entry(R, 1, 1), wp);
        arb_swap(v0, t);
        arb_swap(v1, u);
    }

    /* P_n = u_n / (n! 2^(Bn)), P_(n-1) = u_(n-1) / ((n-1)! 2^(B(n-1))) */
    arb_fac_ui(f, n, prec);
    arb_mul_2exp_si(pn, v0, -B * (slong) n);
    arb_div(pn, pn, f, prec);
    arb_div_ui(f, f, n, prec);
    arb_mul_2exp_si(pn1, v1, -B * ((slong) n - 1));
    arb_div(pn1, pn1, f, prec);

    fmpz_mat_clear(R);
    arb_clear(v0);
    arb_clear(v1);
    arb_clear(t);
    arb_clear(u);
    arb_clear(f);
}

/* The coefficients of a transfer from x = A/2^B by delta = D/2^E:
   alpha_j = twoAD (j+1)^2 2^sh, beta_j = D2 (j (j+1) - lambda),
   gamma_j = sig (j+1) (j+2) 2^(2 sh), sh = E - B. */
typedef struct
{
    fmpz_t twoAD;
    fmpz_t D2;
    fmpz_t sig;
    slong sh;
    ulong lambda;
    slong fft_bits;
}
transfer_ctx_struct;

typedef transfer_ctx_struct transfer_ctx_t[1];

static void
transfer_ctx_init(transfer_ctx_t c, const fmpz_t A, slong B, const fmpz_t D, slong E, ulong n, slong fft_bits)
{
    fmpz_init(c->twoAD);
    fmpz_init(c->D2);
    fmpz_init(c->sig);

    fmpz_mul(c->twoAD, A, D);
    fmpz_mul_2exp(c->twoAD, c->twoAD, 1);
    fmpz_mul(c->D2, D, D);
    fmpz_one(c->sig);
    fmpz_mul_2exp(c->sig, c->sig, 2 * B);
    fmpz_submul(c->sig, A, A);
    c->sh = E - B;
    c->lambda = n * (n + 1);
    c->fft_bits = fft_bits;
}

static void
transfer_ctx_clear(transfer_ctx_t c)
{
    fmpz_clear(c->twoAD);
    fmpz_clear(c->D2);
    fmpz_clear(c->sig);
}

/* s = sum_(j=0)^K a_j, sp = sum_(j=1)^K j a_j from a_0, a_1, in fixed
   point: the a_j are integers in units 2^-W (exact products, truncating
   divisions, the lengths of the integers decreasing with the terms).
   Each division adds less than one unit; with E_j the error of a_j in
   units,
       E_(j+2) <= (|alpha_j| E_(j+1) + |beta_j| E_j) / gamma_j + 1,
   a contracting recurrence since |delta| is small. */
static void
transfer_fixed(arb_t s, arb_t sp, const arb_t a0, const arb_t a1, const transfer_ctx_t c, slong K, slong prec)
{
    fmpz_t f0, f1, f2, fs, fsp, t, al, be, ga;
    mag_t E0, E1, E2, ES, ESP, mal, mbe, msig, u, v;
    slong j, W, sh = c->sh;
    arb_t w;

    fmpz_init(f0);
    fmpz_init(f1);
    fmpz_init(f2);
    fmpz_init(fs);
    fmpz_init(fsp);
    fmpz_init(t);
    fmpz_init(al);
    fmpz_init(be);
    fmpz_init(ga);
    mag_init(E0);
    mag_init(E1);
    mag_init(E2);
    mag_init(ES);
    mag_init(ESP);
    mag_init(mal);
    mag_init(mbe);
    mag_init(msig);
    mag_init(u);
    mag_init(v);
    arb_init(w);

    /* y' = sp / delta is needed to about |a_1| 2^-prec, and s (the new
       residual P_n) only to |P_n'| 2^-prec, which is coarser */
    W = prec + 16 - arf_abs_bound_lt_2exp_si(arb_midref(a1));

    arb_mul_2exp_si(w, a0, W);
    arf_get_fmpz(f0, arb_midref(w), ARF_RND_DOWN);
    mag_add_ui(E0, arb_radref(w), 1);
    arb_mul_2exp_si(w, a1, W);
    arf_get_fmpz(f1, arb_midref(w), ARF_RND_DOWN);
    mag_add_ui(E1, arb_radref(w), 1);

    fmpz_add(fs, f0, f1);
    fmpz_set(fsp, f1);
    mag_add(ES, E0, E1);
    mag_set(ESP, E1);

    mag_set_fmpz(mal, c->twoAD);
    mag_set_fmpz(mbe, c->D2);
    mag_set_fmpz_lower(msig, c->sig);

    for (j = 0; j + 2 <= K; j++)
    {
        /* f2 = trunc(trunc((al f1 2^sh + be f0) / 2^(2 sh)) / ga),
           al = twoAD (j+1)^2, ga = sig (j+1) (j+2) */
        fmpz_mul_ui(al, c->twoAD, (ulong) (j + 1) * (j + 1));
        fmpz_mul_si(be, c->D2, (slong) ((ulong) j * (j + 1)) - (slong) c->lambda);
        fmpz_mul_ui(ga, c->sig, (ulong) (j + 1) * (j + 2));
        fmpz_mul(t, f1, al);
        fmpz_mul_2exp(t, t, sh);
        fmpz_addmul(t, f0, be);
        fmpz_tdiv_q_2exp(t, t, 2 * sh);
        fmpz_tdiv_q(f2, t, ga);

        fmpz_add(fs, fs, f2);
        fmpz_addmul_ui(fsp, f2, j + 2);

        /* E2 = (|al| 2^-sh E1 + |be| 2^(-2 sh) E0) / ga + 1 */
        mag_mul_ui(u, mal, (ulong) (j + 1) * (j + 1));
        mag_mul(u, u, E1);
        mag_mul_2exp_si(u, u, sh);
        mag_mul_ui(v, mbe, FLINT_ABS((slong) ((ulong) j * (j + 1)) - (slong) c->lambda));
        mag_mul(v, v, E0);
        mag_add(u, u, v);
        mag_mul_2exp_si(u, u, -2 * sh);
        mag_mul_ui_lower(v, msig, (ulong) (j + 1) * (j + 2));
        mag_div(E2, u, v);
        mag_add_ui(E2, E2, 1);

        mag_add(ES, ES, E2);
        mag_mul_ui(u, E2, j + 2);
        mag_add(ESP, ESP, u);

        fmpz_swap(f0, f1);
        fmpz_swap(f1, f2);
        mag_swap(E0, E1);
        mag_swap(E1, E2);
    }

    arb_set_fmpz(s, fs);
    mag_add(arb_radref(s), arb_radref(s), ES);
    arb_mul_2exp_si(s, s, -W);
    arb_set_fmpz(sp, fsp);
    mag_add(arb_radref(sp), arb_radref(sp), ESP);
    arb_mul_2exp_si(sp, sp, -W);

    fmpz_clear(f0);
    fmpz_clear(f1);
    fmpz_clear(f2);
    fmpz_clear(fs);
    fmpz_clear(fsp);
    fmpz_clear(t);
    fmpz_clear(al);
    fmpz_clear(be);
    fmpz_clear(ga);
    mag_clear(E0);
    mag_clear(E1);
    mag_clear(E2);
    mag_clear(ES);
    mag_clear(ESP);
    mag_clear(mal);
    mag_clear(mbe);
    mag_clear(msig);
    mag_clear(u);
    mag_clear(v);
    arb_clear(w);
}

/* A block of transfer steps j in [lo, hi): with v_j = (a_(j+1), a_j),
   v_hi = X v_lo / z, and the terms a_(lo+2), ..., a_(hi+1) contribute
   Y_0 v_lo / z to s and Y_1 v_lo / z to sp; XY = [X; Y] (4 x 2).
   A single step is X = [[alpha, beta], [gamma, 0]],
   Y = [[alpha, beta], [(j+2) alpha, (j+2) beta]], z = gamma; blocks
   combine as X = X_H X_L, Y = Y_H X_L + z_H Y_L, z = z_H z_L.  Common
   powers of two are removed from all entries and z. */
static void
transfer_bsplit(fmpz_mat_t XY, fmpz_t z, slong lo, slong hi, const transfer_ctx_t c)
{
    slong i, v;

    if (hi - lo == 1)
    {
        ulong j = lo;
        fmpz * al = fmpz_mat_entry(XY, 0, 0);
        fmpz * be = fmpz_mat_entry(XY, 0, 1);

        fmpz_mul_ui(al, c->twoAD, (j + 1) * (j + 1));
        fmpz_mul_2exp(al, al, c->sh);
        fmpz_mul_si(be, c->D2, (slong) (j * (j + 1)) - (slong) c->lambda);
        fmpz_mul_ui(z, c->sig, (j + 1) * (j + 2));
        fmpz_mul_2exp(z, z, 2 * c->sh);
        fmpz_set(fmpz_mat_entry(XY, 1, 0), z);
        fmpz_zero(fmpz_mat_entry(XY, 1, 1));
        fmpz_set(fmpz_mat_entry(XY, 2, 0), al);
        fmpz_set(fmpz_mat_entry(XY, 2, 1), be);
        fmpz_mul_ui(fmpz_mat_entry(XY, 3, 0), al, j + 2);
        fmpz_mul_ui(fmpz_mat_entry(XY, 3, 1), be, j + 2);
    }
    else
    {
        fmpz_mat_t H, L, X;
        fmpz_t zH, zL;
        slong mid = lo + (hi - lo) / 2;

        fmpz_mat_init(H, 4, 2);
        fmpz_mat_init(L, 4, 2);
        fmpz_init(zH);
        fmpz_init(zL);

        transfer_bsplit(L, zL, lo, mid, c);
        transfer_bsplit(H, zH, mid, hi, c);

        /* [X; Y] = [X_H; Y_H] X_L, then Y += z_H Y_L */
        fmpz_mat_window_init(X, L, 0, 0, 2, 2);
        mat_mul(XY, H, X, c->fft_bits);
        fmpz_mat_window_clear(X);
        for (i = 2; i < 4; i++)
        {
            fmpz_addmul(fmpz_mat_entry(XY, i, 0), zH, fmpz_mat_entry(L, i, 0));
            fmpz_addmul(fmpz_mat_entry(XY, i, 1), zH, fmpz_mat_entry(L, i, 1));
        }
        fmpz_mul(z, zH, zL);

        fmpz_mat_clear(H);
        fmpz_mat_clear(L);
        fmpz_clear(zH);
        fmpz_clear(zL);
    }

    v = fmpz_val2(z);
    for (i = 0; i < 8 && v > 0; i++)
        if (!fmpz_is_zero(fmpz_mat_entry(XY, i / 2, i % 2)))
            v = FLINT_MIN(v, (slong) fmpz_val2(fmpz_mat_entry(XY, i / 2, i % 2)));
    if (v > 0)
    {
        fmpz_tdiv_q_2exp(z, z, v);
        for (i = 0; i < 4; i++)
        {
            fmpz_tdiv_q_2exp(fmpz_mat_entry(XY, i, 0), fmpz_mat_entry(XY, i, 0), v);
            fmpz_tdiv_q_2exp(fmpz_mat_entry(XY, i, 1), fmpz_mat_entry(XY, i, 1), v);
        }
    }
}

/* The same sums as transfer_fixed by exact blocks of steps of about prec
   bits, applied to (a_1, a_0) in ball arithmetic. */
static void
transfer_blocks(arb_t s, arb_t sp, const arb_t a0, const arb_t a1, const transfer_ctx_t c, ulong n, slong K, slong prec)
{
    fmpz_mat_t XY;
    fmpz_t z;
    arb_t b0, b1, t, u, zz;
    slong lo, hi, h, E2;

    E2 = 2 * (c->sh + (slong) fmpz_bits(c->sig) / 2) + 3 * FLINT_BIT_COUNT(n) + 8;
    h = FLINT_MAX(1, (prec + 64) / FLINT_MAX(E2, 1));

    fmpz_mat_init(XY, 4, 2);
    fmpz_init(z);
    arb_init(b0);
    arb_init(b1);
    arb_init(t);
    arb_init(u);
    arb_init(zz);

    arb_set(b0, a0);
    arb_set(b1, a1);
    arb_add(s, a0, a1, prec);
    arb_set(sp, a1);

    for (lo = 0; lo < K - 1; lo = hi)
    {
        hi = FLINT_MIN(K - 1, lo + h);
        transfer_bsplit(XY, z, lo, hi, c);
        arb_set_fmpz(zz, z);

        arb_mul_fmpz(t, b1, fmpz_mat_entry(XY, 2, 0), prec);
        arb_addmul_fmpz(t, b0, fmpz_mat_entry(XY, 2, 1), prec);
        arb_div(t, t, zz, prec);
        arb_add(s, s, t, prec);

        arb_mul_fmpz(t, b1, fmpz_mat_entry(XY, 3, 0), prec);
        arb_addmul_fmpz(t, b0, fmpz_mat_entry(XY, 3, 1), prec);
        arb_div(t, t, zz, prec);
        arb_add(sp, sp, t, prec);

        arb_mul_fmpz(t, b1, fmpz_mat_entry(XY, 0, 0), prec);
        arb_addmul_fmpz(t, b0, fmpz_mat_entry(XY, 0, 1), prec);
        arb_mul_fmpz(u, b1, fmpz_mat_entry(XY, 1, 0), prec);
        arb_addmul_fmpz(u, b0, fmpz_mat_entry(XY, 1, 1), prec);
        arb_div(b1, t, zz, prec);
        arb_div(b0, u, zz, prec);
    }

    fmpz_mat_clear(XY);
    fmpz_clear(z);
    arb_clear(b0);
    arb_clear(b1);
    arb_clear(t);
    arb_clear(u);
    arb_clear(zz);
}

/* (y, yp) = (P_n, P_n') at x = A/2^B -> at x + D/2^E (E >= B), with the
   terms of degree at most K and the tail bound if K <= n; A is replaced
   by A 2^(E-B) + D.  Returns 0 if the tail bound fails. */
static int
transfer(arb_t y, arb_t yp, fmpz_t A, slong B, const fmpz_t D, slong E, ulong n, slong K, int blocks, slong fft_bits, slong prec)
{
    transfer_ctx_t c;
    arb_t x, dl, a1, s, sp;
    mag_t r, tau, T0, T1;
    slong wp;
    int ok = 1;

    transfer_ctx_init(c, A, B, D, E, n, fft_bits);
    arb_init(x);
    arb_init(dl);
    arb_init(a1);
    arb_init(s);
    arb_init(sp);
    mag_init(r);
    mag_init(tau);
    mag_init(T0);
    mag_init(T1);

    arb_set_fmpz(x, A);
    arb_mul_2exp_si(x, x, -B);
    arb_set_fmpz(dl, D);
    arb_mul_2exp_si(dl, dl, -E);

    /* a_0 = y, a_1 = y' delta */
    wp = prec + E - (slong) fmpz_bits(D);
    arb_mul(a1, yp, dl, wp);

    if (blocks)
        transfer_blocks(s, sp, y, a1, c, n, K, wp);
    else
        transfer_fixed(s, sp, y, a1, c, K, prec);

    if (K <= (slong) n)
    {
        _arb_hypgeom_legendre_p_ui_disk_radius(r, x, n);
        arb_get_mag(tau, dl);
        ok = _arb_hypgeom_legendre_p_ui_taylor_tail(T0, T1, n, x, r, tau, K);
        arb_add_error_mag(s, T0);
    }

    /* y' = sp / delta (+ the tail of the derivative) */
    arb_div(yp, sp, dl, prec);
    if (K <= (slong) n)
        arb_add_error_mag(yp, T1);
    arb_set_round(y, s, prec);

    fmpz_mul_2exp(A, A, E - B);
    fmpz_add(A, A, D);

    transfer_ctx_clear(c);
    arb_clear(x);
    arb_clear(dl);
    arb_clear(a1);
    arb_clear(s);
    arb_clear(sp);
    mag_clear(r);
    mag_clear(tau);
    mag_clear(T0);
    mag_clear(T1);

    return ok;
}

/* The stages can at most double the accuracy each */
#define MAX_STAGES (2 * FLINT_BITS)

int
_arb_hypgeom_legendre_p_ui_root_bitburst_param(arb_t res, arb_t weight, ulong n, ulong k, slong prec,
    slong limb_max_prec, double transfer_blocks_min, slong fft_bits)
{
    fmpz_t A, D;
    slong B, E, wp, p0, K, ex, stage, j;
    double acc, lgn, xd, q, lr;
    ulong lambda;
    arb_t X0, y, yp, pn, pn1, dl, t, x, Dl, N, Sd, Spd;
    arb_struct c[3];
    mag_t r, tau, T0, T1, eps;
    int ok = 0, limbs;

    if (n < 2 || k >= n / 2 || n > ARB_HYPGEOM_LEGENDRE_ROOT_MAX_N)
        return 0;

    lambda = n * (n + 1);
    lgn = log2((double) n);
    wp = prec + 64 + 2 * FLINT_BIT_COUNT(n);
    p0 = 64 + 2 * FLINT_BIT_COUNT(n);

    fmpz_init(A);
    fmpz_init(D);
    arb_init(X0);
    arb_init(y);
    arb_init(yp);
    arb_init(pn);
    arb_init(pn1);
    arb_init(dl);
    arb_init(t);
    arb_init(x);
    arb_init(Dl);
    arb_init(N);
    arb_init(Sd);
    arb_init(Spd);
    for (j = 0; j < 3; j++)
        arb_init(c + j);
    mag_init(r);
    mag_init(tau);
    mag_init(T0);
    mag_init(T1);
    mag_init(eps);

    /* an enclosure of x_k containing no other root */
    if (!_arb_hypgeom_legendre_p_ui_root_isolate(X0, n, k, p0))
        goto cleanup;

    /* x_0 = A / 2^B.  The first transfer needs |x_0 - x_k| well below
       the disk radius r_0 of the Taylor bounds, about 1/n in the
       interior and (k+1)^2 / n^2 near 1: 2^-B must be below about r_0 / 64,
       and the double-precision approximation (about 2^-50) suffices for
       r_0 above about 2^-44, else the midpoint of Petras' enclosure
       (accurate to a small fraction of r_0) is used.  The limb
       recurrence needs a word-size numerator, B = FLINT_BITS - 2 -
       bits(2n) with 2B >= FLINT_BITS. */
    xd = root_d(n, k);
    if (!(xd > 0.0 && xd < 1.0))
        goto cleanup;
    lr = -log2(0.99 * FLINT_MIN(sqrt(1.0 - xd * xd) / n, (1.0 - xd) / 4.0));
    if (!(lr > 0.0 && lr < 1e9))
        goto cleanup;

    B = FLINT_BITS - 2 - FLINT_BIT_COUNT(2 * n);
    limbs = (wp < limb_max_prec && 2 * B >= FLINT_BITS && B >= lr + 6 && lr <= 44);
    if (!limbs)
        B = FLINT_MAX(62 - FLINT_BIT_COUNT(2 * n), (slong) lr + 8);

    if (lr <= 44)
    {
        fmpz_set_d(A, ldexp(xd, B));
    }
    else
    {
        arf_mul_2exp_si(arb_midref(t), arb_midref(X0), B);
        arf_get_fmpz(A, arb_midref(t), ARF_RND_NEAR);
    }

    if (fmpz_sgn(A) <= 0 || fmpz_bits(A) > (ulong) B)
        goto cleanup;

    if (limbs)
        stage0_limbs(pn, pn1, n, fmpz_get_ui(A), B, wp);
    else
        stage0_blocks(pn, pn1, n, A, B, wp, fft_bits);
    deriv_from_pair(y, yp, pn, pn1, n, A, B, wp);

    for (stage = 0; ; stage++)
    {
        if (stage >= MAX_STAGES)
            goto cleanup;

        arb_div(dl, y, yp, wp);
        arb_neg(dl, dl);
        if (!arb_is_finite(dl))
            goto cleanup;
        if (arf_is_zero(arb_midref(dl)))
            break;

        acc = -(double) arf_abs_bound_lt_2exp_si(arb_midref(dl));
        if (2 * acc - lgn - 16 >= wp)
            break;

        /* delta = D / 2^E, about 2 acc bits; the point keeps E bits */
        E = (slong) (2 * acc - lgn - 6);
        E = FLINT_MAX(E, B + 16);
        E = FLINT_MIN(E, wp);
        if (E < B)
            goto cleanup;
        arf_mul_2exp_si(arb_midref(t), arb_midref(dl), E);
        arf_get_fmpz(D, arb_midref(t), ARF_RND_NEAR);
        if (fmpz_is_zero(D))
            break;

        /* the terms decrease by about 2^-q: |delta| is about 2^-acc and
           the radius of convergence is about sin(theta) / (2n) */
        xd = fmpz_get_d_2exp(&ex, A);
        xd = ldexp(xd, ex - B);
        q = acc - log2(2.0 * n / sqrt(FLINT_MAX(1.0 - xd * xd, 1e-300))) - 1.0;
        q = FLINT_MAX(q, 1.0);
        K = (slong) (wp / q) + 4;
        K = FLINT_MIN(K, (slong) n + 2);

        if (!transfer(y, yp, A, B, D, E, n, K, (double) K * wp > transfer_blocks_min, fft_bits, wp))
            goto cleanup;

        B = E;
    }

    /* the Taylor polynomial of degree 2 at x = A/2^B:
       c_0 = y, c_1 = y', c_2 = (2 x y' - lambda y) / (2 (1 - x^2)) */
    arb_set_fmpz(x, A);
    arb_mul_2exp_si(x, x, -B);
    arb_set(c + 0, y);
    arb_set(c + 1, yp);
    arb_mul(t, x, x, wp);
    arb_sub_ui(t, t, 1, wp);
    arb_neg(t, t);
    arb_mul(c + 2, x, yp, wp);
    arb_mul_2exp_si(c + 2, c + 2, 1);
    arb_submul_ui(c + 2, y, lambda, wp);
    arb_div(c + 2, c + 2, t, wp);
    arb_mul_2exp_si(c + 2, c + 2, -1);

    /* d = -(c_0 + c_2 d^2) / c_1 from d = -c_0 / c_1 (midpoints) */
    arb_div(dl, c + 0, c + 1, wp);
    arb_neg(dl, dl);
    mag_zero(arb_radref(dl));
    arb_mul(t, dl, dl, wp);
    arb_mul(t, t, c + 2, wp);
    arb_add(t, t, c + 0, wp);
    arb_div(t, t, c + 1, wp);
    arb_neg(t, t);
    arf_set(arb_midref(dl), arb_midref(t));

    /* interval Newton step N = d - S(d) / S'(Dl) within Dl = d +- eps,
       S the Taylor polynomial plus the tail bound for |t| <= |Dl| */
    mag_set_ui_2exp_si(eps, 1, -(prec + 24));
    arb_set(Dl, dl);
    arb_add_error_mag(Dl, eps);
    arb_get_mag(tau, Dl);
    _arb_hypgeom_legendre_p_ui_disk_radius(r, x, n);
    if (!_arb_hypgeom_legendre_p_ui_taylor_tail(T0, T1, n, x, r, tau, 2))
        goto cleanup;

    /* S(d) */
    arb_mul(Sd, c + 2, dl, wp);
    arb_add(Sd, Sd, c + 1, wp);
    arb_mul(Sd, Sd, dl, wp);
    arb_add(Sd, Sd, c + 0, wp);
    arb_add_error_mag(Sd, T0);
    /* S'(Dl) */
    arb_mul(Spd, c + 2, Dl, wp);
    arb_mul_2exp_si(Spd, Spd, 1);
    arb_add(Spd, Spd, c + 1, wp);
    arb_add_error_mag(Spd, T1);

    arb_div(N, Sd, Spd, wp);
    arb_sub(N, dl, N, wp);
    if (!arb_contains(Dl, N))
        goto cleanup;

    arb_add(res, x, N, wp);
    if (!arb_contains(X0, res))
        goto cleanup;

    if (weight != NULL)
    {
        /* P_n'(x_k) = S'(N) +- T1 */
        arb_mul(Spd, c + 2, N, wp);
        arb_mul_2exp_si(Spd, Spd, 1);
        arb_add(Spd, Spd, c + 1, wp);
        arb_add_error_mag(Spd, T1);
        arb_mul(t, res, res, wp);
        arb_sub_ui(t, t, 1, wp);
        arb_neg(t, t);
        arb_mul(Sd, Spd, Spd, wp);
        arb_mul(Sd, Sd, t, wp);
        arb_ui_div(weight, 2, Sd, prec);
    }

    arb_set_round(res, res, prec);
    ok = 1;

cleanup:
    fmpz_clear(A);
    fmpz_clear(D);
    arb_clear(X0);
    arb_clear(y);
    arb_clear(yp);
    arb_clear(pn);
    arb_clear(pn1);
    arb_clear(dl);
    arb_clear(t);
    arb_clear(x);
    arb_clear(Dl);
    arb_clear(N);
    arb_clear(Sd);
    arb_clear(Spd);
    for (j = 0; j < 3; j++)
        arb_clear(c + j);
    mag_clear(r);
    mag_clear(tau);
    mag_clear(T0);
    mag_clear(T1);
    mag_clear(eps);

    return ok;
}

int
_arb_hypgeom_legendre_p_ui_root_bitburst(arb_t res, arb_t weight, ulong n, ulong k, slong prec)
{
    return _arb_hypgeom_legendre_p_ui_root_bitburst_param(res, weight, n, k, prec,
        LIMB_MAX_PREC, TRANSFER_BLOCKS_MIN, MAT_FFT_BITS);
}
