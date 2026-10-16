/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_dft.h"

/* Bound computation for DFT plans executed in nfixed arithmetic
   (see nfixed.c): additions and subtractions are exact, real
   multiplications truncate with error at most 2 ulp, and every value t
   appearing in the computation must satisfy |t| < 1.

   The propagation tracks, per element, a bound B on the complex
   modulus (which also bounds the real and imaginary parts) and a bound
   D on the modulus of the accumulated error, measured in ulp; the
   final D bounds the errors of the real and imaginary parts of the
   output. A running peak records the largest modulus attained by any
   input, intermediate or output value; the transform is safe from
   overflow iff peak < 1. The peak includes a factor sqrt(2) at every
   root multiplication, accounting for the intermediate sum x_r +- x_i
   inside the 2-multiplication rotation classes.

   The rules are (with Dw an ulp bound for the root table entries and
   all updates inflated by a fudge factor covering the rounding of
   these double-precision computations themselves):

       add/sub:   B <- B_a + B_b,  D <- D_a + D_b        (exact)
       rotation:  B <- B,          D <- D + Dw + K

   the constant K bounding the modulus of the rounding of one complex
   multiplication (sqrt(2) times the per-component bound reported by
   _gr_dft_nfixed_cmul_err_ulps for the actual context: 2, 4 or 6 ulp
   per component depending on the limb count and whether the Karatsuba
   product is used), and Dw * B <= Dw absorbing the error of the root
   itself.

   All quantities are affine in the input (B, D), so the recursive
   algorithms are propagated with memoized affine coefficients, giving
   O(log n) cost. */

#define NFB_FUDGE (1.0 + 1e-10)
#define NFB_SQRT2 1.41421356237309515

typedef struct
{
    double B;      /* modulus bound */
    double D;      /* modulus error bound, ulp */
    double peak;   /* running maximum modulus */
    double K;      /* modulus rounding of one complex multiplication */
}
nfb_state;

static void
nfb_peak(nfb_state * st, double val)
{
    st->peak = FLINT_MAX(st->peak, val * NFB_FUDGE);
}

static double
nfb_dw(const gr_dft_pre_t P)
{
    return FLINT_MAX(P->nfixed_root_err, 2.0);
}

/* multiplication by a root of unity */
static void
nfb_rot(nfb_state * st, double Dw)
{
    nfb_peak(st, NFB_SQRT2 * st->B);
    st->D = (st->D + Dw + st->K) * NFB_FUDGE;
}

/* radix-2 Cooley-Tukey, depth levels */
static void
nfb_ct(nfb_state * st, int depth, double Dw)
{
    int lev;

    for (lev = 0; lev < depth; lev++)
    {
        st->B = 2.0 * st->B * NFB_FUDGE;
        st->D = 2.0 * st->D * NFB_FUDGE;
        nfb_peak(st, st->B);
        nfb_rot(st, Dw);
    }
}

/* split radix: affine coefficients per depth, out = (Bm * B0,
   Da * D0 + Db, peak = Pk * B0), memoized */
static void
nfb_split_coeffs(int depth, double Dw, double K,
        double * Bm, double * Da, double * Db, double * Pk)
{
    double bm[FLINT_BITS + 1], da[FLINT_BITS + 1];
    double db[FLINT_BITS + 1], pk[FLINT_BITS + 1];
    int d;

    for (d = 0; d <= depth; d++)
    {
        if (d <= 2)
        {
            /* small basecases: model as plain Cooley-Tukey; extract the
               affine form D_out = Da D_0 + Db by evaluating at
               D_0 = 0 and D_0 = 1 (the propagation is affine) */
            nfb_state st0 = { 1.0, 0.0, 1.0, K };
            nfb_state st1 = { 1.0, 1.0, 1.0, K };
            nfb_ct(&st0, d, Dw);
            nfb_ct(&st1, d, Dw);
            bm[d] = st1.B;
            db[d] = st0.D;
            da[d] = st1.D - st0.D;
            pk[d] = st1.peak;
        }
        else
        {
            /* combine: u1, u2 from the half transform; zk, zpk from
               the quarter transforms; t1 = w^k zk, t2 = w^(3k) zpk;
               t3 = t1 + t2; t1 = t1 - t2; t2 = -i t1 (free);
               outputs u1 +- t3, u2 +- t2 */
            double Bh = bm[d - 1], Bq = bm[d - 2];
            double rotD = da[d - 2];      /* D of t1, t2: quarter + rot */
            double rotDb = db[d - 2] + (Dw + K) * NFB_FUDGE;

            bm[d] = (Bh + 2.0 * Bq) * NFB_FUDGE;
            da[d] = (da[d - 1] + 2.0 * rotD) * NFB_FUDGE;
            db[d] = (db[d - 1] + 2.0 * rotDb) * NFB_FUDGE;

            pk[d] = FLINT_MAX(pk[d - 1], pk[d - 2]);
            pk[d] = FLINT_MAX(pk[d], NFB_SQRT2 * Bq);   /* inside rotations */
            pk[d] = FLINT_MAX(pk[d], 2.0 * Bq);         /* t3, t1 */
            pk[d] = FLINT_MAX(pk[d], bm[d]);            /* outputs */
            pk[d] *= NFB_FUDGE;
        }
    }

    *Bm = bm[depth];
    *Da = da[depth];
    *Db = db[depth];
    *Pk = pk[depth];
}

static void nfb_plan(nfb_state * st, const gr_dft_pre_t P);

/* the direct prime kernel of the mixed-radix algorithm, radix p:
   inputs z_0, ..., z_(p-1) (already twiddled), outputs
   Y_0 = sum z_b and Y_r = (z_0 - z_b*) + sum w^e (z_b - z_b*) */
static void
nfb_prime_kernel(nfb_state * st, ulong p, double Dw)
{
    double B = st->B, D = st->D;
    double Bd = 2.0 * B, Dd = 2.0 * D;
    nfb_state rot = { Bd, Dd, st->peak, st->K };
    double By, Dy;

    /* Y_0 */
    By = p * B;
    Dy = p * D;

    /* rotated differences */
    nfb_rot(&rot, Dw);
    st->peak = rot.peak;
    nfb_peak(st, Bd);

    /* Y_r: accumulation of one difference plus (p - 2) rotated
       differences; the running sums are bounded by the final value */
    By = FLINT_MAX(By, (p - 1) * Bd);
    Dy = FLINT_MAX(Dy, Dd + (p - 2) * rot.D);

    st->B = By * NFB_FUDGE;
    st->D = Dy * NFB_FUDGE;
    nfb_peak(st, st->B);
}

static void
nfb_mixed(nfb_state * st, const gr_dft_pre_t P)
{
    double Dw = nfb_dw(P);
    slong lvl;

    /* propagate from the innermost level outward */
    for (lvl = P->num_radices - 1; lvl >= 0; lvl--)
    {
        ulong p = P->radices[lvl];

        /* twiddle */
        nfb_rot(st, Dw);

        if (P->P1 != NULL && P->P1->n == p)
        {
            nfb_plan(st, P->P1);
        }
        else
        {
            nfb_prime_kernel(st, p, Dw);
        }
    }
}

static void
nfb_naive(nfb_state * st, const gr_dft_pre_t P)
{
    double Dw = nfb_dw(P);
    ulong n = P->n;
    nfb_state term = *st;

    nfb_rot(&term, Dw);
    st->peak = term.peak;

    st->B = n * st->B * NFB_FUDGE;
    st->D = n * term.D * NFB_FUDGE;
    nfb_peak(st, st->B);
}

/* Bluestein over fixed-point arithmetic (modeling the shifted-kernel
   variant used by the fixed-point contexts, where the transformed
   kernel carries the folded scaling 1/(2 conv_len) and the outputs
   are doubled at the end).

   Writing z for the chirped input (modulus at most B_z per entry,
   n nonzero entries), N = conv_len and h for the unit-modulus chirp
   kernel scaled by 1/(2N):

   - the forward transform of z follows the generic model (values up
     to N B_z, rotation intermediates up to sqrt(2) N B_z);

   - the transformed kernel satisfies |h^| <= (2n-1)/(2N) < 1/2
     (2n - 1 nonzero unit entries scaled by 1/(2N); note 2n - 1 < N
     strictly since N is a power of two), so the pointwise products
     are bounded by N B_z (2n-1)/(2N);

   - the values of the unscaled inverse transform -- including all its
     intermediates -- are bounded by N ||z||_1 max|h| <= n B_z / 2:
     each stage-s intermediate is an unscaled length-2^s inverse of a
     subsampled spectrum, which by the aliasing identity equals
     2^s times an alias sum of the true convolution u = z (*) h, and
     since |u_k| <= ||z||_1 max|h| uniformly, any alias sum over
     N / 2^s terms is at most (N / 2^s) ||z||_1 max|h|. This replaces
     the generic stage-doubling bound (which would give N n B_z) and
     is what makes Bluestein usable in fixed point;

   - the final doubling and output chirp then give |X| <= n B_z as for
     any DFT.

   The error propagation uses the generic child model throughout
   (error growth does follow the stage doubling); the kernel error is
   modeled from the chirp table error scaled by 1/(2N) and pushed
   through the child transform. */
static void
nfb_bluestein(nfb_state * st, const gr_dft_pre_t P)
{
    double Dw = nfb_dw(P);
    double n = (double) P->n;
    double N = (double) P->conv_len;
    double Hmax = (2.0 * n - 1.0) / (2.0 * N);
    double Bz, Dh, Bf;
    nfb_state ch;

    /* kernel error: chirp entries (table error Dw), scaled by
       1/(2N) (errors scale down; truncation below 1 ulp per part),
       then the forward child transform */
    ch.B = (1.0 / (2.0 * N)) * NFB_FUDGE;
    ch.D = (Dw / (2.0 * N) + NFB_SQRT2) * NFB_FUDGE;
    ch.peak = ch.B;
    ch.K = st->K;
    nfb_plan(&ch, P->P1);
    Dh = ch.D;

    /* input chirp */
    nfb_rot(st, Dw);
    Bz = st->B;

    /* forward convolution transform */
    nfb_plan(st, P->P1);
    Bf = st->B;

    /* pointwise multiplication by the transformed kernel */
    st->D = (Hmax * st->D + Bf * Dh + st->K) * NFB_FUDGE;
    st->B = (Bf * Hmax) * NFB_FUDGE;
    nfb_peak(st, st->B);

    /* unscaled inverse convolution transform: generic model for the
       errors, aliasing bound for the values */
    {
        nfb_state di = *st;
        nfb_plan(&di, P->P1);
        st->D = di.D;
    }
    st->B = (n * Bz / 2.0) * NFB_FUDGE;
    nfb_peak(st, NFB_SQRT2 * st->B);

    /* doubling (exact) */
    st->B = 2.0 * st->B * NFB_FUDGE;
    st->D = 2.0 * st->D * NFB_FUDGE;
    nfb_peak(st, st->B);

    /* output chirp */
    nfb_rot(st, Dw);
}

static void
nfb_plan(nfb_state * st, const gr_dft_pre_t P)
{
    double Dw = nfb_dw(P);

    nfb_peak(st, st->B);

    switch (P->alg)
    {
        case GR_DFT_ALG_CT:
            nfb_ct(st, P->depth, Dw);
            break;

        case GR_DFT_ALG_SPLIT:
            {
                double Bm, Da, Db, Pk;
                nfb_split_coeffs(P->depth, Dw, st->K, &Bm, &Da, &Db, &Pk);
                nfb_peak(st, Pk * st->B);
                st->D = (Da * st->D + Db) * NFB_FUDGE;
                st->B = Bm * st->B * NFB_FUDGE;
            }
            break;

        case GR_DFT_ALG_BAILEY:
            /* column transforms, twiddles, row transforms; the final
               transposition is exact */
            nfb_plan(st, P->P2);
            nfb_rot(st, Dw);
            nfb_plan(st, P->P1);
            break;

        case GR_DFT_ALG_MIXED:
            nfb_mixed(st, P);
            break;

        case GR_DFT_ALG_PFA:
            /* no twiddles; gather/scatter are exact */
            nfb_plan(st, P->P2);
            nfb_plan(st, P->P1);
            break;

        case GR_DFT_ALG_BLUESTEIN:
            nfb_bluestein(st, P);
            break;

        default:
            nfb_naive(st, P);
            break;
    }

    nfb_peak(st, st->B);
}

/* Computes a bound on the modulus of all input, intermediate and
   output values (written to peak; the computation is free of overflow
   in nfixed arithmetic iff this is < 1), and a bound on the errors of
   the real and imaginary parts of the output measured in ulp (written
   to err_ulps), for a forward or unscaled inverse transform with the
   plan P executed in nfixed arithmetic, given that the input values
   have complex modulus at most in_mag and componentwise errors of at
   most in_err_ulps ulp. */
void
gr_dft_precomp_nfixed_bound(double * peak, double * err_ulps,
        double in_mag, double in_err_ulps, const gr_dft_pre_t P)
{
    nfb_state st;

    st.B = in_mag;
    st.D = NFB_SQRT2 * in_err_ulps;   /* componentwise -> modulus */
    st.peak = in_mag;

    /* rounding of one complex multiplication, as a modulus: the exact
       per-component error of the plan's context when it is a
       fixed-point context, else the worst case of the model */
    if (P->ctx != NULL && _gr_dft_ctx_is_nfixed_complex((gr_ctx_struct *) P->ctx))
        st.K = NFB_SQRT2 * _gr_dft_nfixed_cmul_err_ulps((gr_ctx_struct *) P->ctx);
    else
        st.K = NFB_SQRT2 * 6.0;

    nfb_plan(&st, P);

    *peak = st.peak;
    *err_ulps = st.D;
}
