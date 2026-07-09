/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_DFT_H
#define GR_DFT_H

#include "thread_pool.h"
#include "gr.h"
#include "arf_types.h"
#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Algorithm selection ******************************************************/

#define GR_DFT_ALG_AUTO       0
#define GR_DFT_ALG_NAIVE      1  /* O(n^2) evaluation (any length) */
#define GR_DFT_ALG_CT         2  /* iterative radix-2 Cooley-Tukey (power of two) */
#define GR_DFT_ALG_BAILEY     3  /* cache-friendly four-step (power of two) */
#define GR_DFT_ALG_SPLIT      4  /* recursive split-radix (power of two) */
#define GR_DFT_ALG_PFA        5  /* Good-Thomas prime factor algorithm over
                                    coprime components (composite lengths) */
#define GR_DFT_ALG_MIXED      6  /* recursive mixed-radix Cooley-Tukey with
                                    prime basecases (any length) */
#define GR_DFT_ALG_BLUESTEIN  7  /* chirp-z reduction to a power-of-two
                                    convolution (odd lengths) */

/* Flags ********************************************************************/

/* Forward transform leaves the output in scrambled order (bit-reversed
   for NAIVE/CT, transposed for BAILEY); the inverse transform expects
   its input in the same scrambled order. Only supported for power-of-two
   lengths with the NAIVE, CT and BAILEY algorithms; ignored otherwise. */
#define GR_DFT_SCRAMBLED   1

/* Default threading granularity (elements of work per serial block);
   see gr_dft_precomp_set_serial_block. */
/* default minimum number of elements of per-chunk work for the
   threaded transforms (the number of chunks per phase or pass is at
   most n / serial_block, and threaded recursion stops below this
   size); deliberately small so that every algorithm uses as many
   threads as it structurally can -- users of cheap rings should raise
   it with gr_dft_precomp_set_serial_block */
#define GR_DFT_SERIAL_BLOCK_DEFAULT 256

/* blocks larger than this (in bytes) are processed depth-first by the
   Cooley-Tukey transform, so that cache-sized sub-blocks complete all
   their butterfly levels while resident instead of being streamed
   from memory once per level (machine-tunable; should be at most
   about half the per-core cache) */
#define GR_DFT_CT_BLOCK_BYTES (UWORD(128) << 10)

/* Precomputed plan *********************************************************/

/* Cost classes of the roots in complex mode (internal). Since complex
   mode always uses the standard root w = exp(-2 pi i / n), the classes
   are determined structurally by the exponent: the quarter points are
   free rotations and the odd multiples of n/8 are of the form
   c (1 -+ i), requiring 2 real multiplications. */
#define GR_DFT_ROOT_GENERIC   0  /* 3 real multiplications */
#define GR_DFT_ROOT_DMC_ZERO  1  /* d = c:  w = c (1 + i), 2 multiplications */
#define GR_DFT_ROOT_DPC_ZERO  2  /* d = -c: w = c (1 - i), 2 multiplications */
#define GR_DFT_ROOT_NEG_ONE   3  /* w = -1, free */
#define GR_DFT_ROOT_I         4  /* w = i, free */
#define GR_DFT_ROOT_NEG_I     5  /* w = -i, free */

typedef struct gr_dft_pre_struct
{
    ulong n;                    /* transform length */
    int depth;                  /* n == 2^depth if n is a power of two, else -1 */
    int alg;                    /* resolved algorithm (never AUTO) */
    int flags;
    gr_ctx_struct * ctx;        /* the ring R (or C = R'[i] in complex mode) */
    gr_ctx_struct * real_ctx;   /* complex mode: the real ring R'; else NULL */
    gr_ptr roots;               /* w^0, w^1, ..., w^(n-1), elements of ctx;
                                   NULL for canonical PFA plans (which
                                   perform no root multiplications) and
                                   for plans with packed stage tables */
    gr_ptr stage_tab;           /* CT and split-radix plans (without the
                                   complex Karatsuba tables) store their
                                   twiddles packed per stage instead of
                                   the serial table, so that every table
                                   access during a transform is a
                                   sequential walk; same total size.
                                   CT: stage s at offset n - 2^s with
                                   entries w^(j n / 2^s), j < 2^(s-1).
                                   Split: size-m section at offset n - m
                                   with interleaved pairs
                                   (w^(k n/m), w^(3 k n/m mod n)),
                                   k < m/4. */
    slong stage_len;
    gr_ptr wtab;                /* complex Karatsuba tables (when enabled): 3n
                                   elements of real_ctx, (c, d - c, d + c) for
                                   each root c + d*i; else NULL */
    unsigned char * wclass;     /* complex mode: cost class of each root */
    ulong n1, n2;               /* BAILEY, PFA: n = n1 * n2 */
    ulong pfa_a, pfa_b;         /* PFA: CRT output coefficients */
    ulong * radices;            /* MIXED: prime radix of each recursion level */
    slong num_radices;
    ulong conv_len;             /* BLUESTEIN: power-of-two convolution length */
    gr_ptr bl_kern;             /* BLUESTEIN: transformed convolution kernel,
                                   conv_len elements, scaled by 1/conv_len */
    gr_ptr bl_wtab;             /* BLUESTEIN, with Karatsuba tables enabled:
                                   multiplication tables for bl_kern
                                   (3 conv_len elements) */
    struct gr_dft_pre_struct * P1;  /* BAILEY, PFA: sub-plan of length n1;
                                       MIXED: optional plan for the prime
                                       radix; BLUESTEIN: power-of-two plan
                                       of length conv_len */
    struct gr_dft_pre_struct * P2;  /* BAILEY, PFA: sub-plan of length n2 */
    thread_pool_handle * threads;   /* attached worker threads (borrowed;
                                       NULL: request from the global pool
                                       during transforms) */
    slong num_threads;              /* number of attached workers */
    slong serial_block;             /* threading granularity (0: default) */
    double nfixed_root_err;         /* nfixed contexts: ulp error bound of
                                       the root table (see nfixed.c) */
    int bl_shifted;                 /* BLUESTEIN over nfixed: the kernel
                                       carries an extra 1/2, undone by
                                       doubling the outputs */
}
gr_dft_pre_struct;

typedef gr_dft_pre_struct gr_dft_pre_t[1];

/* Plan initialization ******************************************************/

WARN_UNUSED_RESULT int gr_dft_default_root(gr_ptr w, ulong n, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_dft_precomp_init_root(gr_dft_pre_t P, gr_srcptr w, ulong n, int alg, int flags, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_dft_precomp_init(gr_dft_pre_t P, ulong n, int alg, int flags, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_dft_precomp_init_karatsuba(gr_dft_pre_t P, ulong n, int alg, int flags, gr_ctx_t real_ctx, gr_ctx_t ctx);

void gr_dft_precomp_set_threads(gr_dft_pre_t P, thread_pool_handle * threads, slong num_threads);
void gr_dft_precomp_set_serial_block(gr_dft_pre_t P, slong serial_block);

void gr_dft_precomp_clear(gr_dft_pre_t P);

void gr_dft_precomp_output_perm(ulong * perm, const gr_dft_pre_t P);

/* Transforms ***************************************************************/

WARN_UNUSED_RESULT int gr_dft_precomp(gr_ptr res, gr_srcptr vec, const gr_dft_pre_t P, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_dft_inverse_precomp(gr_ptr res, gr_srcptr vec, const gr_dft_pre_t P, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_dft(gr_ptr res, gr_srcptr vec, ulong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_dft_inverse(gr_ptr res, gr_srcptr vec, ulong n, gr_ctx_t ctx);

/* Fixed-point contexts for numerical DFTs (see nfixed.c) *******************/

int gr_dft_ctx_init_nfixed(gr_ctx_t ctx, slong nlimbs);
int gr_dft_ctx_init_nfixed_complex(gr_ctx_t ctx, slong nlimbs);
int _gr_dft_ctx_is_nfixed(gr_ctx_t ctx);
int _gr_dft_ctx_is_nfixed_complex(gr_ctx_t ctx);

int gr_dft_nfixed_set_arf(gr_ptr res, const arf_t x, gr_ctx_t ctx);
void gr_dft_nfixed_get_arf(arf_t res, gr_srcptr x, gr_ctx_t ctx);

int _gr_dft_nfixed_roots(gr_ptr roots, ulong n, double * err_ulps, gr_ctx_t ctx);
double _gr_dft_nfixed_cmul_err_ulps(gr_ctx_t ctx);
double _gr_dft_nfixed_root_err_bound(ulong n, double kmul);
/* forced complex-multiplication paths, for machine tuning of
   GR_DFT_NFIXED_KARATSUBA_CUTOFF (see profile/p-gr_dft_nfixed) */
int _gr_dft_nfixed_cmul_schoolbook(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);
int _gr_dft_nfixed_cmul_karatsuba(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);

/* Two-phase plan construction: the layout carries the complete
   decomposition without any ring elements, so that cost and error
   bounds can be inspected (e.g. to choose a fixed-point precision)
   before the root tables are computed by the realize step, which
   uses canonical roots. */
int _gr_dft_precomp_init_layout(gr_dft_pre_t P, ulong n, int alg, int flags,
        int complex_mode);
int _gr_dft_precomp_realize(gr_dft_pre_t P, gr_ctx_struct * real_ctx,
        gr_ctx_t ctx);

void gr_dft_precomp_nfixed_bound(double * peak, double * err_ulps,
        double in_mag, double in_err_ulps, const gr_dft_pre_t P);

#define GR_DFT_NFIXED_MAX_NLIMBS 2048

/* Drop-in replacements for acb_dft / acb_dft_inverse (see acb.c) */

void gr_dft_acb(acb_ptr w, acb_srcptr v, slong n, slong prec);
void gr_dft_acb_inverse(acb_ptr w, acb_srcptr v, slong n, slong prec);
int _gr_dft_acb(acb_ptr w, acb_srcptr v, slong n, int inverse, int which,
        slong prec);

/* Precomputed variant, amortizing the plan construction over repeated
   transforms of the same length and precision. which = 1: ball
   arithmetic (arb/acb contexts); which = 2: fixed-point contexts,
   together with the input-independent constants of the scaling and
   error analysis (the input-dependent power-of-two scale is chosen
   per transform). */
typedef struct
{
    slong n;
    slong prec;
    int which;
    slong nl;               /* fixed point: number of limbs */
    slong p1e;              /* fixed point: peak(1) < 2^p1e */
    double in_mag;          /* fixed point: scaled input bound */
    double errulps;         /* fixed point: output roundoff, in ulp */
    gr_ctx_struct rctx[1];
    gr_ctx_struct cctx[1];
    gr_dft_pre_t P;
}
gr_dft_acb_pre_struct;

typedef gr_dft_acb_pre_struct gr_dft_acb_pre_t[1];

WARN_UNUSED_RESULT int gr_dft_acb_precomp_init(gr_dft_acb_pre_t Q, slong n, slong prec);
void gr_dft_acb_precomp_clear(gr_dft_acb_pre_t Q);
void gr_dft_acb_precomp(acb_ptr w, acb_srcptr v, const gr_dft_acb_pre_t Q, slong prec);
void gr_dft_acb_inverse_precomp(acb_ptr w, acb_srcptr v, const gr_dft_acb_pre_t Q, slong prec);
int _gr_dft_acb_precomp(acb_ptr w, acb_srcptr v, int inverse, const gr_dft_acb_pre_t Q, slong prec);

/* Internal functions *******************************************************/

WARN_UNUSED_RESULT int _gr_dft_mul_root(gr_ptr res, gr_srcptr x, ulong e, int inverse, gr_ptr rtmp, const gr_dft_pre_t P);
WARN_UNUSED_RESULT int _gr_dft_mul_const(gr_ptr res, gr_srcptr x, gr_srcptr w, gr_srcptr wtab3, gr_ptr rtmp, const gr_dft_pre_t P);
void _gr_dft_bit_reverse(gr_ptr x, slong stride, int depth, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_dft_ct(gr_ptr x, slong stride, int inverse, int scrambled, const gr_dft_pre_t P, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_dft_ct_threaded(gr_ptr x, slong stride, int inverse, int scrambled, const gr_dft_pre_t P, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_dft_split(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_dft_split_threaded(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx);
#ifdef FLINT_HAVE_FILE
void gr_dft_precomp_fprint(FILE * out, const gr_dft_pre_t P);
#endif
void gr_dft_precomp_print(const gr_dft_pre_t P);
WARN_UNUSED_RESULT int _gr_dft_bailey(gr_ptr x, int inverse, int scrambled, const gr_dft_pre_t P, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_dft_naive(gr_ptr res, gr_srcptr vec, int inverse, int scrambled, const gr_dft_pre_t P, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_dft_mixed(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_dft_pfa(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_dft_bluestein(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_dft_precomp_raw(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
