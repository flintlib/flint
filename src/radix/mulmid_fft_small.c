/*
    Copyright (C) 2022 Daniel Schultz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "radix.h"

/* allow changing this from the test code */
FLINT_DLL int radix_mulmid_force_threading = 0;

#if FLINT_HAVE_FFT_SMALL

#include "thread_pool.h"
#include "thread_support.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "crt_helpers.h"
#include "fft_small.h"

static void _mod_red(
    double* abuf, ulong atrunc,
    const ulong* a, ulong an,
    const sd_fft_ctx_struct* fft,
    nmod_t mod)
{
    double* aI;
    ulong i, j;

    FLINT_ASSERT(atrunc < an);
    FLINT_ASSERT(atrunc%BLK_SZ == 0);

    ulong tt = an%atrunc;

#define UNROLL 8

    for (i = 0; i < atrunc; i += BLK_SZ)
    {
        aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);

        vec8n nn = vec8n_set_n(mod.n);
        vec8d n = vec8d_set_d(fft->p);
        vec8d ninv = vec8d_set_d(fft->pinv);

        for (j = 0; j < BLK_SZ; j += UNROLL)
        {
            if (i+j+UNROLL <= tt || i+j >= tt)
            {
                ulong k = i+j;
                FLINT_ASSERT(k+UNROLL-1 < an);
                vec8n t = vec8n_load_unaligned(a + k);

                if (mod.norm == 0)
                    for (k += atrunc; k < an; k += atrunc)
                        t = vec8n_addmod(t, vec8n_load_unaligned(a + k), nn);
                else
                    for (k += atrunc; k < an; k += atrunc)
                        t = vec8n_addmod_limited(t, vec8n_load_unaligned(a + k), nn);


                vec8d tlo = vec8n_convert_limited_vec8d(vec8n_bit_and(t, vec8n_set_n(n_pow2(32)-1)));
                vec8d thi = vec8n_convert_limited_vec8d(vec8n_bit_shift_right_32(t));
                vec8d_store_aligned(aI + j, vec8d_add(tlo, vec8d_mulmod(thi, vec8d_set_d(n_pow2(32)), n, ninv)));
            }
            else
            {
                for (ulong l = 0; l < UNROLL; l++)
                {
                    ulong k = i+j+l;
                    ulong c = a[k];
                    for (k += atrunc; k < an; k += atrunc)
                        c = nmod_add(c, a[k], mod);

                    aI[j+l] = (slong)(nmod_set_ui(c, fft->mod));
                }
            }
        }
    }
}

static void _mod(
    double* abuf, ulong atrunc,
    const ulong* a, ulong an,
    const sd_fft_ctx_struct* fft,
    nmod_t mod)
{
    double* aI;
    ulong i, j;

    if (atrunc < an)
    {
        _mod_red(abuf, atrunc, a, an, fft, mod);
        return;
    }

    if (mod.n <= fft->mod.n)
    {
        for (i = 0; i + BLK_SZ <= an; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
            for (j = 0; j < BLK_SZ; j += 8)
            {
                vec8n t = vec8n_load_unaligned(a + i + j);
                FLINT_ASSERT(i+j < atrunc);
                vec8d_store_aligned(aI + j, vec8n_convert_limited_vec8d(t));
            }
        }

        aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
        for (j = 0; j < an - i; j++)
            aI[j] = (slong)a[i + j];
    }
    else
    {
        vec8d n = vec8d_set_d(fft->p);
        vec8d ninv = vec8d_set_d(fft->pinv);
        for (i = 0; i + BLK_SZ <= an; i += BLK_SZ)
        {
            aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
#if 1
            for (j = 0; j < BLK_SZ; j += 8)
            {
                vec8n t = vec8n_load_unaligned(a + i + j);
                vec8d tlo = vec8n_convert_limited_vec8d(vec8n_bit_and(t, vec8n_set_n(n_pow2(32)-1)));
                vec8d thi = vec8n_convert_limited_vec8d(vec8n_bit_shift_right_32(t));
                vec8d_store_aligned(aI + j, vec8d_add(tlo, vec8d_mulmod(thi, vec8d_set_d(n_pow2(32)), n, ninv)));
            }
#else
            for (j = 0; j < BLK_SZ; j += 1)
            {
                aI[j] = (slong)nmod_set_ui(a[i+j], fft->mod);
            }
#endif
        }

        aI = sd_fft_ctx_blk_index(abuf, i/BLK_SZ);
        for (j = 0; j < an - i; j++)
            aI[j] = (slong)nmod_set_ui(a[i+j], fft->mod);
    }

    for (i = an; i < atrunc; i++)
        abuf[i] = 0;
}

typedef struct {
    nmod_t mod;
    int write_carry_out;
} _conv_params_struct;

/* the reconstructed coefficients are digits in a chain: each range keeps
   a running carry that is divided out coefficient by coefficient, with
   the final carry exported through `local` for the serial fixup pass */
#define CRT_FN(NP) CAT3(_crt, NP, fn)
#define CRT_Z_TYPE ulong*
#define CRT_HEAD \
    nmod_t mod = ((const _conv_params_struct*) params)->mod; \
    ulong cy[3] = { 0, 0, 0 };
#define CRT_EMIT(zi, r, NP, N, M) \
    if (N == 2) \
    { \
        add_ssaaaa(r[1], r[0], r[1], r[0], cy[1], cy[0]); \
        if (mod.norm == 0) \
            z[zi] = flint_mpn_divrem_2_1_preinv_norm(cy, r, mod.n, mod.ninv); \
        else \
            z[zi] = flint_mpn_divrem_2_1_preinv_unnorm(cy, r, mod.n, mod.ninv, mod.norm); \
    } \
    else \
    { \
        FLINT_ASSERT(N < 4 || r[3] == 0); \
        add_sssaaaaaa(r[2], r[1], r[0], r[2], r[1], r[0], cy[2], cy[1], cy[0]); \
        if (mod.norm == 0) \
            z[zi] = flint_mpn_divrem_3_1_preinv_norm(cy, r, mod.n, mod.ninv); \
        else \
            z[zi] = flint_mpn_divrem_3_1_preinv_unnorm(cy, r, mod.n, mod.ninv, mod.norm); \
    }
#define CRT_TAIL \
    local[0] = cy[0]; \
    local[1] = cy[1]; \
    local[2] = cy[2];
CRT_DEFINE(2, 2, 1)
CRT_DEFINE(3, 3, 2)
CRT_DEFINE(4, 4, 3)
#undef CRT_FN
#undef CRT_Z_TYPE
#undef CRT_HEAD
#undef CRT_EMIT
#undef CRT_TAIL
#undef CRT_DEFINE

static void _crt_1(
    ulong* z, ulong zl, ulong zi_start, ulong zi_stop,
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride,
    crt_data_struct* FLINT_UNUSED(Rcrts),
    ulong * carry_out,
    nmod_t mod)
{
    ulong i, j, jstart, jstop;
    ulong Xs[BLK_SZ*1];
    ulong cy = 0;

    for (i = n_round_down(zi_start, BLK_SZ); i < zi_stop; i += BLK_SZ)
    {
        _convert_block(Xs, Rffts, d, dstride, 1, i/BLK_SZ);

        jstart = (i < zi_start) ? zi_start - i : 0;
        jstop = FLINT_MIN(BLK_SZ, zi_stop - i);

        FLINT_ASSERT(mod.norm != 0);

        for (j = jstart; j < jstop; j += 1)
        {
            z[i+j-zl] = n_divrem_preinv_unnorm(&cy, Xs[j] + cy, mod.n, mod.ninv, mod.norm);
        }
    }

    carry_out[0] = cy;
    carry_out[1] = 0;
    carry_out[2] = 0;
}


/* ------------------------------------------------------------------------ */
/* glue for the generic convolution engine                                  */
/* ------------------------------------------------------------------------ */

static void _mod_fn(double* xbuf, ulong xtrunc,
    const void* x, ulong xn, slong FLINT_UNUSED(xaux),
    const sd_fft_ctx_struct* fft, const void* params)
{
    const _conv_params_struct* par = (const _conv_params_struct*) params;

    _mod(xbuf, xtrunc, (const ulong*) x, xn, fft, par->mod);
}

#define DEFINE_IT(NP) \
static void CAT3(_crt, NP, fn)(void* z, ulong zl, \
    ulong zi_start, ulong zi_stop, \
    sd_fft_ctx_struct* Rffts, double* d, ulong dstride, \
    crt_data_struct* Rcrts, ulong FLINT_UNUSED(min_an_bn), \
    ulong* local, const void* params) \
{ \
    const _conv_params_struct* par = (const _conv_params_struct*) params; \
 \
    CAT(_crt, NP)((ulong*) z, zl, zi_start, zi_stop, Rffts, d, dstride, \
                  Rcrts, local, par->mod); \
}
DEFINE_IT(1)
#undef DEFINE_IT

static fft_small_crt_func _radix_crt_fn(ulong np)
{
    return np == 1 ? _crt_1_fn :
           np == 2 ? _crt_2_fn :
           np == 3 ? _crt_3_fn : _crt_4_fn;
}

static ulong _radix_conv_s1_threads(ulong np, ulong tune_n)
{
    return (tune_n < 1500 || radix_mulmid_force_threading) ? 1 : np;
}

static int _radix_conv_s1_b_worker(ulong FLINT_UNUSED(np), ulong tune_n,
                                   int squaring)
{
    return tune_n > 5000 && !squaring;
}

static int _radix_conv_s2_rethread(ulong np, ulong zn)
{
    return np*zn > 10000 || radix_mulmid_force_threading;
}

static const fft_small_conv_tuning _radix_conv_tuning =
    { _radix_conv_s1_threads, _radix_conv_s1_b_worker,
      _radix_conv_s2_rethread };

/* Serial pass over the stage 2 ranges: reduce each range's carries into
   the following coefficients and write the top carry of the full product
   if requested. */
static void _radix_s2_finish(void* zv, ulong zl, ulong FLINT_UNUSED(zh),
    fft_small_crt_range_struct* ranges, ulong nranges, const void* params)
{
    ulong* z = (ulong*) zv;
    const _conv_params_struct* par = (const _conv_params_struct*) params;
    nmod_t mod = par->mod;
    ulong i;

    for (i = 0; i < nranges; i++)
    {
        slong stop = ranges[i].stop_zi;

        /* Final segment */
        if (i == nranges - 1)
        {
            /* Write carry to the top limb of the full product. */
            if (par->write_carry_out)
            {
                z[stop - zl] = ranges[i].local[0];

                FLINT_ASSERT(ranges[i].local[1] == 0 && ranges[i].local[2] == 0);
            }
            /* Otherwise, the product was truncated. */
        }
        else
        {
            /* Propagate carries from the previous segment */
            slong stop2 = ranges[i + 1].stop_zi;

            ulong cy[3];
            slong j;

            cy[0] = ranges[i].local[0];
            cy[1] = ranges[i].local[1];
            cy[2] = ranges[i].local[2];

            for (j = stop; j < stop2; j++)
            {
                add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0], 0, 0, z[j - zl]);

                if (mod.norm == 0)
                    z[j - zl] = flint_mpn_divrem_3_1_preinv_norm(cy, cy, mod.n, mod.ninv);
                else
                    z[j - zl] = flint_mpn_divrem_3_1_preinv_unnorm(cy, cy, mod.n, mod.ninv, mod.norm);

                if (cy[0] == 0 && cy[1] == 0 && cy[2] == 0)
                    break;
            }

            /* Propagate carry to next segment */
            add_sssaaaaaa(ranges[i + 1].local[2], ranges[i + 1].local[1], ranges[i + 1].local[0],
                          ranges[i + 1].local[2], ranges[i + 1].local[1], ranges[i + 1].local[0],
                          cy[2], cy[1], cy[0]);
        }
    }
}

/* the window [zl, zh) of the convolution of (a, an) and (b, bn), with
   zh <= an + bn - 1, as a window of the full product (wrapped around when
   a power of two allows) */
static void _radix_mul_mid_mpn_ctx(
    ulong * z, ulong zl, ulong zh,
    const ulong* a, ulong an,
    const ulong* b, ulong bn,
    int write_carry_out,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong modbits = FLINT_BITS - mod.norm;
    ulong zn = an + bn - 1;
    ulong atrunc, btrunc;
    int squaring, success;
    fft_small_plan_t P;
    fft_small_conv_arg_struct A;
    _conv_params_struct par;

    FLINT_ASSERT(an > 0);
    FLINT_ASSERT(bn > 0);
    FLINT_ASSERT(zl < zh);
    FLINT_ASSERT(zh <= zn);

    squaring = (a == b) && (an == bn);

    atrunc = n_round_up(an, BLK_SZ);
    btrunc = n_round_up(bn, BLK_SZ);

    P->R = R;
    P->sign = 0;

    _fft_small_plan_set_window(P, zl, zh, zn, n_max(atrunc, btrunc),
                               (ulong) 4);

    /* need prod_of_primes >= min(an, bn) * 4^modbits; the radix chinese
       remaindering handles at most 4 primes */
    success = _fft_small_plan_set_bound(P, FLINT_MIN(an, bn), 2*modbits, 4);
    FLINT_ASSERT(success);
    (void) success;

    _fft_small_plan_set_normalizers(P);

    /* the block rounding above can exceed a short plan's capacity;
       cap at the transform actually chosen (a no-op at block depths,
       where ztrunc dominates both roundings) */
    atrunc = n_min(atrunc, P->ztrunc);
    btrunc = n_min(btrunc, P->ztrunc);

    par.mod = mod;
    par.write_carry_out = write_carry_out;

    A.a = a; A.an = an; A.aaux = 0; A.atrunc = atrunc;
    A.b = b; A.bn = bn; A.baux = 0; A.btrunc = btrunc;
    A.bfft = NULL;
    A.bfft_stride = 0;
    A.squaring = squaring;
    A.transposed = 0;
    A.tune_n = FLINT_MIN(an, bn);
    A.min_an_bn = FLINT_MIN(an, bn);
    A.mod_fn = _mod_fn;
    A.s2_finish = _radix_s2_finish;
    A.crt_fn = _radix_crt_fn(P->np);
    A.params = &par;
    A.tuning = &_radix_conv_tuning;

    _fft_small_conv(z, P, &A);

    fft_small_plan_clear(P);
}

/*
    The transposed product (see _nmod_poly_mul_tr_mpn_ctx): the m
    coefficients z[j] = sum_i r[i] v[i + j], 0 <= j < m, with v read as
    padded with zeros to length rn + m - 1, followed by the radix carry
    chain and, if write_carry_out, its final carry in z[m]. With r the
    reversal of an operand s, this is a window of the convolution of s and
    v at the cost of a product of lengths rn and m.
*/
static void _radix_mul_tr_mpn_ctx(
    ulong * z, ulong m,
    const ulong* r, ulong rn,
    const ulong* v, ulong vn,
    int write_carry_out,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong modbits = FLINT_BITS - mod.norm;
    ulong pn = rn + m - 1;
    ulong rtrunc, mtrunc;
    int success;
    fft_small_plan_t P;
    fft_small_conv_arg_struct A;
    _conv_params_struct par;

    FLINT_ASSERT(rn > 0);
    FLINT_ASSERT(vn > 0);
    FLINT_ASSERT(m > 0);
    FLINT_ASSERT(vn <= pn);

    P->R = R;
    P->sign = 0;

    _fft_small_plan_set_window(P, 0, m, pn,
        n_max(_fft_small_len_trunc(rn), _fft_small_len_trunc(m)), (ulong) 4);

    /* every output coefficient is a sum of at most rn products < B^2 */
    success = _fft_small_plan_set_bound(P, rn, 2*modbits, 4);
    FLINT_ASSERT(success);
    (void) success;

    _fft_small_plan_set_normalizers(P);

    rtrunc = n_min(_fft_small_op_trunc(rn, P), P->ztrunc);
    mtrunc = n_min(_fft_small_op_trunc(m, P), P->ztrunc);

    par.mod = mod;
    par.write_carry_out = write_carry_out;

    /* transposed: the engine's a = v goes through ifft_tr at full length
       and comes back through fft_tr truncated to mtrunc; its b = r is
       transformed as in a multiplication */
    A.a = v; A.an = vn; A.aaux = 0; A.atrunc = mtrunc;
    A.b = r; A.bn = rn; A.baux = 0; A.btrunc = rtrunc;
    A.bfft = NULL;
    A.bfft_stride = 0;
    A.squaring = 0;
    A.transposed = 1;
    A.tune_n = FLINT_MIN(rn, m);
    A.min_an_bn = rn;
    A.mod_fn = _mod_fn;
    A.s2_finish = _radix_s2_finish;
    A.crt_fn = _radix_crt_fn(P->np);
    A.params = &par;
    A.tuning = &_radix_conv_tuning;

    _fft_small_conv(z, P, &A);

    fft_small_plan_clear(P);
}

/*
    Whether the transposed product of lengths sn and m behind the window
    [zl, zh) of the convolution of lengths an and bn has a strictly shorter
    transform than the window of the full product, the transform lengths
    being computed as in _radix_mul_tr_mpn_ctx and _radix_mul_mid_mpn_ctx
    respectively (same window selection, same arguments).
*/
static int
_radix_mul_transposed_shorter(ulong an, ulong bn, ulong zl, ulong zh,
                              ulong sn, ulong m)
{
    fft_small_plan_t Pw, Pt;

    _fft_small_plan_set_window(Pw, zl, zh, an + bn - 1,
        n_max(n_round_up(an, BLK_SZ), n_round_up(bn, BLK_SZ)), (ulong) 4);
    _fft_small_plan_set_window(Pt, 0, m, sn + m - 1,
        n_max(_fft_small_len_trunc(sn), _fft_small_len_trunc(m)), (ulong) 4);

    return Pt->ztrunc < Pw->ztrunc;
}

static void _radix_mul_mpn_ctx(
    ulong * z, ulong zl, ulong zh,
    const ulong* a, ulong an,
    const ulong* b, ulong bn,
    nmod_t mod,
    mpn_ctx_t R)
{
    ulong zn, t, m, sn, ln;
    const ulong * s, * l;
    int write_carry_out;

    an = FLINT_MIN(an, zh);
    bn = FLINT_MIN(bn, zh);

    FLINT_ASSERT(an > 0);
    FLINT_ASSERT(bn > 0);

    if (zl >= zh)
        return;

    /* the convolution of the (truncated) operands has zn coefficients;
       limb zn only receives the final carry */
    zn = an + bn - 1;
    write_carry_out = (zh == an + bn);

    if (zh > zn)
    {
        if (zl >= zn)
        {
            flint_mpn_zero(z, zh - zl);
            return;
        }

        flint_mpn_zero(z + zn - zl, zh - zn);
        zh = zn;
    }

    /* drop the low limbs of each operand that only meet the other below
       zl: they only contribute to the carries from below the window,
       which are omitted. The window shifts with the operands, z does not. */
    if (zl > bn - 1)
    {
        t = zl - (bn - 1);
        a += t;
        an -= t;
        zl -= t;
        zh -= t;
    }

    if (zl > an - 1)
    {
        t = zl - (an - 1);
        b += t;
        bn -= t;
        zl -= t;
        zh -= t;
    }

    FLINT_ASSERT(an > 0 && bn > 0);
    FLINT_ASSERT(zl < zh && zh <= an + bn - 1);

    if (an >= bn)
        s = b, sn = bn, l = a, ln = an;
    else
        s = a, sn = an, l = b, ln = bn;

    m = zh - zl;

    /*
        As for _nmod_poly_mulmid_fft_small: the window [zl, zh) is the
        transpose of the multiplication by s reversed of a polynomial of
        length m, applied to l shifted by sn - 1 - zl, at the cost of a
        product of lengths sn and m instead of sn and ln. The middle
        products of Newton iterations (an operand of length about 2n, a
        window of n limbs starting at n) are the typical case. The carry
        chain starts at zl in both drivers, so the result is the same.
        The window may start before the full sums and end after them, as
        the windows of the Newton iterations of radix_inv_approx and the
        divisions do (a few guard limbs on each side): the operand l is
        then padded with zeros in front (a copy) and at the end (by the
        driver). Not for squarings (two transforms for the window product).
    */
    if (!(a == b && an == bn) && m < ln
        && _radix_mul_transposed_shorter(an, bn, zl, zh, sn, m))
    {
        slong lead = (slong) sn - 1 - (slong) zl;   /* zeros in front of v */
        ulong vn = m + sn - 1;                     /* length of v, zeros included */
        ulong i, * u;
        TMP_INIT;

        TMP_START;

        if (lead <= 0)
        {
            /* v = l[-lead, ..), with the zeros at the end supplied by
               the driver */
            FLINT_ASSERT((slong) ln + lead > 0);
            u = TMP_ARRAY_ALLOC(sn, ulong);
            for (i = 0; i < sn; i++)
                u[i] = s[sn - 1 - i];
            _radix_mul_tr_mpn_ctx(z, m, u, sn, l - lead,
                FLINT_MIN((ulong) ((slong) ln + lead), vn),
                write_carry_out, mod, R);
        }
        else
        {
            ulong glen = FLINT_MIN(ln, vn - (ulong) lead);

            u = TMP_ARRAY_ALLOC(sn + lead + glen, ulong);
            for (i = 0; i < sn; i++)
                u[i] = s[sn - 1 - i];
            flint_mpn_zero(u + sn, lead);
            flint_mpn_copyi(u + sn + lead, l, glen);
            _radix_mul_tr_mpn_ctx(z, m, u, sn, u + sn, lead + glen,
                write_carry_out, mod, R);
        }

        TMP_END;

        return;
    }

    _radix_mul_mid_mpn_ctx(z, zl, zh, a, an, b, bn, write_carry_out, mod, R);
}


void
radix_mulmid_fft_small(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong zlo, slong zhi, const radix_t radix)
{
    if (an < bn)
    {
        FLINT_SWAP(nn_srcptr, a, b);
        FLINT_SWAP(slong, an, bn);
    }

    _radix_mul_mpn_ctx(res, zlo, zhi, a, an, b, bn, radix->B, get_default_mpn_ctx());
}

#else

void
radix_mulmid_fft_small(nn_ptr FLINT_UNUSED(res), nn_srcptr FLINT_UNUSED(a),
    slong FLINT_UNUSED(an), nn_srcptr FLINT_UNUSED(b), slong FLINT_UNUSED(bn),
    slong FLINT_UNUSED(zlo), slong FLINT_UNUSED(zhi),
    const radix_t FLINT_UNUSED(radix))
{
    flint_throw(FLINT_ERROR, "radix_mulmid_fft_small: fft_small not available");
}

#endif

