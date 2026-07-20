/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "ulong_extras.h"

#if FLINT_HAVE_NATIVE_mpn_mulmid_n && FLINT_HAVE_FFT_SMALL
# define MULMID_KARATSUBA_CUTOFF 36
#else
# define MULMID_KARATSUBA_CUTOFF 20
#endif

#if FLINT_HAVE_FFT_SMALL
# define MULMID_FFT_SMALL_CUTOFF 250
#endif

#if !FLINT_HAVE_NATIVE_mpn_mulmid_n && !FLINT_HAVE_FFT_SMALL
/* bare build: above this a full product (GMP/FFT-backed flint_mpn_mul) plus a
   slice beats schoolbook for a wide middle window */
# define MULMID_BARE_MUL_CUTOFF 128
#endif

/* Zeros we are willing to append to a length-len operand, and the test for it
   (also true when target <= len, i.e. when the operand is merely truncated). */
#define MULMID_PAD_LIMIT(len) FLINT_MAX((mp_size_t) 4, (len) / 32)
#define MULMID_PAD_OK(target, len) ((target) - (len) <= MULMID_PAD_LIMIT(len))

void
flint_mpn_mulmid(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn,
                 mp_size_t zlo, mp_size_t zhi)
{
    mp_size_t ac, bc, m, M, zn, top, lim;
#if FLINT_HAVE_NATIVE_mpn_mulmid_n
    mp_size_t vac, vbc, vpa, vqb, La, Lb, zlo2, nb, vhalf;
    int via_ok;
#endif

    FLINT_ASSERT(an >= 1);
    FLINT_ASSERT(bn >= 1);
    FLINT_ASSERT(zlo >= 0);
    FLINT_ASSERT(zhi > zlo);
    FLINT_ASSERT(zhi <= an + bn);

    ac = FLINT_MIN(an, zhi);
    bc = FLINT_MIN(bn, zhi);
    m = FLINT_MIN(ac, bc);
    M = FLINT_MAX(ac, bc);
    zn = zhi - zlo;
    top = an + bn - zhi;            /* limbs above the window */
    lim = MULMID_PAD_LIMIT(m);

    if (m < MULMID_KARATSUBA_CUTOFF || zn * zn <= m)
    {
        flint_mpn_mulmid_classical(z, a, an, b, bn, zlo, zhi);
        return;
    }

    /* (A) full / nearly-full window */
    if (zlo <= lim && top <= lim)
    {
        flint_mpn_mulmid_via_mul(z, a, an, b, bn, zlo, zhi);
        return;
    }

    /* (B) low / nearly-low window */
    if (zlo <= lim && MULMID_PAD_OK(zhi, an) && MULMID_PAD_OK(zhi, bn)
#if FLINT_HAVE_FFT_SMALL
            && (m < 2 * MULMID_FFT_SMALL_CUTOFF || zn < MULMID_FFT_SMALL_CUTOFF)
#endif
        )
    {
        flint_mpn_mulmid_via_mullow_n(z, a, an, b, bn, zlo, zhi);
        return;
    }

    /* (C) high / nearly-high window */
    if (top <= lim && zlo + lim >= M)
    {
        flint_mpn_mulmid_via_mulhigh_n(z, a, an, b, bn, zlo, zhi);
        return;
    }

#if FLINT_HAVE_FFT_SMALL
    /* large operands with a large (non-full, non-low, non-high) window: FFT */
    if (m >= MULMID_FFT_SMALL_CUTOFF && zn >= MULMID_FFT_SMALL_CUTOFF)
    {
        flint_mpn_mulmid_fft_small(z, a, an, b, bn, zlo, zhi);
        return;
    }
#endif

#if FLINT_HAVE_NATIVE_mpn_mulmid_n
    /* Padding flint_mpn_mulmid_via_n_padded would incur: after clamping to zhi
       and slicing the low non-contributing limbs, the longer sliced operand La
       fills a 2n-1 slot and the shorter Lb an n slot. */
    vac = FLINT_MIN(an, zhi);
    vbc = FLINT_MIN(bn, zhi);
    vpa = (zlo > vbc - 1) ? zlo - (vbc - 1) : 0;
    vqb = (zlo > vac - 1) ? zlo - (vac - 1) : 0;
    La = vac - vpa;
    Lb = vbc - vqb;
    via_ok = 0;
    if (La > 0 && Lb > 0)
    {
        zlo2 = zlo - vpa - vqb;
        if (La < Lb) { mp_size_t t2 = La; La = Lb; Lb = t2; }
        nb = zn;
        if (Lb > nb) nb = Lb;
        vhalf = (La + Lb - zlo2 + 1) / 2;
        if (vhalf > nb) nb = vhalf;
        /* via_n_padded reduces to a balanced middle product of size nb; even
           with padding this is worth it whenever nb stays comfortably below the
           M-sized problem the other reductions would build, so a generous
           padding tolerance (max(8, len/8)) is used here. */
        via_ok = (2 * nb - 1) - La <= FLINT_MAX((mp_size_t) 8, La / 8)
              && nb - Lb <= FLINT_MAX((mp_size_t) 8, Lb / 8);
    }

    /* (D) balanced Karatsuba/Toom reduction, when its padding is economical */
    if (via_ok)
    {
        flint_mpn_mulmid_via_n_padded(z, a, an, b, bn, zlo, zhi);
        return;
    }
#endif

#if FLINT_HAVE_NATIVE_mpn_mulmid_n || FLINT_HAVE_FFT_SMALL
    /* medium middle window: schoolbook (large wide windows already handled) */
    flint_mpn_mulmid_classical(z, a, an, b, bn, zlo, zhi);
#else
    /* bare build: a wide middle window is cheaper as a full product plus a slice
       (flint_mpn_mul is GMP/FFT-backed); a narrow one stays schoolbook */
    if (m >= MULMID_BARE_MUL_CUTOFF && zn >= MULMID_BARE_MUL_CUTOFF)
        flint_mpn_mulmid_via_mul(z, a, an, b, bn, zlo, zhi);
    else
        flint_mpn_mulmid_classical(z, a, an, b, bn, zlo, zhi);
#endif
}
