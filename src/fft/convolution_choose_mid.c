/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "fft.h"

/* the number of output points fft_convolution actually processes for a
   given trunc at the given depth (it rounds trunc up internally) */
static slong
_fft_convolution_rounded_trunc(slong trunc, slong depth)
{
    if (depth <= 6)
        return 2 * ((trunc + 1) / 2);
    else
    {
        slong sqrt = WORD(1) << (depth / 2);
        return 2 * sqrt * ((trunc + 2 * sqrt - 1) / (2 * sqrt));
    }
}

/*
    Plan a convolution whose output is only needed at the indices
    [nlo, nhi) of a product of length len_out.

    The standard choice is the truncated transform of depth
    clog2(len_out) - 2 with trunc = len_out. Alternatively, a full
    (untruncated) transform of length M = 4 * 2^depth computes the
    product modulo x^M - 1, where coefficient k receives c_k + c_{k+M};
    the window is uncontaminated as long as M >= nhi (it fits) and
    M >= len_out - nlo (everything that wraps lands below nlo). We take
    the smallest such power of two and use it when it processes fewer
    points than the truncated transform would; for a middle product
    (2m) x m -> m that is 2m points instead of 3m when 2m is a power of
    two. Returns the depth and sets *trunc (the trunc argument for
    fft_convolution, equal to 4 * 2^depth when wrapping).
*/
slong
fft_convolution_choose_mid(slong * trunc, slong len_out, slong nlo, slong nhi)
{
    slong loglen, loglen_wrap, M;

    FLINT_ASSERT(0 <= nlo && nlo < nhi && nhi <= len_out);

    loglen = FLINT_MAX(FLINT_CLOG2(len_out), 2);
    *trunc = len_out;

    loglen_wrap = FLINT_MAX(FLINT_CLOG2(FLINT_MAX(nhi, len_out - nlo)), 2);
    M = WORD(1) << loglen_wrap;

    if (M < _fft_convolution_rounded_trunc(len_out, loglen - 2))
    {
        loglen = loglen_wrap;
        *trunc = M;
    }

    return loglen - 2;
}
