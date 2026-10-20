.. _n-fft:

**n_fft.h** -- FFT modulo word-size primes, integer arithmetic
===============================================================================

This module provides discrete Fourier transforms (and their inverses and
transposes) over `\mathbb{Z}/p\mathbb{Z}` for an odd prime `p < 2^{62}` with
`p - 1` divisible by `8`, using only integer arithmetic on ``ulong`` values
(Shoup's modular multiplication with precomputation, and lazy butterflies in
the style of Harvey). Unlike :ref:`fft_small <fft-small>`, it does not require
SIMD support and it accepts moduli of up to 62 bits.

Transform lengths are powers of two, `2^{depth}`. The evaluation points are
listed in *bit-reversed order*: writing `w` for the chosen primitive
`2^{depth}`-th root of unity, the `i`-th point is `w^{\mathrm{br}(i)}` where
`\mathrm{br}` is the bit-reversal permutation on `depth` bits. Concretely the
list of points starts with `1, -1, I, -I, J, -J, IJ, -IJ, \ldots` where
`I^2 = -1` and `J^2 = I`. Consecutive points come in pairs `\pm w_k`, with
`w_k` stored in the context's table ``tab_w[2*k]``.

Context object
--------------------------------------------------------------------------------

.. type:: n_fft_ctx_struct
          n_fft_ctx_t

    Context object holding the modulus ``mod``, the maximal supported depth
    ``max_depth`` (the `2`-adic valuation of ``mod - 1``), the ``cofactor``
    such that ``mod == 1 + cofactor * 2**max_depth``, the currently supported
    ``depth``, and precomputed tables: powers of a primitive
    `2^{\mathrm{max\_depth}}`-th root of unity ``w`` in bit-reversed order
    (``tab_w``), the same for `1/w` (``tab_iw``), the powers `w^{2^k}`
    (``tab_w2``) and the inverses of `2^k` (``tab_inv2``), each together with
    its precomputed quotient for Shoup multiplication. The tables ``tab_w``
    and ``tab_iw`` have `2^{\mathrm{depth}}` entries each; the total memory
    is about `2^{\mathrm{depth} + 1}` words.

.. function:: void n_fft_ctx_init2(n_fft_ctx_t F, ulong depth, ulong p)
              void n_fft_ctx_init(n_fft_ctx_t F, ulong p)

    Initialize ``F`` for the prime ``p``, finding a primitive root of unity of
    order `2^{\mathrm{max\_depth}}`, and precompute tables supporting
    transforms of length up to `2^{\mathrm{depth}}`. If ``depth`` is less than
    `3` it is taken to be `3`, and if it is more than ``max_depth`` it is
    reduced to ``max_depth``. The second version uses a default depth
    (``N_FFT_CTX_DEFAULT_DEPTH``, currently `12`). Requirements (not checked
    unless assertions are enabled): ``p`` is an odd prime with
    `p < 2^{\mathrm{FLINT\_BITS}-2}` and `8 \mid p - 1`.

.. function:: void n_fft_ctx_init2_root(n_fft_ctx_t F, ulong w, ulong max_depth, ulong cofactor, ulong depth, ulong p)
              void n_fft_ctx_init_root(n_fft_ctx_t F, ulong w, ulong max_depth, ulong cofactor, ulong p)

    Same as above, but with a user-provided primitive root ``w`` of order
    `2^{\mathrm{max\_depth}}`, where ``p == 1 + cofactor * 2**max_depth`` and
    ``max_depth >= 3``.

.. function:: void n_fft_ctx_fit_depth(n_fft_ctx_t F, ulong depth)

    Grow the tables of ``F`` (if needed) so that transforms of length up to
    `2^{\mathrm{depth}}` are supported; ``depth`` is reduced to ``max_depth``
    if larger. This reallocates the tables, so it must not be called
    concurrently with other uses of ``F``.

.. function:: void n_fft_ctx_clear(n_fft_ctx_t F)

    Free the memory allocated by ``F``.

Power of two length transforms
--------------------------------------------------------------------------------

All functions below work in place on an array ``p`` of length
`2^{\mathrm{depth}}`, with entries in `[0, \mathrm{mod})` in input and in
output, and require ``depth <= F->depth`` (not checked). We write `w_i`,
`0 \le i < 2^{\mathrm{depth}}`, for the evaluation points in bit-reversed
order as described above, and `\mathrm{DFT}` for the matrix
`[w_i^j]_{0 \le i, j < 2^{\mathrm{depth}}}`.

.. function:: void n_fft_dft(nn_ptr p, ulong depth, n_fft_ctx_t F)

    Transform the coefficients `p_0, \ldots, p_{\ell-1}` of a polynomial
    (`\ell = 2^{\mathrm{depth}}`) into its evaluations `p(w_0), \ldots,
    p(w_{\ell-1})`, i.e. multiply ``p`` by `\mathrm{DFT}`.

.. function:: void n_fft_idft(nn_ptr p, ulong depth, n_fft_ctx_t F)

    Inverse of :func:`n_fft_dft`: transform the evaluations at
    `w_0, \ldots, w_{\ell-1}` into the coefficients of the unique polynomial
    of degree less than `\ell` taking these values.

.. function:: void n_fft_dft_t(nn_ptr p, ulong depth, n_fft_ctx_t F)

    Transposed transform: multiply the row vector ``p`` by `\mathrm{DFT}`,
    i.e. compute the weighted power sums `\sum_i p_i w_i^j` for
    `0 \le j < \ell`.

.. function:: void n_fft_idft_t(nn_ptr p, ulong depth, n_fft_ctx_t F)

    Transposed inverse transform: multiply the row vector ``p`` by
    `\mathrm{DFT}^{-1}`. This computes the coefficients `q_i` of the partial
    fraction decomposition `\sum_j p_j x^j = \sum_i q_i / (1 - w_i x) + O(x^\ell)`.

Truncated transforms
--------------------------------------------------------------------------------

The truncated variants use the same evaluation points but the input length
and the number of evaluations need not be powers of two. Below, ``len``
denotes the smallest power of two which is at least the number of
evaluations, and ``depth`` its logarithm; the functions require
``depth <= F->depth`` (not checked). The lengths must be multiples of `4`
(except that :func:`n_fft_itft` also accepts lengths at most `2`); the
internal helpers ``n_fft_tft_prepare`` and ``n_fft_itft_prepare`` in
``n_fft/impl.h`` compute suitable parameters from arbitrary lengths together
with the required array length.

.. function:: void n_fft_tft(nn_ptr p, ulong ilen, ulong olen, n_fft_ctx_t F)

    Transform the coefficients `p_0, \ldots, p_{\mathrm{ilen}-1}` of a
    polynomial into its evaluations `p(w_0), \ldots, p(w_{\mathrm{olen}-1})`.
    ``ilen`` and ``olen`` are multiples of `4` (``ilen == 0`` yields zero
    evaluations). The array ``p`` must have space for ``max(ilen, len)``
    entries; the entries beyond ``ilen`` may hold any value and may be
    overwritten, and in output only the first ``olen`` entries are
    meaningful.

.. function:: void n_fft_itft(nn_ptr p, ulong iolen, n_fft_ctx_t F)

    Transform the evaluations `p(w_0), \ldots, p(w_{\mathrm{iolen}-1})` of a
    polynomial of degree less than ``iolen`` into its coefficients. ``iolen``
    is either at most `2` or a multiple of `4`. The array ``p`` must have
    space for ``len`` entries; the entries beyond ``iolen`` may hold any value
    and may be overwritten, and in output only the first ``iolen`` entries
    are meaningful.
