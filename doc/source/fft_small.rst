.. _fft-small:

**fft_small.h** -- FFT modulo word-size primes
===============================================================================

This module currently requires building FLINT with support for
AVX2 or NEON instructions.

The truncated FFT routines used internally by this module follow the
approach described by van der Hoeven in [vdH2004]_.

Integer multiplication
--------------------------------------------------------------------------------

.. type:: mpn_ctx_struct
          mpn_ctx_t

    Context object for multiplications allowing non-FFT moduli.
    The structure contains FFT context objects for multiple FFT primes
    (currently 8) together with tables for Chinese remaindering.

.. function:: void mpn_ctx_init(mpn_ctx_t R, ulong p)

    Initialize multiplication context object with initial prime ``p``.
    Usually ``p`` is a constant provided by :func:`get_default_mpn_ctx`.

.. function:: void mpn_ctx_clear(mpn_ctx_t R)

    Free memory allocated by the context object.

.. function:: mpn_ctx_struct * get_default_mpn_ctx(void)

    Return a pointer to a cached thread-local context object used by default
    for multiplications. Calling :func:`flint_cleanup` or :func:`flint_cleanup_master`
    frees the cache.

.. function:: void * mpn_ctx_fit_buffer(mpn_ctx_t R, ulong n)

    Returns at least *n* bytes of scratch from the context's buffer,
    growing it geometrically when needed. The buffer is a single
    region: a request may reallocate it and invalidate every pointer
    handed out before, so at most one request is outstanding at a time
    and the returned pointer is valid only until the next request.
    Callers never free what they are given; the buffer persists with
    the context, so a warm thread pays no allocation.

.. function:: void * mpn_ctx_fit_buffer_reserve(mpn_ctx_t R, ulong head, ulong tail)
              void mpn_ctx_fit_buffer_release(mpn_ctx_t R)

    Pins the scratch buffer: it is grown once to hold *head* + *tail*
    bytes and the returned tail region keeps its address until the
    release. While the reservation is live the buffer can neither grow
    nor move, so :func:`mpn_ctx_fit_buffer` serves requests of up to
    *head* bytes from it and anything larger from a secondary buffer
    kept in the context under the same contract. Reserving, like any
    growth, invalidates a scratch pointer still outstanding.

    A reservation is intended for a caller that is the only user of
    fft_small on its thread while it holds one; the transformed-mpn
    ring (:func:`gr_ctx_init_transformed_mpn`) reserves at
    construction, performs only ring operations and releases at
    destruction, so in normal use the secondary buffer is never
    touched. It guarantees that an interleaved multiplication -- test
    code checking results with ``fmpz`` arithmetic, or an external
    holder of a ring context -- is served correctly rather than
    aborted. One reservation exists at a time: a second call returns
    ``NULL`` and the caller provides its own storage.

.. function:: void mpn_ctx_mpn_mul(mpn_ctx_t R, ulong * r1, const ulong * i1, ulong n1, const ulong * i2, ulong n2)
              void mpn_mul_default_mpn_ctx(nn_ptr r1, nn_srcptr i1, slong n1, nn_srcptr i2, slong n2)

    Writes to ``r1`` the product of the integers ``(i1, n1)`` and ``(i2, n2)``.
    Assumes that `n_1 \ge n_2 \ge 1`, respectively using a given context
    object ``R`` or the default thread-local object.

Polynomial arithmetic
---------------------------------------------------------------------------------

.. function:: void _nmod_poly_mul_mid_mpn_ctx(ulong * z, ulong zl, ulong zh, const ulong * a, ulong an, const ulong * b, ulong bn, nmod_t mod, mpn_ctx_t R)
              void _nmod_poly_mul_mid_default_mpn_ctx(nn_ptr res, slong zl, slong zh, nn_srcptr a, slong an, nn_srcptr b, slong bn, nmod_t mod)

    Writes to ``z`` the middle product containing coefficients in the
    range `[zl, zh)` of the product of the polynomials  ``(a, an)`` and ``(b, bn)``,
    respectively using a given context object ``R`` or the default thread-local object.
    Assumes that `an, bn \ge 1`. Coefficients of index at least
    `an + bn - 1` are set to zero.

    This uses a cyclic convolution with wraparound. The product `c = ab` of
    length `\ell = an + bn - 1` is computed (modulo each prime of the plan)
    modulo `x^L - 1` with `L` a power of 2. If `L \ge \ell`, nothing wraps
    around. A shorter `L` is used when the coefficients that wrap around only
    affect indices outside the wanted range `[zl, zh)`. Explicitly, this is
    when there is a power of 2 `L` such that `\max(\mathtt{atrunc},
    \mathtt{btrunc}, zh) \le L \le \ell` and `\ell - L \le zl`, where
    ``atrunc`` and ``btrunc`` are the operand lengths rounded up to a multiple
    of ``BLK_SZ`` (currently 256), or to a power of 2 (at least 16) for lengths
    below ``BLK_SZ``. Then `L` is at least the length of each operand, so `\ell
    < 2L` and, for `0 \le k < L`, the coefficient of `x^k` is `c_k + c_{k+L}`.
    For every `k` in `[zl, zh)`, `k + L \ge zl + L \ge \ell`, so `c_{k+L} = 0`
    and the coefficient of `x^k` is exactly `c_k`. The nonzero terms `c_j` with
    `L \le j < \ell` do wrap around, onto the coefficients
    `j - L` in `[0, \ell - L)`, which are below the window and not output.
    Since `\ell < 2 \max(an, bn) \le 2L`, there is at most one such power of 2.
    The cost therefore depends on where the window lies and not only on its
    length. For instance, for `an \ge bn`, the middle product `[bn - 1, an)`
    wraps around whenever there is a power of 2 in `[\mathtt{atrunc}, \ell]`:
    for `an = 319` and `bn = 256` it uses `L = 512` instead of 1024.

.. function:: void _nmod_poly_mul_tr_mpn_ctx(ulong * z, ulong zn, const ulong * a, ulong an, const ulong * v, ulong vn, nmod_t mod, mpn_ctx_t R)
              void _nmod_poly_mul_tr_default_mpn_ctx(nn_ptr z, slong zn, nn_srcptr a, slong an, nn_srcptr v, slong vn, nmod_t mod)

    Transposed product: writes to ``z`` the `\mathtt{zn}` coefficients
    `z_j = \sum_{0 \le i < \mathtt{an}} a_i v_{i+j}` for
    `0 \le j < \mathtt{zn}`, respectively using a given context object ``R``
    or the default thread-local object. Assumes
    `\mathtt{an}, \mathtt{vn}, \mathtt{zn} \ge 1` and
    `\mathtt{vn} \le \mathtt{an} + \mathtt{zn} - 1`, that is, ``v`` is read
    as padded with zeros up to the length `\mathtt{an} + \mathtt{zn} - 1`.

    This is the transpose of the multiplication by ``(a, an)`` of
    polynomials of length `\mathtt{zn}` [BLS2003]_, and it is computed that
    way, with the transposed transforms :func:`sd_fft_trunc_tr` and
    :func:`sd_ifft_trunc_tr`. It therefore costs the same as one product of
    lengths `\mathtt{an}` and `\mathtt{zn}`. Moreover, ``a`` is transformed
    exactly as in a multiplication by ``a``, so the same transform of ``a``
    could serve both operations.

    Link with the middle product: if ``a`` is the reversal of ``(f,
    fn)`` (so `\mathtt{an} = \mathtt{fn}`), the output is the range
    `[\mathtt{fn} - 1, \mathtt{fn} - 1 + \mathtt{zn})` of the product of ``(f,
    fn)`` and ``(v, vn)``. Computing it this way costs one product of lengths
    `\mathtt{fn}` and `\mathtt{zn}`, whereas :func:`_nmod_poly_mul_mid_mpn_ctx`
    computes a product of lengths `\mathtt{fn}` and `\mathtt{vn}`, up to its
    wraparound. In the balanced case `\mathtt{fn} = \mathtt{zn}`, that is a
    transform of length about `2\mathtt{fn}` instead of `3\mathtt{fn}`. The
    dispatcher :func:`_nmod_poly_mulmid_fft_small` uses this function for the
    windows where it gives the shorter transform. No dispatch to a classical
    algorithm is done here, so short operands are better served by
    :func:`_nmod_poly_mulmid`: the same output is ``_nmod_poly_mulmid(z, v, vn,
    ar, an, an - 1, an - 1 + zn, mod)`` with ``ar`` the reversal of ``(a,
    an)``.

.. function:: int _fmpz_poly_mul_mid_mpn_ctx(fmpz * z, ulong zl, ulong zh, const fmpz * a, ulong an, const fmpz * b, ulong bn, mpn_ctx_t R)
              int _fmpz_poly_mul_mid_default_mpn_ctx(fmpz * z, ulong zl, ulong zh, const fmpz * a, ulong an, const fmpz * b, ulong bn)

    Like the ``nmod`` functions. Performs the multiplication and returns 1
    if there are sufficiently many primes ``R`` to compute the result;
    otherwise returns 0 without touching the output.

.. function:: void _nmod_poly_divrem_mpn_ctx(ulong * q, ulong * r, const ulong * a, ulong an, const ulong * b, ulong bn, nmod_t mod, mpn_ctx_t R)

    Polynomial division with remainder.

Preconditioned polynomial arithmetic
---------------------------------------------------------------------------------

.. type:: mul_precomp_struct

.. function:: void _mul_precomp_init(mul_precomp_struct * M, const ulong * b, ulong bn, ulong btrunc, ulong depth, nmod_t mod, mpn_ctx_t R)
              void _mul_precomp_clear(mul_precomp_struct * M)

    Represents ``(b, bn)`` in transformed form for preconditioned multiplication.

.. function:: int _nmod_poly_mul_mid_precomp(ulong * z, ulong zl, ulong zh, const ulong * a, ulong an, mul_precomp_struct * M, nmod_t mod, mpn_ctx_t R)

    Polynomial multiplication given a precomputed transform ``M``.
    Returns 1 if successful, 0 if the precomputed transform is too short.

.. type:: nmod_poly_divrem_precomp_struct

.. function:: void _nmod_poly_divrem_precomp_init(nmod_poly_divrem_precomp_struct * M, const ulong * b, ulong bn, ulong Bn, nmod_t mod, mpn_ctx_t R)
              void _nmod_poly_divrem_precomp_clear(nmod_poly_divrem_precomp_struct * M)

    Represents ``(b, bn)`` and its inverse in transformed form for preconditioned multiplication.

.. function:: int _nmod_poly_divrem_precomp(ulong * q, ulong * r, const ulong * a, ulong an, nmod_poly_divrem_precomp_struct * M, nmod_t mod, mpn_ctx_t R)

    Division with remainder of ``(a, an)`` by ``(b, bn)`` given the
    precomputation ``M``, setting ``(q, an - bn + 1)`` and ``(r, bn - 1)``.
    Returns 1 if successful, 0 if the quotient length ``an - bn + 1``
    exceeds the precision ``Bn`` of the stored inverse or the precomputed
    transform is too short.


Transform plans and operands
--------------------------------------------------------------------------------

.. type:: fft_small_plan_struct
          fft_small_plan_t

    Describes one convolution shape: the primes used, transform depth,
    truncation lengths and coefficient windows. A plan is computed once
    and reused for any number of transforms of that shape.

.. type:: fft_small_op_struct
          fft_small_op_t

    One operand in transformed representation: a buffer of evaluations
    for each prime, together with the plan parameters it was built for
    and its current domain (primal data or pointwise products).

.. function:: void fft_small_op_init(fft_small_op_t X, const fft_small_plan_t P)
              void fft_small_op_init_borrowed(fft_small_op_t X, const fft_small_plan_t P, double * data)
              ulong fft_small_op_sizeof_data(const fft_small_plan_t P)
              void fft_small_op_clear(fft_small_op_t X)

    Operand storage management. The *borrowed* variant places the
    operand on caller-provided storage of ``fft_small_op_sizeof_data``
    bytes, 4096-aligned, which must outlive the operand and is not freed
    by *clear*.

.. function:: void fft_small_fft_mpn(fft_small_op_t X, const ulong * a, ulong an, const fft_small_plan_t P)
              void _fft_small_fft_mpn_trunc(fft_small_op_t X, const ulong * a, ulong an, ulong otrunc, const fft_small_plan_t P)

    Sets *X* to the forward transforms of the nonnegative integer
    ``(a, an)``, split into chunks of the plan's bit size, modulo each of
    the plan's primes. The first computes the evaluations up to the
    plan's output truncation, for products of any length the plan
    admits; the second only the first *otrunc* evaluations (at most the
    plan's truncation, a multiple of the block size), enough for
    products of at most *otrunc* coefficients: a partner shorter than
    the plan's longest product then costs only its own length.

.. function:: int fft_small_plan_init_mpn(fft_small_plan_t P, mpn_ctx_t R, ulong an_max, ulong bn_max, ulong len_bound, int is_signed)

    Plan for products of operands of at most *an_max* and *bn_max*
    limbs through the op interface, selecting the pipeline (two-prime
    or wide) automatically. The bound *len_bound* counts accumulated
    elementary products; slot values must stay below
    `\mathrm{len\_bound} \cdot 2^{2 \cdot \mathrm{bits}}` for
    exactness. When *is_signed* is nonzero the spare factor of two
    required by the centered signed exports is applied to the
    internal capacity, after the pipeline choice has been made from
    the plain accumulation count, so that short signed dot products
    (complex multiplication, small matrix multiplication) still reach
    the two-prime geometry; callers must not double the bound
    themselves. Returns 1 on success.

.. function:: int fft_small_plan_init_mpn_np(fft_small_plan_t P, mpn_ctx_t R, ulong an_max, ulong bn_max, ulong len_bound, int is_signed, ulong np)

    As :func:`fft_small_plan_init_mpn` with the prime regime
    selectable, for profiling and test code: ``np == 0`` makes the
    automatic choice, ``np == 2`` forces the two-prime geometry
    (failing when its exactness bound cannot be met), and any other
    value forces the wide selection.

.. function:: int fft_small_plan_init_mpn_cyclic(fft_small_plan_t P, mpn_ctx_t R, ulong * nn, ulong len_bound)

    Plan for products of nonnegative integers of at most ``*nn`` limbs
    reduced modulo `B^{nn} - 1` where `B = 2^{\mathrm{FLINT\_BITS}}`,
    rounding ``*nn`` up to the nearest length admitting an exact chunk
    splitting `\mathrm{bits} \cdot 2^{\mathrm{depth}} =
    \mathrm{FLINT\_BITS} \cdot nn` (the chosen value is written back).
    Returns 0 if no admissible geometry exists. Operands longer than
    ``nn`` limbs are not folded by the packing and must be reduced by
    the caller. As in the linear plans, ``len_bound`` accounts for
    bilinear accumulation performed in transform space.

    Exports are destructive on the transform data, the default
    semantics of output conversion throughout fft_small; a caller
    needing the transform afterwards copies it first, as the
    non-destructive ring conversions do.
    :func:`fft_small_export_mpn` yields a nonnegative representative of
    the wrapped product; see the source for its size bound and the
    caller-side fold. The geometry selection was tuned with
    ``profile/p-mpn_cyclic_tuning``.

.. function:: void fft_small_export_mpn_signed(ulong * z, ulong zn, int * sign, fft_small_op_t X, ulong nslots, const fft_small_plan_t P)
              void fft_small_export_mpn_signed_trunc(ulong * z, ulong zn, int * sign, fft_small_op_t X, ulong nslots, ulong lo_limbs, const fft_small_plan_t P)

    Chinese-remainder reconstruction of a signed integer result: the
    value is assembled in two's complement and the sign resolved from
    the top limb, with the magnitude written to *z*. The *trunc* variant
    writes the limbs of the value starting at limb *lo_limbs*, with an error against the
    exact value within `(-1.5, +0.5)` ulp of the lowest returned limb --
    equivalently, at most 1 from the floor-truncated value: every slot
    whose coefficient span can reach the first returned limb is included
    with carries propagated, so the error is the truncation itself
    (one-sided, below 1 ulp) plus the total mass of the wholly dropped
    slots, under half an ulp either way. With ``lo_limbs = 0`` the
    export is exact.

.. function:: void fft_small_export_mpn_signed_low(ulong * z, ulong zn, fft_small_op_t X, ulong nslots, const fft_small_plan_t P)

    The two's-complement residue of the signed reconstruction modulo
    `2^{\mathtt{zn} \cdot \mathtt{FLINT\_BITS}}`, written to *z* with no sign
    resolution. Only the slots starting below limb *zn* are read: carries
    and borrows propagate upward only, so the slots above cannot affect
    the residue, and unlike the other signed exports the window need not
    hold the value.

.. function:: void fft_small_export_mpn_signed_window(ulong * z, ulong zn, fft_small_op_t X, ulong nslots, ulong lo_limbs, const fft_small_plan_t P)

    The limbs `[\mathtt{lo\_limbs}, \mathtt{lo\_limbs} + \mathtt{zn})` of the
    two's-complement sum of the slots that can reach them, with no sign
    resolution, reading only the slots starting below limb
    ``lo_limbs + zn``; the value need not fit. With ``lo_limbs = 0`` this
    is :func:`fft_small_export_mpn_signed_low`; otherwise the slots below
    the first one that can reach limb ``lo_limbs`` are dropped as in the
    truncated export, and since every such slot lies below limb
    ``lo_limbs`` whatever value it holds (values are lifted to the
    centered range of the prime product), the window is within one unit
    of the exact one even when those slots hold garbage, such as the
    aliased coefficients of a cyclic transform.

.. function:: int fft_small_plan_init_mpn_window(fft_small_plan_t P, mpn_ctx_t R, ulong an_max, ulong bn_max, ulong len_bound, int is_signed, ulong lo, ulong hi)
              ulong fft_small_plan_mpn_nslots(const fft_small_plan_t P, ulong an_max, ulong bn_max)
              int fft_small_plan_is_cyclic(const fft_small_plan_t P)

    :func:`fft_small_plan_init_mpn` for products of which only the limbs
    `[lo, hi)` are needed. Only the chunk coefficients that can reach
    those limbs are reconstructed, and when a power of two `2^d` at least
    the operand lengths and the highest such coefficient exists with the
    coefficients past `2^d` aliasing onto ones below the lowest, the
    transforms are cyclic of length `2^d`; exports must then read windows
    inside `[lo, hi)` (:func:`fft_small_export_mpn_signed_window`), which
    stay within one unit of the exact window. ``lo = 0`` with *hi* at
    least the product length gives the ordinary plan. *nslots* is the
    number of slots the inverse transform computes for operands of at
    most *an_max* and *bn_max* limbs; *is_cyclic* tells whether slots
    alias. Shared by the complex products and the gr transformed ring.

.. function:: ulong fft_small_export_mpn_signed_limbs(const fft_small_plan_t P, ulong nslots)
              ulong fft_small_export_mpn_limbs(const fft_small_plan_t P, ulong len_bound, ulong nslots)

    The number of limbs the signed, respectively unsigned (with up to
    *len_bound* accumulated products per slot), reconstruction of
    *nslots* slots writes.

.. function:: void fft_small_ifft_trunc(fft_small_op_t X, ulong itrunc, ulong clear_to, const ulong * scale, const fft_small_plan_t P)

    Inverse transform of the first *itrunc* slots in place, threaded
    over the primes. When *scale* is not ``NULL`` it holds per-prime
    normalizers deferred past the pointwise stage, applied after the
    inverse transform (the gr ring defers them). The slots
    `[\mathrm{itrunc}, \mathrm{clear\_to})` are zeroed, for unsigned
    exports reading whole chunks past the transformed range.

.. function:: void fft_small_ifft_tr_nmod(fft_small_op_t X, const ulong * a, ulong an, nmod_t mod, const fft_small_plan_t P)
              void fft_small_fft_tr(fft_small_op_t X, ulong itrunc, const fft_small_plan_t P)

    The transposed transforms, for transposed convolutions.
    :func:`fft_small_ifft_tr_nmod` converts ``(a, an)`` and applies
    :func:`sd_ifft_trunc_tr`, leaving a *primal* operand for the pointwise
    stage; :func:`fft_small_fft_tr` applies :func:`sd_fft_trunc_tr` for an
    operand of truncation ``itrunc``, consuming the normalization the
    pointwise stage applied, exactly as :func:`fft_small_ifft` does. The
    transpose, with respect to ``a``, of the multiplication by ``b`` given
    by ``fft_small_fft_nmod`` / ``fft_small_op_mul`` / ``fft_small_ifft`` /
    ``fft_small_export_nmod`` is therefore

    .. code ::

        fft_small_fft_nmod(B, b, bn, btrunc, mod, P);
        fft_small_ifft_tr_nmod(X, v, vn, mod, P);
        fft_small_op_mul(X, X, B, P);
        fft_small_fft_tr(X, atrunc, P);
        fft_small_export_nmod_range(z, X, 0, an, mod, P);

    which performs the same operations as the multiplication and so costs
    the same. This is how :func:`_nmod_poly_mul_tr_mpn_ctx` is built.

.. function:: void fft_small_op_pool_init(fft_small_op_pool_t pool)
              int fft_small_op_pool_fit(fft_small_op_pool_t pool, mpn_ctx_t R, ulong slab_size, slong count, ulong stage_bytes, int flags)
              void fft_small_op_pool_clear(fft_small_op_pool_t pool)
              double * fft_small_op_pool_take(fft_small_op_pool_t pool)
              void fft_small_op_pool_give(fft_small_op_pool_t pool, double * d)
              int fft_small_op_pool_owns(const fft_small_op_pool_t pool, const void * d)
              ulong * fft_small_op_pool_scratch_take(fft_small_op_pool_t pool, ulong limbs, int use_slabs)
              void fft_small_op_pool_scratch_give(fft_small_op_pool_t pool, ulong * t)

    Operand storage pools: *count* slabs of *slab_size* bytes (a
    multiple of :macro:`FLINT_FFT_SMALL_ALIGNMENT`, normally
    :func:`fft_small_op_sizeof_data`) for :func:`fft_small_op_init_borrowed`,
    plus a staging region of *stage_bytes* for conversions. With
    ``FFT_SMALL_OP_POOL_RESERVE`` in *flags* the block is a reservation
    in *R*'s buffer (:func:`mpn_ctx_fit_buffer_reserve`), warm memory;
    with ``FFT_SMALL_OP_POOL_ALLOC`` it is allocated when no reservation
    is possible. *fit* returns 0 if neither applies, leaving an empty
    pool. *take* returns ``NULL`` when no slab is free. Take and give
    are thread safe; a reserved pool must be cleared on the thread
    owning *R*, after every slab has been given back. Scratch comes from
    the staging region when free and large enough, else from a free
    slab if *use_slabs*, else from a plain allocation, so concurrent
    conversions never share it. Used by the complex products and the
    gr transformed ring.

.. function:: ulong * fft_small_stage_begin(ulong * z, ulong zn, ulong need, fft_small_op_pool_struct * pool, int use_slabs)
              int fft_small_stage_end(ulong * z, ulong zn, ulong * buf, ulong len, fft_small_op_pool_struct * pool)

    Staged conversion out of a value whose reconstruction needs *need*
    limbs into a window of *zn* limbs: *begin* returns *z* when it is
    wide enough and pool scratch otherwise (*pool* may be ``NULL``:
    plain allocation). After exporting there, and possibly adjusting
    the value, *end* with its normalized length *len* copies it out,
    zero pads the window and returns 0 if it does not fit.

.. function:: int _fft_small_mpn_mulmid_complex(nn_ptr zr, nn_ptr zi, slong lo, slong hi, nn_srcptr ar, slong arn, int ar_sgn, nn_srcptr ai, slong ain, int ai_sgn, nn_srcptr br, slong brn, int br_sgn, nn_srcptr bi, slong bin, int bi_sgn)
              int _fft_small_mpn_sqrmid_complex(nn_ptr zr, nn_ptr zi, slong lo, slong hi, nn_srcptr ar, slong arn, int ar_sgn, nn_srcptr ai, slong ain, int ai_sgn)
              int _fft_small_mpn_mul_complex(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, slong lo, slong hi, nn_srcptr ar, slong arn, int ar_sgn, nn_srcptr ai, slong ain, int ai_sgn, nn_srcptr br, slong brn, int br_sgn, nn_srcptr bi, slong bin, int bi_sgn)
              int _fft_small_mpn_sqr_complex(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, slong lo, slong hi, nn_srcptr ar, slong arn, int ar_sgn, nn_srcptr ai, slong ain, int ai_sgn)

    The transformed complex products behind all the
    ``flint_mpn_*_complex`` functions (full, low, high and windowed
    products and squares), one engine. *mulmid* and *sqrmid* write the
    limbs `[lo, hi)` of the real and imaginary parts as two's-complement
    windows, within one unit of `\lfloor x / B^{lo} \rfloor` and exact
    for `lo = 0` (the residues mod `B^{hi}`); the window plan is used,
    cyclic when the window allows. *mul* and *sqr* write magnitudes and
    signed lengths of the limbs `[lo, hi)` (which must hold the parts
    shifted right by *lo* limbs), exact for `lo = 0` and otherwise within
    `(-1.5, +0.5)` units of the lowest limb, through the staged signed
    exports. Parts may be zero: only the nonzero parts are transformed,
    so a real or imaginary operand is transformed once for both output
    parts. The signs are applied in transform space for signed exports;
    an output part that is a single product with `lo = 0` takes the
    cheaper unsigned reconstruction with its known sign attached, and
    when every part does, the plan is unsigned, packing more bits per
    chunk. Each output is formed over an operand dead after its first
    product (five transforms' storage for a product, three for a
    square), from a pool whose freed operand slabs serve the staging.
    A full product (`lo = 0`) whose output parts are single products
    sharing a part -- a real or imaginary operand times a complex one
    -- is computed by :func:`_fft_small_mpn_mul_shared`, whose fused
    per-product pipeline is faster on that task than the engine.
    Return 0 if no plan exists.

.. function:: void fft_small_export_nmod_range(ulong * z, const fft_small_op_t X, ulong zl, ulong zh, nmod_t mod, const fft_small_plan_t P)

    Chinese remaindering restricted to the coefficient window
    `[zl, zh)`, which must lie inside the plan's window; *z* receives
    `zh - zl` reduced coefficients. Used for polynomial middle products.

.. macro:: FLINT_FFT_SMALL_ALIGNMENT

    Alignment, in bytes, of transform data buffers.

.. var:: ulong flint_fft_small_max_transformed_ring_size

    Upper bound, in bytes, on the transform storage a single transformed
    ring or context may plan for: its expected number of simultaneously
    live elements times the per-element data size. Constructors decline
    above the bound, so drivers fall back to slower algorithms rather
    than exhaust memory; the matrix multiplication drivers respond by
    multiplying in blocks instead. Mutable for tuning; the default is
    4 GiB on 64-bit machines.


Low-level transforms
-------------------------------------------------------------------------------

These are the truncated transforms on which the rest of the module is built.
They operate in place on an array of ``double`` of size at least `2^L`,
allocated with alignment ``FLINT_FFT_SMALL_ALIGNMENT``, holding residues modulo
the prime `n` of the context object as ``double`` values.  :func:`sd_fft_trunc`
and :func:`sd_ifft_trunc_tr` accept and produce values in `(-3n, 3n)`;
:func:`sd_fft_trunc_tr` accepts and produces values in `(-2n, 2n)`;
:func:`sd_ifft_trunc` accepts values in `(-2n, 2n)` and produces values in
`(-3n, 3n)` (see the end of ``mulmod_satisfies_bounds.c``).
The inverse transforms need the tighter input range because their butterflies
sum four entries before any reduction, and `8n < 2^{53}` for the primes `n <
2^{50}` of this module. The pointwise products leave values in `(-3n/2, 3n/2)`,
so that a forward transform, pointwise product and inverse transform can be
chained without reductions in between.
The four transforms below grow the context to depth `L` themselves, with
``sd_fft_ctx_fit_depth``.

.. function:: void sd_fft_trunc(sd_fft_ctx_t Q, double * d, ulong L, ulong itrunc, ulong otrunc)

    Truncated forward transform of length `2^L`. Reads the first ``itrunc``
    entries of ``d``, treating the remaining entries as zero, and writes the
    first ``otrunc`` entries. With
    `w_i = \mathtt{sd\_fft\_ctx\_w}(Q, \mathtt{sd\_fft\_ctx\_trunc\_index}(L, i))`
    the evaluation point at position `i`, the output is the evaluation
    `\mathtt{out}[i] = \sum_{0 \le t < \mathtt{itrunc}} \mathtt{in}[t] \, w_i^t`
    for `0 \le i < \mathtt{otrunc}`. Any
    `1 \le \mathtt{itrunc}, \mathtt{otrunc} \le 2^L` is allowed. The
    evaluation points do not depend on `L`: if
    `\max(\mathtt{itrunc}, \mathtt{otrunc}) \le 2^{L-1}`, the result is that
    of depth `L - 1`, at about the same cost, so it usually only makes sense to
    have `\max(\mathtt{itrunc}, \mathtt{otrunc}) \ge 2^{L-1}`. Recovering the
    input with :func:`sd_ifft_trunc` requires `\mathtt{otrunc} \ge
    \mathtt{itrunc}`. For `L \le \mathtt{LG\_BLK\_SZ}` the entries from
    ``itrunc`` on are set to zero and all `2^L` outputs are computed; for
    larger `L` both truncations are rounded up to multiples of
    `\mathtt{BLK\_SZ}`.

.. function:: void sd_ifft_trunc(sd_fft_ctx_t Q, double * d, ulong L, ulong trunc)

    Truncated inverse transform, the inverse of :func:`sd_fft_trunc` with
    ``itrunc`` `=` ``otrunc`` `=` ``trunc`` up to the factor `2^L`, which is
    left to the caller. Both are linear maps on the first
    `\mathtt{BLK\_SZ} \lceil \mathtt{trunc}/\mathtt{BLK\_SZ} \rceil` entries of
    ``d`` (on the first `2^L` entries when `L \le \mathtt{LG\_BLK\_SZ}`).
    In particular, for `L \le \mathtt{LG\_BLK\_SZ}` the parameter ``trunc``
    is ignored and the full inverse of the `2^L` entries is computed: this is
    the truncated inverse only if all `2^L` entries hold evaluations, as
    :func:`sd_fft_trunc` writes them at such depths. For larger `L`, the
    entries up to the next multiple of `\mathtt{BLK\_SZ}` after ``trunc``
    are read as evaluations. See [vdH2004]_.

.. function:: void sd_fft_trunc_tr(sd_fft_ctx_t Q, double * d, ulong L, ulong itrunc, ulong otrunc)
              void sd_ifft_trunc_tr(sd_fft_ctx_t Q, double * d, ulong L, ulong trunc)

    Transposes of the two transforms above, in the sense of the transposition
    principle (see [BLS2003]_): with the same parameters,
    `\langle \mathtt{sd\_fft\_trunc}(a), v \rangle = \langle a, \mathtt{sd\_fft\_trunc\_tr}(v) \rangle`
    and
    `\langle \mathtt{sd\_ifft\_trunc}(a), v \rangle = \langle a, \mathtt{sd\_ifft\_trunc\_tr}(v) \rangle`.
    In particular :func:`sd_fft_trunc_tr` reads the first ``otrunc`` entries and
    writes the first ``itrunc`` entries, the other way round from
    :func:`sd_fft_trunc`; the parameters are named after the transform being
    transposed so that transposing a call means adding the suffix and nothing
    else. With `w_i` as above, its output is the power projection
    `\mathtt{out}[t] = \sum_{0 \le i < \mathtt{otrunc}} \mathtt{in}[i] \, w_i^t`
    for `0 \le t < \mathtt{itrunc}`, so that the input is indexed exactly as
    the output of :func:`sd_fft_trunc`. The two functions perform the same
    operations as the transforms they transpose, in reverse order, and so have
    the same cost. Composed as
    `v \mapsto \mathtt{fft\_tr}(\mathtt{fft}(b) \cdot \mathtt{ifft\_tr}(v))`,
    they give the transpose of multiplication by a fixed `b`.


Products with a shared or precomputed transform
--------------------------------------------------------------------------------

The functions in this section multiply nonnegative integers given as
limb arrays, reusing the transform of one operand across several
products: two products with a common operand, or any number of
products by an operand fixed in advance (the power tables of binary
splitting). Like :func:`mpn_ctx_mpn_mul` they have a plan for all
lengths: the two-prime geometry up to ``FLINT_MPN_MUL_NP2_MAX_AN`` limbs,
where each product runs fused (pack, forward transforms, pointwise
product, inverse transforms and reconstruction while its data is hot),
and the wide geometry above, threaded over the primes, with each
partner's transforms truncated to its own product. They use the
calling thread's default context (:func:`get_default_mpn_ctx`) for
scratch. The program ``profile/p-mpn_mul_shared`` compares them with
independent products.

.. function:: void _fft_small_mpn_mul_shared(nn_ptr z1, nn_srcptr a, ulong an, nn_ptr z2, nn_srcptr b, ulong bn, nn_srcptr x, ulong xn)

    Sets ``(z1, an + xn)`` to `a x` and ``(z2, bn + xn)`` to `b x`, with
    the transform of `x` computed once: three forward and two inverse
    transforms against four and two for two independent products. All
    lengths must be positive; the operands may have any relative
    lengths. The output *z1* may overlap *a* and *x* but not *b*; *z2*
    may overlap any input but not *z1*. Above the FFT threshold of
    :func:`flint_mpn_mul` this is faster than two independent products
    (by about 10% for products of similar lengths); see
    :func:`flint_mpn_mul_shared`, which chooses between the two.

.. type:: fft_small_mpn_mul_precomp_struct
          fft_small_mpn_mul_precomp_t

    The transform of a fixed operand `x`, for products by partners of
    up to a planned length.

.. function:: ulong _fft_small_mpn_mul_precomp_size(ulong xn, ulong yn_max)
              void _fft_small_mpn_mul_precomp_init(fft_small_mpn_mul_precomp_t P, nn_srcptr x, ulong xn, ulong yn_max, double * storage)
              void _fft_small_mpn_mul_precomp_clear(fft_small_mpn_mul_precomp_t P)

    Prepares *P* for products of ``(x, xn)`` by partners of 1 to
    *yn_max* limbs, computing the transform of *x* once. The transform
    is stored in *storage*, which must hold the number of doubles given
    by *_size* (64-byte aligned), or in an allocation of *init*'s own
    when *storage* is ``NULL``; only the evaluations the products read
    are kept, up to half of a full transform. The array *x* is not
    referenced after *init*. Lengths must be positive.

.. function:: void _fft_small_mpn_mul_precomp(nn_ptr z, const fft_small_mpn_mul_precomp_t P, nn_srcptr y, ulong yn)

    Sets ``(z, xn + yn)`` to `x y`, with `x` the operand of *P* and
    `1 \le yn \le` *yn_max*; *z* may overlap *y*. Each product costs
    one conversion and forward transform less than an independent
    product (26-31% faster in the two-prime range, 16-22% in the wide
    range, for partners of about twice the length of `x`). *P* is read
    only, so products may run concurrently on several threads.


Negacyclic pointwise multiplication
-------------------------------------------------------------------------------

These functions multiply residues modulo `2^N + 1` by splitting the
operands into `m` digits of `b = N/m` bits, so that the product is a
length-`m` negacyclic convolution of the digit sequences, computed via
the small-prime FFT machinery of this module and Chinese remaindering.
They serve as the pointwise stage of Schoenhage-Strassen convolutions:
see :func:`fft_convolution_sd_fft` and the ``fmpz_poly`` multiplication
functions, which engage this engine for large coefficient sizes.

Operands and results occupy `n + 1 = N/\mathtt{FLINT\_BITS} + 1` limbs.
An *operand* may have top limb `0` or `1` with arbitrary remaining
limbs; its value is interpreted modulo `2^N + 1` (in particular the
representation need not be reduced, and both representatives of a
residue with top limb set are accepted). The *result* is returned
fully reduced: the unique representative in `[0, 2^N]`, with top limb
`1` exactly when the value is `2^N` and the remaining limbs then
zero.

.. function:: void sd_fft_mpn_mulmod_2expp1_ctx_init(sd_fft_mpn_mulmod_2expp1_ctx_struct * C, mpn_ctx_struct * R, slong N, slong m)

    Initializes a context for multiplication modulo `2^N + 1` with
    digit count *m*, requiring `N = mb` for an integer digit size
    *b*, *m* a power of two at least 16, and *N* a multiple of
    ``FLINT_BITS``. Passing `m = 0` chooses the digit count
    internally by an estimated-cost comparison over the instantiated
    digit sizes dividing *N*. The number of FFT primes is chosen as
    the smallest count, between 3 and 8, whose exact prime product
    exceeds the coefficient bound `2^{2b + \log_2 m + 1}` and for
    which a digit-splitting input pass is instantiated. The primes
    and transform contexts are taken from *R*, normally
    ``get_default_mpn_ctx()``.

    An exception is raised when *N* is not a positive multiple of
    ``FLINT_BITS`` or no supported configuration exists.
    With `m = 0` and the default context, initialization succeeds
    precisely when `N = 2^k b` for some `k \geq 4` and an
    instantiated digit size
    `b \in \{64, 84, 88, 92, 112, 116, 120, 126, 128, 136, 140,
    144, 160, 164, 168, 184, 188, 192\}`, subject to prime
    capacity: unlimited in practice (`k \geq 40` supported) for
    `b \leq 168`, and `k \leq 28`, `20`, `12` for
    `b = 184`, `188`, `192` respectively, since those digit sizes
    admit no larger prime count to grow into. In particular every
    power of two `N \geq 1024` is supported through `b = 64` and
    `b = 128`.

.. function:: void sd_fft_mpn_mulmod_2expp1(sd_fft_mpn_mulmod_2expp1_ctx_struct * C, nn_ptr z, nn_srcptr x, nn_srcptr y, sd_fft_mpn_mulmod_2expp1_scratch_struct * S)
              void sd_fft_mpn_mulmod_2expp1_scratch_init(sd_fft_mpn_mulmod_2expp1_scratch_struct * S, const sd_fft_mpn_mulmod_2expp1_ctx_struct * C)
              void sd_fft_mpn_mulmod_2expp1_scratch_clear(sd_fft_mpn_mulmod_2expp1_scratch_struct * S)

    The context is read-only during multiplication; all mutable state
    lives in a scratch object, so a single context may serve several
    threads concurrently, each with its own scratch.
    The pointwise stage of the
    ``fmpz_poly`` Schoenhage-Strassen convolutions multiplies its
    coefficients in parallel over the available threads this way.

.. function:: void sd_fft_mpn_mulmod_2expp1_ctx_clear(sd_fft_mpn_mulmod_2expp1_ctx_struct * C)

    Releases the memory held by the context.

