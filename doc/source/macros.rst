.. _macros:

**Macros**
===============================================================================

Flint Macros
-------------------------------------------------------------------------------

The file ``flint.h`` contains various useful macros.

The macro constant ``FLINT_BITS`` is set at compile time to be the
number of bits per limb on the machine.  FLINT requires it to be either
32 or 64 bits.  Other architectures are not currently supported.

The macro constant ``FLINT_D_BITS`` is set at compile time to be the
number of bits per double on the machine or one less than the number of
bits per limb, whichever is smaller.  This will have the value `53` or `31`
on currently supported architectures.  Numerous internal functions using
precomputed inverses only support operands up to ``FLINT_D_BITS`` bits,
hence the macro.

The macro ``FLINT_ABS(x)`` returns the absolute value of `x`
for primitive signed numerical types.  It might fail for least negative
values such as ``INT_MIN`` and ``WORD_MIN``.

The macro ``FLINT_MIN(x, y)`` returns the minimum of `x` and
`y` for primitive signed or unsigned numerical types.  This macro
is only safe to use when `x` and `y` are of the same type,
to avoid problems with integer promotion.

Similar to the previous macro, ``FLINT_MAX(x, y)`` returns the
maximum of `x` and `y`.

The function ``FLINT_BIT_COUNT(x)`` returns the number of binary bits
required to represent an ``ulong x``.  If `x` is zero, returns `0`.

Derived from this there are the two macros ``FLINT_FLOG2(x)`` and
``FLINT_CLOG2(x)`` which, for any `x \geq 1`, compute `\lfloor \log_2 x  \rfloor`
and `\lceil \log_2 x \rceil`.

To determine the current FLINT version a number of macros are available.
For example, if the current FLINT version is ``2.4.0`` then
``__FLINT_VERSION`` will have the value `2`, ``__FLINT_MINOR``
will have the value `4` and ``__FLINT_PATCHLEVEL`` will have the value
`0`.

The ``__FLINT_RELEASE`` macro gives a single number representing the FLINT
version. For example, it will have the value ``20400`` for version ``2.4.0``.

The ``FLINT_VERSION`` macro is a static text string giving the version
number, e.g. "2.4" or "2.4.1". Note that if the final digit is a zero
it is suppressed.


