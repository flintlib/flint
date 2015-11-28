.. _ulong_extras:

**ulong_extras** -- machine word arithmetic
:::::::::::::::::::::::::::::::::::::::::::

The *ulong_extras* module provides fast primitives for unsigned machine word
arithmetic:

- bit manipulation
- division/modular arithmetic with precomputed inverses
- greatest common divisor
- integer square/cube/n-th roots
- integer logarithms
- jacobi symbols
- probable prime testing
- primality proving
- integer factorisation
- number theoretic functions (Jacobi, Euler phi, Mobius)
- random generation

Also see *longlong.h* for assembly optimised low level word operations.

Types, macros and constants
---------------------------

The following types, macros and constants relevant to *ulong_extras* are
defined in *flint.h* and available throughout the whole of Flint.

.. type:: ulong

    This is the basic type that *ulong_extras* functions operate on. It is
    defined to be the same as GMP's *mp_limb_t*, which is to say a 32 or 64 bit
    unsigned integer depending on the ABI of the machine. It's usually
    implemented as either an *unsigned long* or *unsigned long long*.

    This type is intended as a portable replacement for *unsigned long* which
    is not portable, especially to Windows and MIPS.

.. type:: slong

    As per *ulong*, but signed. This type is more often used for loop counters
    and lengths of arrays than for arithmetic. It is defined to be the same as
    GMP's *mp_limb_signed_t* and is usually implemented as either a *long* or
    *long long*.

    This type is intended as a portable replacement for *long* which is not
    portable, especially to Windows and MIPS.

.. macro:: UWORD_MAX

    The largest integer able to be stored in a *ulong*. Usually either
    `2^{32} - 1` or `2^{64} - 1`.

.. macro:: UWORD_MIN

    The smallest integer able to be stored in a *ulong*. Usually `0`.

.. macro:: WORD_MAX

    The largest integer able to be stored in an *slong*. Usually either
    `2^{31} - 1` or `2^{63} - 1`.

.. macro:: WORD_MIN

    The smallest integer able to be stored in an *slong*. Usually either
    `-2^{31}` or `-2^{63}`.

Input/output
------------

Flint provides its own printing/reading functions which can deal with the
*ulong* and *slong* types. This makes printing, reading, file and stream
I/O portable between unixes and other operating systems such as Windows.

.. function:: int flint_printf(char * fmt, ...)

    As per the standard C *printf* function, but supports the following
    additional format specifiers:

    - *%wd* : print an *slong* in decimal format
    - *%wu* : print a *ulong* in decimal format
    - *%wx* : print a *ulong* in hexadecimal format
    - *%\*wd* : print an *slong*, space padded in a field of the given width

    Returns the number of characters printed.

.. function:: int flint_sprintf(char * str, const char * fmt, ...)

    As per the standard C *sprintf* function, but with the additional format
    specifiers provided by *flint_printf*.

    Returns the number of characters printed.

.. function:: int flint_fprintf(FILE * stream, const char * fmt, ...)

    As per the standard C *fprintf* function, but with the additional format
    specifiers provided by *flint_printf*.

    Returns the number of characters printed.

.. function:: int flint_scanf(const char * fmt, ...)

    As per the standard C *scanf* function, but supports the following
    additional format specifiers:

    - *%wd* : read an *slong* in decimal format
    - *%wu* : read a *ulong* in decimal format
    - *%wx* : read a *ulong* in hexadecimal format

    Returns the number of items in the argument list successfully filled.

.. function:: int flint_sscanf(const char * str, const char * fmt, ...)

    As per the standard C *sscanf* function, but supports the additional format
    specifiers provided by *flint_scanf*.

    Returns the number of items in the argument list successfully filled.

.. function:: int flint_fscanf(FILE * stream, const char * fmt, ...)

    As per the standard C *fscanf* function, but supports the additional format
    specifiers provided by *flint_scanf*.

    Returns the number of items in the argument list successfully filled.

Bit manipulation
----------------

.. function:: ulong n_revbin(ulong n, ulong b)

    Considering `n` to be a binary number with `b` bits, return the number with
    the reverse binary bit representation. For example *n_revbin(3, 4)*
    would consider `3` as the `4` bit binary number `0011`, which it would
    reverse to give the return value of `12`, with binary representation `1100`.

    **Conditions:** We require `b \leq` *FLINT_BITS*. Only the lower `b` bits of
    `n` are read and the remaining bits can be arbitrary.

    **Algorithm:** If `b \leq 8` we swap the bits in the least significant byte
    using a lookup table, then shift to the required number of bits, otherwise
    we swap all bits in the word and then shift to the required number of bits.

    To swap all the bits in a word we mask and swap alternate bits, then
    alternate pairs of bits, then alternate nibbles and finally swap all bytes
    in the word using the *byte_swap* macro from *longlong.h*.

Arithmetic with precomputed inverse
-----------------------------------

.. function:: ulong n_preinvert_limb(ulong n)

    Return a Moller-Granlund precomputed inverse of the word `n`.

    **Conditions:** We require `n > 0`.

    **Algorithm:** This function normalises `n` (shifts it so that its most
    significant bit is set) and then computes the precomputed inverse
    `(\beta^2 - 1)/n - \beta` where `\beta = 2^w` with `w =` *FLINT_BITS*,
    using the *invert_limb* macro from longlong.h.

    For details of the precomputed inverse see [MolGra2011]_. 

.. function:: ulong n_div2_preinv(ulong a, ulong n, ulong ninv)

    Return the Euclidean quotient of `a` by `n` given a precomputed inverse
    *ninv* provided by *n_preinvert_limb*.

    **Conditions:** We require `n > 0`.

    **Algorithm:** Both `a` and `n` are shifted left by the same number of
    bits `b` so that `n` is normalised. In general `a` now occupies two limbs.
    We then compute the quotient using Algorithm 4 of [MolGra2011]_.

.. function:: ulong n_mod2_preinv(ulong a, ulong n, ulong ninv)

    Return the Euclidean remainder of `a` divided by `n` given a precomputed
    inverse *ninv* provided by *n_preinvert_limb*.

    **Conditions:** We require `n > 0`.

    **Algorithm:** Both `a` and `n` are shifted left by the same number of
    bits `b` so that `n` is normalised. In general `a` now occupies two limbs.
    We then reduce `a` modulo `n` using Algorithm 4 of [MolGra2011]_.
    Finally the resulting remainder is shifted right by `b` bits.

.. function:: ulong n_divrem2_preinv(ulong * q, ulong a, ulong n, ulong ninv)

    Return the Euclidean remainder and set `q` to the Euclidean quotient of `a`
    by `n`, given a precomputed inverse *ninv* provided by *n_preinvert_limb*.

    **Conditions:** We require `n > 0`.

    **Algorithm:** Both `a` and `n` are shifted left by the same number of
    bits `b` so that `n` is normalised. In general `a` now occupies two limbs.
    We then reduce `a` modulo `n` using Algorithm 4 of [MolGra2011]_.
    Finally the resulting remainder is shifted right by `b` bits.

.. function:: ulong n_ll_mod_preinv(ulong a1, ulong a0, ulong n, ulong ninv)

    Given a two word input `a = \langle a_1, a_0 \rangle`, return the
    remainder upon division of `a` by `n`, given a precomputed inverse *ninv*
    provided by *n_preinvert_limb*.

    This function is useful for delayed reduction and reduction of products
    that have accumulated in no more than two words.

    **Conditions:** We require `n > 0`.

    **Algorithm:** If the word `a_1` is not reduced modulo `n` we reduce it
    as per *n_mod2_preinv*. The remainder of `a` divided by `n` is then
    computed using Algorithm 4 of [MolGra2011]_.


