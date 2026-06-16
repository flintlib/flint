.. _padic-nmod:

**padic_nmod.h** -- floating point p-adic numbers
===============================================================================


Introduction
--------------------------------------------------------------------------------

The ``padic_nmod_t`` data type represents elements of `\mathbf{Q}_p` to
precision `N`, stored in the form `x = p^v u` with `u \in \mathbf{Z}` and
`v \in \mathbf{Z} / p^N \mathbf{Z}`.
Arithmetic operations can be carried out with respect to a context
containing the prime number `p` and various pieces of pre-computed data.
Multiple precisions `N` are possible.

Independent of the context, we consider a `p`-adic number
`x = u p^v` to be in canonical form whenever either
`p \nmid u` or `u = v = 0`.

We briefly describe the interface:

The functions in this module expect arguments of type ``padic_nmod_t``, and
each variable carries its own precision.  The functions have an interface that
is similar to the MPFR functions.  In particular, they have the same semantics,
specified as follows:  Compute the requested operation exactly and then reduce
the result to the precision of the output variable.

Data structures
--------------------------------------------------------------------------------

A `p`-adic number of type ``padic_nmod_t`` comprises a mantissa `man`, and a
valuation `val`.


Context
--------------------------------------------------------------------------------

A context object for `p`-adic arithmetic contains data pertinent to 
`p`-adic computations, but which we choose not to store with each 
element individually.
Currently, this includes the prime number `p`, its ``double`` 
inverse in case of word-sized primes, precomputed powers of `p` 
in the range between `1` and `N`, and the integers modulo `p` 
and `p^N`.

.. function:: int padic_nmod_ctx_init(gr_ctx_t ctx, ulong p, slong n)

    Initialises the context ``ctx`` with the given data.

    Assumes that `p` is a prime.  This is not verified but the subsequent 
    behaviour is undefined if `p` is a composite number.

.. function:: void padic_nmod_ctx_clear(gr_ctx_t ctx)


Memory management
--------------------------------------------------------------------------------


.. function:: void padic_nmod_init(padic_nmod_t res, gr_ctx_t ctx)

    Initialises the `p`-adic number.

.. function:: void _padic_nmod_canonicalise(padic_nmod_t x, gr_ctx_t ctx)

    Brings the `p`-adic number ``x`` into canonical form.

    That is to say, ensures that either `u = v = 0` or 
    `p \nmid u`.  There is no reduction modulo a power 
    of `p`.


Randomisation
--------------------------------------------------------------------------------


.. function:: int padic_nmod_randtest(padic_nmod_t rop, flint_rand_t state, gr_ctx_t ctx)

    Sets ``rop`` to a random `p`-adic number.

.. function:: int padic_nmod_randtest_not_zero(padic_nmod_t rop, flint_rand_t state, gr_ctx_t ctx)

    Sets ``rop`` to a random non-zero `p`-adic number.


Assignments and conversions
--------------------------------------------------------------------------------

All assignment functions set the value of ``res`` from ``x``, 
reduced to the precision of ``ctx``.

.. function:: int padic_nmod_set(padic_nmod_t res, const padic_nmod_t x, gr_ctx_t ctx)

    Sets ``res`` to the `p`-adic number ``x``.

.. function:: int padic_nmod_set_ui(padic_nmod_t res, ulong x, gr_ctx_t ctx)

    Sets the `p`-adic number ``res`` to the ``ulong`` 
    integer ``x``.

.. function:: int padic_nmod_zero(padic_nmod_t res, gr_ctx_t ctx)

    Sets the `p`-adic number ``res`` to zero.

.. function:: int padic_nmod_one(padic_nmod_t res, gr_ctx_t ctx)

    Sets the `p`-adic number ``res`` to one, reduced modulo the 
    precision of ``res``.


Comparison
--------------------------------------------------------------------------------


.. function:: truth_t padic_nmod_is_zero(const padic_nmod_t x, gr_ctx_t ctx)

    Returns whether ``x`` is equal to zero.

.. function:: truth_t padic_nmod_is_one(const padic_nmod_t x, gr_ctx_t ctx)

    Returns whether ``x`` is equal to one, that is, whether 
    `u = 1` and `v = 0`.


Arithmetic operations
--------------------------------------------------------------------------------


.. function:: int padic_nmod_add(padic_nmod_t res, const padic_nmod_t a, const padic_nmod_t b, gr_ctx_t ctx)

    Sets ``res`` to the sum of ``a`` and ``b``.

.. function:: int padic_nmod_div(padic_nmod_t res, const padic_nmod_t a, const padic_nmod_t b, gr_ctx_t ctx)

    Sets ``res`` to the quotient of ``a`` and ``b``.

.. function:: int padic_nmod_inv(padic_nmod_t res, const padic_nmod_t a, gr_ctx_t ctx)

    Sets ``res`` to the inverse of ``a``.

.. function:: int padic_nmod_mul(padic_nmod_t res, const padic_nmod_t a, const padic_nmod_t b, gr_ctx_t ctx)

    Sets ``res`` to the product of ``a`` and ``b``.

.. function:: int padic_nmod_neg(padic_nmod_t res, const padic_nmod_t a, gr_ctx_t ctx)

    Sets ``res`` to the additive inverse of ``a``.

.. function:: int padic_nmod_sub(padic_nmod_t res, const padic_nmod_t a, const padic_nmod_t b, gr_ctx_t ctx)

    Sets ``res`` to the difference of ``a`` and ``b``.


Input and output
--------------------------------------------------------------------------------


.. function:: void padic_nmod_println(const padic_nmod_t x, gr_ctx_t ctx)

    Prints the string representation of the `p`-adic number ``x`` 
    to the stream ``stdout``.

