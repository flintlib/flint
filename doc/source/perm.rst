.. _perm:

**perm.h** -- permutations
===============================================================================


Memory management
--------------------------------------------------------------------------------


.. function:: slong * _perm_init(slong n)

    Initialises the permutation for use.

.. function:: void _perm_clear(slong * vec)

    Clears the permutation.


Assignment
--------------------------------------------------------------------------------


.. function:: void _perm_set(slong * res, const slong * vec, slong n)

    Sets the permutation ``res`` to the same as the permutation ``vec``.

.. function:: void _perm_one(slong * vec, slong n)

    Sets the permutation to the identity permutation.

.. function:: void _perm_inv(slong * res, const slong * vec, slong n)

    Sets ``res`` to the inverse permutation of ``vec``.
    Allows aliasing of ``res`` and ``vec``.


Composition
--------------------------------------------------------------------------------


.. function:: void _perm_compose(slong * res, const slong * vec1, const slong * vec2, slong n)

    Forms the composition `\pi_1 \circ \pi_2` of two permutations 
    `\pi_1` and `\pi_2`.  Here, `\pi_2` is applied first, that is, 
    `(\pi_1 \circ \pi_2)(i) = \pi_1(\pi_2(i))`.

    Allows aliasing of ``res``, ``vec1`` and ``vec2``.


Parity
--------------------------------------------------------------------------------


.. function:: int _perm_parity(const slong * vec, slong n)

    Returns the parity of ``vec``, 0 if the permutation is even and 1 if
    the permutation is odd.


Randomisation
--------------------------------------------------------------------------------


.. function:: int _perm_randtest(slong * vec, slong n, flint_rand_t state)

    Generates a random permutation vector of length `n` and returns
    its parity, 0 or 1.

    This function uses the Knuth shuffle algorithm to generate a uniformly 
    random permutation without retries.


Iteration
--------------------------------------------------------------------------------

.. function:: int _perm_next_lex(slong * res, slong n)

    Sets ``res`` to the next permutation after the input ``res``, with
    respect to the lexicographic ordering, if possible. Returns 1 on success,
    i.e. if the input ``res`` was not the lexicographically last permutation,
    and returns 0 otherwise (in that case ``res`` is not modified).

.. function:: int _perm_next_heap(slong * res, slong n, slong * stack)

    Sets ``res`` to the next permutation after the input ``res``, according to
    an iterative version of B.R. Heap's algorithm, if possible. The
    implementation uses a vector ``stack`` of size ``n`` to keep track of
    iteration counters. The standard way to use this function is to call it
    first with ``res`` set to the identity permutation and ``stack`` set to a
    vector of zeros; each call will try to update ``res`` and ``stack`` in
    place. Returns 1 on success, i.e. if ``res`` was updated, and 0 otherwise.
