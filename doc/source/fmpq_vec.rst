.. _fmpq-vec:

**fmpq_vec.h** -- vectors over rational numbers
===============================================================================

Memory management
--------------------------------------------------------------------------------


.. function:: fmpq * _fmpq_vec_init(slong n)

    Initialises a vector of ``fmpq`` values of length `n` and sets
    all values to 0. This is equivalent to generating a ``fmpz`` vector
    of length `2n` with ``_fmpz_vec_init`` and setting all denominators
    to 1.

.. function:: void _fmpq_vec_clear(fmpq * vec, slong n)

    Frees an ``fmpq`` vector.


Randomisation
--------------------------------------------------------------------------------


.. function:: void _fmpq_vec_randtest(fmpq * f, flint_rand_t state, slong len, flint_bitcnt_t bits)

    Sets the entries of a vector of the given length to random rationals with
    numerator and denominator having up to the given number of bits per entry.

.. function:: void _fmpq_vec_randtest_uniq_sorted(fmpq * vec, flint_rand_t state, slong len, flint_bitcnt_t bits)

    Sets the entries of a vector of the given length to random distinct
    rationals with numerator and denominator having up to the given number
    of bits per entry. The entries in the vector are sorted.


Bit sizes and heights
--------------------------------------------------------------------------------

.. function:: void _fmpq_vec_max_height(fmpz_t height, const fmpq * vec, slong len);

    Computes the maximum of the height of any coefficient of ``(vec, len)``,
    each height being computed by :func:`fmpq_height`.

.. function:: flint_bitcnt_t _fmpq_vec_max_height_bits(const fmpq * vec, slong len);

    Computes the maximum number of bits of the height of any coefficient
    of ``(vec, len)``, each one being computed by :func:`fmpq_height_bits`.


Comparison
--------------------------------------------------------------------------------


.. function:: int _fmpq_vec_equal(const fmpq * vec1, const fmpq * vec2, slong len)

    Compares two vectors of the given length and returns `1` if they are
    equal, otherwise returns `0`.


Sorting
--------------------------------------------------------------------------------


.. function:: void _fmpq_vec_sort(fmpq * vec, slong len)

    Sorts the entries of ``(vec, len)``.


Conversions
--------------------------------------------------------------------------------


.. function:: void _fmpq_vec_set_fmpz_vec(fmpq * res, const fmpz * vec, slong len)

    Sets ``(res, len)`` to ``(vec, len)``.

.. function:: void _fmpq_vec_get_fmpz_vec_fmpz(fmpz * num, fmpz_t den, const fmpq * a, slong len)

    Find a common denominator ``den`` of the entries of ``a`` and set ``(num, len)`` to the corresponding numerators.


Dot product
--------------------------------------------------------------------------------


.. function:: void _fmpq_vec_dot(fmpq_t res, const fmpq * vec1, const fmpq * vec2, slong len)

    Sets ``res`` to the dot product of the vectors ``(vec1, len)`` and
    ``(vec2, len)``.


Input and output
--------------------------------------------------------------------------------


.. function:: int _fmpq_vec_fprint(FILE * file, const fmpq * vec, slong len)

    Prints the vector of given length to the stream ``file``. The
    format is the length followed by two spaces, then a space separated
    list of coefficients. If the length is zero, only `0` is printed.

    In case of success, returns a positive value. In case of failure,
    returns a non-positive value.

.. function:: int _fmpq_vec_print(const fmpq * vec, slong len)

    Prints the vector of given length to ``stdout``.

    For further details, see :func:`_fmpq_vec_fprint()`.
