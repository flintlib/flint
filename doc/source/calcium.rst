.. _calcium:

**calcium.h** -- global definitions
===============================================================================

Version
-------------------------------------------------------------------------------

.. function:: const char * calcium_version(void)

    Returns a pointer to the version of the library as a string ``X.Y.Z``.

Test code
-------------------------------------------------------------------------------

.. function:: double calcium_test_multiplier(void)

    Multiplier for the number of iterations to run in each unit test.
    The value can be changed by setting the environment variable
    ``CALCIUM_TEST_MULTIPLIER``. The default value is 1.0.

Triple-valued logic
-------------------------------------------------------------------------------

This library uses two kinds of predicate functions:

* Predicates with signature ``int foo_is_X(const foo_t x)`` 
  return the usual C boolean values ``1`` for true and  ``0`` for false,
  unless otherwise documented. Some functions may return ``0`` also when
  truth cannot be certified (this will be documented explicitly).

* Predicates with signature ``truth_t foo_check_is_X(const foo_t x)`` check a
  mathematical property that may not be decidable (or may be too costly to
  decide). The return value is a :type:`truth_t` (``T_TRUE``,
  ``T_FALSE`` or ``T_UNKNOWN``).

.. enum:: truth_t

    Represents one of the following truth values:

    .. macro:: T_TRUE

    .. macro:: T_FALSE

    .. macro:: T_UNKNOWN

    Warning: the constants ``T_TRUE`` and ``T_FALSE`` do not correspond to 1 and 0.
    It is erroneous to write, for example ``!t`` if ``t`` is a 
    :type:`truth_t`. One should instead write ``t != T_TRUE``, ``t == T_FALSE``,
    etc. depending on whether the unknown case should be included
    or excluded.

Flint, Arb and Antic types
-------------------------------------------------------------------------------

The following types from Flint, Arb and Antic are used throughout Calcium.
Although not included automatically by ``calcium.h``, we document them
here for convenience.

.. type:: slong

    Signed full-word integer (64 bits on a 64-bit system).

.. type:: ulong

    Unsigned full-word integer (64 bits on a 64-bit system).

.. type:: fmpz_t

    Flint integer.

.. type:: fmpq_t

    Flint rational number.

.. type:: fmpz_poly_t

    Flint dense univariate polynomial over the integers.

.. type:: fmpq_poly_t

    Flint dense univariate polynomial over the rational numbers.

.. type:: fmpz_mpoly_t

    Flint sparse multivariate integer polynomial.

.. type:: fmpz_mpoly_ctx_t

    Context for Flint sparse multivariate integer polynomial (defining the
    number of variables and monomial order).

.. type:: fmpz_mat_t

    Flint dense matrix over the integers.

.. type:: fmpq_mat_t

    Flint dense matrix over the rational numbers.

.. type:: arb_t

    Arb real number.

.. type:: acb_t

    Arb complex number.

.. type:: nf_t

    Antic number field.

.. type:: nf_elem_t

    Antic number field element.

Flint, Arb and Antic extras
-------------------------------------------------------------------------------

Here we collect various utility methods for Flint, Arb and Antic
types that are missing in those libraries. Some of these functions
may be migrated upstream in the future.

.. function:: ulong calcium_fmpz_hash(const fmpz_t x)

    Hash function for integers. The algorithm may change;
    presently, this simply extracts the low word (with sign).

Input and output
-------------------------------------------------------------------------------

.. type:: calcium_stream_struct

.. type:: calcium_stream_t

    A stream object which can hold either a file pointer or a
    string (with automatic resizing).

.. function:: void calcium_stream_init_file(calcium_stream_t out, FILE * fp)

    Initializes the stream *out* for writing to the file *fp*.
    The file can be *stdout*, *stderr*, or any file opened for writing
    by the user.

.. function:: void calcium_stream_init_str(calcium_stream_t out)

    Initializes the stream *out* for writing to a string in memory.
    When finished, the user should free the string (the *s* member
    of *out* with ``flint_free()``).

.. function:: void calcium_write(calcium_stream_t out, const char * s)

    Writes the string *s* to *out*.

.. function:: void calcium_write_free(calcium_stream_t out, char * s)

    Writes *s* to *out* and then frees *s* by calling ``flint_free()``.

.. function:: void calcium_write_si(calcium_stream_t out, slong x)

    Writes the integer *x* to *out*.

.. function:: void calcium_write_acb(calcium_stream_t out, const acb_t z, slong digits, ulong flags)

    Writes the Arb complex number *z* to *out*, showing *digits*
    digits and with the display style specified by *flags*
    (``ARB_STR_NO_RADIUS``, etc.).



.. raw:: latex

    \newpage

