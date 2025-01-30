.. _calcium:

**calcium.h** -- global definitions
===============================================================================

Version
-------------------------------------------------------------------------------

.. function:: const char * calcium_version(void)

    Returns a pointer to the version of the library as a string ``X.Y.Z``.

Triple-valued logic
-------------------------------------------------------------------------------

The Calcium modules use two kinds of predicate functions:

* Predicates with signature ``int foo_is_X(const foo_t x)`` 
  return the usual C boolean values ``1`` for true and  ``0`` for false,
  unless otherwise documented. Some functions may return ``0`` also when
  truth cannot be certified (this will be documented explicitly).

* Predicates with signature ``truth_t foo_check_is_X(const foo_t x)`` check a
  mathematical property that may not be decidable (or may be too costly to
  decide). The return value is a :type:`truth_t` (``T_TRUE``,
  ``T_FALSE`` or ``T_UNKNOWN``).

FLINT, Arb and Antic extras
-------------------------------------------------------------------------------

Here we collect various utility methods for FLINT, Arb and Antic
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
              void calcium_write_fmpz(calcium_stream_t out, const fmpz_t x)

    Writes the integer *x* to *out*.

.. function:: void calcium_write_arb(calcium_stream_t out, const arb_t z, slong digits, ulong flags)
              void calcium_write_acb(calcium_stream_t out, const acb_t z, slong digits, ulong flags)

    Writes the Arb number *z* to *out*, showing *digits*
    digits and with the display style specified by *flags*
    (``ARB_STR_NO_RADIUS``, etc.).



.. raw:: latex

    \newpage

