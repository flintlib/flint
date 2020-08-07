.. calcium:

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

.. type:: arb_t

    Arb real number.

.. type:: acb_t

    Arb complex number.

.. type:: nf_t

    Antic number field.

.. type:: nf_elem_t

    Antic number field element.


.. raw:: latex

    \newpage

