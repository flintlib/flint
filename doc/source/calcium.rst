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

