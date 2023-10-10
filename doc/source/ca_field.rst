.. _ca-field:

**ca_field.h** -- extension fields
===============================================================================

A :type:`ca_field_t` represents the parent field
`K = \mathbb{Q}(a_1,\ldots,a_n)` of a :type:`ca_t` element.
A :type:`ca_field_t` contains a list of pointers to
:type:`ca_ext_t` objects as well as a reduction ideal.

The user does not normally need to create :type:`ca_field_t` objects
manually: a :type:`ca_ctx_t` context object manages a cache of
fields automatically.

Internally, three types of field representation are used:

* The trivial field `\mathbb{Q}`.

* An Antic number field `\mathbb{Q}(a)` where *a* is defined by a :type:`qqbar_t`

* A generic field `\mathbb{Q}(a_1,\ldots,a_n)` where `n \ge 1`,
  and `a_1` is not defined by a :type:`qqbar_t` if `n = 1`.

The field type mainly affects the internal storage of the field
elements; the distinction is mostly transparent to the external interface.

Type and macros
-------------------------------------------------------------------------------

For all types, a *type_t* is defined as an array of length one of type
*type_struct*, permitting a *type_t* to be passed by reference.

.. type:: ca_field_struct

.. type:: ca_field_t

    Represents a formal field.

.. type:: ca_field_ptr

   Alias for ``ca_field_struct *``.

.. type:: ca_field_srcptr

   Alias for ``const ca_field_struct *``.

.. macro:: CA_FIELD_LENGTH(K)

    Accesses the number *n* of extension numbers of *K*. This is 0 if
    `K = \mathbb{Q}`.

.. macro:: CA_FIELD_EXT(K)

    Accesses the array of extension numbers as a :type:`ca_ext_ptr`.

.. macro:: CA_FIELD_EXT_ELEM(K, i)

    Accesses the extension number at position *i* (indexed from zero)
    as a :type:`ca_ext_t`.

.. macro:: CA_FIELD_HASH(K)

    Accesses the hash value of *K*.

.. macro:: CA_FIELD_IS_QQ(K)

    Returns whether *K* is the trivial field `\mathbb{Q}`.

.. macro:: CA_FIELD_IS_NF(K)

    Returns whether *K* represents an Antic number field
    `K = \mathbb{Q}(a)` where *a* is represented by a :type:`qqbar_t`.

.. macro:: CA_FIELD_IS_GENERIC(K)

    Returns whether *K* represents a generic field.

.. macro:: CA_FIELD_NF(K)

    Assuming that *K* represents an Antic number field `K = \mathbb{Q}(a)`,
    accesses the :type:`nf_t` object representing this field.

.. macro:: CA_FIELD_NF_QQBAR(K)

    Assuming that *K* represents an Antic number field `K = \mathbb{Q}(a)`,
    accesses the :type:`qqbar_t` object representing *a*.

.. macro:: CA_FIELD_IDEAL(K)

    Assuming that *K* represents a multivariate field, accesses the
    reduction ideal as a :type:`fmpz_mpoly_t` array.

.. macro:: CA_FIELD_IDEAL_ELEM(K, i)

    Assuming that *K* represents a multivariate field, accesses element *i*
    (indexed from zero) of the reduction ideal as a ``::fmpz_mpoly_t``.

.. macro:: CA_FIELD_IDEAL_LENGTH(K)

    Assuming that *K* represents a multivariate field, accesses the number
    of polynomials in the reduction ideal.

.. doxygendefine:: CA_FIELD_MCTX

Memory management
-------------------------------------------------------------------------------

.. doxygenfunction:: ca_field_init_qq

.. doxygenfunction:: ca_field_init_nf

.. doxygenfunction:: ca_field_init_const

.. doxygenfunction:: ca_field_init_fx

.. doxygenfunction:: ca_field_init_fxy

.. doxygenfunction:: ca_field_init_multi

.. doxygenfunction:: ca_field_set_ext

.. doxygenfunction:: ca_field_clear

Input and output
-------------------------------------------------------------------------------

.. doxygenfunction:: ca_field_print

Ideal
-------------------------------------------------------------------------------

.. doxygenfunction:: ca_field_build_ideal

.. doxygenfunction:: ca_field_build_ideal_erf


Structure operations
-------------------------------------------------------------------------------

.. doxygenfunction:: ca_field_cmp

Cache
-------------------------------------------------------------------------------

.. doxygenstruct:: ca_field_cache_struct

.. doxygentypedef:: ca_field_cache_t

.. doxygenfunction:: ca_field_cache_init

.. doxygenfunction:: ca_field_cache_clear

.. doxygenfunction:: ca_field_cache_insert_ext

.. raw:: latex

    \newpage

.. doxygenstruct:: ca_field_struct
    :members:
