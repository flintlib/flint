.. ca-field:

**ca_field.h** -- extension fields
===============================================================================

A :type:`ca_field_t` represents the parent field
`K = \mathbb{Q}(a_1,\ldots,a_n)` of a :type:`ca_t` element.
A :type:`ca_field_t` contains a list of pointers to
:type:`ca_ext_t` objects as well as a reduction ideal.

The user does not normally need to create :type:`ca_field_t` objects
manually: a :type:`ca_ctx_t` context object manages a cache of
fields automatically.

Type and macros
-------------------------------------------------------------------------------

For all types, a *type_t* is defined as an array of length one of type
*type_struct*, permitting a *type_t* to be passed by reference.

.. type:: ca_field_struct

.. type:: ca_field_t

    Represents a formal field.

.. macro:: CA_FIELD_LENGTH(K)

    Accesses the number of extension elements *n* of *K*. This is 0 if
    `K = \mathbb{Q}`.

.. macro:: CA_FIELD_EXT(K)

    Accesses the array of extension elements as a :type:`ca_ext_ptr`.

.. macro:: CA_FIELD_EXT_ELEM(K, i)

    Accesses the extension element at position *i* (indexed from zero)
    as a :type:`ca_ext_t`.

.. macro:: CA_FIELD_HASH(K)

    Accesses the hash value of *K*.

.. macro:: CA_FIELD_IS_NF(K)

    Returns whether *K* represents an Antic number field
    `K = \mathbb{Q}(a)` where *a* is represented by a :type:`qqbar_t`.

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
    (indexed from zero) of the reduction ideal as a :type:`fmpz_mpoly_t`.

.. macro:: CA_FIELD_IDEAL_LENGTH(K)

    Assuming that *K* represents a multivariate field, accesses the number
    of polynomials in the reduction ideal.

.. macro:: CA_FIELD_MCTX(K, ctx)

    Assuming that *K* represents a multivariate field, accesses the
    :type:`fmpz_mpoly_ctx_t` context object for multivariate polynomial
    arithmetic on the internal representation of elements in this field.

Memory management
-------------------------------------------------------------------------------

.. function:: void ca_field_init_qq(ca_field_t K, ca_ctx_t ctx)

    Initializes *K* to represent the trivial field `\mathbb{Q}`.

.. function:: void ca_field_init_nf(ca_field_t K, const qqbar_t x, ca_ctx_t ctx)

    Initializes *K* to represent the algebraic number field `\mathbb{Q}(x)`.

.. function:: void ca_field_init_const(ca_field_t K, ulong func, ca_ctx_t ctx)

    Initializes *K* to represent the field
    `\mathbb{Q}(x)` where *x* is a builtin constant defined by
    *func* (example: *func* = *CA_Pi* for `x = \pi`).

.. function:: void ca_field_init_fx(ca_field_t K, ulong func, const ca_t x, ca_ctx_t ctx)

    Initializes *K* to represent the field
    `\mathbb{Q}(a)` where `a = f(x)`, given a number *x* and a builtin
    univariate function *func* (example: *func* = *CA_Exp* for `e^x`).

.. function:: void ca_field_init_multi(ca_field_t K, slong len, ca_ctx_t ctx)

    Initializes *K* to represent a multivariate field
    `\mathbb{Q}(a_1, \ldots, a_n)` in *n*
    extension numbers. The extension numbers must subsequently be
    assigned one by one using :func:`ca_field_set_ext`.

.. function:: void ca_field_set_ext(ca_field_t K, slong i, slong x_index, ca_ctx_t ctx)

    Sets the extension number at position *i* (here indexed from 0) of *K*
    to the generator of the field with index *x_index* in *ctx*.
    (It is assumed that the generating field is a univariate field.)

    This only inserts a shallow reference: the field at index *x_index* must
    be kept alive until *K* has been cleared.

.. function:: void ca_field_clear(ca_field_t K, ca_ctx_t ctx)

    Clears the field *K*. This does not clear the individual extension
    elements, which are only held as references.

Input and output
-------------------------------------------------------------------------------

.. function:: void ca_field_print(const ca_field_t K, const ca_ctx_t ctx)

    Prints a description of the field *K* to standard output.

Structure operations
-------------------------------------------------------------------------------

.. function:: int ca_field_cmp(const ca_field_t K1, const ca_field_t K2, ca_ctx_t ctx)

    Compares the field objects *K1* and *K2* in a canonical sort order,
    returning -1, 0 or 1. This only performs a lexicographic comparison
    of the representations of *K1* and *K2*; the return value does not say
    anything meaningful about the relative structures of *K1* and *K2*
    as mathematical fields.


.. raw:: latex

    \newpage

