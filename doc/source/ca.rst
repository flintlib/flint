.. ca:

**ca.h** -- exact real and complex numbers
===============================================================================

A :type:`ca_t` represents a real or complex number in a form suitable
for exact field arithmetic or comparison.
Exceptionally, a :type:`ca_t` may also hold a special nonnumerical value,
such as an infinity.

Introduction: numbers
-------------------------------------------------------------------------------

A *Calcium number* is a real or complex number
represented as an element of a formal field `K = \mathbb{Q}(a_1, \ldots, a_n)`
where the symbols `a_k` denote fixed algebraic or transcendental numbers.
For example, `e^{-2 \pi} - 3 i` may be represented as
`(1 - 3 a_2^2 a_1) / a_2^2` in the field `\mathbb{Q}(a_1,a_2)` with
`a_1 = i, a_2 = e^{\pi}`.
Extension elements `a_k` can be algebraic numbers represented
in canonical form by :type:`qqbar_t` instances (example: `i`),
symbolic constants (example: `\pi`),
or symbolic functions applied to :type:`ca_t` arguments
(example: `e^{\pi} = \exp(x)`, where `x = 1 \cdot \pi`
is a :type:`ca_t` representing an element of `\mathbb{Q}(\pi)`).
The user does not need to construct formal fields explicitly:
operations on Calcium numbers generate and cache fields automatically as needed
to express the results.

This representation is not canonical (in general). There can be
many representations for the same number, depending on
the choice of formal field *K*. Even within a fixed field *K*,
a number can have many representations if there are algebraic
relations between the extension elements.
Two numbers *x* and *y* can be tested for inequality using numerical
evaluation; to test for equality, it may be necessary to
eliminate redundant extension elements.
One of the central goals of Calcium will be to
implement heuristics for such elimination.

Introduction: special values
-------------------------------------------------------------------------------

In order to provide a closed arithmetic system and express limiting
cases of operators and special functions, a :type:`ca_t` can hold
any of the following special values besides ordinary numbers:

* *Unsigned infinity*, a formal object `{\tilde \infty}` representing
  the value of `1 / 0`. More generally, this is the
  value of meromorphic functions at poles.

* *Signed infinity*, a formal object `a \cdot \infty` where the sign `a`
  is a Calcium number with `|a| = 1`.
  The most common values are `+\infty, -\infty, +i \infty, -i \infty`.
  Signed infinities are used to denote directional limits and logarithmic
  singularities (for example, `\log(0) = -\infty`).

* *Undefined*, a formal object representing the value of indeterminate
  forms such as `0 / 0` and essential singularities such as 
  `\exp(\tilde \infty)`, where a number or infinity would not make sense
  as an answer.

* *Unknown*, a meta-value used to signal that the actual desired value
  could not be computed, either because Calcium does not (yet) have a data
  structure or algorithm for that case, or because doing so would be
  unreasonably expensive. This occurs, for example, if Calcium performs
  a division and is unable to decide whether the result is a number,
  unsigned infinity or undefined (because testing for zero fails).
  Wrappers may want to check output variables for
  *Unknown* and throw an exception (e.g. *NotImplementedError* in Python).

The distinction between *Calcium numbers* (which must represent
elements of `\mathbb{C}`) and the different kinds of nonnumerical values
(infinities, Undefined or Unknown) is essential. Nonnumerical values may
not be used as field extension elements `a_k`, and the denominator of
a formal field element must always represent a nonzero complex number.
Accordingly, for any given Calcium value *x* that is not *Unknown*,
it is exactly known whether *x* represents A) a number, B) unsigned infinity,
C) a signed infinity, or D) Undefined.


Number objects
-------------------------------------------------------------------------------

For all types, a *type_t* is defined as an array of length one of type
*type_struct*, permitting a *type_t* to be passed by reference.

.. type:: ca_struct

.. type:: ca_t

    A :type:`ca_t` contains an index to a field *K*, and
    data representing an element *x* of *K*.
    The data is either an inline rational number (:type:`fmpq_t`), an
    inline Antic number field element (:type:`nf_elem_t`)
    when *K* is an absolute algebraic number field `\mathbb{Q}(a)`,
    or a pointer to a heap-allocated :type:`fmpz_mpoly_q_t` representing
    an element of a generic field `\mathbb{Q}(a_1,\ldots,a_n)`.
    Special values are encoded using magic bits in the field index.

.. macro:: CA_FMPQ(x)

.. macro:: CA_FMPQ_NUMREF(x)

.. macro:: CA_FMPQ_DENREF(x)

    Assuming that *x* holds an element of the trivial field `\mathbb{Q}`,
    returns a pointer to this field which can be used as an :type:`fmpq_t`,
    or respectively to the numerator or denominator as an :type:`fmpz_t`.

.. macro:: CA_MPOLY_Q(x)

    Assuming that *x* holds a generic field element as data,
    this macro returns a pointer to this data which can be used as
    an :type:`fmpz_mpoly_q_t`.

.. macro:: CA_NF_ELEM(x)

    Assuming that *x* holds an Antic number field element as data,
    this macro returns a pointer to this data which can be used as
    an :type:`nf_elem_t`.

.. function:: ca_init(ca_t x, ca_ctx_t ctx)

    Initializes the variable *x* for use, associating it with the
    context object *ctx*. The value of *x* is set to the rational number 0.

.. function:: ca_clear(ca_t x, ca_ctx_t ctx)

    Clears the variable *x*.

Context objects
-------------------------------------------------------------------------------

.. type:: ca_ctx_struct

.. type:: ca_ctx_t

    A :type:`ca_ctx_t` context object holds a cache of fields *K* and
    constituent extension elements `a_k`.
    The field index in an individual :type:`ca_t` instance represents
    a shallow reference to the object defining the field *K* within the
    context object, so creating many elements of the same field is cheap.

    Using a :type:`ca_t` as an argument to an operation that creates a new
    field extension causes an immutable copy of that value to be placed in
    the context object.
    Since context objects are mutable, they must not be accessed
    simultaneously by different threads: in multithreaded environments, the
    user must use a separate context object for each thread.

.. function:: ca_ctx_init(ca_ctx_t ctx)

    Initializes the context object *ctx* for use.
    Any evaluation options stored in the context object
    are set to default values.

.. function:: ca_ctx_clear(ca_ctx_t ctx)

    Clears the context object *ctx*, freeing any memory allocated internally.
    This function should only be called after all :type:`ca_t` instances
    referring to this context have been cleared.

Extension and field objects
-------------------------------------------------------------------------------

.. type:: ca_extension_struct

.. type:: ca_extension_t

    Represents an extension element.

.. type:: ca_field_struct

.. type:: ca_field_t

    Represents a formal field.

.. function:: void ca_extension_init_qqbar(ca_extension_t ext, const qqbar_t x)

    Initializes *ext* to an extension element representing the
    algebraic number *x*.

.. function:: void ca_extension_init_const(ca_extension_t ext, ulong func)

    Initializes *ext* to an extension element representing the
    builtin constant *func* (example: *func* = *CA_Pi* for `\pi`).

.. function:: void ca_extension_init_fx(ca_extension_t ext, ulong func, const ca_t x)

    Initializes *ext* to an extension element representing the value `f(x)`
    where *func* is a builtin univariate function (example: *func* = *CA_Exp*
    for `e^x`).

.. function:: void ca_extension_clear(ca_extension_t ext)

    Clears the extension element *ext*.

.. function:: void ca_field_init_qq(ca_field_t K)

    Initializes *K* to represent the trivial field `\mathbb{Q}`.

.. function:: void ca_field_init_mpoly_q(ca_field_t K, slong n)

    Initializes *K* to represent a field `\mathbb{Q}(a_1, \ldots, a_n)` in *n*
    extension elements.

.. function:: void ca_field_set_ext(ca_field_t K, slong i, ca_extension_struct * ext)

    Sets the extension element at position *i* (here indexed from 0) of *K*
    to *ext*. This only inserts a shallow reference: the original object
    *ext* must be kept alive until *K* has been cleared.

.. function:: void ca_field_clear(ca_field_t K)

    Clears the field *K*.

