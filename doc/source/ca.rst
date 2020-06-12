.. ca:

**ca.h** -- exact real and complex numbers
===============================================================================

A :type:`ca_t` represents a real or complex number in a form suitable
for exact field arithmetic or comparison.
Exceptionally, a :type:`ca_t` may represent a special nonnumerical value,
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

Context objects
-------------------------------------------------------------------------------

.. type:: ca_ctx_struct

.. type:: ca_ctx_t

    A :type:`ca_ctx_t` context object holds a cache of fields *K* and
    constituent extension elements `a_k`.
    The field index in an individual :type:`ca_t` instance represents
    a shallow reference to the object defining the field *K* within the
    context object, so creating many elements of the same field is cheap.

    Since context objects are mutable (and may be mutated even when
    performing read-only operations on :type:`ca_t` instances), they must not
    be accessed simultaneously by different threads: in multithreaded
    environments, the user must use a separate context object for each thread.

.. function:: void ca_ctx_init(ca_ctx_t ctx)

    Initializes the context object *ctx* for use.
    Any evaluation options stored in the context object
    are set to default values.

.. function:: void ca_ctx_clear(ca_ctx_t ctx)

    Clears the context object *ctx*, freeing any memory allocated internally.
    This function should only be called after all :type:`ca_t` instances
    referring to this context have been cleared.

Memory management for numbers
-------------------------------------------------------------------------------

.. function:: void ca_init(ca_t x, ca_ctx_t ctx)

    Initializes the variable *x* for use, associating it with the
    context object *ctx*. The value of *x* is set to the rational number 0.

.. function:: void ca_clear(ca_t x, ca_ctx_t ctx)

    Clears the variable *x*.

.. function:: void ca_swap(ca_t x, ca_t y, ca_ctx_t ctx)

    Efficiently swaps the variables *x* and *y*.

Input and output
-------------------------------------------------------------------------------

.. function:: void ca_ctx_print(const ca_ctx_t ctx)

    Prints a description of the context *ctx* to standard output.
    This will give a complete listing of the cached fields in *ctx*.

.. function:: void ca_print(ca_t x, ca_ctx_t ctx)

    Prints a description of *x* to standard output.

Assignment and specific values
-------------------------------------------------------------------------------

.. function:: void ca_set(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to a copy of *x*.

.. function:: void ca_zero(ca_t x, ca_ctx_t ctx)
              void ca_one(ca_t x, ca_ctx_t ctx)

    Sets *x* to the integer 0 or 1. This creates a canonical representation
    of this number as an element of the trivial field `\mathbb{Q}`.

.. function:: void ca_set_si(ca_t x, slong v, ca_ctx_t ctx)
              void ca_set_ui(ca_t x, ulong v, ca_ctx_t ctx)
              void ca_set_fmpz(ca_t x, const fmpz_t v, ca_ctx_t ctx)
              void ca_set_fmpq(ca_t x, const fmpq_t v, ca_ctx_t ctx)

    Sets *x* to the integer or rational number *v*. This creates a canonical
    representation of this number as an element of the trivial field
    `\mathbb{Q}`.

.. function:: void ca_i(ca_t x, ca_ctx_t ctx)

    Sets *x* to the imaginary unit `i = \sqrt{-1}`. This creates a canonical
    representation of `i` as the generator of the algebraic number field
    `\mathbb{Q}(i)`.

.. function:: void ca_unknown(ca_t x, ca_ctx_t ctx)

    Sets *x* to the meta-value *Unknown*.

.. function:: void ca_undefined(ca_t x, ca_ctx_t ctx)

    Sets *x* to *Undefined*.

.. function:: void ca_uinf(ca_t x, ca_ctx_t ctx)

    Sets *x* to unsigned infinity `{\tilde \infty}`.

.. function:: void ca_pos_inf(ca_t x, ca_ctx_t ctx)
              void ca_neg_inf(ca_t x, ca_ctx_t ctx)
              void ca_pos_i_inf(ca_t x, ca_ctx_t ctx)
              void ca_neg_i_inf(ca_t x, ca_ctx_t ctx)

    Sets *x* to the signed infinity `+\infty`, `-\infty`, `+i \infty` or `-i \infty`.

Assignment of algebraic numbers
-------------------------------------------------------------------------------

.. function:: void ca_set_qqbar(ca_t res, const qqbar_t x, ca_ctx_t ctx)

    Sets *res* to the algebraic number *x*.

    If *x* is rational, *res* is set to the canonical representation as
    an element in the trivial field `\mathbb{Q}`.

    If *x* is irrational, this function always sets *res* to an element of
    a univariate number field `\mathbb{Q}(a)`. It will not, for example,
    identify `\sqrt{2} + \sqrt{3}`
    as an element of `\mathbb{Q}(\sqrt{2}, \sqrt{3})`. However, it may
    attempt to find a simpler number field than that generated by *x*
    itself. For example:

    * If *x* is quadratic, it will be expressed as an element of
      `\mathbb{Q}(\sqrt{N})` where *N* has no small repeated factors
      (obtained by performing a smooth factorisation of the discriminant).

    * TODO: if possible, coerce *x* to a low-degree cyclotomic field.


Representation properties
-------------------------------------------------------------------------------

The following predicates deal with the representation of a :type:`ca_t` and
hence can always be decided quickly. The return value is 0 for false
and 1 for true.

.. function:: int ca_is_unknown(const ca_t x, ca_ctx_t ctx)

    Returns 1 if *x* is Unknown, and 0 otherwise.

Value predicates
-------------------------------------------------------------------------------

The following predicates check a mathematical property which might
not be effectively decidable. The result is a :type:`truth_t` to allow
representing an unknown outcome.

.. function:: truth_t ca_check_is_number(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is a number. The result is ``T_TRUE`` is *x* is
    a field element (and hence a complex number), ``T_FALSE`` if *x* is
    an infinity or *Undefined*, and ``T_UNKNOWN`` if *x* is *Unknown*.

.. function:: truth_t ca_check_is_zero(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_one(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_neg_one(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_i(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_neg_i(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is equal to the number `0`, `1`, `-1`, `i`, or `-i`.

.. function:: truth_t ca_check_is_algebraic(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_rational(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_integer(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is respectively an algebraic number, a rational number,
    or an integer.

.. function:: truth_t ca_check_is_real(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is a real number. Warning: this returns ``T_FALSE`` if *x* is an
    infinity with real sign.

.. function:: truth_t ca_check_is_imaginary(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is an imaginary number. Warning: this returns ``T_FALSE`` if
    *x* is an infinity with imaginary sign.

.. function:: truth_t ca_check_is_nonreal(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is a nonreal complex number.

.. function:: truth_t ca_check_is_undefined(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is the special value *Undefined*.

.. function:: truth_t ca_check_is_infinity(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is any infinity (unsigned or signed).

.. function:: truth_t ca_check_is_uinf(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is unsigned infinity `{\tilde \infty}`.

.. function:: truth_t ca_check_is_signed_inf(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is any signed infinity.

.. function:: truth_t ca_check_is_pos_inf(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_neg_inf(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_pos_i_inf(const ca_t x, ca_ctx_t ctx)
              truth_t ca_check_is_neg_i_inf(const ca_t x, ca_ctx_t ctx)

    Tests if *x* is equal to the signed infinity `+\infty`, `-\infty`,
    `+i \infty`, `-i \infty`, respectively.

Comparisons
-------------------------------------------------------------------------------

.. function:: truth_t ca_check_equal(const ca_t x, const ca_t y, ca_ctx_t ctx)

    Tests `x = y`.
    The result is ``T_UNKNOWN`` if either operand is *Unknown*.
    The result may also be ``T_UNKNOWN`` if *x* and *y* are numerically
    indistinguishable and cannot be proved equal or unequal by
    an exact computation.

.. function:: truth_t ca_check_lt(const ca_t x, const ca_t y, ca_ctx_t ctx)
              truth_t ca_check_le(const ca_t x, const ca_t y, ca_ctx_t ctx)
              truth_t ca_check_gt(const ca_t x, const ca_t y, ca_ctx_t ctx)
              truth_t ca_check_ge(const ca_t x, const ca_t y, ca_ctx_t ctx)

    Compares *x* and *y*, implementing the respective operations
    `x < y`, `x \le y`, `x > y`, `x \ge y`.
    Only real numbers and `-\infty` and `+\infty` are considered comparable.
    The result is ``T_FALSE`` (not ``T_UNKNOWN``) if either operand is not
    comparable (being a nonreal complex number, unsigned infinity, or
    undefined).

Field structure operations
-------------------------------------------------------------------------------

.. function:: void ca_condense_field(ca_t res, ca_ctx_t ctx)

    Attempts to demote the value of *res* to a trivial subfield of its
    current field. In particular, this demotes any obviously rational
    value to the trivial field `\mathbb{Q}`.

    This function is applied automatically in most operations
    (arithmetic operations, etc.).

Arithmetic
-------------------------------------------------------------------------------

.. function:: void ca_neg(ca_t res, const ca_t x, ca_ctx_t ctx)

    Sets *res* to the negation of *x*.

.. function:: void ca_add_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
              void ca_add_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
              void ca_add_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
              void ca_add_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
              void ca_add(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Sets *res* to the sum of *x* and *y*. For special values, the following
    rules apply (`c \infty` denotes a signed infinity, `|c| = 1`):

    * `c \infty + d \infty = c \infty` if `c = d`

    * `c \infty + d \infty = \text{Undefined}` if `c \ne d`

    * `\tilde \infty + c \infty = \tilde \infty + \tilde \infty = \text{Undefined}`

    * `c \infty + z = c \infty` if `z \in \mathbb{C}`

    * `\tilde \infty + z = \tilde \infty` if `z \in \mathbb{C}`

    * `z + \text{Undefined} = \text{Undefined}` for any value *z* (including *Unknown*)

    In any other case, or if the correct case cannot be distinguished,
    the result is *Unknown*.

.. function:: void ca_sub_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
              void ca_sub_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
              void ca_sub_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
              void ca_sub_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
              void ca_sub(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Sets *res* to the difference of *x* and *y*.

.. function:: void ca_mul_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx)
              void ca_mul_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx)
              void ca_mul_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx)
              void ca_mul_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx)
              void ca_mul(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Sets *res* to the product of *x* and *y*. For special values, the following
    rules apply (`c \infty` denotes a signed infinity, `|c| = 1`):

    * `c \infty \cdot d \infty = c d \infty`

    * `c \infty \cdot \tilde \infty = \tilde \infty`

    * `\tilde \infty \cdot \tilde \infty = \tilde \infty`

    * `c \infty \cdot z = \operatorname{sgn}(z) c \infty` if `z \in \mathbb{C} \setminus \{0\}`

    * `c \infty \cdot 0 = \text{Undefined}`

    * `\tilde \infty \cdot 0 = \text{Undefined}`

    * `z \cdot  \text{Undefined} = \text{Undefined}` for any value *z* (including *Unknown*)

    In any other case, or if the correct case cannot be distinguished,
    the result is *Unknown*.

Internal representation
-------------------------------------------------------------------------------

.. macro:: CA_FMPQ(x)

.. macro:: CA_FMPQ_NUMREF(x)

.. macro:: CA_FMPQ_DENREF(x)

    Assuming that *x* holds an element of the trivial field `\mathbb{Q}`,
    this macro returns a pointer which can be used as an :type:`fmpq_t`,
    or respectively to the numerator or denominator as an :type:`fmpz_t`.

.. macro:: CA_MPOLY_Q(x)

    Assuming that *x* holds a generic field element as data,
    this macro returns a pointer which can be used as
    an :type:`fmpz_mpoly_q_t`.

.. macro:: CA_NF_ELEM(x)

    Assuming that *x* holds an Antic number field element as data,
    this macro returns a pointer which can be used as
    an :type:`nf_elem_t`.


.. function:: void _ca_make_field_element(ca_t x, slong new_index, ca_ctx_t ctx)

    Changes the internal representation of *x* to that of an element
    of the field with index *new_index* in the context object *ctx*.
    This may destroy the value of *x*.

.. function:: void _ca_make_fmpq(ca_t x, ca_ctx_t ctx)

    Changes the internal representation of *x* to that of an element of
    the trivial field `\mathbb{Q}`. This may destroy the value of *x*.

Extension field objects
-------------------------------------------------------------------------------

The following methods are intended for internal use.
The user should normally only manipulate :type:`ca_t` instances,
leaving the construction of field objects to the context object.

.. type:: ca_field_struct

.. type:: ca_field_t

    Represents a formal field.

.. function:: void ca_field_init_qq(ca_field_t K)

    Initializes *K* to represent the trivial field `\mathbb{Q}`.

.. function:: void ca_field_init_nf(ca_field_t K, const qqbar_t x)

    Initializes *K* to represent the algebraic number field `\mathbb{Q}(x)`.

.. function:: void ca_field_init_const(ca_field_t K, ulong func)

    Initializes *K* to represent the field
    `\mathbb{Q}(x)` where *x* is a builtin constant defined by
    *func* (example: *func* = *CA_Pi* for `x = \pi`).

.. function:: void ca_field_init_fx(ca_field_t K, ulong func, const ca_t x, ca_ctx_t ctx)

    Initializes *K* to represent the field
    `\mathbb{Q}(a)` where `a = f(x)`, given a number *x* and a builtin
    univariate function *func* (example: *func* = *CA_Exp* for `e^x`).

.. function:: void ca_field_init_multi(ca_field_t K, slong len)

    Initializes *K* to represent a multivariate field
    `\mathbb{Q}(a_1, \ldots, a_n)` in *n*
    extension elements. The extension elements must subsequently be
    assigned one by one using :func:`ca_field_set_ext`.

.. function:: void ca_field_set_ext(ca_field_t K, slong i, slong x_index, ca_ctx_t ctx)

    Sets the extension element at position *i* (here indexed from 0) of *K*
    to the generator of the field with index *x_index* in *ctx*.
    (It is assumed that the generating field is a univariate field.)

    This only inserts a shallow reference: the field at index *x_index* must
    be kept alive until *K* has been cleared.

.. function:: void ca_field_clear(ca_field_t K)

    Clears the field *K*.

.. function:: void ca_field_print(const ca_field_t K)

    Prints a description of the field *K* to standard output.

