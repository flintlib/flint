.. _introduction-calcium:

Introduction
===============================================================================

Exact numbers in Calcium
-----------------------------------------------------------------------

The core idea behind Calcium is to represent real and complex
numbers as elements of extension fields

.. math::

    \mathbb{Q}(a_1, \ldots, a_n)

of the rational numbers, where the extension numbers `a_k` are
described by symbolic expressions (which may depend on other fields
recursively).
The system constructs such fields automatically as needed
to represent the results of computations.
Any extension field is isomorphic to a formal field

.. math::

    \mathbb{Q}(a_1,\ldots,a_n) \;\; \cong \;\; K_{\text{formal}} := \operatorname{Frac}(\mathbb{Q}[X_1,\ldots,X_n] / I)

where *I* is the ideal of algebraic relations among the
extension numbers.
The relations may involve algebraic numbers
(for example: `i^2 + 1 = 0`), transcendental numbers
(for example: `e^{-\pi} \cdot e^{\pi} = 1`),
or combinations thereof.

Computation in the formal field depends (in general) on
multivariate polynomial arithmetic together with
use of a Gröbner basis for the ideal.
The map from the formal field to the true complex field
is maintained using arbitrary-precision ball arithmetic where
necessary.

As an important special case, Calcium can be used for
arithmetic in algebraic number fields (embedded explicitly in
`\mathbb{C}`)

.. math::

    \mathbb{Q}(a) \cong \mathbb{Q}[X] / \langle f(X) \rangle

with excellent performance thanks to internal use of the Antic library.

It will not always work perfectly: although
Calcium by design should never give a mathematically erroneous
answer, it may be unable to simplify a result as much as expected
and it may be unable to decide a predicate
(in which case it can return "Unknown").
Equality is at least decidable over the algebraic numbers
`\overline{\mathbb{Q}}` (for practical
degrees and bit sizes of the numbers!), and in certain
cases involving transcendentals.
We hope to improve Calcium's capabilities gradually
through enhancements to its built-in algorithms
and through customization options.

Usage details
.......................................................................

To understand how Calcium works more concretely, see
:ref:`examples-calcium` and the documentation for the main Calcium
number type (:type:`ca_t`):

* :ref:`ca`

Implementation details for
extension numbers and formal fields
can be found in the documentation of the corresponding modules:

* :ref:`ca-ext`
* :ref:`ca-field`

The following modules are used internally for arithmetic
in transcendental number fields (rational function fields)
`\mathbb{Q}(x_1,\ldots,x_n)` and over the field of algebraic
numbers `\overline{\mathbb{Q}}`, respectively. They may
be of independent interest:

* :ref:`fmpz-mpoly-q`
* :ref:`qqbar`


FAQ
-----------------------------------------------------------------------

**Isn't x = 0 undecidable?**

In general, yes: equality over the reals is undecidable.
In practice, much of calculus
and elementary number theory can be done with numbers that are
simple algebraic combinations of well-known elementary
and special functions, and there are heuristics that work quite well
for deciding predicates about such numbers.
Calcium will be able to give a definitive answer at least in 
simple cases (for example, proving 
`16 \operatorname{atan}(\tfrac{1}{5}) - 4 \operatorname{atan}(\tfrac{1}{239}) = \pi`
or `\sqrt{5+2\sqrt{6}} = \sqrt{2}+\sqrt{3}`),
and will simply answer "Unknown" when its heuristics are not powerful enough.

**How does Calcium compare to ordinary numerical computing?**

Calcium is far too slow to replace floating-point numbers
for 99.93% of scientific computing. The target is symbolic and
algebraic computation.
Nevertheless, Calcium may well be useful as a tool to test
and enhance the capabilities of numerical programs.

**How does Calcium compare to Arb arithmetic?**

The main advantage of Calcium over ball arithmetic alone is the ability
to do exact comparisons. The automatic precision management in Calcium
can also be convenient.

Calcium will usually be slower than Arb arithmetic.
If a computation is mostly numerical, it is probably better to try
using Arb first, and fall back on an exact calculation with Calcium
only if that fails because an exact comparison is needed.

**How does Calcium compare to symbolic computation systems (Mathematica, SymPy, etc.)?**

Calculating with constant values is only a small part of what
such systems have to do, but it is one of the most complex parts.
Existing computer algebra systems sometimes manage this very well,
and sometimes fail horribly. The most common problems are
1) getting numerical error bounds or branch cuts wrong, and
2) slowing down too much
when the expressions get large.
Calcium is intended to address both problems (through rigorous
numerical evaluation and use of fast polynomial arithmetic).

Ultimately, Calcium will no doubt handle some problems better
and others worse, and it should be considered a complement
to existing computer algebra systems rather than a replacement.
A symbolic expression simplifier may use Calcium evaluation
as one of its tools, but this probably needs to be done selectively
and in combination with many other heuristics.

**Why is Calcium written in C?**

The main advantage of developing Calcium as a C library is that it
will not be tied to a particular programming language
ecosystem: C is uniquely easy to interface from
almost any other language.
The second most important reason is familiarity: Calcium follows
the design of Flint and Arb
(coding style, naming, module layout, memory management,
test code, etc.) which has proved to work quite well for
libraries of this type.

There is also the performance argument. Some core functions will
benefit from
optimizations that are natural in C such as in-place operations
and fine-grained manual memory management. However, the performance
aspect should not be overemphasized: Calcium will spend
most of its time in Flint and Arb kernel functions
and this would probably still be true even if it were written
in a slower language.

There are certainly types of mathematical functionality that will be too
inconvenient to implement in C. Our intention is indeed to leave such
functionality to projects written in Python, Julia, etc. which may
then opt to depend on Calcium for basic operations.

**What is the development status of Calcium?**

Calcium is presently in early development and should be considered
experimental software.
The interfaces are subject to change and many important
functions and optimizations have not been implemented.
A more stable and functional release can be expected in late 2021.


.. raw:: latex

    \newpage

