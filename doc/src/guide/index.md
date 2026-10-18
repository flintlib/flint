---
title: User guide
summary: Installing and using FLINT: building, conventions, correctness and its limits, choosing types and generic programming, ball arithmetic, exact numbers, interfaces, and how to contribute.
---

The user guide explains how to install FLINT, how to compile programs
against it, and the conventions that all modules share: the word types,
the naming of functions, aliasing, memory management, threads and
portability. The module documentation assumes these conventions.

The guide then explains what FLINT guarantees and what it does not: which
algorithms are proved correct, how FLINT is tested, which inputs are
checked, the limits on sizes, and why FLINT should not be used for
cryptography or on untrusted input. It helps choose among FLINT's many
types for a task and describes the different ways in which FLINT's code
is generic, from templates to the `gr` generic rings. It introduces
rigorous ball arithmetic, which has its own way of working, gives a short
overview of the exact real and complex numbers, describes the interfaces
to other languages, and explains how to report bugs and contribute.

New users should start with [](guide/getting-started/introduction) and
[](guide/getting-started/building), and then follow the first tutorial,
[](tutorials/getting-started/first-program). The release notes of all
versions are in [](history).
