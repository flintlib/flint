---
title: Reporting bugs
summary: Where to report bugs in FLINT and what information to include.
order: 2
---

The maintainers want to know about all bugs: wrong results, crashes,
memory leaks, build failures, poor performance, and errors in the
documentation. Please open an issue at the GitHub repository,
https://github.com/flintlib/flint/issues, or write to the `flint-devel`
mailing list, https://groups.google.com/g/flint-devel.

A good bug report contains:

- the FLINT version (`flint_version`, or the release or git commit), the
  versions of GMP and MPFR, the compiler and its version, and the
  operating system and processor;
- how FLINT was configured and built (the `configure` command line, or
  the package that was installed);
- a minimal, complete C program that shows the problem, with its output
  and the expected output. Programs that need only FLINT's public
  functions and run in a fraction of a second are easiest to work with.

If a computation gives a wrong result, it helps to check whether the bug
is also present in a build with assertions enabled (`--enable-assert`):
the failing assertion often locates the problem, or shows that a
precondition of a function was violated. Crashes and memory errors are
best investigated with Valgrind or the address sanitizer (`--enable-asan`),
calling `flint_cleanup_master` at the end of the program so that FLINT's
caches are not reported as leaks.

FLINT requires GMP 6.2.1 or later and MPFR 4.1.0 or later; problems with
older versions are not bugs.
