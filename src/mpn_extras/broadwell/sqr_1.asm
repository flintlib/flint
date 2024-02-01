#
#   Copyright (C) 2023 Albin Ahlbäck
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#

include(`config.m4')dnl
dnl
.text

.global	FUNC(flint_mpn_sqr_1)
.p2align	4, 0x90
TYPE(flint_mpn_sqr_1)

FUNC(flint_mpn_sqr_1):
	.cfi_startproc
	mov	0*8(%rsi), %rdx
	mulx	%rdx, %rcx, %rax
	mov	%rcx, 0*8(%rdi)
	mov	%rax, 1*8(%rdi)

	ret
.flint_mpn_sqr_1_end:
SIZE(flint_mpn_sqr_1, .flint_mpn_sqr_1_end)
.cfi_endproc
