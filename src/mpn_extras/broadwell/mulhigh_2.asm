#
#   Copyright (C) 2024 Albin Ahlbäck
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

.global	FUNC(flint_mpn_mulhigh_2)
.p2align	4, 0x90
TYPE(flint_mpn_mulhigh_2)

FUNC(flint_mpn_mulhigh_2):
	.cfi_startproc
	mov	1*8(%rdx), %r11
	mov	0*8(%rdx), %rdx
	xor	%r10d, %r10d
	mulx	0*8(%rsi), %r9, %rax
	mulx	1*8(%rsi), %r9, %rcx
	adcx	%r9, %rax
	adcx	%r10, %rcx
	mov	%r11, %rdx
	mulx	0*8(%rsi), %r9, %r8
	adcx	%r9, %rax
	adcx	%r8, %rcx
	mulx	1*8(%rsi), %r9, %r8
	adox	%r9, %rcx
	adox	%r10, %r8
	adcx	%r10, %r8
	mov	%rcx, 0*8(%rdi)
	mov	%r8, 1*8(%rdi)

	ret
.flint_mpn_mulhigh_2_end:
SIZE(flint_mpn_mulhigh_2, .flint_mpn_mulhigh_2_end)
.cfi_endproc
