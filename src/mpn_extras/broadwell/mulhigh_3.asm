#
#   Copyright (C) 2024 Albin Ahlbäck
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 2.1 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#

include(`config.m4')dnl
dnl
.text

.global	FUNC(flint_mpn_mulhigh_3)
.p2align	4, 0x90
TYPE(flint_mpn_mulhigh_3)

FUNC(flint_mpn_mulhigh_3):
	.cfi_startproc
	push	%rbx

	mov	%rdx, %rcx
	mov	0*8(%rdx), %rdx
	xor	%r9d, %r9d

	mulx	1*8(%rsi), %r8, %rax
	mulx	2*8(%rsi), %r8, %r10
	adcx	%r8, %rax
	adcx	%r9, %r10

	mov	1*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rbx
	mulx	1*8(%rsi), %r8, %r11
	adcx	%rbx, %rax
	adox	%r8, %rax
	adcx	%r11, %r10
	mulx	2*8(%rsi), %r8, %r11
	adox	%r8, %r10
	adcx	%r9, %r11
	adox	%r9, %r11

	mov	2*8(%rcx), %rdx
	mulx	0*8(%rsi), %r8, %rbx
	adcx	%r8, %rax
	adcx	%rbx, %r10
	mulx	1*8(%rsi), %r8, %rbx
	adox	%r8, %r10
	adox	%rbx, %r11
	mulx	2*8(%rsi), %r8, %rbx
	adcx	%r8, %r11
	adcx	%r9, %rbx
	adox	%r9, %rbx

	mov	%r10, 0*8(%rdi)
	mov	%r11, 1*8(%rdi)
	mov	%rbx, 2*8(%rdi)

	pop	%rbx

	ret
.flint_mpn_mulhigh_3_end:
SIZE(flint_mpn_mulhigh_3, .flint_mpn_mulhigh_3_end)
.cfi_endproc
