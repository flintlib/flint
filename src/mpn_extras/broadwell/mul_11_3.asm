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

.macro	m3 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	mov	\scr1, \res_offset*8(\res)
	adcx	\scr2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	adcx	\scr1, \r1
	adcx	\zero, \r2
.endm

.macro	am3 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, scr, zero
	mulx	(0+\ap_offset)*8(\ap), \scr, \r3
	adcx	\scr, \r0
	mov	\r0, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r3, \r1
	adox	\r0, \r1
	mulx	(2+\ap_offset)*8(\ap), \r0, \r3
	adcx	\scr, \r2
	adox	\r0, \r2
	adcx	\zero, \r3
	adox	\zero, \r3
.endm

.global	FUNC(flint_mpn_mul_11_3)
.p2align	4, 0x90
TYPE(flint_mpn_mul_11_3)

FUNC(flint_mpn_mul_11_3):
	.cfi_startproc
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx

	xor	%ebx, %ebx

	m3	%rdi, 0, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx

	mov	1*8(%rsi), %rdx
	am3	%rdi, 1, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx
	mov	2*8(%rsi), %rdx
	am3	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %rax, %r11, %rbx
	mov	3*8(%rsi), %rdx
	am3	%rdi, 3, %rcx, 0, %r9, %r10, %rax, %r8, %r11, %rbx
	mov	4*8(%rsi), %rdx
	am3	%rdi, 4, %rcx, 0, %r10, %rax, %r8, %r9, %r11, %rbx
	mov	5*8(%rsi), %rdx
	am3	%rdi, 5, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx
	mov	6*8(%rsi), %rdx
	am3	%rdi, 6, %rcx, 0, %r8, %r9, %r10, %rax, %r11, %rbx
	mov	7*8(%rsi), %rdx
	am3	%rdi, 7, %rcx, 0, %r9, %r10, %rax, %r8, %r11, %rbx
	mov	8*8(%rsi), %rdx
	am3	%rdi, 8, %rcx, 0, %r10, %rax, %r8, %r9, %r11, %rbx
	mov	9*8(%rsi), %rdx
	am3	%rdi, 9, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx
	mov	10*8(%rsi), %rdx
	am3	%rdi, 10, %rcx, 0, %r8, %r9, %r10, %rax, %r11, %rbx

	mov	%r9, 11*8(%rdi)
	mov	%r10, 12*8(%rdi)
	mov	%rax, 13*8(%rdi)

	pop	%rbx

	ret
.flint_mpn_mul_11_3_end:
SIZE(flint_mpn_mul_11_3, .flint_mpn_mul_11_3_end)
.cfi_endproc
