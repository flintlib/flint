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

.macro	m5 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	mov	\scr1, \res_offset*8(\res)
	adcx	\scr2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	mulx	(3+\ap_offset)*8(\ap), \scr2, \r3
	adcx	\scr1, \r1
	adcx	\scr2, \r2
	mulx	(4+\ap_offset)*8(\ap), \scr1, \r4
	adcx	\scr1, \r3
	adcx	\zero, \r4
.endm

.macro	am5 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5, scr, zero
	mulx	(0+\ap_offset)*8(\ap), \scr, \r5
	adcx	\scr, \r0
	mov	\r0, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r5, \r1
	adox	\r0, \r1
	mulx	(2+\ap_offset)*8(\ap), \r0, \r5
	adcx	\scr, \r2
	adox	\r0, \r2
	mulx	(3+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r5, \r3
	adox	\r0, \r3
	mulx	(4+\ap_offset)*8(\ap), \r0, \r5
	adcx	\scr, \r4
	adox	\r0, \r4
	adcx	\zero, \r5
	adox	\zero, \r5
.endm

.global	FUNC(flint_mpn_mul_10_5)
.p2align	4, 0x90
TYPE(flint_mpn_mul_10_5)

FUNC(flint_mpn_mul_10_5):
	.cfi_startproc
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12

	xor	%r12d, %r12d

	m5	%rdi, 0, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12

	mov	1*8(%rsi), %rdx
	am5	%rdi, 1, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12
	mov	2*8(%rsi), %rdx
	am5	%rdi, 2, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %r8, %rbp, %r12
	mov	3*8(%rsi), %rdx
	am5	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %r8, %rax, %rbp, %r12
	mov	4*8(%rsi), %rdx
	am5	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %r8, %rax, %r9, %rbp, %r12
	mov	5*8(%rsi), %rdx
	am5	%rdi, 5, %rcx, 0, %r11, %rbx, %r8, %rax, %r9, %r10, %rbp, %r12
	mov	6*8(%rsi), %rdx
	am5	%rdi, 6, %rcx, 0, %rbx, %r8, %rax, %r9, %r10, %r11, %rbp, %r12
	mov	7*8(%rsi), %rdx
	am5	%rdi, 7, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12
	mov	8*8(%rsi), %rdx
	am5	%rdi, 8, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %r8, %rbp, %r12
	mov	9*8(%rsi), %rdx
	am5	%rdi, 9, %rcx, 0, %r9, %r10, %r11, %rbx, %r8, %rax, %rbp, %r12

	mov	%r10, 10*8(%rdi)
	mov	%r11, 11*8(%rdi)
	mov	%rbx, 12*8(%rdi)
	mov	%r8, 13*8(%rdi)
	mov	%rax, 14*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
.flint_mpn_mul_10_5_end:
SIZE(flint_mpn_mul_10_5, .flint_mpn_mul_10_5_end)
.cfi_endproc
