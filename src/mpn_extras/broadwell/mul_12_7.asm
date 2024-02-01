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

.macro	m7 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5, r6, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	mov	\scr1, \res_offset*8(\res)
	adcx	\scr2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	mulx	(3+\ap_offset)*8(\ap), \scr2, \r3
	adcx	\scr1, \r1
	adcx	\scr2, \r2
	mulx	(4+\ap_offset)*8(\ap), \scr1, \r4
	mulx	(5+\ap_offset)*8(\ap), \scr2, \r5
	adcx	\scr1, \r3
	adcx	\scr2, \r4
	mulx	(6+\ap_offset)*8(\ap), \scr1, \r6
	adcx	\scr1, \r5
	adcx	\zero, \r6
.endm

.macro	am7 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5, r6, r7, scr, zero
	mulx	(0+\ap_offset)*8(\ap), \scr, \r7
	adcx	\scr, \r0
	mov	\r0, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r7, \r1
	adox	\r0, \r1
	mulx	(2+\ap_offset)*8(\ap), \r0, \r7
	adcx	\scr, \r2
	adox	\r0, \r2
	mulx	(3+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r7, \r3
	adox	\r0, \r3
	mulx	(4+\ap_offset)*8(\ap), \r0, \r7
	adcx	\scr, \r4
	adox	\r0, \r4
	mulx	(5+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r7, \r5
	adox	\r0, \r5
	mulx	(6+\ap_offset)*8(\ap), \r0, \r7
	adcx	\scr, \r6
	adox	\r0, \r6
	adcx	\zero, \r7
	adox	\zero, \r7
.endm

.global	FUNC(flint_mpn_mul_12_7)
.p2align	4, 0x90
TYPE(flint_mpn_mul_12_7)

FUNC(flint_mpn_mul_12_7):
	.cfi_startproc
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14

	mov	1*8(%rsi), %rdx
	am7	%rdi, 1, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14
	mov	2*8(%rsi), %rdx
	am7	%rdi, 2, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r8, %r13, %r14
	mov	3*8(%rsi), %rdx
	am7	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rbp, %r12, %r8, %rax, %r13, %r14
	mov	4*8(%rsi), %rdx
	am7	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rbp, %r12, %r8, %rax, %r9, %r13, %r14
	mov	5*8(%rsi), %rdx
	am7	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %r12, %r8, %rax, %r9, %r10, %r13, %r14
	mov	6*8(%rsi), %rdx
	am7	%rdi, 6, %rcx, 0, %rbx, %rbp, %r12, %r8, %rax, %r9, %r10, %r11, %r13, %r14
	mov	7*8(%rsi), %rdx
	am7	%rdi, 7, %rcx, 0, %rbp, %r12, %r8, %rax, %r9, %r10, %r11, %rbx, %r13, %r14
	mov	8*8(%rsi), %rdx
	am7	%rdi, 8, %rcx, 0, %r12, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r13, %r14
	mov	9*8(%rsi), %rdx
	am7	%rdi, 9, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14
	mov	10*8(%rsi), %rdx
	am7	%rdi, 10, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r8, %r13, %r14
	mov	11*8(%rsi), %rdx
	am7	%rdi, 11, %rcx, 0, %r9, %r10, %r11, %rbx, %rbp, %r12, %r8, %rax, %r13, %r14

	mov	%r10, 12*8(%rdi)
	mov	%r11, 13*8(%rdi)
	mov	%rbx, 14*8(%rdi)
	mov	%rbp, 15*8(%rdi)
	mov	%r12, 16*8(%rdi)
	mov	%r8, 17*8(%rdi)
	mov	%rax, 18*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
.flint_mpn_mul_12_7_end:
SIZE(flint_mpn_mul_12_7, .flint_mpn_mul_12_7_end)
.cfi_endproc
