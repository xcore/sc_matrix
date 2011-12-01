/*
 * Matrix op worker threads
 *
 * Copyright (C) 2011 Steve Kerrison <github@stevekerrison.com>
 *
 * This software is freely distributable under a derivative of the
 * University of Illinois/NCSA Open Source License posted in
 * LICENSE.txt and at <http://github.xcore.com/>
 */

#ifndef MATRIX_WORKER_H_
#define MATRIX_WORKER_H_

#include "matrix.h"
	
void matrix_mul_worker(int ptA, int ptDimA, int ptB, int ptDimB, int ptC,
	short startC, short lenC, int ptOps);
	
/* Piecewise (array) ops */
	
void matrix_arr_worker_add(int ptA, int ptB, int ptC,
	int ptOps, short offset, short len);

void matrix_arr_worker_sub(int ptA, int ptB, int ptC,
	int ptOps, short offset, short len);

void matrix_arr_worker_mul(int ptA, int ptB, int ptC,
	int ptOps, short offset, short len);

void matrix_arr_worker_div(int ptA, int ptB, int ptC,
	int ptOps, short offset, short len);

void matrix_arr_worker_udiv(int ptA, int ptB, int ptC,
	int ptOps, short offset, short len);
	
/* Scalar ops */
	
void matrix_sca_worker_add(int ptA, int S, int ptC,
	int ptOps, short offset, short len);

void matrix_sca_worker_sub(int ptA, int S, int ptC,
	int ptOps, short offset, short len);

void matrix_sca_worker_mul(int ptA, int S, int ptC,
	int ptOps, short offset, short len);

void matrix_sca_worker_div(int ptA, int S, int ptC,
	int ptOps, short offset, short len);

void matrix_sca_worker_udiv(int ptA, int S, int ptC,
	int ptOps, short offset, short len);

void matrix_sca_worker_rand(int ptC, int ptOps, short offset, short len);

void matrix_sca_worker_set(int S, int ptC, int ptOps, short offset, short len);

void matrix_sca_worker_shr(int ptA, int S, int ptC, int ptOps, short offset, short len);

void matrix_sca_worker_ashr(int ptA, int S, int ptC, int ptOps, short offset, short len);

void matrix_sca_worker_shl(int ptA, int S, int ptC, int ptOps, short offset, short len);

#endif /* MATRIX_WORKER_H_ */
