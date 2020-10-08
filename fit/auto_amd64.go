package fit

/*
#include <gsl/gsl_fit.h>
#include <stdlib.h>
#cgo LDFLAGS: -lgsl -lgslcblas
*/
import "C"

import (
	"unsafe"
)

// gsl_fit_linear
func Linear(x []float64, xstride int, y []float64, ystride int, n int) (c0 float64, c1 float64, cov00 float64, cov01 float64, cov11 float64, sumsq float64, ret int32) {
	_x := (*C.double)(nil)
	if len(x) > 0 {
		_x = (*C.double)(unsafe.Pointer(&x[0]))
	}
	_xstride := C.size_t(xstride)
	_y := (*C.double)(nil)
	if len(y) > 0 {
		_y = (*C.double)(unsafe.Pointer(&y[0]))
	}
	_ystride := C.size_t(ystride)
	_n := C.size_t(n)
	_c0 := (*C.double)(unsafe.Pointer(&c0))
	_c1 := (*C.double)(unsafe.Pointer(&c1))
	_cov00 := (*C.double)(unsafe.Pointer(&cov00))
	_cov01 := (*C.double)(unsafe.Pointer(&cov01))
	_cov11 := (*C.double)(unsafe.Pointer(&cov11))
	_sumsq := (*C.double)(unsafe.Pointer(&sumsq))
	_ret := C.gsl_fit_linear(_x, _xstride, _y, _ystride, _n, _c0, _c1, _cov00, _cov01, _cov11, _sumsq)
	ret = int32(_ret)
	return
}

// gsl_fit_linear_est
func LinearEst(x float64, c0 float64, c1 float64, cov00 float64, cov01 float64, cov11 float64) (y float64, yErr float64, ret int32) {
	_x := C.double(x)
	_c0 := C.double(c0)
	_c1 := C.double(c1)
	_cov00 := C.double(cov00)
	_cov01 := C.double(cov01)
	_cov11 := C.double(cov11)
	_y := (*C.double)(unsafe.Pointer(&y))
	_yErr := (*C.double)(unsafe.Pointer(&yErr))
	_ret := C.gsl_fit_linear_est(_x, _c0, _c1, _cov00, _cov01, _cov11, _y, _yErr)
	ret = int32(_ret)
	return
}

// gsl_fit_mul
func Mul(x []float64, xstride int, y []float64, ystride int, n int) (c1 float64, cov11 float64, sumsq float64, ret int32) {
	_x := (*C.double)(nil)
	if len(x) > 0 {
		_x = (*C.double)(unsafe.Pointer(&x[0]))
	}
	_xstride := C.size_t(xstride)
	_y := (*C.double)(nil)
	if len(y) > 0 {
		_y = (*C.double)(unsafe.Pointer(&y[0]))
	}
	_ystride := C.size_t(ystride)
	_n := C.size_t(n)
	_c1 := (*C.double)(unsafe.Pointer(&c1))
	_cov11 := (*C.double)(unsafe.Pointer(&cov11))
	_sumsq := (*C.double)(unsafe.Pointer(&sumsq))
	_ret := C.gsl_fit_mul(_x, _xstride, _y, _ystride, _n, _c1, _cov11, _sumsq)
	ret = int32(_ret)
	return
}

// gsl_fit_mul_est
func MulEst(x float64, c1 float64, cov11 float64) (y float64, yErr float64, ret int32) {
	_x := C.double(x)
	_c1 := C.double(c1)
	_cov11 := C.double(cov11)
	_y := (*C.double)(unsafe.Pointer(&y))
	_yErr := (*C.double)(unsafe.Pointer(&yErr))
	_ret := C.gsl_fit_mul_est(_x, _c1, _cov11, _y, _yErr)
	ret = int32(_ret)
	return
}

// gsl_fit_wlinear
func Wlinear(x []float64, xstride int, w []float64, wstride int, y []float64, ystride int, n int) (c0 float64, c1 float64, cov00 float64, cov01 float64, cov11 float64, chisq float64, ret int32) {
	_x := (*C.double)(nil)
	if len(x) > 0 {
		_x = (*C.double)(unsafe.Pointer(&x[0]))
	}
	_xstride := C.size_t(xstride)
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_y := (*C.double)(nil)
	if len(y) > 0 {
		_y = (*C.double)(unsafe.Pointer(&y[0]))
	}
	_ystride := C.size_t(ystride)
	_n := C.size_t(n)
	_c0 := (*C.double)(unsafe.Pointer(&c0))
	_c1 := (*C.double)(unsafe.Pointer(&c1))
	_cov00 := (*C.double)(unsafe.Pointer(&cov00))
	_cov01 := (*C.double)(unsafe.Pointer(&cov01))
	_cov11 := (*C.double)(unsafe.Pointer(&cov11))
	_chisq := (*C.double)(unsafe.Pointer(&chisq))
	_ret := C.gsl_fit_wlinear(_x, _xstride, _w, _wstride, _y, _ystride, _n, _c0, _c1, _cov00, _cov01, _cov11, _chisq)
	ret = int32(_ret)
	return
}

// gsl_fit_wmul
func Wmul(x []float64, xstride int, w []float64, wstride int, y []float64, ystride int, n int) (c1 float64, cov11 float64, sumsq float64, ret int32) {
	_x := (*C.double)(nil)
	if len(x) > 0 {
		_x = (*C.double)(unsafe.Pointer(&x[0]))
	}
	_xstride := C.size_t(xstride)
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_y := (*C.double)(nil)
	if len(y) > 0 {
		_y = (*C.double)(unsafe.Pointer(&y[0]))
	}
	_ystride := C.size_t(ystride)
	_n := C.size_t(n)
	_c1 := (*C.double)(unsafe.Pointer(&c1))
	_cov11 := (*C.double)(unsafe.Pointer(&cov11))
	_sumsq := (*C.double)(unsafe.Pointer(&sumsq))
	_ret := C.gsl_fit_wmul(_x, _xstride, _w, _wstride, _y, _ystride, _n, _c1, _cov11, _sumsq)
	ret = int32(_ret)
	return
}
