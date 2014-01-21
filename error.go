// Copyright 2014, Hǎiliàng Wáng. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package gsl

/*
#include <gsl/gsl_errno.h>
#cgo LDFLAGS: -lgsl -lgslcblas
*/
import "C"

func init() {
	C.gsl_set_error_handler_off()
}

type Error C.int

func toError(e C.int) error {
	if e != 0 {
		return Error(e)
	}
	return nil
}

func (e Error) Error() string {
	return C.GoString(C.gsl_strerror(C.int(e)))
}
