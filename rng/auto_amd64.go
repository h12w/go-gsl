package rng

/*
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#cgo LDFLAGS: -lgsl -lgslcblas
*/
import "C"

import (
	"unsafe"
)
