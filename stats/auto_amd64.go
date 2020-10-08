package stats

/*
#include <gsl/gsl_statistics_double.h>
#include <stdlib.h>
#cgo LDFLAGS: -lgsl -lgslcblas
*/
import "C"

import (
	"unsafe"
)

// gsl_stats_Qn0_from_sorted_data
func Qn0FromSortedData(sortedData []float64, stride int, n int) (work float64, workInt int32, ret float64) {
	_sortedData := (*C.double)(nil)
	if len(sortedData) > 0 {
		_sortedData = (*C.double)(unsafe.Pointer(&sortedData[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_work := (*C.double)(unsafe.Pointer(&work))
	_workInt := (*C.int)(unsafe.Pointer(&workInt))
	_ret := C.gsl_stats_Qn0_from_sorted_data(_sortedData, _stride, _n, _work, _workInt)
	ret = float64(_ret)
	return
}

// gsl_stats_Qn_from_sorted_data
func QnFromSortedData(sortedData []float64, stride int, n int) (work float64, workInt int32, ret float64) {
	_sortedData := (*C.double)(nil)
	if len(sortedData) > 0 {
		_sortedData = (*C.double)(unsafe.Pointer(&sortedData[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_work := (*C.double)(unsafe.Pointer(&work))
	_workInt := (*C.int)(unsafe.Pointer(&workInt))
	_ret := C.gsl_stats_Qn_from_sorted_data(_sortedData, _stride, _n, _work, _workInt)
	ret = float64(_ret)
	return
}

// gsl_stats_Sn0_from_sorted_data
func Sn0FromSortedData(sortedData []float64, stride int, n int) (work float64, ret float64) {
	_sortedData := (*C.double)(nil)
	if len(sortedData) > 0 {
		_sortedData = (*C.double)(unsafe.Pointer(&sortedData[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_work := (*C.double)(unsafe.Pointer(&work))
	_ret := C.gsl_stats_Sn0_from_sorted_data(_sortedData, _stride, _n, _work)
	ret = float64(_ret)
	return
}

// gsl_stats_Sn_from_sorted_data
func SnFromSortedData(sortedData []float64, stride int, n int) (work float64, ret float64) {
	_sortedData := (*C.double)(nil)
	if len(sortedData) > 0 {
		_sortedData = (*C.double)(unsafe.Pointer(&sortedData[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_work := (*C.double)(unsafe.Pointer(&work))
	_ret := C.gsl_stats_Sn_from_sorted_data(_sortedData, _stride, _n, _work)
	ret = float64(_ret)
	return
}

// gsl_stats_absdev
func Absdev(data []float64, stride int, n int) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_absdev(_data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_absdev_m
func AbsdevM(data []float64, stride int, n int, mean float64) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_mean := C.double(mean)
	_ret := C.gsl_stats_absdev_m(_data, _stride, _n, _mean)
	ret = float64(_ret)
	return
}

// gsl_stats_correlation
func Correlation(data1 []float64, stride1 int, data2 []float64, stride2 int, n int) (ret float64) {
	_data1 := (*C.double)(nil)
	if len(data1) > 0 {
		_data1 = (*C.double)(unsafe.Pointer(&data1[0]))
	}
	_stride1 := C.size_t(stride1)
	_data2 := (*C.double)(nil)
	if len(data2) > 0 {
		_data2 = (*C.double)(unsafe.Pointer(&data2[0]))
	}
	_stride2 := C.size_t(stride2)
	_n := C.size_t(n)
	_ret := C.gsl_stats_correlation(_data1, _stride1, _data2, _stride2, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_covariance
func Covariance(data1 []float64, stride1 int, data2 []float64, stride2 int, n int) (ret float64) {
	_data1 := (*C.double)(nil)
	if len(data1) > 0 {
		_data1 = (*C.double)(unsafe.Pointer(&data1[0]))
	}
	_stride1 := C.size_t(stride1)
	_data2 := (*C.double)(nil)
	if len(data2) > 0 {
		_data2 = (*C.double)(unsafe.Pointer(&data2[0]))
	}
	_stride2 := C.size_t(stride2)
	_n := C.size_t(n)
	_ret := C.gsl_stats_covariance(_data1, _stride1, _data2, _stride2, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_covariance_m
func CovarianceM(data1 []float64, stride1 int, data2 []float64, stride2 int, n int, mean1 float64, mean2 float64) (ret float64) {
	_data1 := (*C.double)(nil)
	if len(data1) > 0 {
		_data1 = (*C.double)(unsafe.Pointer(&data1[0]))
	}
	_stride1 := C.size_t(stride1)
	_data2 := (*C.double)(nil)
	if len(data2) > 0 {
		_data2 = (*C.double)(unsafe.Pointer(&data2[0]))
	}
	_stride2 := C.size_t(stride2)
	_n := C.size_t(n)
	_mean1 := C.double(mean1)
	_mean2 := C.double(mean2)
	_ret := C.gsl_stats_covariance_m(_data1, _stride1, _data2, _stride2, _n, _mean1, _mean2)
	ret = float64(_ret)
	return
}

// gsl_stats_gastwirth_from_sorted_data
func GastwirthFromSortedData(sortedData []float64, stride int, n int) (ret float64) {
	_sortedData := (*C.double)(nil)
	if len(sortedData) > 0 {
		_sortedData = (*C.double)(unsafe.Pointer(&sortedData[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_gastwirth_from_sorted_data(_sortedData, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_kurtosis
func Kurtosis(data []float64, stride int, n int) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_kurtosis(_data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_kurtosis_m_sd
func KurtosisMSd(data []float64, stride int, n int, mean float64, sd float64) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_mean := C.double(mean)
	_sd := C.double(sd)
	_ret := C.gsl_stats_kurtosis_m_sd(_data, _stride, _n, _mean, _sd)
	ret = float64(_ret)
	return
}

// gsl_stats_lag1_autocorrelation
func Lag1Autocorrelation(data []float64, stride int, n int) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_lag1_autocorrelation(_data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_lag1_autocorrelation_m
func Lag1AutocorrelationM(data []float64, stride int, n int, mean float64) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_mean := C.double(mean)
	_ret := C.gsl_stats_lag1_autocorrelation_m(_data, _stride, _n, _mean)
	ret = float64(_ret)
	return
}

// gsl_stats_mad
func Mad(data []float64, stride int, n int) (work float64, ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_work := (*C.double)(unsafe.Pointer(&work))
	_ret := C.gsl_stats_mad(_data, _stride, _n, _work)
	ret = float64(_ret)
	return
}

// gsl_stats_mad0
func Mad0(data []float64, stride int, n int) (work float64, ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_work := (*C.double)(unsafe.Pointer(&work))
	_ret := C.gsl_stats_mad0(_data, _stride, _n, _work)
	ret = float64(_ret)
	return
}

// gsl_stats_max
func Max(data []float64, stride int, n int) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_max(_data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_max_index
func MaxIndex(data []float64, stride int, n int) (ret int) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_max_index(_data, _stride, _n)
	ret = int(_ret)
	return
}

// gsl_stats_mean
func Mean(data []float64, stride int, n int) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_mean(_data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_median
func Median(stride int, n int) (sortedData float64, ret float64) {
	_sortedData := (*C.double)(unsafe.Pointer(&sortedData))
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_median(_sortedData, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_median_from_sorted_data
func MedianFromSortedData(sortedData []float64, stride int, n int) (ret float64) {
	_sortedData := (*C.double)(nil)
	if len(sortedData) > 0 {
		_sortedData = (*C.double)(unsafe.Pointer(&sortedData[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_median_from_sorted_data(_sortedData, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_min
func Min(data []float64, stride int, n int) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_min(_data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_min_index
func MinIndex(data []float64, stride int, n int) (ret int) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_min_index(_data, _stride, _n)
	ret = int(_ret)
	return
}

// gsl_stats_minmax
func Minmax(data []float64, stride int, n int) (min float64, max float64) {
	_min := (*C.double)(unsafe.Pointer(&min))
	_max := (*C.double)(unsafe.Pointer(&max))
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	C.gsl_stats_minmax(_min, _max, _data, _stride, _n)
	return
}

// gsl_stats_minmax_index
func MinmaxIndex(data []float64, stride int, n int) (minIndex int, maxIndex int) {
	_minIndex := (*C.size_t)(unsafe.Pointer(&minIndex))
	_maxIndex := (*C.size_t)(unsafe.Pointer(&maxIndex))
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	C.gsl_stats_minmax_index(_minIndex, _maxIndex, _data, _stride, _n)
	return
}

// gsl_stats_pvariance
func Pvariance(data1 []float64, stride1 int, n1 int, data2 []float64, stride2 int, n2 int) (ret float64) {
	_data1 := (*C.double)(nil)
	if len(data1) > 0 {
		_data1 = (*C.double)(unsafe.Pointer(&data1[0]))
	}
	_stride1 := C.size_t(stride1)
	_n1 := C.size_t(n1)
	_data2 := (*C.double)(nil)
	if len(data2) > 0 {
		_data2 = (*C.double)(unsafe.Pointer(&data2[0]))
	}
	_stride2 := C.size_t(stride2)
	_n2 := C.size_t(n2)
	_ret := C.gsl_stats_pvariance(_data1, _stride1, _n1, _data2, _stride2, _n2)
	ret = float64(_ret)
	return
}

// gsl_stats_quantile_from_sorted_data
func QuantileFromSortedData(sortedData []float64, stride int, n int, f float64) (ret float64) {
	_sortedData := (*C.double)(nil)
	if len(sortedData) > 0 {
		_sortedData = (*C.double)(unsafe.Pointer(&sortedData[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_f := C.double(f)
	_ret := C.gsl_stats_quantile_from_sorted_data(_sortedData, _stride, _n, _f)
	ret = float64(_ret)
	return
}

// gsl_stats_sd
func Sd(data []float64, stride int, n int) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_sd(_data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_sd_m
func SdM(data []float64, stride int, n int, mean float64) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_mean := C.double(mean)
	_ret := C.gsl_stats_sd_m(_data, _stride, _n, _mean)
	ret = float64(_ret)
	return
}

// gsl_stats_sd_with_fixed_mean
func SdWithFixedMean(data []float64, stride int, n int, mean float64) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_mean := C.double(mean)
	_ret := C.gsl_stats_sd_with_fixed_mean(_data, _stride, _n, _mean)
	ret = float64(_ret)
	return
}

// gsl_stats_select
func Select(stride int, n int, k int) (data float64, ret float64) {
	_data := (*C.double)(unsafe.Pointer(&data))
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_k := C.size_t(k)
	_ret := C.gsl_stats_select(_data, _stride, _n, _k)
	ret = float64(_ret)
	return
}

// gsl_stats_skew
func Skew(data []float64, stride int, n int) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_skew(_data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_skew_m_sd
func SkewMSd(data []float64, stride int, n int, mean float64, sd float64) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_mean := C.double(mean)
	_sd := C.double(sd)
	_ret := C.gsl_stats_skew_m_sd(_data, _stride, _n, _mean, _sd)
	ret = float64(_ret)
	return
}

// gsl_stats_spearman
func Spearman(data1 []float64, stride1 int, data2 []float64, stride2 int, n int) (work float64, ret float64) {
	_data1 := (*C.double)(nil)
	if len(data1) > 0 {
		_data1 = (*C.double)(unsafe.Pointer(&data1[0]))
	}
	_stride1 := C.size_t(stride1)
	_data2 := (*C.double)(nil)
	if len(data2) > 0 {
		_data2 = (*C.double)(unsafe.Pointer(&data2[0]))
	}
	_stride2 := C.size_t(stride2)
	_n := C.size_t(n)
	_work := (*C.double)(unsafe.Pointer(&work))
	_ret := C.gsl_stats_spearman(_data1, _stride1, _data2, _stride2, _n, _work)
	ret = float64(_ret)
	return
}

// gsl_stats_trmean_from_sorted_data
func TrmeanFromSortedData(trim float64, sortedData []float64, stride int, n int) (ret float64) {
	_trim := C.double(trim)
	_sortedData := (*C.double)(nil)
	if len(sortedData) > 0 {
		_sortedData = (*C.double)(unsafe.Pointer(&sortedData[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_trmean_from_sorted_data(_trim, _sortedData, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_tss
func Tss(data []float64, stride int, n int) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_tss(_data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_tss_m
func TssM(data []float64, stride int, n int, mean float64) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_mean := C.double(mean)
	_ret := C.gsl_stats_tss_m(_data, _stride, _n, _mean)
	ret = float64(_ret)
	return
}

// gsl_stats_ttest
func Ttest(data1 []float64, stride1 int, n1 int, data2 []float64, stride2 int, n2 int) (ret float64) {
	_data1 := (*C.double)(nil)
	if len(data1) > 0 {
		_data1 = (*C.double)(unsafe.Pointer(&data1[0]))
	}
	_stride1 := C.size_t(stride1)
	_n1 := C.size_t(n1)
	_data2 := (*C.double)(nil)
	if len(data2) > 0 {
		_data2 = (*C.double)(unsafe.Pointer(&data2[0]))
	}
	_stride2 := C.size_t(stride2)
	_n2 := C.size_t(n2)
	_ret := C.gsl_stats_ttest(_data1, _stride1, _n1, _data2, _stride2, _n2)
	ret = float64(_ret)
	return
}

// gsl_stats_variance
func Variance(data []float64, stride int, n int) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_variance(_data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_variance_m
func VarianceM(data []float64, stride int, n int, mean float64) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_mean := C.double(mean)
	_ret := C.gsl_stats_variance_m(_data, _stride, _n, _mean)
	ret = float64(_ret)
	return
}

// gsl_stats_variance_with_fixed_mean
func VarianceWithFixedMean(data []float64, stride int, n int, mean float64) (ret float64) {
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_mean := C.double(mean)
	_ret := C.gsl_stats_variance_with_fixed_mean(_data, _stride, _n, _mean)
	ret = float64(_ret)
	return
}

// gsl_stats_wabsdev
func Wabsdev(w []float64, wstride int, data []float64, stride int, n int) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_wabsdev(_w, _wstride, _data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_wabsdev_m
func WabsdevM(w []float64, wstride int, data []float64, stride int, n int, wmean float64) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_wmean := C.double(wmean)
	_ret := C.gsl_stats_wabsdev_m(_w, _wstride, _data, _stride, _n, _wmean)
	ret = float64(_ret)
	return
}

// gsl_stats_wkurtosis
func Wkurtosis(w []float64, wstride int, data []float64, stride int, n int) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_wkurtosis(_w, _wstride, _data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_wkurtosis_m_sd
func WkurtosisMSd(w []float64, wstride int, data []float64, stride int, n int, wmean float64, wsd float64) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_wmean := C.double(wmean)
	_wsd := C.double(wsd)
	_ret := C.gsl_stats_wkurtosis_m_sd(_w, _wstride, _data, _stride, _n, _wmean, _wsd)
	ret = float64(_ret)
	return
}

// gsl_stats_wmean
func Wmean(w []float64, wstride int, data []float64, stride int, n int) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_wmean(_w, _wstride, _data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_wsd
func Wsd(w []float64, wstride int, data []float64, stride int, n int) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_wsd(_w, _wstride, _data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_wsd_m
func WsdM(w []float64, wstride int, data []float64, stride int, n int, wmean float64) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_wmean := C.double(wmean)
	_ret := C.gsl_stats_wsd_m(_w, _wstride, _data, _stride, _n, _wmean)
	ret = float64(_ret)
	return
}

// gsl_stats_wsd_with_fixed_mean
func WsdWithFixedMean(w []float64, wstride int, data []float64, stride int, n int, mean float64) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_mean := C.double(mean)
	_ret := C.gsl_stats_wsd_with_fixed_mean(_w, _wstride, _data, _stride, _n, _mean)
	ret = float64(_ret)
	return
}

// gsl_stats_wskew
func Wskew(w []float64, wstride int, data []float64, stride int, n int) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_wskew(_w, _wstride, _data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_wskew_m_sd
func WskewMSd(w []float64, wstride int, data []float64, stride int, n int, wmean float64, wsd float64) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_wmean := C.double(wmean)
	_wsd := C.double(wsd)
	_ret := C.gsl_stats_wskew_m_sd(_w, _wstride, _data, _stride, _n, _wmean, _wsd)
	ret = float64(_ret)
	return
}

// gsl_stats_wtss
func Wtss(w []float64, wstride int, data []float64, stride int, n int) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_wtss(_w, _wstride, _data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_wtss_m
func WtssM(w []float64, wstride int, data []float64, stride int, n int, wmean float64) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_wmean := C.double(wmean)
	_ret := C.gsl_stats_wtss_m(_w, _wstride, _data, _stride, _n, _wmean)
	ret = float64(_ret)
	return
}

// gsl_stats_wvariance
func Wvariance(w []float64, wstride int, data []float64, stride int, n int) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_ret := C.gsl_stats_wvariance(_w, _wstride, _data, _stride, _n)
	ret = float64(_ret)
	return
}

// gsl_stats_wvariance_m
func WvarianceM(w []float64, wstride int, data []float64, stride int, n int, wmean float64) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_wmean := C.double(wmean)
	_ret := C.gsl_stats_wvariance_m(_w, _wstride, _data, _stride, _n, _wmean)
	ret = float64(_ret)
	return
}

// gsl_stats_wvariance_with_fixed_mean
func WvarianceWithFixedMean(w []float64, wstride int, data []float64, stride int, n int, mean float64) (ret float64) {
	_w := (*C.double)(nil)
	if len(w) > 0 {
		_w = (*C.double)(unsafe.Pointer(&w[0]))
	}
	_wstride := C.size_t(wstride)
	_data := (*C.double)(nil)
	if len(data) > 0 {
		_data = (*C.double)(unsafe.Pointer(&data[0]))
	}
	_stride := C.size_t(stride)
	_n := C.size_t(n)
	_mean := C.double(mean)
	_ret := C.gsl_stats_wvariance_with_fixed_mean(_w, _wstride, _data, _stride, _n, _mean)
	ret = float64(_ret)
	return
}
