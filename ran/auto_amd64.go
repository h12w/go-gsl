package ran

/*
#include <gsl/gsl_randist.h>
#include <stdlib.h>
#cgo LDFLAGS: -lgsl -lgslcblas
*/
import "C"

import (
	"go-gsl/rng"
	"unsafe"
)

// gsl_ran_discrete_t
type Discrete struct {
	K uint64
	A *uint64
	F *float64
}

// gsl_ran_discrete_free
func (g *Discrete) Free() {
	_g := (*C.gsl_ran_discrete_t)(unsafe.Pointer(g))
	C.gsl_ran_discrete_free(_g)
}

// gsl_ran_bernoulli
func Bernoulli(r *rng.GslRng, p float64) (ret uint32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_p := C.double(p)
	_ret := C.gsl_ran_bernoulli(_r, _p)
	ret = uint32(_ret)
	return
}

// gsl_ran_bernoulli_pdf
func BernoulliPdf(k uint32, p float64) (ret float64) {
	_k := C.uint(k)
	_p := C.double(p)
	_ret := C.gsl_ran_bernoulli_pdf(_k, _p)
	ret = float64(_ret)
	return
}

// gsl_ran_beta
func Beta(r *rng.GslRng, a float64, b float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_beta(_r, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_beta_pdf
func BetaPdf(x float64, a float64, b float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_beta_pdf(_x, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_binomial
func Binomial(r *rng.GslRng, p float64, n uint32) (ret uint32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_p := C.double(p)
	_n := C.uint(n)
	_ret := C.gsl_ran_binomial(_r, _p, _n)
	ret = uint32(_ret)
	return
}

// gsl_ran_binomial_knuth
func BinomialKnuth(r *rng.GslRng, p float64, n uint32) (ret uint32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_p := C.double(p)
	_n := C.uint(n)
	_ret := C.gsl_ran_binomial_knuth(_r, _p, _n)
	ret = uint32(_ret)
	return
}

// gsl_ran_binomial_pdf
func BinomialPdf(k uint32, p float64, n uint32) (ret float64) {
	_k := C.uint(k)
	_p := C.double(p)
	_n := C.uint(n)
	_ret := C.gsl_ran_binomial_pdf(_k, _p, _n)
	ret = float64(_ret)
	return
}

// gsl_ran_binomial_tpe
func BinomialTpe(r *rng.GslRng, p float64, n uint32) (ret uint32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_p := C.double(p)
	_n := C.uint(n)
	_ret := C.gsl_ran_binomial_tpe(_r, _p, _n)
	ret = uint32(_ret)
	return
}

// gsl_ran_bivariate_gaussian
func BivariateGaussian(r *rng.GslRng, sigmaX float64, sigmaY float64, rho float64) (x float64, y float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_sigmaX := C.double(sigmaX)
	_sigmaY := C.double(sigmaY)
	_rho := C.double(rho)
	_x := (*C.double)(unsafe.Pointer(&x))
	_y := (*C.double)(unsafe.Pointer(&y))
	C.gsl_ran_bivariate_gaussian(_r, _sigmaX, _sigmaY, _rho, _x, _y)
	return
}

// gsl_ran_bivariate_gaussian_pdf
func BivariateGaussianPdf(x float64, y float64, sigmaX float64, sigmaY float64, rho float64) (ret float64) {
	_x := C.double(x)
	_y := C.double(y)
	_sigmaX := C.double(sigmaX)
	_sigmaY := C.double(sigmaY)
	_rho := C.double(rho)
	_ret := C.gsl_ran_bivariate_gaussian_pdf(_x, _y, _sigmaX, _sigmaY, _rho)
	ret = float64(_ret)
	return
}

// gsl_ran_cauchy
func Cauchy(r *rng.GslRng, a float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_ret := C.gsl_ran_cauchy(_r, _a)
	ret = float64(_ret)
	return
}

// gsl_ran_cauchy_pdf
func CauchyPdf(x float64, a float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_ret := C.gsl_ran_cauchy_pdf(_x, _a)
	ret = float64(_ret)
	return
}

// gsl_ran_chisq
func Chisq(r *rng.GslRng, nu float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_nu := C.double(nu)
	_ret := C.gsl_ran_chisq(_r, _nu)
	ret = float64(_ret)
	return
}

// gsl_ran_chisq_pdf
func ChisqPdf(x float64, nu float64) (ret float64) {
	_x := C.double(x)
	_nu := C.double(nu)
	_ret := C.gsl_ran_chisq_pdf(_x, _nu)
	ret = float64(_ret)
	return
}

// gsl_ran_choose
func Choose(r *rng.GslRng, dest uintptr, k int, src uintptr, n int, size int) (ret int32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_dest := unsafe.Pointer(dest)
	_k := C.size_t(k)
	_src := unsafe.Pointer(src)
	_n := C.size_t(n)
	_size := C.size_t(size)
	_ret := C.gsl_ran_choose(_r, _dest, _k, _src, _n, _size)
	ret = int32(_ret)
	return
}

// gsl_ran_dir_2d
func Dir2d(r *rng.GslRng) (x float64, y float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_x := (*C.double)(unsafe.Pointer(&x))
	_y := (*C.double)(unsafe.Pointer(&y))
	C.gsl_ran_dir_2d(_r, _x, _y)
	return
}

// gsl_ran_dir_2d_trig_method
func Dir2dTrigMethod(r *rng.GslRng) (x float64, y float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_x := (*C.double)(unsafe.Pointer(&x))
	_y := (*C.double)(unsafe.Pointer(&y))
	C.gsl_ran_dir_2d_trig_method(_r, _x, _y)
	return
}

// gsl_ran_dir_3d
func Dir3d(r *rng.GslRng) (x float64, y float64, z float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_x := (*C.double)(unsafe.Pointer(&x))
	_y := (*C.double)(unsafe.Pointer(&y))
	_z := (*C.double)(unsafe.Pointer(&z))
	C.gsl_ran_dir_3d(_r, _x, _y, _z)
	return
}

// gsl_ran_dir_nd
func DirNd(r *rng.GslRng, n int) (x float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_n := C.size_t(n)
	_x := (*C.double)(unsafe.Pointer(&x))
	C.gsl_ran_dir_nd(_r, _n, _x)
	return
}

// gsl_ran_dirichlet
func Dirichlet(r *rng.GslRng, K int, alpha []float64) (theta float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_K := C.size_t(K)
	_alpha := (*C.double)(nil)
	if len(alpha) > 0 {
		_alpha = (*C.double)(unsafe.Pointer(&alpha[0]))
	}
	_theta := (*C.double)(unsafe.Pointer(&theta))
	C.gsl_ran_dirichlet(_r, _K, _alpha, _theta)
	return
}

// gsl_ran_dirichlet_lnpdf
func DirichletLnpdf(K int, alpha []float64, theta []float64) (ret float64) {
	_K := C.size_t(K)
	_alpha := (*C.double)(nil)
	if len(alpha) > 0 {
		_alpha = (*C.double)(unsafe.Pointer(&alpha[0]))
	}
	_theta := (*C.double)(nil)
	if len(theta) > 0 {
		_theta = (*C.double)(unsafe.Pointer(&theta[0]))
	}
	_ret := C.gsl_ran_dirichlet_lnpdf(_K, _alpha, _theta)
	ret = float64(_ret)
	return
}

// gsl_ran_dirichlet_pdf
func DirichletPdf(K int, alpha []float64, theta []float64) (ret float64) {
	_K := C.size_t(K)
	_alpha := (*C.double)(nil)
	if len(alpha) > 0 {
		_alpha = (*C.double)(unsafe.Pointer(&alpha[0]))
	}
	_theta := (*C.double)(nil)
	if len(theta) > 0 {
		_theta = (*C.double)(unsafe.Pointer(&theta[0]))
	}
	_ret := C.gsl_ran_dirichlet_pdf(_K, _alpha, _theta)
	ret = float64(_ret)
	return
}

// gsl_ran_discrete
func Discrete_(r *rng.GslRng, g *Discrete) (ret int) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_g := (*C.gsl_ran_discrete_t)(unsafe.Pointer(g))
	_ret := C.gsl_ran_discrete(_r, _g)
	ret = int(_ret)
	return
}

// gsl_ran_discrete_pdf
func DiscretePdf(k int, g *Discrete) (ret float64) {
	_k := C.size_t(k)
	_g := (*C.gsl_ran_discrete_t)(unsafe.Pointer(g))
	_ret := C.gsl_ran_discrete_pdf(_k, _g)
	ret = float64(_ret)
	return
}

// gsl_ran_discrete_preproc
func DiscretePreproc(K int, P []float64) (ret *Discrete) {
	_K := C.size_t(K)
	_P := (*C.double)(nil)
	if len(P) > 0 {
		_P = (*C.double)(unsafe.Pointer(&P[0]))
	}
	_ret := C.gsl_ran_discrete_preproc(_K, _P)
	ret = (*Discrete)(unsafe.Pointer(_ret))
	return
}

// gsl_ran_erlang
func Erlang(r *rng.GslRng, a float64, n float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_n := C.double(n)
	_ret := C.gsl_ran_erlang(_r, _a, _n)
	ret = float64(_ret)
	return
}

// gsl_ran_erlang_pdf
func ErlangPdf(x float64, a float64, n float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_n := C.double(n)
	_ret := C.gsl_ran_erlang_pdf(_x, _a, _n)
	ret = float64(_ret)
	return
}

// gsl_ran_exponential
func Exponential(r *rng.GslRng, mu float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_mu := C.double(mu)
	_ret := C.gsl_ran_exponential(_r, _mu)
	ret = float64(_ret)
	return
}

// gsl_ran_exponential_pdf
func ExponentialPdf(x float64, mu float64) (ret float64) {
	_x := C.double(x)
	_mu := C.double(mu)
	_ret := C.gsl_ran_exponential_pdf(_x, _mu)
	ret = float64(_ret)
	return
}

// gsl_ran_exppow
func Exppow(r *rng.GslRng, a float64, b float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_exppow(_r, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_exppow_pdf
func ExppowPdf(x float64, a float64, b float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_exppow_pdf(_x, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_fdist
func Fdist(r *rng.GslRng, nu1 float64, nu2 float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_nu1 := C.double(nu1)
	_nu2 := C.double(nu2)
	_ret := C.gsl_ran_fdist(_r, _nu1, _nu2)
	ret = float64(_ret)
	return
}

// gsl_ran_fdist_pdf
func FdistPdf(x float64, nu1 float64, nu2 float64) (ret float64) {
	_x := C.double(x)
	_nu1 := C.double(nu1)
	_nu2 := C.double(nu2)
	_ret := C.gsl_ran_fdist_pdf(_x, _nu1, _nu2)
	ret = float64(_ret)
	return
}

// gsl_ran_flat
func Flat(r *rng.GslRng, a float64, b float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_flat(_r, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_flat_pdf
func FlatPdf(x float64, a float64, b float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_flat_pdf(_x, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_gamma
func Gamma(r *rng.GslRng, a float64, b float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_gamma(_r, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_gamma_int
func GammaInt(r *rng.GslRng, a uint32) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.uint(a)
	_ret := C.gsl_ran_gamma_int(_r, _a)
	ret = float64(_ret)
	return
}

// gsl_ran_gamma_knuth
func GammaKnuth(r *rng.GslRng, a float64, b float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_gamma_knuth(_r, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_gamma_mt
func GammaMt(r *rng.GslRng, a float64, b float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_gamma_mt(_r, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_gamma_pdf
func GammaPdf(x float64, a float64, b float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_gamma_pdf(_x, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_gaussian
func Gaussian(r *rng.GslRng, sigma float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_sigma := C.double(sigma)
	_ret := C.gsl_ran_gaussian(_r, _sigma)
	ret = float64(_ret)
	return
}

// gsl_ran_gaussian_pdf
func GaussianPdf(x float64, sigma float64) (ret float64) {
	_x := C.double(x)
	_sigma := C.double(sigma)
	_ret := C.gsl_ran_gaussian_pdf(_x, _sigma)
	ret = float64(_ret)
	return
}

// gsl_ran_gaussian_ratio_method
func GaussianRatioMethod(r *rng.GslRng, sigma float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_sigma := C.double(sigma)
	_ret := C.gsl_ran_gaussian_ratio_method(_r, _sigma)
	ret = float64(_ret)
	return
}

// gsl_ran_gaussian_tail
func GaussianTail(r *rng.GslRng, a float64, sigma float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_sigma := C.double(sigma)
	_ret := C.gsl_ran_gaussian_tail(_r, _a, _sigma)
	ret = float64(_ret)
	return
}

// gsl_ran_gaussian_tail_pdf
func GaussianTailPdf(x float64, a float64, sigma float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_sigma := C.double(sigma)
	_ret := C.gsl_ran_gaussian_tail_pdf(_x, _a, _sigma)
	ret = float64(_ret)
	return
}

// gsl_ran_gaussian_ziggurat
func GaussianZiggurat(r *rng.GslRng, sigma float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_sigma := C.double(sigma)
	_ret := C.gsl_ran_gaussian_ziggurat(_r, _sigma)
	ret = float64(_ret)
	return
}

// gsl_ran_geometric
func Geometric(r *rng.GslRng, p float64) (ret uint32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_p := C.double(p)
	_ret := C.gsl_ran_geometric(_r, _p)
	ret = uint32(_ret)
	return
}

// gsl_ran_geometric_pdf
func GeometricPdf(k uint32, p float64) (ret float64) {
	_k := C.uint(k)
	_p := C.double(p)
	_ret := C.gsl_ran_geometric_pdf(_k, _p)
	ret = float64(_ret)
	return
}

// gsl_ran_gumbel1
func Gumbel1(r *rng.GslRng, a float64, b float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_gumbel1(_r, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_gumbel1_pdf
func Gumbel1Pdf(x float64, a float64, b float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_gumbel1_pdf(_x, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_gumbel2
func Gumbel2(r *rng.GslRng, a float64, b float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_gumbel2(_r, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_gumbel2_pdf
func Gumbel2Pdf(x float64, a float64, b float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_gumbel2_pdf(_x, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_hypergeometric
func Hypergeometric(r *rng.GslRng, n1 uint32, n2 uint32, t uint32) (ret uint32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_n1 := C.uint(n1)
	_n2 := C.uint(n2)
	_t := C.uint(t)
	_ret := C.gsl_ran_hypergeometric(_r, _n1, _n2, _t)
	ret = uint32(_ret)
	return
}

// gsl_ran_hypergeometric_pdf
func HypergeometricPdf(k uint32, n1 uint32, n2 uint32, t uint32) (ret float64) {
	_k := C.uint(k)
	_n1 := C.uint(n1)
	_n2 := C.uint(n2)
	_t := C.uint(t)
	_ret := C.gsl_ran_hypergeometric_pdf(_k, _n1, _n2, _t)
	ret = float64(_ret)
	return
}

// gsl_ran_landau
func Landau(r *rng.GslRng) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_ret := C.gsl_ran_landau(_r)
	ret = float64(_ret)
	return
}

// gsl_ran_landau_pdf
func LandauPdf(x float64) (ret float64) {
	_x := C.double(x)
	_ret := C.gsl_ran_landau_pdf(_x)
	ret = float64(_ret)
	return
}

// gsl_ran_laplace
func Laplace(r *rng.GslRng, a float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_ret := C.gsl_ran_laplace(_r, _a)
	ret = float64(_ret)
	return
}

// gsl_ran_laplace_pdf
func LaplacePdf(x float64, a float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_ret := C.gsl_ran_laplace_pdf(_x, _a)
	ret = float64(_ret)
	return
}

// gsl_ran_levy
func Levy(r *rng.GslRng, c float64, alpha float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_c := C.double(c)
	_alpha := C.double(alpha)
	_ret := C.gsl_ran_levy(_r, _c, _alpha)
	ret = float64(_ret)
	return
}

// gsl_ran_levy_skew
func LevySkew(r *rng.GslRng, c float64, alpha float64, beta float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_c := C.double(c)
	_alpha := C.double(alpha)
	_beta := C.double(beta)
	_ret := C.gsl_ran_levy_skew(_r, _c, _alpha, _beta)
	ret = float64(_ret)
	return
}

// gsl_ran_logarithmic
func Logarithmic(r *rng.GslRng, p float64) (ret uint32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_p := C.double(p)
	_ret := C.gsl_ran_logarithmic(_r, _p)
	ret = uint32(_ret)
	return
}

// gsl_ran_logarithmic_pdf
func LogarithmicPdf(k uint32, p float64) (ret float64) {
	_k := C.uint(k)
	_p := C.double(p)
	_ret := C.gsl_ran_logarithmic_pdf(_k, _p)
	ret = float64(_ret)
	return
}

// gsl_ran_logistic
func Logistic(r *rng.GslRng, a float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_ret := C.gsl_ran_logistic(_r, _a)
	ret = float64(_ret)
	return
}

// gsl_ran_logistic_pdf
func LogisticPdf(x float64, a float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_ret := C.gsl_ran_logistic_pdf(_x, _a)
	ret = float64(_ret)
	return
}

// gsl_ran_lognormal
func Lognormal(r *rng.GslRng, zeta float64, sigma float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_zeta := C.double(zeta)
	_sigma := C.double(sigma)
	_ret := C.gsl_ran_lognormal(_r, _zeta, _sigma)
	ret = float64(_ret)
	return
}

// gsl_ran_lognormal_pdf
func LognormalPdf(x float64, zeta float64, sigma float64) (ret float64) {
	_x := C.double(x)
	_zeta := C.double(zeta)
	_sigma := C.double(sigma)
	_ret := C.gsl_ran_lognormal_pdf(_x, _zeta, _sigma)
	ret = float64(_ret)
	return
}

// gsl_ran_multinomial
func Multinomial(r *rng.GslRng, K int, N uint32, p []float64) (n uint32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_K := C.size_t(K)
	_N := C.uint(N)
	_p := (*C.double)(nil)
	if len(p) > 0 {
		_p = (*C.double)(unsafe.Pointer(&p[0]))
	}
	_n := (*C.uint)(unsafe.Pointer(&n))
	C.gsl_ran_multinomial(_r, _K, _N, _p, _n)
	return
}

// gsl_ran_multinomial_lnpdf
func MultinomialLnpdf(K int, p []float64, n []uint32) (ret float64) {
	_K := C.size_t(K)
	_p := (*C.double)(nil)
	if len(p) > 0 {
		_p = (*C.double)(unsafe.Pointer(&p[0]))
	}
	_n := (*C.uint)(nil)
	if len(n) > 0 {
		_n = (*C.uint)(unsafe.Pointer(&n[0]))
	}
	_ret := C.gsl_ran_multinomial_lnpdf(_K, _p, _n)
	ret = float64(_ret)
	return
}

// gsl_ran_multinomial_pdf
func MultinomialPdf(K int, p []float64, n []uint32) (ret float64) {
	_K := C.size_t(K)
	_p := (*C.double)(nil)
	if len(p) > 0 {
		_p = (*C.double)(unsafe.Pointer(&p[0]))
	}
	_n := (*C.uint)(nil)
	if len(n) > 0 {
		_n = (*C.uint)(unsafe.Pointer(&n[0]))
	}
	_ret := C.gsl_ran_multinomial_pdf(_K, _p, _n)
	ret = float64(_ret)
	return
}

// gsl_ran_multivariate_gaussian
func MultivariateGaussian(r *rng.GslRng, mu *[40]byte, L *[48]byte, result *[40]byte) (ret int32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_mu := (*C.gsl_vector)(unsafe.Pointer(mu))
	_L := (*C.gsl_matrix)(unsafe.Pointer(L))
	_result := (*C.gsl_vector)(unsafe.Pointer(result))
	_ret := C.gsl_ran_multivariate_gaussian(_r, _mu, _L, _result)
	ret = int32(_ret)
	return
}

// gsl_ran_multivariate_gaussian_log_pdf
func MultivariateGaussianLogPdf(x *[40]byte, mu *[40]byte, L *[48]byte, work *[40]byte) (result float64, ret int32) {
	_x := (*C.gsl_vector)(unsafe.Pointer(x))
	_mu := (*C.gsl_vector)(unsafe.Pointer(mu))
	_L := (*C.gsl_matrix)(unsafe.Pointer(L))
	_result := (*C.double)(unsafe.Pointer(&result))
	_work := (*C.gsl_vector)(unsafe.Pointer(work))
	_ret := C.gsl_ran_multivariate_gaussian_log_pdf(_x, _mu, _L, _result, _work)
	ret = int32(_ret)
	return
}

// gsl_ran_multivariate_gaussian_mean
func MultivariateGaussianMean(X *[48]byte, muHat *[40]byte) (ret int32) {
	_X := (*C.gsl_matrix)(unsafe.Pointer(X))
	_muHat := (*C.gsl_vector)(unsafe.Pointer(muHat))
	_ret := C.gsl_ran_multivariate_gaussian_mean(_X, _muHat)
	ret = int32(_ret)
	return
}

// gsl_ran_multivariate_gaussian_pdf
func MultivariateGaussianPdf(x *[40]byte, mu *[40]byte, L *[48]byte, work *[40]byte) (result float64, ret int32) {
	_x := (*C.gsl_vector)(unsafe.Pointer(x))
	_mu := (*C.gsl_vector)(unsafe.Pointer(mu))
	_L := (*C.gsl_matrix)(unsafe.Pointer(L))
	_result := (*C.double)(unsafe.Pointer(&result))
	_work := (*C.gsl_vector)(unsafe.Pointer(work))
	_ret := C.gsl_ran_multivariate_gaussian_pdf(_x, _mu, _L, _result, _work)
	ret = int32(_ret)
	return
}

// gsl_ran_multivariate_gaussian_vcov
func MultivariateGaussianVcov(X *[48]byte, sigmaHat *[48]byte) (ret int32) {
	_X := (*C.gsl_matrix)(unsafe.Pointer(X))
	_sigmaHat := (*C.gsl_matrix)(unsafe.Pointer(sigmaHat))
	_ret := C.gsl_ran_multivariate_gaussian_vcov(_X, _sigmaHat)
	ret = int32(_ret)
	return
}

// gsl_ran_negative_binomial
func NegativeBinomial(r *rng.GslRng, p float64, n float64) (ret uint32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_p := C.double(p)
	_n := C.double(n)
	_ret := C.gsl_ran_negative_binomial(_r, _p, _n)
	ret = uint32(_ret)
	return
}

// gsl_ran_negative_binomial_pdf
func NegativeBinomialPdf(k uint32, p float64, n float64) (ret float64) {
	_k := C.uint(k)
	_p := C.double(p)
	_n := C.double(n)
	_ret := C.gsl_ran_negative_binomial_pdf(_k, _p, _n)
	ret = float64(_ret)
	return
}

// gsl_ran_pareto
func Pareto(r *rng.GslRng, a float64, b float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_pareto(_r, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_pareto_pdf
func ParetoPdf(x float64, a float64, b float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_pareto_pdf(_x, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_pascal
func Pascal(r *rng.GslRng, p float64, n uint32) (ret uint32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_p := C.double(p)
	_n := C.uint(n)
	_ret := C.gsl_ran_pascal(_r, _p, _n)
	ret = uint32(_ret)
	return
}

// gsl_ran_pascal_pdf
func PascalPdf(k uint32, p float64, n uint32) (ret float64) {
	_k := C.uint(k)
	_p := C.double(p)
	_n := C.uint(n)
	_ret := C.gsl_ran_pascal_pdf(_k, _p, _n)
	ret = float64(_ret)
	return
}

// gsl_ran_poisson
func Poisson(r *rng.GslRng, mu float64) (ret uint32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_mu := C.double(mu)
	_ret := C.gsl_ran_poisson(_r, _mu)
	ret = uint32(_ret)
	return
}

// gsl_ran_poisson_array
func PoissonArray(r *rng.GslRng, n int, mu float64) (array uint32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_n := C.size_t(n)
	_array := (*C.uint)(unsafe.Pointer(&array))
	_mu := C.double(mu)
	C.gsl_ran_poisson_array(_r, _n, _array, _mu)
	return
}

// gsl_ran_poisson_pdf
func PoissonPdf(k uint32, mu float64) (ret float64) {
	_k := C.uint(k)
	_mu := C.double(mu)
	_ret := C.gsl_ran_poisson_pdf(_k, _mu)
	ret = float64(_ret)
	return
}

// gsl_ran_rayleigh
func Rayleigh(r *rng.GslRng, sigma float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_sigma := C.double(sigma)
	_ret := C.gsl_ran_rayleigh(_r, _sigma)
	ret = float64(_ret)
	return
}

// gsl_ran_rayleigh_pdf
func RayleighPdf(x float64, sigma float64) (ret float64) {
	_x := C.double(x)
	_sigma := C.double(sigma)
	_ret := C.gsl_ran_rayleigh_pdf(_x, _sigma)
	ret = float64(_ret)
	return
}

// gsl_ran_rayleigh_tail
func RayleighTail(r *rng.GslRng, a float64, sigma float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_sigma := C.double(sigma)
	_ret := C.gsl_ran_rayleigh_tail(_r, _a, _sigma)
	ret = float64(_ret)
	return
}

// gsl_ran_rayleigh_tail_pdf
func RayleighTailPdf(x float64, a float64, sigma float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_sigma := C.double(sigma)
	_ret := C.gsl_ran_rayleigh_tail_pdf(_x, _a, _sigma)
	ret = float64(_ret)
	return
}

// gsl_ran_sample
func Sample(r *rng.GslRng, dest uintptr, k int, src uintptr, n int, size int) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_dest := unsafe.Pointer(dest)
	_k := C.size_t(k)
	_src := unsafe.Pointer(src)
	_n := C.size_t(n)
	_size := C.size_t(size)
	C.gsl_ran_sample(_r, _dest, _k, _src, _n, _size)
}

// gsl_ran_shuffle
func Shuffle(r *rng.GslRng, base uintptr, nmembm int, size int) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_base := unsafe.Pointer(base)
	_nmembm := C.size_t(nmembm)
	_size := C.size_t(size)
	C.gsl_ran_shuffle(_r, _base, _nmembm, _size)
}

// gsl_ran_tdist
func Tdist(r *rng.GslRng, nu float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_nu := C.double(nu)
	_ret := C.gsl_ran_tdist(_r, _nu)
	ret = float64(_ret)
	return
}

// gsl_ran_tdist_pdf
func TdistPdf(x float64, nu float64) (ret float64) {
	_x := C.double(x)
	_nu := C.double(nu)
	_ret := C.gsl_ran_tdist_pdf(_x, _nu)
	ret = float64(_ret)
	return
}

// gsl_ran_ugaussian
func Ugaussian(r *rng.GslRng) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_ret := C.gsl_ran_ugaussian(_r)
	ret = float64(_ret)
	return
}

// gsl_ran_ugaussian_pdf
func UgaussianPdf(x float64) (ret float64) {
	_x := C.double(x)
	_ret := C.gsl_ran_ugaussian_pdf(_x)
	ret = float64(_ret)
	return
}

// gsl_ran_ugaussian_ratio_method
func UgaussianRatioMethod(r *rng.GslRng) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_ret := C.gsl_ran_ugaussian_ratio_method(_r)
	ret = float64(_ret)
	return
}

// gsl_ran_ugaussian_tail
func UgaussianTail(r *rng.GslRng, a float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_ret := C.gsl_ran_ugaussian_tail(_r, _a)
	ret = float64(_ret)
	return
}

// gsl_ran_ugaussian_tail_pdf
func UgaussianTailPdf(x float64, a float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_ret := C.gsl_ran_ugaussian_tail_pdf(_x, _a)
	ret = float64(_ret)
	return
}

// gsl_ran_weibull
func Weibull(r *rng.GslRng, a float64, b float64) (ret float64) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_weibull(_r, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_weibull_pdf
func WeibullPdf(x float64, a float64, b float64) (ret float64) {
	_x := C.double(x)
	_a := C.double(a)
	_b := C.double(b)
	_ret := C.gsl_ran_weibull_pdf(_x, _a, _b)
	ret = float64(_ret)
	return
}

// gsl_ran_wishart
func Wishart(r *rng.GslRng, df float64, L *[48]byte, result *[48]byte, work *[48]byte) (ret int32) {
	_r := (*C.gsl_rng)(unsafe.Pointer(r))
	_df := C.double(df)
	_L := (*C.gsl_matrix)(unsafe.Pointer(L))
	_result := (*C.gsl_matrix)(unsafe.Pointer(result))
	_work := (*C.gsl_matrix)(unsafe.Pointer(work))
	_ret := C.gsl_ran_wishart(_r, _df, _L, _result, _work)
	ret = int32(_ret)
	return
}

// gsl_ran_wishart_log_pdf
func WishartLogPdf(X *[48]byte, lX *[48]byte, df float64, L *[48]byte, work *[48]byte) (result float64, ret int32) {
	_X := (*C.gsl_matrix)(unsafe.Pointer(X))
	_lX := (*C.gsl_matrix)(unsafe.Pointer(lX))
	_df := C.double(df)
	_L := (*C.gsl_matrix)(unsafe.Pointer(L))
	_result := (*C.double)(unsafe.Pointer(&result))
	_work := (*C.gsl_matrix)(unsafe.Pointer(work))
	_ret := C.gsl_ran_wishart_log_pdf(_X, _lX, _df, _L, _result, _work)
	ret = int32(_ret)
	return
}

// gsl_ran_wishart_pdf
func WishartPdf(X *[48]byte, lX *[48]byte, df float64, L *[48]byte, work *[48]byte) (result float64, ret int32) {
	_X := (*C.gsl_matrix)(unsafe.Pointer(X))
	_lX := (*C.gsl_matrix)(unsafe.Pointer(lX))
	_df := C.double(df)
	_L := (*C.gsl_matrix)(unsafe.Pointer(L))
	_result := (*C.double)(unsafe.Pointer(&result))
	_work := (*C.gsl_matrix)(unsafe.Pointer(work))
	_ret := C.gsl_ran_wishart_pdf(_X, _lX, _df, _L, _result, _work)
	ret = int32(_ret)
	return
}
