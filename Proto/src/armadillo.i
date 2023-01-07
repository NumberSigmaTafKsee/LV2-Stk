%module arma
%{
#include<armadillo>
#include<vector>
#include<complex>
#include<string>

using namespace arma;
%}
%include "std_complex.i"
%include "std_vector.i"
%include "std_string.i"
%include "arma.i"

%template(fcomplex) std::complex<float>;
%template(xcomplex) std::complex<double>;

%template(fmat) arma::Mat<float>;
%template(dmat) arma::Mat<double>;

%template(cfmat) arma::Mat<std::complex<float>>;
%template(cxmat) arma::Mat<std::complex<double>>;
%template(imat) arma::Mat<int>;
%template(umat) arma::Mat<unsigned>;
%template(lmat) arma::Mat<long>;
%template(ulmat) arma::Mat<unsigned long>;
%template(llmat) arma::Mat<long long>;
%template(ullmat) arma::Mat<unsigned long long>;

%template(fvec) arma::Col<float>;
%template(dvec) arma::Col<double>;
%template(cfvec) arma::Col<std::complex<float>>;
%template(cxvec) arma::Col<std::complex<double>>;

// i'm not sure what armadillo calls row atm
%template(frow) arma::Row<float>;
%template(drow) arma::Row<double>;
%template(cfrow) arma::Row<std::complex<float>>;
%template(cxrow) arma::Row<std::complex<double>>;




// they only work for double
%template(linspace) ops::linspace<double>;
%template(logpace) ops::logspace<double>;
%template(regpace) ops::logspace<double>;
//%template(randperm) ops::randperm<double>;
%template(eye) ops::eye<double>;
%template(col_zero) ops::zero<arma::Col<double>>;
%template(col_ones) ops::zero<arma::Col<double>>;
%template(col_random) ops::zero<arma::Col<double>>;
%template(col_randu) ops::randu<arma::Col<double>>;
%template(random_vector) ops::random_vector<double>;
%template(toeplitz) ops::toeplitz<double>;
%template(circ_toeplitz) ops::circ_toeplitz<double>;
%template(accu) ops::accu<double>;
%template(arg) ops::arg<double>;
%template(sum) ops::sum<double>;
%template(row_sum) ops::row_sum<double>;
%template(col_sum) ops::col_sum<double>;
%template(affmul) ops::affmul<double>;
%template(cdot) ops::cdot<double>;
%template(approx_equal) ops::approx_equal<double>;
%template(as_scalar) ops::as_scalar<double>;
%template(clamp) ops::clamp<double>;
%template(cumsum) ops::cumsum<double>;
%template(cumprod) ops::cumprod<double>;
%template(det) ops::det<double>;
%template(diagmat) ops::diagmat<double>;
%template(diagvec) ops::diagvec<double>;
%template(diff) ops::diff<double>;
%template(dot) ops::dot<double>;
%template(norm_dot) ops::norm_dot<double>;
%template(inplace_trans) ops::inplace_trans<double>;
%template(max) ops::max<double>;
%template(col_max) ops::col_max<double>;
%template(row_max) ops::row_max<double>;
%template(min) ops::min<double>;
%template(col_min) ops::col_min<double>;
%template(col_max) ops::col_max<double>;
%template(normalize) ops::normalize<double>;
%template(prod) ops::prod<double>;
%template(col_prod) ops::col_prod<double>;
%template(row_prod) ops::row_prod<double>;
%template(powmat) ops::powmat<double>;
%template(repelem) ops::repelem<double>;
%template(repmat) ops::repmat<double>;
%template(reshape) ops::reshape<double>;
%template(reverse) ops::reverse<double>;
%template(shift) ops::shift<double>;
%template(shuffle) ops::shuffle<double>;
%template(size) ops::size<double>;
%template(sort) ops::sort<double>;
%template(symmatu) ops::symmatu<double>;
%template(trace) ops::trace<double>;
%template(trans) ops::trans<double>;
%template(trimatu) ops::trimatu<double>;
%template(vectorize) ops::vectorize<double>;

%template(abs) ops::abs<double>;
%template(exp) ops::exp<double>;
%template(exp2) ops::exp2<double>;
%template(exp10) ops::exp10<double>;
%template(expm1) ops::expm1<double>;
%template(trunc_exp) ops::trunc_exp<double>;
%template(log) ops::log<double>;
%template(log10) ops::log10<double>;
%template(log2) ops::log2<double>;
%template(log1p) ops::log1p<double>;
%template(trunc_log) ops::trunc_log<double>;
%template(pow) ops::pow<double>;
%template(square) ops::square<double>;
%template(sqrt) ops::sqrt<double>;
%template(floor) ops::floor<double>;
%template(ceil) ops::ceil<double>;
%template(round) ops::round<double>;
%template(trunc) ops::trunc<double>;
%template(erf) ops::erf<double>;
%template(erfc) ops::erfc<double>;
%template(tgamma) ops::tgamma<double>;
%template(lgamma) ops::lgamma<double>;
%template(sign) ops::sign<double>;
%template(sin) ops::sin<double>;
%template(cos) ops::cos<double>;
%template(tan) ops::tan<double>;
%template(asin) ops::asin<double>;
%template(acos) ops::acos<double>;
%template(atan) ops::atan<double>;
%template(sinh) ops::sinh<double>;
%template(cosh) ops::cosh<double>;
%template(tanh) ops::tanh<double>;
%template(asinh) ops::asinh<double>;
%template(acosh) ops::acosh<double>;
%template(atanh) ops::atanh<double>;
%template(atan2) ops::atan2<double>;
%template(hypot) ops::hypot<double>;
%template(sinc) ops::sinc<double>;

//%template(addToEachRow) ops::addToEachRow<double>;
//%template(addToEachCol) ops::addToEachCol<double>;

%template(set_real) ops::set_real<double>;
%template(set_imag) ops::set_imag<double>;

%template(st) ops::st<double>;
%template(is_hermitian) ops::is_hermitian<double>;


%template(real) ops::real<double>;
%template(imag) ops::imag<double>;
%template(log_det) ops::log_det<double>;
%template(log_det_sympd) ops::log_det_sympd<double>;
%template(logmat) ops::logmat<double>;
%template(logmat_sympd) ops::logmat_sympd<double>;
%template(sqrtmat) ops::sqrtmat<double>;
%template(sqrtmat_sympd) ops::sqrtmat_sympd<double>;

%template(chol) ops::chol<double>;
%template(eig_sym) ops::eig_sym<double>;
%template(eig_gen) ops::eig_gen<double>;
%template(eig_pair) ops::eig_pair<double>;
%template(hess) ops::hess<double>;
%template(inv) ops::inv<double>;
%template(inv_sympd) ops::inv_sympd<double>;
%template(lu) ops::lu<double>;
%template(null) ops::null<double>;
%template(orth) ops::orth<double>;
%template(pinv) ops::pinv<double>;
%template(qr) ops::qr<double>;
%template(qr_econ) ops::qr_econ<double>;
%template(schur) ops::schur<double>;
%template(solve) ops::solve<double>;
%template(svd) ops::svd<double>;
%template(conv) ops::conv<double>;
%template(conv2) ops::conv2<double>;
%template(fft) ops::fft<double>;
%template(ifft) ops::ifft<double>;
%template(fft2) ops::fft2<double>;
%template(ifft2) ops::ifft2<double>;
%template(interp1) ops::interp1<double>;
%template(interp2) ops::interp2<double>;
%template(polyfit) ops::polyfit<double>;
%template(polyval) ops::polyval<double>;
%template(mean) ops::mean<double>;
%template(median) ops::median<double>;
%template(stddev) ops::stddev<double>;
%template(var) ops::var<double>;
%template(range) ops::range<double>;
%template(cov) ops::cov<double>;
%template(cor) ops::cor<double>;
%template(hist) ops::hist<double>;
//%template(histc) ops::histc<double>;
%template(quantile) ops::quantile<double>;
%template(princomp) ops::princomp<double>;

    
