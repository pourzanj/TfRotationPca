#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(rstan)]]

#include "stan/math.hpp"
#include "stan/math/fwd/mat.hpp"
#include "linear_regression_model.hpp"
#include <stan/io/dump.hpp>
#include <iostream>
#include <fstream>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

stan::io::dump readData(std::string filename) {
  std::ifstream dumpFile(filename.c_str());
  
  if(!dumpFile.good()) {
    dumpFile.close();
    throw std::domain_error("Error opening dump file");
  }
  
  stan::io::dump dFile(dumpFile);
  
  dumpFile.close();
  return dFile;
}

static linear_regression_model_namespace::linear_regression_model *model = NULL;

// [[Rcpp::export]]
void set_data(std::string filename) {
  stan::io::dump dfile = readData(filename);
  if(model != NULL) {
    delete model;
  }
  model = new linear_regression_model_namespace::linear_regression_model(dfile, &Rcpp::Rcout);
}

// [[Rcpp::export]]
Rcpp::List jacobian(std::vector<double> params) {
  using namespace Rcpp;
  using stan::math::var;
  using stan::math::fvar;
  
  if(model == NULL) {
    throw std::invalid_argument("Must call set_data before jacobian");
  }
  
  NumericVector jac(params.size());
  
  std::vector<var> params_r;
  std::vector<int> params_i({});
  
  params_r.insert(params_r.begin(), params.begin(), params.end());
  
  var lp = model->log_prob<true, true, var>(params_r, params_i, &Rcpp::Rcout);
  
  lp.grad();
  
  for(size_t i = 0; i < params_r.size(); i++) {
    jac(i) = params_r[i].adj();
  }
  
  List out;
  
  out["u"] = lp.val();
  out["jac"] = jac;
  
  stan::math::recover_memory();
  
  return out;
}

// [[Rcpp::export]]
Rcpp::List hessian(std::vector<double> params) {
  using namespace Rcpp;
  using stan::math::var;
  using stan::math::fvar;
  
  if(model == NULL) {
    throw std::invalid_argument("Must call set_data before jacobian");
  }
  
  NumericVector jac(params.size());
  NumericMatrix hess(params.size(), params.size());
  
  std::vector<int> params_i({});
  
  double lp_ = 0.0;
  
  for(size_t i = 0; i < params.size(); i++) {
    std::vector<fvar<var> > params_r;
    for(auto v : params)
      params_r.push_back(v);
    
    params_r[i].d_ = 1.0;
    fvar<var> lp = model->log_prob<true, true, fvar<var> >(params_r, params_i, &Rcpp::Rcout);
    
    jac(i) = lp.tangent().val();
    
    lp.d_.grad();
    for(size_t j = 0; j < params_r.size(); j++) {
      hess(i, j) = params_r[j].val().adj();
    }
    
    lp_ = lp.val().val(); // Same every time
  }
  
  List out;
  
  out["u"] = lp_;
  out["jac"] = jac;
  out["hess"] = hess;
  
  stan::math::recover_memory();
  
  return out;
}

// [[Rcpp::export]]
Rcpp::List hessian_vector(std::vector<double> params, std::vector<double> vector) {
  using namespace Rcpp;
  using stan::math::var;
  using stan::math::fvar;
  
  if(model == NULL) {
    throw std::invalid_argument("Must call set_data before jacobian");
  }
  
  double jacv;
  NumericVector hessv(params.size());
  
  std::vector<int> params_i({});
  
  double lp_ = 0.0;
  
  std::vector<fvar<var> > params_r;
  for(int i = 0; i < params.size(); ++i) {
    params_r.emplace_back(params[i], vector[i]);
  }
  
  fvar<var> lp = model->log_prob<true, true, fvar<var> >(params_r, params_i, &Rcpp::Rcout);
  
  lp_ = lp.val().val(); // Same every time
  jacv = lp.tangent().val();
  
  lp.d_.grad();
  for(size_t j = 0; j < params_r.size(); j++) {
    hessv(j) = params_r[j].val().adj();
  }
  
  List out;
  
  out["u"] = lp_;
  out["jacv"] = jacv;
  out["hessv"] = hessv;
  
  stan::math::recover_memory();
  
  return out;
}

// The code for this was taken from https://eigen.tuxfamily.org/dox/group__MatrixfreeSolverExample.html
class HessianMatrixReplacement;

namespace Eigen {
namespace internal {
// HessianMatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
template<>
struct traits<HessianMatrixReplacement> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >{};
}
}

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class HessianMatrixReplacement : public Eigen::EigenBase<HessianMatrixReplacement> {
public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };
  Index rows() const { return Index(params_.size()); }
  Index cols() const { return Index(params_.size()); }
  
  template<typename Rhs>
  Eigen::Product<HessianMatrixReplacement,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<HessianMatrixReplacement,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }
  
  // Custom API:
  HessianMatrixReplacement(std::vector<double> params, std::vector<double> M, double h) : params_(params), M_(M), h_(h) {
  }
  
  const std::vector<double> &get_params() const { return params_; }
  const std::vector<double> &get_M() const { return M_; }
  double get_h() const { return h_; }
private:
  std::vector<double> params_;
  std::vector<double> M_;
  double h_;
};

// Implementation of HessianMatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {
template<typename Rhs>
struct generic_product_impl<HessianMatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<HessianMatrixReplacement,Rhs,generic_product_impl<HessianMatrixReplacement,Rhs> >
{
  typedef typename Product<HessianMatrixReplacement,Rhs>::Scalar Scalar;
  template<typename Dest>
  static void scaleAndAddTo(Dest& dst, const HessianMatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha)
  {
    // This method should implement "dst += alpha * lhs * rhs" inplace
    auto M = lhs.get_M();
    auto h = lhs.get_h();
    std::vector<double> rhs_vector(rhs.size());
    for(int i = 0; i < rhs.size(); i++)
      rhs_vector[i] = rhs(i);
    Rcpp::NumericVector out = hessian_vector(lhs.get_params(), rhs_vector)["hessv"];
    
    // hessian should be negative since we're working in term of U(q) = -log(p(q))
    for(int i = 0; i < out.size(); i++) {
      out(i) = -out(i);
    }
    
    // std::cout << "lhs.get_params():" << std::endl; for(int i = 0; i < lhs.get_params().size(); i++) std::cout << lhs.get_params()[i] << std::endl;
    // std::cout << "rhs_vector:" << std::endl; for(int i = 0; i < rhs.size(); i++) std::cout << rhs_vector[i] << std::endl;
    // std::cout << "out:" << std::endl; for(int i = 0; i < rhs.size(); i++) std::cout << out[i] << std::endl;
    
    for(int i = 0; i < out.size(); i++) {
      out(i) = rhs(i) + h * h * out(i) / (4.0 * M[i]);
    }
    // std::cout << "out:" << std::endl; for(int i = 0; i < rhs.size(); i++) std::cout << out[i] << std::endl;
    // std::cout << "~~~~~~~~~~~~~" << std::endl;
    for(int i = 0; i < dst.size(); i++)
      dst(i) += alpha * out(i);
  }
};
}
}

// [[Rcpp::export]]
Rcpp::List hessian_solve(std::vector<double> params, std::vector<double> rhs, std::vector<double> M, double h, std::vector<double> guess, double tolerance) {
  HessianMatrixReplacement A(params, M, h);
  Eigen::VectorXd rhs_eigen(rhs.size()),
  guess_eigen(guess.size());
  Rcpp::NumericVector x(rhs.size());
  Rcpp::List out;
  
  for(int i = 0; i < rhs.size(); i++)
    rhs_eigen(i) = rhs[i];
  
  for(int i = 0; i < rhs.size(); i++)
    guess_eigen(i) = guess[i];
  
  Eigen::GMRES<HessianMatrixReplacement, Eigen::IdentityPreconditioner> solver;
  
  solver.setTolerance(tolerance);
  solver.compute(A);
  
  Eigen::VectorXd x_eigen = solver.solveWithGuess(rhs_eigen, guess_eigen);
  
  for (int i = 0; i < x_eigen.size(); ++i)
    x(i) = x_eigen(i);
  
  out["x"] = x;
  out["iterations"] = solver.iterations();
  out["error"] = solver.error();
  
  return out;
}