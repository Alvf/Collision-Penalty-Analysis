#include "PCG.h"
#include <iostream>

using namespace std;

PCG::PCG(const SPARSE_MATRIX& A, const DIAGONAL& diagonal) : 
  _A(A),
  _diagonal(diagonal),
  //_tolerance(1e-9),
  _tolerance(1e-6),
  //_maxIterations(10000),
  _maxIterations(1000),
  _x0(A.rows())
{
  _x0.setZero();
}

PCG::~PCG()
{
}

VECTOR PCG::solve(const VECTOR& rhs)
{
  //return solveEigenStyle(rhs);
  return solveAMGCL(rhs);
}

VECTOR PCG::solveAMGCL(const VECTOR& rhs)
{
  const REAL tol = _tolerance;
  const int maxiter = _maxIterations;

  //static const coef_type one  = math::identity<coef_type>();
  //static const coef_type zero = math::zero<coef_type>();
  const REAL one = 1.0;
  const REAL zero = 0.0;

  //ios_saver ss(std::cout);

  //scalar_type norm_rhs = norm(rhs);
  REAL norm_rhs = rhs.norm();
  /*
  if (norm_rhs < amgcl::detail::eps<scalar_type>(1)) {
    if (prm.ns_search) {
      norm_rhs = math::identity<scalar_type>();
    } else {
      backend::clear(x);
      return std::make_tuple(0, norm_rhs);
    }
  }
  */

  //scalar_type eps = std::max(prm.tol * norm_rhs, prm.abstol);
  REAL eps = tol * norm_rhs;

  //coef_type rho1 = 2 * eps * one;
  REAL rho1 = 2 * eps * one;
  //coef_type rho2 = zero;
  REAL rho2 = zero;

  VECTOR x(rhs.size());
  x.setZero();

  //backend::residual(rhs, A, x, *r);
  VECTOR r = rhs - _A * x;
  REAL res_norm = r.norm();

  VECTOR p;

  size_t iter = 0;
  //for(; iter < prm.maxiter && math::norm(res_norm) > eps; ++iter) 
  for(; iter < maxiter && res_norm > eps; ++iter) 
  {
    //P.apply(*r, *s);
    VECTOR s = r;

    rho2 = rho1;
    //rho1 = inner_product(*r, *s);
    rho1 = r.dot(s);

    if (iter)
      //backend::axpby(one, *s, rho1 / rho2, *p);
      p = s + (rho1 / rho2) * p;
    else
      //backend::copy(*s, *p);
      p = s;

    //backend::spmv(one, A, *p, zero, *q);
    VECTOR q = _A * p;

    //coef_type alpha = rho1 / inner_product(*q, *p);
    REAL alpha = rho1 / q.dot(p);

    //backend::axpby( alpha, *p, one,  x);
    x += alpha * p;
    //backend::axpby(-alpha, *q, one, *r);
    r -= alpha * q;

    //res_norm = norm(*r);
    res_norm = r.norm();
    //if (prm.verbose && iter % 5 == 0)
    //if (prm.verbose)
    //std::cout << "iteration: " << iter << "\t" << std::scientific << res_norm / norm_rhs << std::endl;
    std::cout << "iteration: " << iter << "\t" << std::scientific << res_norm << std::endl;
  }
  _iterations = iter;
  _error = res_norm;

  //return std::make_tuple(iter, res_norm / norm_rhs);
  return x;
}

VECTOR PCG::solveEigenStyle(const VECTOR& rhs)
{
  VECTOR x = _x0;
  const REAL rhsNorm2 = rhs.squaredNorm();
  _iterations = 0;

#if 0
  VECTOR p = _diagonal.apply(rhs);
  VECTOR Ap = _A * p;
  VECTOR residual = rhs - Ap;
  REAL residualNorm2 = residual.squaredNorm();
  const REAL considerAsZero = (std::numeric_limits<REAL>::min)();
  const REAL threshold = Eigen::numext::maxi(_tolerance * _tolerance * rhsNorm2, considerAsZero);
  if (residualNorm2 < threshold)
  {
    _iterations = 0;
    return p;
  }
#else
  VECTOR residual = rhs - _A * x;
  
  const REAL considerAsZero = (std::numeric_limits<REAL>::min)();
  const REAL threshold = Eigen::numext::maxi(_tolerance * _tolerance * rhsNorm2, considerAsZero);
  REAL residualNorm2 = residual.squaredNorm();

  if (residualNorm2 < threshold)
  {
    _iterations = 0;
    return x;
  }

  VECTOR p = _diagonal.apply(residual);
#endif

  bool verbose = true;
  if (verbose)
  {
    /*
    std::cout << " Initial guess norm: " << x.norm() << std::endl;
    std::cout << " A norm: " << _A.norm() << std::endl;
    std::cout << " Residual 0: " << residual.norm() << std::endl;
    std::cout << " rhs norm: " << rhs.norm() << endl;
    std::cout << " p norm: " << p.norm() << std::endl;
    */

    /*
    std::cout << " Residual after: " << (rhs - _A * p).norm() << endl;

    //Eigen::SimplicialLDLT<SPARSE_MATRIX> solver(_A);
    Eigen::SimplicialLDLT<SPARSE_MATRIX> solver;
    solver.compute(_A);
    VECTOR solution = solver.solve(rhs);

    //std::cout << " Residual here: " << (rhs - _A * solution).norm() << endl;
    //cout << " Solution here: " << solution.transpose() << endl;
    //exit(0);
    */
  }

  VECTOR z = p;
  VECTOR tmp(p.size());
  REAL rdotzNew = Eigen::numext::real(residual.dot(p));  // the square of the absolute value of r scaled by invM
  while (_iterations < _maxIterations)
  {
    tmp.noalias() = _A * p;

    REAL alpha = rdotzNew/p.dot(tmp);

    x += alpha * p;

    residual -= alpha * tmp; 
    //residual = rhs - _A * x;

    residualNorm2 = residual.squaredNorm();

    if (verbose)
    {
      //std::cout << " Ap norm: " << tmp.norm();
      //std::cout << " alpha: " << alpha;
      std::cout << " PCG residual " << _iterations << ": " << residual.norm() / rhs.norm() << std::endl;
      //std::cout << " residual " << _iterations << ": " << residualNorm2 << std::endl;
      //std::cout << " residual: " << residual.transpose() << std::endl;
    }

    if (residualNorm2 < threshold)
      break;

    z = _diagonal.apply(residual);
    //z = residual;
    REAL rdotzOld = rdotzNew;
    rdotzNew = Eigen::numext::real(residual.dot(z));
    REAL beta = rdotzNew/rdotzOld;

    p = z + beta * p;

    _iterations++;
  }

  _error = sqrt(residualNorm2 / rhsNorm2);

  return x;
}

VECTOR PCG::solvePCR(const VECTOR& rhs)
{
  bool verbose = false;
  VECTOR x = _x0;
  const REAL rhsNorm2 = rhs.squaredNorm();
  _iterations = 0;

  VECTOR Ax = _A * x;
  VECTOR r = _diagonal.apply(rhs - Ax);
  
  const REAL considerAsZero = (std::numeric_limits<REAL>::min)();
  const REAL threshold = Eigen::numext::maxi(_tolerance * _tolerance * rhs.norm(), considerAsZero);

  //VECTOR p = _diagonal.apply(r);
  VECTOR p = r;

  VECTOR z(p.size());
  //REAL absNew = Eigen::numext::real(residual.dot(p));  // the square of the absolute value of r scaled by invM
  VECTOR Ar = _A * r;
  VECTOR Ap = _A * p;
  REAL rAr = r.dot(Ar);
  while (_iterations < _maxIterations)
  {
    VECTOR MAp = _diagonal.apply(_A * p);
    REAL alpha = rAr/Ap.dot(MAp);

    x += p * alpha;

    r -= alpha * MAp;

    if (verbose)
    {
      //std::cout << " Ap norm: " << Ap.norm() << std::endl;
      //std::cout << " alpha: " << alpha << std::endl;
      std::cout << " residual " << _iterations << ": " << r.norm() / rhs.norm() << std::endl;
      //std::cout << " residual " << _iterations << ": " << residualNorm2 << std::endl;
    }

    //if (r.norm() / rhs.norm() < threshold)
    if (r.squaredNorm() < threshold)
      break;

    REAL rArOld= rAr;
    Ar = _A * r;

    rAr = r.dot(Ar);
    REAL beta = rAr/rArOld;

    p = r + beta * p;

    Ap = Ar + beta * Ap;

    _iterations++;
  }

  _error = r.norm() / rhs.norm();

  return x;
}

VECTOR PCG::solveCR(const VECTOR& rhs)
{
  bool verbose = false;
  VECTOR x = _x0;
  const REAL rhsNorm2 = rhs.squaredNorm();
  _iterations = 0;

  VECTOR Ax = _A * x;
  VECTOR r = rhs - Ax;
  
  const REAL considerAsZero = (std::numeric_limits<REAL>::min)();
  const REAL threshold = Eigen::numext::maxi(_tolerance * _tolerance * rhs.norm(), considerAsZero);

  //VECTOR p = _diagonal.apply(r);
  VECTOR p = r;

  VECTOR z(p.size());
  //REAL absNew = Eigen::numext::real(residual.dot(p));  // the square of the absolute value of r scaled by invM
  VECTOR Ar = _A * r;
  VECTOR Ap = _A * p;
  REAL rAr = r.dot(Ar);
  while (_iterations < _maxIterations)
  {
    REAL alpha = rAr/Ap.dot(Ap);

    x += p * alpha;

    r -= alpha * Ap; 

    if (verbose)
    {
      //std::cout << " Ap norm: " << Ap.norm() << std::endl;
      //std::cout << " alpha: " << alpha << std::endl;
      std::cout << " residual " << _iterations << ": " << r.norm() / rhs.norm() << std::endl;
      //std::cout << " residual " << _iterations << ": " << residualNorm2 << std::endl;
    }

    if (r.norm() / rhs.norm() < threshold)
      break;

    REAL rArOld= rAr;
    Ar = _A * r;

    rAr = r.dot(Ar);
    REAL beta = rAr/rArOld;

    p = r + beta * p;

    Ap = Ar + beta * Ap;

    _iterations++;
  }

  _error = r.norm() / rhs.norm();

  return x;
}

