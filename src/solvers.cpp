// solvers.cpp -- compiled inner solver for ppsynth, mirroring fast_core.py
// (the numba engine behind the Stata command): simplex projection, FISTA for
// the penalized simplex-constrained quadratic, block-coordinate descent over
// outcomes, and the grid sweep with exact corner solutions and warm starts.
//
// Solver settings are passed in from R and are the same constants as the
// pure-R engine (R/engine.R) and the Python engines: corner solves 8000
// iterations / tol 1e-12, outer BCD loop 300 / 1e-11, inner FISTA 800 / 1e-9.
// No linear-algebra library is used: with N donors (tens) the matrix-vector
// products are explicit loops, which keeps the arithmetic order fixed and
// the dependency footprint at Rcpp alone.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

// Euclidean projection of v onto the probability simplex (in place result).
static void project_simplex(const std::vector<double>& v, std::vector<double>& out,
                            std::vector<double>& u) {
  const int n = (int) v.size();
  u.assign(v.begin(), v.end());
  std::sort(u.begin(), u.end(), std::greater<double>());
  double css = 0.0, theta = 0.0;
  for (int j = 0; j < n; ++j) {
    css += u[j];
    const double t = (css - 1.0) / (j + 1);
    if (u[j] - t > 0.0) theta = t;
  }
  for (int i = 0; i < n; ++i) {
    const double d = v[i] - theta;
    out[i] = d > 0.0 ? d : 0.0;
  }
}

// FISTA for  min_x x'Ax - 2 b'x + lam ||x - center||^2  on the simplex.
// A is N x N (column-major, as R stores it). Returns the iterate in x.
static void fista(const double* A, const double* b, int N, double lam,
                  const double* center, const double* x0, double Lmax,
                  int n_iter, double tol, std::vector<double>& x) {
  const double step = 1.0 / (2.0 * Lmax + 2.0 * lam + 1e-12);
  std::vector<double> y(x0, x0 + N), xn(N), g(N), z(N), u(N);
  x.assign(x0, x0 + N);
  double t = 1.0;
  for (int it = 0; it < n_iter; ++it) {
    for (int i = 0; i < N; ++i) {
      double s = 0.0;
      for (int j = 0; j < N; ++j) s += A[i + (size_t) j * N] * y[j];
      g[i] = 2.0 * (s - b[i]);
    }
    if (lam > 0.0)
      for (int i = 0; i < N; ++i) g[i] += 2.0 * lam * (y[i] - center[i]);
    for (int i = 0; i < N; ++i) z[i] = y[i] - step * g[i];
    project_simplex(z, xn, u);
    const double tn = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * t * t));
    double mv = 0.0;
    for (int i = 0; i < N; ++i) {
      y[i] = xn[i] + ((t - 1.0) / tn) * (xn[i] - x[i]);
      const double d = std::fabs(xn[i] - x[i]);
      if (d > mv) mv = d;
    }
    x.swap(xn);
    t = tn;
    if (mv < tol) break;
  }
}

// Block-coordinate descent over K outcomes at a fixed interior lambda.
// gam is K x N (row k = outcome k), updated in place from its warm start.
static void bcd(const std::vector<const double*>& A, const std::vector<const double*>& b,
                const std::vector<double>& L, int K, int N, double lam,
                std::vector< std::vector<double> >& gam, int n_outer, int n_inner, double tol) {
  std::vector<double> gbar(N), gn(N);
  for (int i = 0; i < N; ++i) {
    double s = 0.0;
    for (int k = 0; k < K; ++k) s += gam[k][i];
    gbar[i] = s / K;
  }
  for (int it = 0; it < n_outer; ++it) {
    double move = 0.0;
    for (int k = 0; k < K; ++k) {
      fista(A[k], b[k], N, lam, gbar.data(), gam[k].data(), L[k], n_inner, 1e-9, gn);
      for (int i = 0; i < N; ++i) {
        const double d = std::fabs(gn[i] - gam[k][i]);
        if (d > move) move = d;
      }
      gam[k] = gn;
    }
    for (int i = 0; i < N; ++i) {
      double s = 0.0;
      for (int k = 0; k < K; ++k) s += gam[k][i];
      gbar[i] = s / K;
    }
    if (move < tol) break;
  }
}

//' @keywords internal
//' @noRd
// [[Rcpp::export]]
NumericVector fit_all_lams_cpp(List A_list, List b_list, NumericVector L,
                               NumericMatrix Acat, NumericVector bcat, double Lcat,
                               NumericVector lam_grid,
                               int corner_iter = 8000, double corner_tol = 1e-12,
                               int n_outer = 300, double outer_tol = 1e-11,
                               int n_inner = 800) {
  const int K = A_list.size();
  const int M = lam_grid.size();
  std::vector<NumericMatrix> Amats(K);
  std::vector<NumericVector> bvecs(K);
  std::vector<const double*> A(K), b(K);
  for (int k = 0; k < K; ++k) {
    Amats[k] = as<NumericMatrix>(A_list[k]);
    bvecs[k] = as<NumericVector>(b_list[k]);
    A[k] = Amats[k].begin();
    b[k] = bvecs[k].begin();
  }
  const int N = Amats[0].nrow();
  std::vector<double> Lv(L.begin(), L.end());
  std::vector<double> unif(N, 1.0 / N), zero(N, 0.0);

  // exact corners: separate fits (lam = 0) and concatenated fit (lam = Inf)
  std::vector< std::vector<double> > sep(K, std::vector<double>(N));
  for (int k = 0; k < K; ++k)
    fista(A[k], b[k], N, 0.0, zero.data(), unif.data(), Lv[k], corner_iter, corner_tol, sep[k]);
  std::vector<double> gcat(N);
  fista(Acat.begin(), bcat.begin(), N, 0.0, zero.data(), unif.data(), Lcat,
        corner_iter, corner_tol, gcat);

  NumericVector out(Dimension(M, K, N));
  std::vector< std::vector<double> > warm = sep;
  for (int j = 0; j < M; ++j) {
    const double lam = lam_grid[j];
    if (lam == 0.0) {
      for (int k = 0; k < K; ++k)
        for (int i = 0; i < N; ++i) out[j + (size_t) M * (k + (size_t) K * i)] = sep[k][i];
    } else if (std::isinf(lam)) {
      for (int k = 0; k < K; ++k)
        for (int i = 0; i < N; ++i) out[j + (size_t) M * (k + (size_t) K * i)] = gcat[i];
    } else {
      bcd(A, b, Lv, K, N, lam, warm, n_outer, n_inner, outer_tol);
      for (int k = 0; k < K; ++k)
        for (int i = 0; i < N; ++i) out[j + (size_t) M * (k + (size_t) K * i)] = warm[k][i];
    }
  }
  return out;
}
