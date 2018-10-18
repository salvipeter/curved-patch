#include "curved-mean.hh"

#include <algorithm>

#include <gsl/gsl_integration.h>

#include "curved-domain.hh"

CurvedMean::~CurvedMean() {
}

namespace {

  struct IntegralData {
    const std::vector<BSCurve> &curves;
    const Point2D &p;
    size_t i;
  };

  double integrand(double theta, void *data) {
    auto id = reinterpret_cast<IntegralData *>(data);
    const auto &curves = id->curves;
    const auto &p = id->p;
    size_t i = id->i;
    double result = 0.0;
    size_t n = curves.size();
    Vector2D d(cos(theta), sin(theta));
    for (size_t j = 0; j < curves.size(); ++j) {
      if (i == j)
        continue;
      auto intersections = curves[j].intersectWithPlane({p[0], p[1], 0.0}, {d[1], -d[0], 0.0});
      std::vector<std::pair<double, double>> pairs;
      for (double s : intersections) {
        auto q = curves[j].eval(s);
        auto deviation = Point2D(q[0], q[1]) - p;
        if (deviation * d < 0.0)
          continue;
        pairs.emplace_back(deviation.norm(), s);
      }
      std::sort(pairs.begin(), pairs.end());
      for (size_t k = 0; k < pairs.size(); ++k) {
        if (pairs[k].first == 0.0)
          return i == n ? 1.0 : 0.0;
        double sign = k % 2 ? -1.0 : 1.0;
        double fval = 1.0;
        if (i < n && j == (i + 1) % n)
          fval = pairs[k].second;
        else if (i < n && (j + 1) % n == i)
          fval = 1.0 - pairs[k].second;
        result += sign / pairs[k].first * fval;
      }
    }
    return result;
  }

}

Point2D
CurvedMean::mapToRibbon(size_t i, const Point2D &uv) const {
  // Ch. Dyken, M. S. Floater, Transfinite mean value interpolation. CAGD 26(1), pp. 117-134, 2009.
  /*
\[\frac{1}{\phi(\mathbf{x})}\int_0^{2\pi}\sum_{j=1}^{n(\mathbf{x},\theta)}
\frac{(-1)^{j-1}}{\rho_j(\mathbf{x},\theta)}f(p_j(\mathbf{x},\theta))\,d\theta\]
where
\begin{itemize}
\item $\mathbf{x}$ is the point to evaluate at,
\item $p_j(\mathbf{x},\theta)$ is the $j$-th intersection with the boundary
      at search direction $\theta$,
\item $n(\mathbf{x},\theta)$ is the number of such intersections,
\item $\rho_j(\mathbf{x},\theta)$ is the distance of $\mathbf{x}$ to $p_j(\mathbf{x},\theta)$,
\item $f$ is the boundary function (0 on the current side, 1 on the distant sides, 
      linear on the adjacent sides) and
\end{itemize}
\[\phi(\mathbf{x})=\int_0^{2\pi}\sum_{j=1}^{n(\mathbf{x},\theta)}
\frac{(-1)^{j-1}}{\rho_j(\mathbf{x},\theta)}\,d\theta.\]
  */
  const auto &curves = dynamic_cast<CurvedDomain *>(domain_.get())->boundaries();
  gsl_function numerator, denominator;
  IntegralData id = {curves, uv, i};
  numerator.function   = &integrand; numerator.params   = &id;
  denominator.function = &integrand; denominator.params = &id;
  double num_int, denom_int, err;
  size_t neval;
  gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(100);
  gsl_integration_cquad(&numerator, 0.0, 2.0 * M_PI, 1.0e-5, 1.0e-3, w, &num_int, &err, &neval);
  id.i = curves.size();
  gsl_integration_cquad(&denominator, 0.0, 2.0 * M_PI, 1.0e-5, 1.0e-3, w, &denom_int, &err, &neval);
  gsl_integration_cquad_workspace_free(w);
  return Point2D(0.0, num_int / denom_int);         // dummy s parameter (!)
}

void
CurvedMean::update() {
  n_ = dynamic_cast<CurvedDomain *>(domain_.get())->boundaries().size();
}
