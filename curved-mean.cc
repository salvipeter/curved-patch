#include "curved-mean.hh"

#include <gsl/gsl_integration.h>

CurvedMean::~CurvedMean() {
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
  double d = 42;
  return Point2D(0, d);         // dummy s parameter (!)
}
