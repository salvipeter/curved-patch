#include "curved-domain.hh"

#include <sstream>

#define ANSI_DECLARATORS
#define REAL double
#define VOID void
extern "C" {
#include <triangle.h>
}

#include "lsq-plane.hh"

CurvedDomain::CurvedDomain() : mesh_updated(false) {
}

CurvedDomain::~CurvedDomain() {
}

namespace {

  inline Point3D to3D(const Point2D &p) {
    return { p[0], p[1], 0.0 };
  }
  
  void scalePoints(Point2DVector &points) {
    Point2D min = points[0], max = points[0];
    for (size_t i = 1; i < points.size(); ++i) {
      for (size_t j = 0; j < 2; ++j) {
        min[j] = std::min(min[j], points[i][j]);
        max[j] = std::max(max[j], points[i][j]);
      }
    }
    auto d = max - min;
    double margin = 0.025;
    double len = std::max(d[0], d[1]) * (1.0 + 2.0 * margin);
    for (auto &p : points) {
      p[0] = (p[0] - min[0]) / len + margin;
      p[1] = (p[1] - min[1]) / len + margin;
    }
  }
  
}

bool
CurvedDomain::update() {
  PointVector pv;
  for (const auto &c : curves_) {
    const auto &cp = c->controlPoints();
    for (size_t i = 1; i < cp.size(); ++i)
      pv.push_back(cp[i]);
  }
  auto projected = LSQPlane::projectToBestFitPlane(pv);
  scalePoints(projected);

  plane_curves_.clear();
  int index = -1;
  for (const auto &c : curves_) {
    plane_curves_.push_back(*c);
    auto &cp = plane_curves_.back().controlPoints();
    cp[0] = to3D(index >= 0 ? projected[index] : projected.back());
    for (size_t i = 1; i < cp.size(); ++i)
      cp[i] = to3D(projected[++index]);
  }
  n_ = plane_curves_.size();

  mesh_updated = false;

  return true;
}

void
CurvedDomain::updateMesh(size_t resolution) {
  DoubleVector points;
  for (const auto &c : plane_curves_) {
    for (size_t i = 0; i < resolution; ++i) {
      double u = (double)i / resolution;
      auto p = c.eval(u);
      points.push_back(p[0]);
      points.push_back(p[1]);
    }
  }
  std::vector<int> segments;
  int n = points.size() / 2;
  for (int i = 0; i < n; ++i) {
    segments.push_back(i);
    segments.push_back(i + 1);
  }
  segments.back() = 0;

  // Setup output data structure
  struct triangulateio in, out;
  in.pointlist = &points[0];
  in.numberofpoints = n;
  in.numberofpointattributes = 0;
  in.pointmarkerlist = nullptr;
  in.segmentlist = &segments[0];
  in.numberofsegments = n;
  in.segmentmarkerlist = nullptr;
  in.numberofholes = 0;
  in.numberofregions = 0;

  // Setup output data structure
  out.pointlist = nullptr;
  out.pointattributelist = nullptr;
  out.pointmarkerlist = nullptr;
  out.trianglelist = nullptr;
  out.triangleattributelist = nullptr;
  out.segmentlist = nullptr;
  out.segmentmarkerlist = nullptr;

  double max_area = 0.0;
  for (const auto &c : plane_curves_)
    max_area = std::max(max_area, c.arcLength(0.0, 1.0));
  max_area /= resolution * 2;
  max_area *= max_area * std::sqrt(3.0) / 4.0;
  std::stringstream cmd;
  cmd << "pqa" << std::fixed << max_area << "DBPzQ";
  triangulate(const_cast<char *>(cmd.str().c_str()), &in, &out, (struct triangulateio *)nullptr);

  for (int i = 0; i < out.numberofpoints; ++i)
    parameters_.emplace_back(out.pointlist[2*i], out.pointlist[2*i+1]);
  mesh_.resizePoints(parameters_.size());
  for (int i = 0; i < out.numberoftriangles; ++i)
    mesh_.addTriangle(out.trianglelist[3*i+2],
                      out.trianglelist[3*i+1],
                      out.trianglelist[3*i+0]);

  mesh_updated = true;

  // Domain test:
  if (false) {
    PointVector pv;
    for (const auto &p : parameters_)
      pv.emplace_back(p[0], p[1], 0.0);
    mesh_.setPoints(pv);
    mesh_.writeOBJ("/tmp/domain.obj");
  }
}

const Point2DVector &
CurvedDomain::parameters(size_t resolution) const {
  if (!mesh_updated)
    const_cast<CurvedDomain *>(this)->updateMesh(resolution);
  return parameters_;
}

TriMesh
CurvedDomain::meshTopology(size_t resolution) const {
  if (!mesh_updated)
    const_cast<CurvedDomain *>(this)->updateMesh(resolution);
  return mesh_;
}

const std::vector<BSCurve> &
CurvedDomain::boundaries() const {
  return plane_curves_;
}

