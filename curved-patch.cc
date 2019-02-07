#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <domain.hh>
#include <parameterization.hh>
#include <ribbon.hh>
#include <surface-corner-based.hh>
#include <surface-generalized-coons.hh>

#include "curved-cb.hh"
#include "curved-cr.hh"
#include "curved-gc.hh"
#include "perpendicular-cb.hh"

CurveVector readLOP(std::string filename) {
  size_t n, deg, nk, nc;
  DoubleVector knots;
  PointVector cpts;
  CurveVector result;
  std::ifstream f(filename);
  f >> n;
  result.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    f >> deg;
    f >> nk;
    knots.resize(nk);
    for (size_t j = 0; j < nk; ++j)
      f >> knots[j];
    f >> nc;
    cpts.resize(nc);
    for (size_t j = 0; j < nc; ++j)
      f >> cpts[j][0] >> cpts[j][1] >> cpts[j][2];
    result.push_back(std::make_shared<BSCurve>(deg, knots, cpts));
  }
  if (!f)
    return CurveVector();
  return result;
}

void ribbonTest(const std::shared_ptr<Surface> &surf, size_t resolution, std::string filename) {
  double ribbon_length = 0.25;

  size_t n = surf->domain()->size();
  TriMesh ribbon_mesh;
  PointVector pv; pv.reserve(n * (resolution + 1) * 2);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j <= resolution; ++j) {
      double u = (double)j / resolution;
      pv.push_back(surf->ribbon(i)->curve()->eval(u));
      pv.push_back(surf->ribbon(i)->eval(Point2D(u, ribbon_length)));
    }
  }
  ribbon_mesh.setPoints(pv);
  size_t index = 0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < resolution; ++j) {
      ribbon_mesh.addTriangle(index, index+1, index+2);
      ++index;
      ribbon_mesh.addTriangle(index, index+2, index+1);
      ++index;
    }
    index += 2;
  }
  ribbon_mesh.writeOBJ(filename);  
}

void fixMesh(TriMesh &mesh, const CurveVector &cv, size_t resolution) {
  size_t index = 0;
  for (const auto &c : cv)
    for (size_t i = 0; i < resolution; ++i) {
      double u = (double)i / resolution;
      mesh[index++] = c->eval(u);
    }
}

std::vector<std::pair<Point3D, Point3D>>
slicer(const TriMesh &mesh, const std::vector<Point2DVector> &params) {
  const double density = 0.1;
  const size_t nr_lines = 5;
  std::vector<std::pair<Point3D, Point3D>> result;
  for (size_t param = 0; param < params[0].size(); ++param) {
    auto slice = [&](size_t i, size_t j, Point3D &p) {
                   double x = params[i][param][1], y = params[j][param][1]; // d-coordinates
                   if (y > x) {
                     std::swap(x, y);
                     std::swap(i, j);
                   }
                   size_t q1 = x / density, q2 = y / density;
                   if (q1 - q2 == 1 && q1 <= nr_lines) {
                     double alpha = (density * q1 - y) / (x - y);
                     p = mesh[j] * (1 - alpha) + mesh[i] * alpha;
                     return true;
                   }
                   return false;
                 };
    for (const auto &tri : mesh.triangles()) {
      Point3D p;
      PointVector ps;
      if (slice(tri[0], tri[1], p))
        ps.push_back(p);
      if (slice(tri[0], tri[2], p))
        ps.push_back(p);
      if (slice(tri[1], tri[2], p))
        ps.push_back(p);
      if (ps.size() == 2)
        result.push_back({ ps[0], ps[1] });
    }
  }
  return result;
}

void
domainEval3D(const std::shared_ptr<Surface> &surf, size_t resolution, std::string filename) {
  auto mesh = surf->domain()->meshTopology(resolution);
  auto uvs = surf->domain()->parameters(resolution);
  std::vector<Point2DVector> params; params.reserve(uvs.size());
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(params),
                 [&](const Point2D &uv) { return surf->parameterization()->mapToRibbons(uv); });
  PointVector points; points.reserve(uvs.size());
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(points),
                 [&](const Point2D &uv) { return surf->eval(uv); });
  mesh.setPoints(points);
  auto segments = slicer(mesh, params);
  std::ofstream f(filename);
  if (!f.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return;
  }
  for (const auto &s : segments) {
    f << "v " << s.first[0] << ' ' << s.first[1] << ' ' << s.first[2] << std::endl;
    f << "v " << s.second[0] << ' ' << s.second[1] << ' ' << s.second[2] << std::endl;
  }
  for (size_t i = 1; i <= segments.size(); ++i)
    f << "l " << 2*i-1 << ' ' << 2*i << std::endl;
  f.close();
}

void
domainEval(const std::shared_ptr<Surface> &surf, size_t resolution, std::string filename) {
  auto mesh = surf->domain()->meshTopology(resolution);
  auto uvs = surf->domain()->parameters(resolution);
  std::vector<Point2DVector> points; points.reserve(uvs.size());
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(points),
                 [&](const Point2D &uv) {
                   auto sds = surf->parameterization()->mapToRibbons(uv);
                   sds.push_back(uv);
                   return sds;
                 });
  std::ofstream f(filename);
  if (!f.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return;
  }
  for (const auto &p : points) {
    f << 'v';
    for (auto coord : p)
      f << ' ' << coord[0] << ' ' << coord[1];
    f << std::endl;
  }
  for (const auto &t : mesh.triangles())
    f << "f " << t[0] + 1 << ' ' << t[1] + 1 << ' ' << t[2] + 1 << std::endl;
  f.close();
}

void surfaceTest(std::string name, std::shared_ptr<Surface> &&surf, const CurveVector &cv,
                 std::string filename, size_t resolution, bool fix_mesh = false) {
  std::cout << name << ":" << std::endl;
  std::chrono::steady_clock::time_point begin, end;

  begin = std::chrono::steady_clock::now();
  surf->setCurves(cv);
  surf->setupLoop();
  surf->update();
  end = std::chrono::steady_clock::now();
  std::cout << "  Setup time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << "ms" << std::endl;

  if (name == "CCB") {
    begin = std::chrono::steady_clock::now();
    ribbonTest(surf, resolution, filename + "-ribbons.obj");
    end = std::chrono::steady_clock::now();
    std::cout << "  Ribbon output time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
              << "ms" << std::endl;

    begin = std::chrono::steady_clock::now();
    domainEval(surf, resolution, filename + "-domain.obj");
    domainEval3D(surf, resolution, filename + "-domain3D.obj");
    end = std::chrono::steady_clock::now();
    std::cout << "  Domain output time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
              << "ms" << std::endl;
  }
  
  begin = std::chrono::steady_clock::now();
  auto mesh = surf->eval(resolution);
  end = std::chrono::steady_clock::now();
  std::cout << "  Evaluation time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << "ms" << std::endl;

  if (fix_mesh)
    fixMesh(mesh, cv, resolution); // computes exact boundaries
  mesh.writeOBJ(filename + "-" + name + ".obj");
}

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0]
              << " basename [resolution]" << std::endl;
    return 1;
  }
  std::string fname(argv[1]);

  CurveVector cv = readLOP(fname + ".lop");
  if (cv.empty()) {
    std::cerr << "Cannot read file: " << argv[1] << std::endl;
    return 2;
  }
  
  size_t resolution = 30;
  if (argc == 3)
    resolution = std::atoi(argv[2]);

  // surfaceTest("CGC", std::make_shared<CurvedGC>(), cv, fname, resolution, true);
  surfaceTest("CCB", std::make_shared<CurvedCB>(), cv, fname, resolution, true);
  // surfaceTest("CCR", std::make_shared<CurvedCR>(), cv, fname, resolution, false);
  // surfaceTest("GC", std::make_shared<Transfinite::SurfaceGeneralizedCoons>(),
  //             cv, fname, resolution);
  // surfaceTest("CB", std::make_shared<Transfinite::SurfaceCornerBased>(),
  //             cv, fname, resolution);
  surfaceTest("PCB", std::make_shared<PerpCB>(), cv, fname, resolution);

  return 0;
}
