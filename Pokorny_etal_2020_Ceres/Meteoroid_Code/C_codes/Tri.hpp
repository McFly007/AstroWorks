#ifndef __FAST_BVH_TRI_HPP__
#define __FAST_BVH_TRI_HPP__

#include "Object.h"

constexpr double TOL = 1e-6;

struct Tri: public Object {
  Vector3 v0, v1, v2, n, p;
  float r;
  int index;

  Tri(Vector3 const & v0, Vector3 const & v1, Vector3 const & v2,
      Vector3 const & n, int index):
    v0(v0),
    v1(v1),
    v2(v2),
    n(n),
    p((v0 + v1 + v2)/3),
    r(std::max(length(v0 - p), std::max(length(v1 - p), length(v2 - p)))),
    index(index)
  {}

  virtual ~Tri() {}
  
  bool getIntersection(Ray const & ray, IntersectionInfo * info) const {
    auto const e1 = v1 - v0;
    auto const e2 = v2 - v0;
    auto const pvec = ray.d^e2;
    auto const det = e1*pvec;
    if (-TOL < det && det < TOL) {
      return false;
    }
    auto const inv_det = 1.0/det;
    auto const tvec = ray.o - v0;
    auto const u = inv_det*(tvec*pvec);
    if (u < 0 || u > 1) {
      return false;
    }
    auto const qvec = tvec^e1;
    auto const v = inv_det*(ray.d*qvec);
    if (v < 0 || u + v > 1) {
      return false;
    }
    auto const t = inv_det*(e2*qvec);
    if (t < 0) {
      return false;
    }
    info->object = this;
    info->t = t;
    return true;
  }

  // TODO: we could modify this so that triangles either interpolate
  // the normal based on where they're struck or whether they have a
  // single face normal supplied
  Vector3 getNormal(IntersectionInfo const &) const {
    return n;
  }

  void getFrenetFrame(Vector3 & t, Vector3 & n, Vector3 & b) const {
    t = normalize(v1 - v0);
    n = this->n;
    b = t^n;
  }

  float getBoundingRadius() const {
    return r;
  }

  BBox getBBox() const {
    return BBox {min(v0, min(v1, v2)), max(v0, max(v1, v2))};
  }

  Vector3 getCentroid() const {
    return p;
  }

  float getArea() const {
    return length((v1 - v0)^(v2 - v0))/2;
  }

  friend std::ostream &
  operator<<(std::ostream & o, Tri const & t) {
    return o << "Tri {" << std::endl
             << "  v0 = " << t.v0 << "," << std::endl
             << "  v1 = " << t.v1 << "," << std::endl
             << "  v2 = " << t.v2 << "," << std::endl
             << "  n = " << t.n << "," << std::endl
             << "  index = " << t.index << std::endl
             << "}" << std::endl;
  }
};

#endif // __FAST_BVH_TRI_HPP__
