
#include "../core/algorithms.hpp"
#include "../core/containers.hpp"

#include "solve_connection.hpp"

namespace panoramix {
    namespace experimental {


        //std::vector<double> ProjectiveObject::depths(const Vec3 & direction) const{
        //    auto coeffs = coefficients(direction);
        //    std::vector<double> dps(coeffs.size(), 0.0);
        //    for (int i = 0; i < dps.size(); i++){
        //        assert(coeffs[i].size() == _params.size());
        //        for (int j = 0; j < _params.size(); j++){
        //            dps[i] += coeffs[i][j] * _params[j];
        //        }
        //        dps[i] = 1.0 / dps[i];
        //    }
        //    return dps;
        //}


        LineDoF1::LineDoF1(const Line3 & l) : ProjectiveComponent(1), line(l) { }


        std::vector<std::vector<double>> LineDoF1::coefficients(const Vec3 & direction) const {
            Ray3 infLine = (line / norm(line.center())).infiniteLine();
            // variable is 1.0/centerDepth
            // corresponding coeff is 1.0/depthRatio
            // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
            double depthRatio = norm(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), direction), infLine).second.first);
            return{ std::vector<double>{1.0 / depthRatio} };
        }

        LineDoF2::LineDoF2(const Line3 & l) : ProjectiveComponent(2), line(l){ }


        std::vector<std::vector<double>> LineDoF2::coefficients(const Vec3 & direction) const {
            double theta = AngleBetweenDirections(normalize(line.first), normalize(line.second));
            double phi = AngleBetweenDirections(normalize(line.first), direction);
            /*           | sin(theta) | | p | | q |
                len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
                | p sin(phi) - q sin(phi - theta) |
                */
            // variables[0] -> 1/p
            // variables[1] -> 1/q
            double coeffFor1_p = -sin(phi - theta) / sin(theta);
            double coeffFor1_q = sin(phi) / sin(theta);
            assert(!IsInfOrNaN(coeffFor1_p) && !IsInfOrNaN(coeffFor1_q));
            return{ std::vector<double>{coeffFor1_p, coeffFor1_q} };
        }

        RegionDoF1::RegionDoF1(const Plane3 & p) : ProjectiveComponent(1), plane(p) {}


        std::vector<std::vector<double>> RegionDoF1::coefficients(const Vec3 & direction) const {
            Plane3 plane(normalize(plane.anchor), plane.normal);
            // variable is 1.0/centerDepth
            // corresponding coeff is 1.0/depthRatio
            // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
            double depthRatio = norm(IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), direction), plane).position);
            return{ std::vector<double>{1.0 / depthRatio} };
        }


        RegionDoF2::RegionDoF2(const Vec3 & ad) : ProjectiveComponent(2), axisDirection(ad) { }


        int SwappedComponent(const Vec3 & orientation) {
            for (int i = 0; i < 2; i++){
                if (abs(orientation[i]) >= 1e-8){
                    return i;
                }
            }
            return 2;
        }


        std::vector<std::vector<double>> RegionDoF2::coefficients(const Vec3 & direction) const {
            // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
            // -> 1.0/depth = ax + by + cz
            // -> 1.0/depth = ax + by + (- o1a - o2b)/o3 * z = a(x-o1*z/o3) + b(y-o2*z/o3)
            auto orientation = normalize(axisDirection);
            int sc = SwappedComponent(orientation);
            Vec3 forientation = orientation;
            std::swap(forientation[sc], forientation[2]);
            Vec3 fdirection = direction;
            std::swap(fdirection[sc], fdirection[2]);
            return{ std::vector<double>{
                fdirection[0] - forientation[0] * fdirection[2] / forientation[2],
                    fdirection[1] - forientation[1] * fdirection[2] / forientation[2]
            } };
        }

        RegionDoF3::RegionDoF3() : ProjectiveComponent(3){ }

        std::vector<std::vector<double>> RegionDoF3::coefficients(const Vec3 & direction) const {
            // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
            // -> 1.0/depth = ax + by + cz
            return{ std::vector<double>{direction[0], direction[1], direction[2]} };
        }



        using Mesh3 = Mesh<Point3, Line3, Polygon3>;

        class MeshDoF1Data {
        public:
            explicit MeshDoF1Data(const Mesh<Point3> & m) {
                // build rtree
                mesh.internalVertices().reserve(m.internalVertices().size());
                mesh.internalHalfEdges().reserve(m.internalHalfEdges().size());
                mesh.internalFaces().reserve(m.internalFaces().size());
                for (auto & v : m.internalVertices()){
                    vertRTree.insert(BoundingBox(v.data).expand(0.1), mesh.addVertex(v.data));
                }
                for (auto & h : m.internalHalfEdges()){
                    Mesh3::VertHandle vh1 = h.topo.from().id;
                    Mesh3::VertHandle vh2 = h.topo.to().id;
                    Line3 line(m.data(h.topo.from()), m.data(h.topo.to()));
                    halfRTree.insert(BoundingBox(line), mesh.addEdge(vh1, vh2, line, line.reversed(), false));
                }
                for (auto & f : m.internalFaces()){
                    std::vector<Mesh3::HalfHandle> hhs;
                    hhs.reserve(f.topo.halfedges.size());
                    std::vector<Point3> points;
                    points.reserve(f.topo.halfedges.size());
                    for (auto hh : f.topo.halfedges){
                        hhs.emplace_back(hh.id);
                        points.push_back(mesh.data(mesh.topo(hhs.back()).to()));
                    }
                    assert(points.size() >= 3);
                    auto normal = normalize((points.at(0) - points.at(1)).cross(points.at(1) - points.at(2)));
                    // trianglulate
                    std::vector<std::array<Point3, 3>> triangles;
                    triangles.reserve(points.size() - 2);
                    Vec3 x, y;
                    std::tie(x, y) = ProposeXYDirectionsFromZDirection(normal);
                    TriangulatePolygon(points.begin(), points.end(), [&x, &y, &points](const Point3 & p){
                        auto pp = p - points.front();
                        return Point2(pp.dot(x), pp.dot(y));
                    }, [&triangles](const Point3 & a, const Point3 & b, const Point3 & c){
                        triangles.push_back(std::array<Point3, 3>{{ a, b, c }});
                    });
                    // polygon
                    Polygon3 poly(std::move(points), normal);
                    Mesh3::FaceHandle fh;
                    faceRTree.insert(BoundingBox(poly), fh = mesh.addFace(std::move(hhs), std::move(poly)));
                    faceTriangles[fh] = std::move(triangles);
                }
            }

            Mesh3 mesh;
            RTree<Box3, Mesh3::VertHandle> vertRTree;
            RTree<Box3, Mesh3::HalfHandle> halfRTree;
            RTree<Box3, Mesh3::FaceHandle> faceRTree;
            HandledTable<Mesh3::FaceHandle, std::vector<std::array<Point3, 3>>> faceTriangles;
        };

        MeshDoF1::MeshDoF1(const Mesh<Point3> & m) : ProjectiveComponent(1), _data(std::make_unique<MeshDoF1Data>(m)) {}


        template <class TT>
        bool TriangleIntersection(const Vec<TT, 3> &  V1,  // Triangle vertices
            const Vec<TT, 3> &  V2,
            const Vec<TT, 3> &  V3,
            const Vec<TT, 3> &   O,  //Ray origin
            const Vec<TT, 3> &   D,  //Ray direction

            TT* out,
            TT epsilon) {

            Vec<TT, 3> e1, e2;  //Edge1, Edge2
            Vec<TT, 3> P, Q, T;
            TT det, inv_det, u, v;
            TT t;

            //Find vectors for two edges sharing V1
            e1 = V2 - V1;
            e2 = V3 - V1;
            //Begin calculating determinant - also used to calculate u parameter
            P = D.cross(e2);
            //if determinant is near zero, ray lies in plane of triangle
            det = e1.dot(P);
            //NOT CULLING
            if (det > -epsilon && det < epsilon) return false;
            inv_det = 1.f / det;

            //calculate distance from V1 to ray origin
            T = O - V1;

            //Calculate u parameter and test bound
            u = T.dot(P) * inv_det;
            //The intersection lies outside of the triangle
            if (u < 0.f || u > 1.f) return false;

            //Prepare to test v parameter
            Q = T.cross(e1);

            //Calculate V parameter and test bound
            v = D.dot(Q) * inv_det;
            //The intersection lies outside of the triangle
            if (v < 0.f || u + v  > 1.f) return false;

            t = e2.dot(Q) * inv_det;

            if (t > epsilon) { //ray intersection
                *out = t;
                return true;
            }

            // No hit, no win
            return false;
        }

        std::vector<Point3> MeshDoF1::intersectionsOnPlacedMesh(const Vec3 & direction) const{
            Box3 bbox;
            for (auto & v : _data->mesh.vertices()){
                bbox |= BoundingBox(v.data);
            }
            auto bballAll = bbox.outerSphere();
            Ray3 centerRay(Point3(0, 0, 0), normalize(direction));

            // discretize the ray
            double startLen = std::max(0.0, norm(bballAll.center) - bballAll.radius);
            double stopLen = std::max(0.0, norm(bballAll.center) + bballAll.radius);
            double stepNum = 1000.0;
            double stepLen = (stopLen - startLen) / stepNum;

            double epsilon = 1e-20;

            std::vector<Point3> hitPoints;
            for (int i = 0; i <= stepNum; i++){
                double centerLen = startLen + stepLen * i;
                double nextCenterLen = centerLen + stepLen;

                Point3 centerP = centerRay.anchor + centerLen * centerRay.direction;
                Point3 nextCenterP = centerRay.anchor + nextCenterLen * centerRay.direction;
                Line3 lineSegment(centerP, nextCenterP);

                Box3 detectionBox = BoundingBox(centerP) | BoundingBox(nextCenterP);
                detectionBox.expand(stepLen);

                _data->faceRTree.search(detectionBox,
                    [&hitPoints, this, &centerRay, epsilon](const Mesh3::FaceHandle & fh){
                    const auto & triangles = _data->faceTriangles[fh];
                    for (auto & t : triangles){
                        double out = 0.0;
                        bool intersected = TriangleIntersection(
                            t[0], t[1], t[2],
                            centerRay.anchor, centerRay.direction,
                            &out, epsilon);
                        if (intersected){
                            hitPoints.push_back(centerRay.anchor + centerRay.direction * out);
                            return true;
                        }
                    }
                    return true;
                });
            }

            return hitPoints;
        }


        std::vector<std::vector<double>> MeshDoF1::coefficients(const Vec3 & direction) const {
            auto pts = intersectionsOnPlacedMesh(direction);
            std::vector<std::vector<double>> coeffs(pts.size(), std::vector<double>(1, 0.0));
            for (int i = 0; i < pts.size(); i++){
                coeffs[i][0] = 1.0 / norm(pts[i]);
            }
            return coeffs;
        }



    }
}