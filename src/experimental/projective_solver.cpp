
#include "../misc/matlab_api.hpp"

#include "../core/algorithms.hpp"
#include "../core/containers.hpp"

#include "projective_solver.hpp"

namespace pano {
    namespace experimental {

        using CoeffVec = ProjectiveComponent::CoeffVec;

        // [1.0/norm(line.center())]
        class LineDoF1 : public ProjectiveComponent {
        public:
            LineDoF1(Line3 & l);
            CoeffVec coefficients(const Vec3 & direction) const override;
            void updateUsingParams(const double * params) const override;
            Line3 & line;
        };

        LineDoF1::LineDoF1(Line3 & l) : ProjectiveComponent(1), line(l) { }
        CoeffVec LineDoF1::coefficients(const Vec3 & direction) const {
            Ray3 infLine = (line / norm(line.center())).ray();
            // variable is 1.0/centerDepth
            // corresponding coeff is 1.0/depthRatio
            // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
            double depthRatio = norm(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), direction), infLine).second.first);
            return{ 1.0 / depthRatio };
        }
        void LineDoF1::updateUsingParams(const double * params) const {
            double centerDepth = 1.0 / params[0];
            line = line / norm(line.center()) * centerDepth;
        }




        // [1.0/norm(line.first), 1.0/norm(line.second)]
        class LineDoF2 : public ProjectiveComponent {
        public:
            LineDoF2(Line3 & l);
            CoeffVec coefficients(const Vec3 & direction) const override;
            void updateUsingParams(const double * params) const override;
            Line3 & line;
        };
        LineDoF2::LineDoF2(Line3 & l) : ProjectiveComponent(2), line(l){ }
        CoeffVec LineDoF2::coefficients(const Vec3 & direction) const {
            double theta = AngleBetweenDirections(normalize(line.first), normalize(line.second));
            double phi = AngleBetweenDirections(normalize(line.first), direction);
            /*    | sin(theta) | | p | | q |
            len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
               | p sin(phi) - q sin(phi - theta) |
            */
            // variables[0] -> 1/p
            // variables[1] -> 1/q
            double coeffFor1_p = -sin(phi - theta) / sin(theta);
            double coeffFor1_q = sin(phi) / sin(theta);
            assert(!IsInfOrNaN(coeffFor1_p) && !IsInfOrNaN(coeffFor1_q));
            return{ coeffFor1_p, coeffFor1_q };
        }
        void LineDoF2::updateUsingParams(const double * params) const {
            double p = 1.0 / params[0];
            double q = 1.0 / params[1];
            line.first = normalize(line.first) * p;
            line.second = normalize(line.second) * q;
        }



        // [1.0/norm(center)]
        class PlaneDoF1 : public ProjectiveComponent {
        public:
            PlaneDoF1(Plane3 & p);
            CoeffVec coefficients(const Vec3 & direction) const override;
            void updateUsingParams(const double * params) const override;
            Plane3 & plane;
        };
        PlaneDoF1::PlaneDoF1(Plane3 & p) : ProjectiveComponent(1), plane(p) {}
        CoeffVec PlaneDoF1::coefficients(const Vec3 & direction) const {
            Plane3 plane(normalize(plane.anchor), plane.normal);
            // variable is 1.0/centerDepth
            // corresponding coeff is 1.0/depthRatio
            // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
            double depthRatio = norm(IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), direction), plane).position);
            return{ 1.0 / depthRatio };
        }
        void PlaneDoF1::updateUsingParams(const double * params) const {
            double anchorDepth = 1.0 / params[0];
            plane = Plane3(normalize(plane.anchor) * anchorDepth, plane.normal);
        }



        // [a, b] in ax+by+cz=1
        class PlaneDoF2 : public ProjectiveComponent {
        public:
            PlaneDoF2(Plane3 & p, const Vec3 & ad);
            CoeffVec coefficients(const Vec3 & direction) const override;
            void updateUsingParams(const double * params) const override;
            Plane3 & plane;
            Vec3 axisDirection;
        };

        PlaneDoF2::PlaneDoF2(Plane3 & p, const Vec3 & ad) : ProjectiveComponent(2), plane(p), axisDirection(ad) { }
        int SwappedComponent(const Vec3 & orientation) {
            for (int i = 0; i < 2; i++){
                if (abs(orientation[i]) >= 1e-8){
                    return i;
                }
            }
            return 2;
        }
        CoeffVec PlaneDoF2::coefficients(const Vec3 & direction) const {
            // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
            // -> 1.0/depth = ax + by + cz
            // -> 1.0/depth = ax + by + (- o1a - o2b)/o3 * z = a(x-o1*z/o3) + b(y-o2*z/o3)
            auto orientation = normalize(axisDirection);
            int sc = SwappedComponent(orientation);
            Vec3 forientation = orientation;
            std::swap(forientation[sc], forientation[2]);
            Vec3 fdirection = direction;
            std::swap(fdirection[sc], fdirection[2]);
            return{ fdirection[0] - forientation[0] * fdirection[2] / forientation[2],
                fdirection[1] - forientation[1] * fdirection[2] / forientation[2] };
        }
        void PlaneDoF2::updateUsingParams(const double * params) const {
            double vs[] = { params[0], params[1], 0.0 }; // fake vs
            // v1 * o1 + v2 * o2 + v3 * o3 = 0, since region.orientation is orthogonal to normal
            auto orientation = axisDirection;
            int c = SwappedComponent(orientation);
            std::swap(orientation[c], orientation[2]); // now fake orientation
            vs[2] = (-vs[0] * orientation[0] - vs[1] * orientation[1])
                / orientation[2];
            std::swap(vs[c], vs[2]); // now real vs
            plane = Plane3FromEquation(vs[0], vs[1], vs[2]);
        }




        // [a, b, c] in ax+by+cz=1
        class PlaneDoF3 : public ProjectiveComponent {
        public:
            PlaneDoF3(Plane3 & p);
            CoeffVec coefficients(const Vec3 & direction) const override;
            void updateUsingParams(const double * params) const override;
            Plane3 & plane;
        };
        PlaneDoF3::PlaneDoF3(Plane3 & p) : ProjectiveComponent(3), plane(p) { }
        CoeffVec PlaneDoF3::coefficients(const Vec3 & direction) const {
            // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
            // -> 1.0/depth = ax + by + cz
            return{ direction[0], direction[1], direction[2] };
        }
        void PlaneDoF3::updateUsingParams(const double * params) const {
            plane = Plane3FromEquation(params[0], params[1], params[2]);
        }




        // [scale of placed mesh / scale of real mesh]
        struct MeshDoF1Data;
        class MeshDoF1 : public ProjectiveComponent{
        public:
            MeshDoF1(Mesh<Point3> & m);
            std::vector<Point3> intersectionsOnPlacedMesh(const Vec3 & direction) const;
            CoeffVec coefficients(const Vec3 & direction) const override;
            CoeffVec coefficients(const ProjectiveAnchor::Ptr & anchor) const override;
            void updateUsingParams(const double * params) const override;
            const MeshDoF1Data & data() const { return *_data; }
            Mesh<Point3> & mesh;
        private:
            std::unique_ptr<MeshDoF1Data> _data;
        };
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

        MeshDoF1::MeshDoF1(Mesh<Point3> & m) : ProjectiveComponent(1), mesh(m), _data(std::make_unique<MeshDoF1Data>(m)) {}

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

        // scalar relative to current mesh
        // 1.0/depthofcurp * [1.0/scalar] = 1.0/depthofrealp
        CoeffVec MeshDoF1::coefficients(const Vec3 & direction) const {
            auto pts = intersectionsOnPlacedMesh(direction);
            assert(!pts.empty());
            return{ 1.0 / norm(pts.front()) };
        }

        class NormalAnchor : public ProjectiveAnchor {
        public:
            NormalAnchor(const Vec3 & dd) : d(dd) {}
            Vec3 direction() const override { return d; }
            Vec3 d;
        };

        class MeshVertexAnchor : public ProjectiveAnchor {
        public:
            MeshVertexAnchor(const ProjectiveComponent * mcomp, Mesh<Point3>::VertHandle v)
                : meshComp(dynamic_cast<const MeshDoF1*>(mcomp)), vh(v) {
                assert(meshComp);
            }
            Vec3 direction() const override { return normalize(meshComp->mesh.data(vh)); }
            bool isGeneratedFrom(const ProjectiveComponent * c) const override {
                return meshComp == c;
            }
            const MeshDoF1 * meshComp;
            Mesh<Point3>::VertHandle vh;
        };

        CoeffVec MeshDoF1::coefficients(const ProjectiveAnchor::Ptr & anchor) const {
            if (anchor->isGeneratedFrom(this)){
                const MeshVertexAnchor * ma = static_cast<const MeshVertexAnchor*>(anchor.get());
                return{ 1.0 / norm(_data->mesh.data(ma->vh)) };
            }
            else{
                return coefficients(anchor->direction());
            }
        }

        void MeshDoF1::updateUsingParams(const double * params) const {
            for (auto & v : mesh.vertices()){
                v.data /= params[0];
            }
        }





        int ProjectiveSolver::bindLineDoF1(Line3 & l) { return append(std::make_unique<LineDoF1>(l)); }
        int ProjectiveSolver::bindLineDoF2(Line3 & l) { return append(std::make_unique<LineDoF2>(l)); }
        int ProjectiveSolver::bindPlaneDoF1(Plane3 & p) { return append(std::make_unique<PlaneDoF1>(p)); }
        int ProjectiveSolver::bindPlaneDoF2(Plane3 & p, const Vec3 & axisDirection) { return append(std::make_unique<PlaneDoF2>(p, axisDirection)); }
        int ProjectiveSolver::bindPlaneDoF3(Plane3 & p) { return append(std::make_unique<PlaneDoF3>(p)); }
        int ProjectiveSolver::bindDoF1(Mesh<Point3> & mesh) { return append(std::make_unique<MeshDoF1>(mesh)); }

        int ProjectiveSolver::makeNormalAnchor(const Vec3 & d) { return append(std::make_unique<NormalAnchor>(d)); }
        int ProjectiveSolver::makeMeshVertexAnchor(int meshComponentId, Mesh<Point3>::VertHandle vh) {
            return append(std::make_unique<MeshVertexAnchor>(_components.at(meshComponentId).get(), vh));
        }


        int ProjectiveSolver::makeAEqualToBAt(int a, int b, int anchorId){
            auto & anchor = _anchors.at(anchorId);
            auto coeffA = _components.at(a)->coefficients(anchor);
            auto coeffB = _components.at(b)->coefficients(anchor);
            int eid = _neqs;
            for (int i = 0; i < _components.at(a)->nparams(); i++){
                _Amat.emplace_back(eid, i + _compStartPositionsInX.at(a), coeffA.at(i));
            }
            for (int i = 0; i < _components.at(b)->nparams(); i++){
                _Amat.emplace_back(eid, i + _compStartPositionsInX.at(b), - coeffB.at(i));
            }
            _Bmat.push_back(0.0);
            _ops.push_back(Equal);
            _neqs++;
            return eid;
        }

        // invd(A) > invd(B) -> .. B - .. A < 0
        int ProjectiveSolver::makeACloserThanBAt(int a, int b, int anchorId){
            auto & anchor = _anchors.at(anchorId);
            auto coeffA = _components.at(a)->coefficients(anchor);
            auto coeffB = _components.at(b)->coefficients(anchor);
            int eid = _neqs;
            for (int i = 0; i < _components.at(a)->nparams(); i++){
                _Amat.emplace_back(eid, i + _compStartPositionsInX.at(a), - coeffA.at(i));
            }
            for (int i = 0; i < _components.at(b)->nparams(); i++){
                _Amat.emplace_back(eid, i + _compStartPositionsInX.at(b), coeffB.at(i));
            }
            _Bmat.push_back(0.0);
            _ops.push_back(LowerThan);
            _neqs++;
            return eid;
        }

        int ProjectiveSolver::makeAEqualToDepthAt(int a, double d, int anchorId){
            auto & anchor = _anchors.at(anchorId);
            auto coeffA = _components.at(a)->coefficients(anchor);
            int eid = _neqs;
            for (int i = 0; i < _components.at(a)->nparams(); i++){
                _Amat.emplace_back(eid, i + _compStartPositionsInX.at(a), coeffA.at(i));
            }
            _Bmat.push_back(1.0 / d);
            _ops.push_back(Equal);
            _neqs++;
            return eid;
        }

        int ProjectiveSolver::makeACloserThanDepthAt(int a, double d, int anchorId){
            auto & anchor = _anchors.at(anchorId);
            auto coeffA = _components.at(a)->coefficients(anchor);
            int eid = _neqs;
            for (int i = 0; i < _components.at(a)->nparams(); i++){
                _Amat.emplace_back(eid, i + _compStartPositionsInX.at(a), - coeffA.at(i));
            }
            _Bmat.push_back(- 1.0 / d);
            _ops.push_back(LowerThan);
            _neqs++;
            return eid;
        }

        int ProjectiveSolver::makeAFartherThanDepthAt(int a, double d, int anchorId){
            auto & anchor = _anchors.at(anchorId);
            auto coeffA = _components.at(a)->coefficients(anchor);
            int eid = _neqs;
            for (int i = 0; i < _components.at(a)->nparams(); i++){
                _Amat.emplace_back(eid, i + _compStartPositionsInX.at(a), coeffA.at(i));
            }
            _Bmat.push_back(1.0 / d);
            _ops.push_back(LowerThan);
            _neqs++;
            return eid;
        }


        int ProjectiveSolver::append(ProjectiveComponent::Ptr && p) {
            _components.push_back(std::move(p));
            _compStartPositionsInX.push_back(_nvars);
            _nvars += _components.back()->nparams();
            return _components.size() - 1;
        }
        int ProjectiveSolver::append(ProjectiveAnchor::Ptr && a) {
            _anchors.push_back(std::move(a));
            return _anchors.size() - 1;
        }



        bool ProjectiveSolver::solve(double * nanOrInfRatio, bool quiet) const {
            SparseMatd A = MakeSparseMatFromElements(_neqs, _nvars, _Amat.begin(), _Amat.end());

            assert(_Bmat.size() == _neqs && _ops.size() == _neqs);

            bool succeed = true;
            misc::Matlab matlab;
            matlab << "clear;";
            succeed = matlab.setVar("A", A) &&
                matlab.setVar("B", cv::Mat(_Bmat));
            assert(succeed);
            matlab << "B = B';";

            std::vector<double> opTypes(_ops.size());
            for (int i = 0; i < _ops.size(); i++){
                opTypes[i] = _ops[i] == Equal ? 0.0 : 1.0;
            }
            succeed = matlab.setVar("Op", cv::Mat(opTypes));
            assert(succeed);
            matlab << "Op = Op';";

            bool result = matlab.run("neq = size(A, 1);") &&
                matlab.run("nvar = size(A, 2);") &&
                matlab.run("equalpart = Op(:) == 0;") &&
                matlab.run("ltpart = Op(:) == 1;") &&

                matlab.run(quiet ? "cvx_begin quiet" : "cvx_begin") &&
                matlab.run("variable X(nvar);") &&
                matlab.run("minimize norm(A(equalpart, :) * X - B(equalpart, :))") &&
                matlab.run("subject to ") &&
                matlab.run("    A(ltpart, :) * X <= B(ltpart, :);") &&
                matlab.run("cvx_end") &&

                matlab.run("X = X(:)' / norm(X);");
            if (!result)
                return false;

            std::vector<double> X = matlab.var("X").toCVMat();

            if (nanOrInfRatio) {
                *nanOrInfRatio = 
                    std::count_if(X.begin(), X.end(), [](double x) {return IsInfOrNaN(x); }) / double(X.size());
            }
            if (std::any_of(X.begin(), X.end(), [](double x) {return IsInfOrNaN(x); }) ||
                std::all_of(X.begin(), X.end(), [](double x) {return x == 0.0; }))
                return false;

            for (int i = 0; i < _components.size(); i++){
                _components[i]->updateUsingParams(X.data() + _compStartPositionsInX[i]);
            }

            return true;
        }



    }
}