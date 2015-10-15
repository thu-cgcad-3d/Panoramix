#include "../core/containers.hpp"
#include "pi_graph_cg.hpp"

namespace pano {
    namespace experimental {

        namespace {
            inline double AngleDistanceBetweenPointAndRay(const Vec3 & p, const Ray3 & ray) {
                auto p1 = ray.anchor;
                auto p2 = ray.direction;
                auto normal = p1.cross(p2);
                return M_PI_2 - AngleBetweenUndirectedVectors(p, normal);
            }
            template <class T, int N>
            bool IsFuzzyNormalized(const Vec<T, N> & v, T epsilon = 1e-3) {
                return abs(norm(v) - 1.0) < epsilon;
            }
            int SwappedComponent(const Vec3 & orientation) {
                //for (int i = 0; i < 2; i++) {
                //    if (abs(orientation[i]) >= 1e-8) {
                //        return i;
                //    }
                //}
                //return 2;
                int maxComp = std::max_element(orientation.val, orientation.val + 3, 
                    [](double a, double b) {
                    return abs(a) < abs(b); 
                }) - orientation.val;
                assert(abs(orientation[maxComp]) > 1e-8);
                return maxComp;
            }
        }


        DenseMatd PIConstraintGraph::Entity::SupportingPlane::matFromVarsToPlaneCoeffs() const {
            assert(dof == 1 || dof == 2 || dof == 3);
            if (dof == 1) {
                assert(IsFuzzyNormalized(toward));
                auto planeNormal = normalize(toward);
                // variable is 1.0/centerDepth
                //  let the normalized planeNormal = [nx, ny, nz], let normalized center = [cx, cy, cz]
                //  then the equation should be k nx x + k ny y + k nz z = 1
                //  if the depth of center is Dc, then
                //      k cx Dc nx + k cy Dc ny + k cz Dc nz = 1
                //  and then k = 1.0 / (Dc*dot(n, c))
                //  the equation coefficients thus are
                //      k n* = n* / (Dc*dot(n, c))
                //  therefore the equation coefficients corresponding to the variable (inverse planeDistanceToOrigin, 1.0/Dc) are
                //      n*/dot(n, c)
                DenseMatd coeffs(3, 1, 0.0);
                for (int i = 0; i < 3; i++) {
                    coeffs(i, 0) = planeNormal[i] / planeNormal.dot(center);
                }
                return coeffs;
            } else if (dof == 2) {
                assert(IsFuzzyNormalized(along));
                auto m = normalize(along);
                int sc = SwappedComponent(m); // find a non zero component
                assert(sc == 0 || sc == 1 || sc == 2);

                // original along orientation = [mx, my, mz]
                // let the plane equation be ax + by + cz = 1
                // there should be a * mx + b * my + c * mz = 0

                // if sc = 0, meaning that mx != 0
                // then we can write as: a = b * (- my / mx) + c * (- mz / mx)
                // then [a]   [-my/mx -mz/mx] 
                //      [b] = [   1      0  ]*[b]
                //      [c]   [   0      1  ] [c]
                if (sc == 0) {
                    DenseMatd coeffs(3, 2, 0.0);
                    coeffs(0, 0) = -m[1] / m[0];    coeffs(0, 1) = -m[2] / m[0];
                    coeffs(1, 0) = 1;               coeffs(1, 1) = 0;
                    coeffs(2, 0) = 0;               coeffs(2, 1) = 1;
                    return coeffs;
                }

                // if sc = 1, meaning that my != 0
                // then we can write as: b = a * (- mx / my) + c * (- mz / my)
                // then [a]   [    1     0  ] [a]
                //      [b] = [-mx/my -mz/my]*
                //      [c]   [    0     1  ] [c]
                if (sc == 1) {
                    DenseMatd coeffs(3, 2, 0.0);
                    coeffs(0, 0) = 1;               coeffs(0, 1) = 0;
                    coeffs(1, 0) = -m[0] / m[1];    coeffs(1, 1) = -m[2] / m[1];
                    coeffs(2, 0) = 0;               coeffs(2, 1) = 1;
                    return coeffs;
                }

                // if sc = 2, meaning that mz != 0
                // then we can write as: c = a * (- mx / mz) + b * (- my / mz)
                // then [a]   [    1     0  ] [a]
                //      [b] = [    0     1  ]*[b]
                //      [c]   [-mx/mz -my/mz]
                /*if (sc == 2)*/ {
                    DenseMatd coeffs(3, 2, 0.0);
                    coeffs(0, 0) = 1;               coeffs(0, 1) = 0;
                    coeffs(1, 0) = 0;               coeffs(1, 1) = 1;
                    coeffs(2, 0) = -m[0] / m[2];    coeffs(2, 1) = -m[1] / m[2];
                    return coeffs;
                }
            } else /*if (dof == 3)*/ {
                return DenseMatd::eye(3, 3);
            }
        }


        // the staibility info
        struct Stability {
            virtual int dof() const = 0;
            virtual void addAnchor(const Vec3 & anchor, double angleThres) = 0;
            virtual Stability * clone() const = 0;
            virtual void stablize() = 0;
        };
        bool supportersAreValid(const std::vector<Vec3> & supporters, double angleThres) {
            for (int i = 0; i < supporters.size(); i++) {
                for (int j = i + 1; j < supporters.size(); j++) {
                    if (AngleBetweenDirections(supporters[i], supporters[j]) <= angleThres) {
                        return false;
                    }
                }
            }
            return true;
        }
        struct SegDof3Stability : Stability {
            std::vector<Vec3> supporters;
            SegDof3Stability() {}
            SegDof3Stability(const SegDof3Stability & s) = default;
            virtual int dof() const override {
                return 3 - supporters.size();
            }
            virtual void addAnchor(const Vec3 & anchor, double angleThres) override {
                assert(supporters.size() <= 3);
                if (supporters.size() == 3) {
                    return;
                }
                if (supporters.empty()) {
                    supporters.push_back(anchor);
                } else if (supporters.size() == 1) {
                    if (AngleBetweenDirections(supporters[0], anchor) > angleThres) {
                        supporters.push_back(anchor);
                    }
                } else if (supporters.size() == 2) {
                    if (AngleDistanceBetweenPointAndRay(anchor, Line3(supporters[0], supporters[1]).ray()) > angleThres) {
                        supporters.push_back(anchor);
                    }
                }
                assert(supportersAreValid(supporters, angleThres));
            }
            virtual Stability * clone() const override {
                return new SegDof3Stability(*this);
            }
            virtual void stablize() override {
                supporters.resize(3);
            }
        };
        struct SegDof2Stability : Stability {
            Vec3 axis;
            std::vector<Vec3> supporters;
            explicit SegDof2Stability(const Vec3 & a) : axis(a) {}
            SegDof2Stability(const SegDof2Stability & s) = default;
            virtual int dof() const override {
                return 2 - supporters.size();
            }
            virtual void addAnchor(const Vec3 & anchor, double angleThres) override {
                assert(supporters.size() <= 2);
                if (supporters.size() == 2) {
                    return;
                }
                if (supporters.empty()) {
                    supporters.push_back(anchor);
                } else if (supporters.size() == 1) {
                    if (AngleDistanceBetweenPointAndRay(anchor, Ray3(supporters[0], axis)) > angleThres) {
                        supporters.push_back(anchor);
                    }
                }
                assert(supportersAreValid(supporters, angleThres));
            }
            virtual Stability * clone() const override {
                return new SegDof2Stability(*this);
            }
            virtual void stablize() override {
                supporters.resize(2);
            }
        };
        struct SegDof1Stability : Stability {
            bool stable;
            SegDof1Stability() : stable(false) {}
            SegDof1Stability(const SegDof1Stability & s) = default;
            virtual int dof() const override { return stable ? 0 : 1; }
            virtual void addAnchor(const Vec3 & anchor, double angleThres) override {
                stable = true;
            }
            virtual Stability * clone() const override {
                return new SegDof1Stability(*this);
            }
            virtual void stablize() override {
                stable = true;
            }
        };

        struct LineDof2Stability : Stability {
            std::vector<Vec3> supporters;
            LineDof2Stability() {}
            LineDof2Stability(const LineDof2Stability & s) = default;
            virtual int dof() const override {
                return 2 - supporters.size();
            }
            virtual void addAnchor(const Vec3 & anchor, double angleThres) override {
                assert(supporters.size() <= 2);
                if (supporters.empty()) {
                    supporters.push_back(anchor);
                } else if (supporters.size() == 1) {
                    if (AngleBetweenDirections(supporters[0], anchor) > angleThres) {
                        supporters.push_back(anchor);
                    }
                }
            }
            virtual Stability * clone() const override {
                return new LineDof2Stability(*this);
            }
            virtual void stablize() override {
                supporters.resize(2);
            }
        };
        struct LineDof1Stability : Stability {
            bool stable;
            LineDof1Stability() : stable(false) {}
            LineDof1Stability(const LineDof1Stability & s) = default;
            virtual int dof() const override { return stable ? 0 : 1; }
            virtual void addAnchor(const Vec3 & anchor, double angleThres) override {
                stable = true;
            }
            virtual Stability * clone() const override {
                return new LineDof1Stability(*this);
            }
            virtual void stablize() override {
                stable = true;
            }
        };





        PIConstraintGraph BuildPIConstraintGraph(const PIGraph & mg, double minAngleThresForAWideEdge, double angleThresForColinearity) {
            std::vector<std::unique_ptr<Stability>> ent2stab;

            PIConstraintGraph cg;
            cg.seg2ent.resize(mg.nsegs, -1);
            cg.line2ent.resize(mg.nlines(), -1);

            auto & seg2ent = cg.seg2ent;
            auto & line2ent = cg.line2ent;
            auto & entities = cg.entities;
            auto & constraints = cg.constraints;

            using Entity = PIConstraintGraph::Entity;
            using Constraint = PIConstraintGraph::Constraint;

            // add entities
            // add segs
            for (int i = 0; i < mg.nsegs; i++) {
                Entity e;
                e.type = Entity::IsSeg;
                e.id = i;
                e.size = mg.seg2areaRatio[i] * 100;

                auto & control = mg.seg2control[i];
                e.supportingPlane.dof = mg.seg2control[i].dof();
                e.supportingPlane.center = normalize(mg.seg2center[i]);
                std::unique_ptr<Stability> stability;
                if (e.supportingPlane.dof == 1) {
                    e.supportingPlane.toward = normalize(mg.vps[control.orientationClaz]);
                    if (e.supportingPlane.toward.dot(e.supportingPlane.center) < 0) {
                        e.supportingPlane.toward = -e.supportingPlane.toward;
                    }
                    stability = std::make_unique<SegDof1Stability>();
                } else if (e.supportingPlane.dof == 2) {
                    e.supportingPlane.along = normalize(mg.vps[control.orientationNotClaz]);
                    stability = std::make_unique<SegDof2Stability>(mg.vps[mg.seg2control[i].orientationNotClaz]);
                } else {
                    stability = std::make_unique<SegDof3Stability>();
                }

                entities.push_back(e);
                ent2stab.push_back(std::move(stability));

                int ent = entities.size() - 1;
                seg2ent[i] = ent;              
            }
            // add lines
            for (int i = 0; i < mg.nlines(); i++) {
                Entity e;
                e.type = Entity::IsLine;
                e.id = i;
                e.size = AngleBetweenDirections(mg.lines[i].component.first, mg.lines[i].component.second);

                e.supportingPlane.dof = mg.lines[i].claz != -1 ? 1 : 2;
                e.supportingPlane.center = normalize(mg.lines[i].component.center());
                std::unique_ptr<Stability> stability;
                if (e.supportingPlane.dof == 1) {
                    const Vec3 & lineDirection = mg.vps[mg.lines[i].claz];
                    auto linePerp = normalize(mg.lines[i].component.center().cross(lineDirection));
                    assert(abs(norm(linePerp) - 1.0) < 1e-2);
                    auto linePlaneNormal = lineDirection.cross(linePerp);                    
                    assert(IsFuzzyPerpendicular(linePlaneNormal, mg.vps[mg.lines[i].claz]));
                    e.supportingPlane.toward = normalize(linePlaneNormal);
                    stability = std::make_unique<LineDof1Stability>();
                } else {
                    const Vec3 & maybeLineDirection = mg.lines[i].component.direction();
                    auto linePerp = mg.lines[i].component.center().cross(maybeLineDirection);
                    assert(IsFuzzyPerpendicular(linePerp, maybeLineDirection));
                    e.supportingPlane.along = normalize(linePerp);
                    stability = std::make_unique<LineDof2Stability>();
                }

                entities.push_back(e);
                ent2stab.push_back(std::move(stability));

                int ent = entities.size() - 1;
                line2ent[i] = ent;
            }


            auto & ent2cons = cg.ent2cons;
            ent2cons.resize(entities.size()); // enttity -> constraints


            // add constraints
            // bndpieces, seg-seg
            for (int i = 0; i < mg.nbndPieces(); i++) {
                if (mg.bndPiece2segRelation[i] != SegRelation::Connected) {
                    continue;
                }
                Constraint connect;
                connect.type = Constraint::Connection;
                auto & segPair = mg.bnd2segs[mg.bndPiece2bnd[i]];
                connect.ent1 = seg2ent[segPair.first];
                connect.ent2 = seg2ent[segPair.second];
                double len = mg.bndPiece2length[i];
                connect.weight = len;
                if (len >= minAngleThresForAWideEdge) {
                    connect.anchors = { mg.bndPiece2dirs[i].front(), mg.bndPiece2dirs[i].back() };
                } else {
                    connect.anchors = { normalize(mg.bndPiece2dirs[i].front() + mg.bndPiece2dirs[i].back()) };
                }
                constraints.push_back(connect);
                int connectCon = constraints.size() - 1;
                ent2cons[connect.ent1].push_back(connectCon);
                ent2cons[connect.ent2].push_back(connectCon);

                // check whether this bp lies on a line
                if (!mg.bndPiece2linePieces[i].empty() || len < minAngleThresForAWideEdge) {
                    // if so, no planarity will be assigned
                    continue;
                }

                // the coplanarity constraint if there are no lines on the bndpiece
                if (mg.bndPiece2linePieces[i].empty()) {
                    Constraint coplanar;
                    coplanar.type = Constraint::Coplanarity;
                    coplanar.ent1 = connect.ent1;
                    coplanar.ent2 = connect.ent2;
                    coplanar.weight = len;
                    //coplanar.anchorOrientation = -1;
                    constraints.push_back(coplanar);
                    int coplanarCon = constraints.size() - 1;
                    ent2cons[coplanar.ent1].push_back(coplanarCon);
                    ent2cons[coplanar.ent2].push_back(coplanarCon);
                }
            }

            // linepieces, seg-line
            for (int i = 0; i < mg.nlinePieces(); i++) {
                if (mg.linePiece2segLineRelation[i] == SegLineRelation::Detached) {
                    continue;
                }
                int bndPiece = mg.linePiece2bndPiece[i];
                int line = mg.linePiece2line[i];
                if (bndPiece == -1) {
                    int seg = mg.linePiece2seg[i];
                    assert(seg != -1);
                    Constraint c;
                    c.type = Constraint::Connection;
                    c.ent1 = seg2ent[seg];
                    c.ent2 = line2ent[line];
                    double len = mg.linePiece2length[i];
                    c.weight = len;
                    if (len >= minAngleThresForAWideEdge) {
                        c.anchors = { mg.linePiece2samples[i].front(), mg.linePiece2samples[i].back() };
                    } else {
                        c.anchors = { normalize(mg.linePiece2samples[i].front() + mg.linePiece2samples[i].back()) };
                    }
                    constraints.push_back(c);
                    int con = constraints.size() - 1;
                    ent2cons[c.ent1].push_back(con);
                    ent2cons[c.ent2].push_back(con);
                } else {
                    int segPair[] = { -1, -1 };
                    std::tie(segPair[0], segPair[1]) = mg.bnd2segs[mg.bndPiece2bnd[bndPiece]];
                    auto segRelation = mg.bndPiece2segRelation[bndPiece];
                    bool connected[] = {
                        segRelation == SegRelation::Connected || segRelation == SegRelation::LeftIsFront,
                        segRelation == SegRelation::Connected || segRelation == SegRelation::RightIsFront
                    };
                    for (int k = 0; k < 2; k++) {
                        int seg = segPair[k];
                        if (!connected[k]) {
                            continue;
                        }
                        Constraint c;
                        c.type = Constraint::Connection;
                        c.ent1 = seg2ent[seg];
                        c.ent2 = line2ent[line];
                        double len = mg.linePiece2length[i];
                        c.weight = len;
                        if (len >= minAngleThresForAWideEdge) {
                            c.anchors = { mg.linePiece2samples[i].front(), mg.linePiece2samples[i].back() };
                        } else {
                            c.anchors = { normalize(mg.linePiece2samples[i].front() + mg.linePiece2samples[i].back()) };
                        }
                        constraints.push_back(c);
                        int con = constraints.size() - 1;
                        ent2cons[c.ent1].push_back(con);
                        ent2cons[c.ent2].push_back(con);
                    }
                }
            }

            // linerelations, line-line
            for (int i = 0; i < mg.nlineRelations(); i++) {
                if (mg.lineRelations[i] == LineRelation::Detached) {
                    continue;
                }
                Constraint c;
                c.type = Constraint::Connection;
                auto & linePair = mg.lineRelation2lines[i];
                c.ent1 = line2ent[linePair.first];
                c.ent2 = line2ent[linePair.second];
                c.weight = mg.lineRelation2weight[i];
                c.anchors = { mg.lineRelation2anchor[i] };
                //c.anchorOrientation = -1;
                constraints.push_back(c);
                int con = constraints.size() - 1;
                ent2cons[c.ent1].push_back(con);
                ent2cons[c.ent2].push_back(con);
            }






            // select the largest dof1 ent as the root!
            std::vector<int> rootCands;
            for (int i = 0; i < entities.size(); i++) {
                int sz = entities[i].size;
                if (ent2stab[i]->dof() != 1) {
                    continue;
                }
                rootCands.push_back(i);
            }
            std::sort(rootCands.begin(), rootCands.end(), [&entities](int a, int b) {return entities[a].size > entities[b].size; });
            if (rootCands.empty()) {
                std::cout << "we can't find any entity whose dof is 1 !!!!!" << std::endl;
                return cg;
            }

            for (int root : rootCands) {
                std::cout << "root: " << root << std::endl;
                std::set<int> entsCollected;

                // initialize stabilities of entities
                std::vector<std::unique_ptr<Stability>> ent2stabHere(entities.size());
                for (int ent = 0; ent < entities.size(); ent++) {
                    ent2stabHere[ent] = std::unique_ptr<Stability>(ent2stab[ent]->clone());
                }
                ent2stabHere[root]->stablize();
                assert(ent2stabHere[root]->dof() == 0);

                MaxHeap<int, int, std::greater<int>> Q; // a min heap recording dofs
                Q.set(root, ent2stabHere[root]->dof());
                while (!Q.empty()) {
                    int curEnt = Q.top();
                    int curEntDoF = ent2stabHere[curEnt]->dof();
                    if (curEntDoF != 0) { // all remaining adjacent entities are not stable, stop the search
                        break;
                    }
                    assert(std::all_of(Q.begin(), Q.end(), [curEntDoF](const Scored<int, int> & ent) {return ent.score >= curEntDoF; }));
                    Q.pop();
                    entsCollected.insert(curEnt);
                    for (int con : ent2cons[curEnt]) {
                        auto & c = constraints[con];
                        if (c.weight == 0.0) {
                            continue;
                        }
                        int adjEnt = c.ent1 == curEnt ? c.ent2 : c.ent1;
                        if (Contains(entsCollected, adjEnt)) {
                            continue;
                        }
                        if (c.isConnection()) {
                            for (auto & anchor : c.anchors) {
                                ent2stabHere[adjEnt]->addAnchor(anchor, angleThresForColinearity);
                            }
                        } else if (c.isCoplanarity()) {
                            assert(entities[adjEnt].isSeg());
                            ent2stabHere[adjEnt]->stablize();
                        }
                        Q.set(adjEnt, ent2stabHere[adjEnt]->dof());
                    }
                }
                if (entsCollected.size() > entities.size() / 2) {
                    cg.rootEnt = root;
                    cg.determinableEnts = std::move(entsCollected);
                    break;
                }

            }

            cg.consBetweenDeterminableEnts.clear();
            for (int i = 0; i < cg.constraints.size(); i++) {
                auto & c = cg.constraints[i];
                if (Contains(cg.determinableEnts, c.ent1) && Contains(cg.determinableEnts, c.ent2)) {
                    cg.consBetweenDeterminableEnts.insert(i);
                }
            }

            return cg;
        }

        PIConstraintGraph BuildPIConstraintGraph(const PILayoutAnnotation & anno, double minAngleThresForAWideEdge) {

            PIConstraintGraph cg;
            cg.seg2ent.resize(anno.nfaces(), -1);

            auto & seg2ent = cg.seg2ent;
            auto & line2ent = cg.line2ent;
            auto & entities = cg.entities;
            auto & constraints = cg.constraints;

            using Entity = PIConstraintGraph::Entity;
            using Constraint = PIConstraintGraph::Constraint;

            // add entities
            // add segs
            for (int i = 0; i < anno.nfaces(); i++) {
                Entity e;
                e.type = Entity::IsSeg;
                e.id = i;
                e.size = 1;

                auto & control = anno.face2control[i];
                e.supportingPlane.dof = control.dof();
                Vec3 center = Origin();
                for (int c : anno.face2corners[i]) {
                    center += anno.corners[c];
                }
                e.supportingPlane.center = normalize(center);
                if (e.supportingPlane.dof == 1) {
                    e.supportingPlane.toward = normalize(anno.vps[control.orientationClaz]);
                    if (e.supportingPlane.toward.dot(e.supportingPlane.center) < 0) {
                        e.supportingPlane.toward = -e.supportingPlane.toward;
                    }
                } else if (e.supportingPlane.dof == 2) {
                    e.supportingPlane.along = normalize(anno.vps[control.orientationNotClaz]);
                } else {
                }

                entities.push_back(e);

                int ent = entities.size() - 1;
                seg2ent[i] = ent;
            }


            auto & ent2cons = cg.ent2cons;
            ent2cons.resize(entities.size()); // enttity -> constraints


            // add constraints
            // face connect face
            // get corners2border and border2face
            std::map<std::pair<int, int>, int> corners2border;
            for (int i = 0; i < anno.nborders(); i++) {
                auto & cs = anno.border2corners[i];
                corners2border[cs] = i;
            }
            std::vector<std::pair<int, int>> border2face(anno.nborders(), std::make_pair(-1, -1));
            for (int i = 0; i < anno.nfaces(); i++) {
                auto & cs = anno.face2corners[i];
                for (int j = 0; j < cs.size(); j++) {
                    int c1 = cs[j];
                    int c2 = cs[(j + 1) % cs.size()];
                    if (Contains(corners2border, std::make_pair(c1, c2))) {
                        int b = corners2border.at(std::make_pair(c1, c2));
                        border2face[b].first = i;
                    } else if (Contains(corners2border, std::make_pair(c2, c1))) {
                        int b = corners2border.at(std::make_pair(c2, c1));
                        border2face[b].second = i;
                    } else {
                        SHOULD_NEVER_BE_CALLED();
                    }
                }
            }

            // connect
            for (int i = 0; i < anno.nborders(); i++) {
                if (!anno.border2connected[i]) {
                    continue;
                }
                Constraint connect;
                connect.type = Constraint::Connection;
                
                auto & segPair = border2face[i];
                connect.ent1 = seg2ent[segPair.first];
                connect.ent2 = seg2ent[segPair.second];
                double len = AngleBetweenDirections(anno.corners[anno.border2corners[i].first],
                    anno.corners[anno.border2corners[i].second]);
                connect.weight = len;
                if (len >= minAngleThresForAWideEdge) {
                    connect.anchors = { anno.corners[anno.border2corners[i].first], anno.corners[anno.border2corners[i].second] };
                } else {
                    connect.anchors = { normalize(anno.corners[anno.border2corners[i].first] + anno.corners[anno.border2corners[i].second]) };
                }
                constraints.push_back(connect);
                int connectCon = constraints.size() - 1;
                ent2cons[connect.ent1].push_back(connectCon);
                ent2cons[connect.ent2].push_back(connectCon);
            }

            // coplanar
            for (auto & fp : anno.coplanarFacePairs) {
                int f1 = fp.first;
                int f2 = fp.second;
                Constraint coplanar;
                coplanar.type = Constraint::Coplanarity;
                coplanar.ent1 = f1;
                coplanar.ent2 = f2;
                coplanar.weight = 1.0;
                constraints.push_back(coplanar);
                int coplanarCon = constraints.size() - 1;
                ent2cons[coplanar.ent1].push_back(coplanarCon);
                ent2cons[coplanar.ent2].push_back(coplanarCon);
            }


            cg.rootEnt = -1;
            for (int i = 0; i < cg.entities.size(); i++) {
                if (cg.entities[i].supportingPlane.dof == 1) {
                    cg.rootEnt = i;
                }
                cg.determinableEnts.insert(i);
            }
            for (int i = 0; i < cg.constraints.size(); i++) {
                cg.consBetweenDeterminableEnts.insert(i);
            }

            return cg;

        }

        PIConstraintGraph BuildPIConstraintGraphWithLines(const PILayoutAnnotation & anno, double minAngleThresForAWideEdge) {
            PIConstraintGraph cg;
            cg.seg2ent.resize(anno.nfaces(), -1);
            cg.line2ent.resize(anno.nborders(), -1);

            auto & seg2ent = cg.seg2ent;
            auto & line2ent = cg.line2ent;
            auto & entities = cg.entities;
            auto & constraints = cg.constraints;

            using Entity = PIConstraintGraph::Entity;
            using Constraint = PIConstraintGraph::Constraint;

            // add entities
            // add segs
            for (int i = 0; i < anno.nfaces(); i++) {
                Entity e;
                e.type = Entity::IsSeg;
                e.id = i;
                e.size = 1;

                auto & control = anno.face2control[i];
                e.supportingPlane.dof = control.dof();
                Vec3 center = Origin();
                for (int c : anno.face2corners[i]) {
                    center += anno.corners[c];
                }
                e.supportingPlane.center = normalize(center);
                if (e.supportingPlane.dof == 1) {
                    e.supportingPlane.toward = normalize(anno.vps[control.orientationClaz]);
                    if (e.supportingPlane.toward.dot(e.supportingPlane.center) < 0) {
                        e.supportingPlane.toward = -e.supportingPlane.toward;
                    }
                } else if (e.supportingPlane.dof == 2) {
                    e.supportingPlane.along = normalize(anno.vps[control.orientationNotClaz]);
                } else {
                }

                entities.push_back(e);

                int ent = entities.size() - 1;
                seg2ent[i] = ent;
            }

            // add lines
            for (int i = 0; i < anno.nborders(); i++) {
                if (!anno.border2connected[i]) {
                    continue;
                }

                Entity e;
                e.type = Entity::IsLine;
                e.id = i;
                e.size = 1;

                Line3 line(normalize(anno.corners[anno.border2corners[i].first]),
                    normalize(anno.corners[anno.border2corners[i].second]));
                e.supportingPlane.dof = 2;
                e.supportingPlane.center = normalize(line.center());
                if (e.supportingPlane.dof == 1) {
                } else {
                    const Vec3 & maybeLineDirection = line.direction();
                    auto linePerp = line.center().cross(maybeLineDirection);
                    assert(IsFuzzyPerpendicular(linePerp, maybeLineDirection));
                    e.supportingPlane.along = normalize(linePerp);
                }

                entities.push_back(e);

                int ent = entities.size() - 1;
                line2ent[i] = ent;
            }


            auto & ent2cons = cg.ent2cons;
            ent2cons.resize(entities.size()); // enttity -> constraints


            // add constraints
            // face connect face
            // get corners2border and border2face
            std::map<std::pair<int, int>, int> corners2border;
            for (int i = 0; i < anno.nborders(); i++) {
                auto & cs = anno.border2corners[i];
                corners2border[cs] = i;
            }
            std::vector<std::pair<int, int>> border2face(anno.nborders(), std::make_pair(-1, -1));
            for (int i = 0; i < anno.nfaces(); i++) {
                auto & cs = anno.face2corners[i];
                for (int j = 0; j < cs.size(); j++) {
                    int c1 = cs[j];
                    int c2 = cs[(j + 1) % cs.size()];
                    if (Contains(corners2border, std::make_pair(c1, c2))) {
                        int b = corners2border.at(std::make_pair(c1, c2));
                        border2face[b].first = i;
                    } else if (Contains(corners2border, std::make_pair(c2, c1))) {
                        int b = corners2border.at(std::make_pair(c2, c1));
                        border2face[b].second = i;
                    } else {
                        SHOULD_NEVER_BE_CALLED();
                    }
                }
            }

            // border-face connect
            for (int i = 0; i < anno.nborders(); i++) {
                if (!anno.border2connected[i]) {
                    continue;
                }
                {
                    int segPair[] = { -1, -1 };
                    std::tie(segPair[0], segPair[1]) = border2face[i];
                    for (int k = 0; k < 2; k++) {
                        int seg = segPair[k];
                        Constraint c;
                        c.type = Constraint::Connection;
                        c.ent1 = seg2ent[seg];
                        c.ent2 = line2ent[i];
                        double len = 1.0;
                        c.weight = len;
                        if (len >= minAngleThresForAWideEdge) {
                            c.anchors = { anno.corners[anno.border2corners[i].first], anno.corners[anno.border2corners[i].second] };
                        } else {
                            c.anchors = { anno.corners[anno.border2corners[i].first] + anno.corners[anno.border2corners[i].second] };
                        }
                        constraints.push_back(c);
                        int con = constraints.size() - 1;
                        ent2cons[c.ent1].push_back(con);
                        ent2cons[c.ent2].push_back(con);
                    }
                }
            }


            // coplanar
            for (auto & fp : anno.coplanarFacePairs) {
                int f1 = fp.first;
                int f2 = fp.second;
                Constraint coplanar;
                coplanar.type = Constraint::Coplanarity;
                coplanar.ent1 = f1;
                coplanar.ent2 = f2;
                coplanar.weight = 1.0;
                constraints.push_back(coplanar);
                int coplanarCon = constraints.size() - 1;
                ent2cons[coplanar.ent1].push_back(coplanarCon);
                ent2cons[coplanar.ent2].push_back(coplanarCon);
            }

          

            cg.rootEnt = -1;
            for (int i = 0; i < cg.entities.size(); i++) {
                if (cg.entities[i].supportingPlane.dof == 1) {
                    cg.rootEnt = i;
                }
                cg.determinableEnts.insert(i);
            }
            for (int i = 0; i < cg.constraints.size(); i++) {
                cg.consBetweenDeterminableEnts.insert(i);
            }

            return cg;
        }


        

    }
}