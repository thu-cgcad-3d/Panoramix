#include "../gui/scene.hpp"
#include "../core/algorithms.hpp"

#include "pi_graph_solve.hpp"

namespace pano {
    namespace experimental {

#if 0

        namespace {
            template <class T>
            inline std::pair<T, T> MakeOrderedPair(const T & a, const T & b) {
                return a < b ? std::make_pair(a, b) : std::make_pair(b, a);
            }
        }

        void BuildConstraintGraph(PIGraph & mg) {
            
            /// build the constraint graph
            std::cout << "building constraint graph" << std::endl;
            
            // verts
            mg.seg2vert.resize(mg.nsegs, -1);
            mg.line2vert.resize(mg.nlines(), -1);

            auto & verts = mg.verts;
            verts.clear();
            verts.reserve(mg.nsegs + mg.nlines()); 
            auto & seg2vert = mg.seg2vert;
            
            for (int seg = 0; seg < mg.nsegs; seg++) { 
                if (mg.seg2control[seg].used) {
                    verts.push_back({ true, seg });
                    seg2vert[seg] = verts.size() - 1;
                }
            }

            auto & line2vert = mg.line2vert;
            for (int line = 0; line < mg.nlines(); line++) { 
                verts.push_back({ false, line }); 
                line2vert[line] = verts.size() - 1;
            }

            // edges
            std::map<std::pair<int, int>, Edge> edgesMap; // vertids -> edge
            // seg - bnd - seg
            for (int bnd = 0; bnd < mg.nbnds(); bnd++) {
                int seg1 = mg.bnd2segs[bnd].first;
                int seg2 = mg.bnd2segs[bnd].second;
                if (seg1 == seg2) {
                    continue;
                }
                if (seg2vert[seg1] == -1 || seg2vert[seg2] == -1) {
                    continue;
                }
                auto & edge = edgesMap[MakeOrderedPair(seg2vert[seg1], seg2vert[seg2])];
                for (int bp : mg.bnd2bndPieces[bnd]) {
                    if (mg.bndPiece2segRelation[bp] != SegRelation::Connected ||
                        mg.bndPiece2segRelation[bp] != SegRelation::Coplanar) {
                        continue;
                    }
                    auto & anchors = mg.bndPiece2dirs[bp];
                    //edge.anchors.insert(edge.anchors.end(), anchors.begin(), anchors.end());
                    edge.anchors.push_back(anchors.front());
                    edge.anchors.push_back(anchors.back());
                    edge.weight += mg.bndPiece2length[bp];
                }                
            }
            // seg - bndpiece - linepiece - line
            for (int lp = 0; lp < mg.nlinePieces(); lp++) {
                int bp = mg.linePiece2bndPiece[lp];
                if (bp == -1) {
                    continue;
                }
                int line = mg.linePiece2line[lp];
                int bnd = mg.bndPiece2bnd[bp];
                int seg1 = mg.bnd2segs[bnd].first;
                int seg2 = mg.bnd2segs[bnd].second;
                SegRelation occ = mg.bndPiece2segRelation[bp];
                if ((occ == SegRelation::Connected || occ == SegRelation::LeftIsFront) && seg2vert[seg1] != -1) {
                    // connect seg1 and line
                    auto & edge = edgesMap[MakeOrderedPair(seg2vert[seg1], line2vert[line])];
                    auto & anchors = mg.bndPiece2dirs[bp];
                    //edge.anchors.insert(edge.anchors.end(), anchors.begin(), anchors.end());
                    edge.anchors.push_back(anchors.front());
                    edge.anchors.push_back(anchors.back());
                    edge.weight += mg.bndPiece2length[bp] * 30;
                }
                if ((occ == SegRelation::Connected || occ == SegRelation::RightIsFront) && seg2vert[seg2] != -1) {
                    // connect seg2 and line
                    auto & edge = edgesMap[MakeOrderedPair(seg2vert[seg2], line2vert[line])];
                    auto & anchors = mg.bndPiece2dirs[bp];
                    //edge.anchors.insert(edge.anchors.end(), anchors.begin(), anchors.end());
                    edge.anchors.push_back(anchors.front());
                    edge.anchors.push_back(anchors.back());
                    edge.weight += mg.bndPiece2length[bp] * 30;
                }
            }
            // seg - linepiece - line
            for (int lp = 0; lp < mg.nlinePieces(); lp++) {
                if (mg.linePiece2bndPiece[lp] != -1) {
                    continue;
                }
                int line = mg.linePiece2line[lp];
                int seg = mg.linePiece2seg[lp];
                if (seg2vert[seg] == -1) {
                    continue;
                }
                SegLineRelation att = mg.linePiece2segLineRelation[lp];
                if (att == SegLineRelation::Attached) {
                    // connect seg and line
                    auto & edge = edgesMap[MakeOrderedPair(seg2vert[seg], line2vert[line])];
                    auto & anchors = mg.linePiece2samples[lp];
                    //edge.anchors.insert(edge.anchors.end(), anchors.begin(), anchors.end());
                    edge.anchors.push_back(anchors.front());
                    edge.anchors.push_back(anchors.back());
                    edge.weight += mg.linePiece2length[lp] * 30;
                }
            }
            // line - linerelation - line
            for (int lr = 0; lr < mg.nlineRelations(); lr++) {
                if (mg.lineRelation2weight[lr] == 0.0) {
                    continue;
                }
                int line1 = mg.lineRelation2lines[lr].first;
                int line2 = mg.lineRelation2lines[lr].second;
                auto & edge = edgesMap[MakeOrderedPair(line2vert[line1], line2vert[line2])];
                edge.anchors.push_back(mg.lineRelation2anchor[lr]);
                edge.weight += mg.lineRelation2weight[lr] * 5;
            }

            // remove empty edges
            auto & edges = mg.edges;
            edges.clear();
            edges.reserve(edgesMap.size());
            for (auto & ep : edgesMap) {
                auto & e = ep.second;
                if (e.anchors.empty() || e.weight == 0.0 || ep.first.first == ep.first.second) {
                    continue;
                }
                edges.push_back(std::move(e));
                edges.back().vert1 = ep.first.first;
                edges.back().vert2 = ep.first.second;
            }


            /// cc decompose on constriant graph
            std::cout << "cc decomposition" << std::endl;
            auto & vert2cc = mg.vert2cc;
            vert2cc.clear();
            vert2cc.resize(verts.size(), -1);
            auto & ccid2verts = mg.cc2verts;
            ccid2verts.clear();

            std::vector<int> vertids(verts.size());
            std::iota(vertids.begin(), vertids.end(), 0);
            mg.nccs = core::ConnectedComponents(vertids.begin(), vertids.end(), [&verts, &edges](int vertid) {
                std::vector<int> neighbors;
                for (auto & e : edges) {
                    if (e.vert1 == vertid) {
                        neighbors.push_back(e.vert2);
                    } else if (e.vert2 == vertid) {
                        neighbors.push_back(e.vert1);
                    }
                }
                return neighbors;
            }, [&vert2cc, &ccid2verts](int vertid, int ccid) {
                vert2cc[vertid] = ccid;
                ccid2verts[ccid].push_back(vertid);
            });

            std::cout << "nccs = " << mg.nccs << std::endl;


            // order ccs from large to small
            auto & ccidsBigToSmall = mg.ccidsBigToSmall;
            ccidsBigToSmall.resize(mg.nccs);
            std::iota(ccidsBigToSmall.begin(), ccidsBigToSmall.end(), 0);
            std::sort(ccidsBigToSmall.begin(), ccidsBigToSmall.end(), [&ccid2verts](int a, int b) {
                return ccid2verts.at(a).size() > ccid2verts.at(b).size();
            });

        }

        namespace {
            int SwappedComponent(const Vec3 & orientation) {
                for (int i = 0; i < 2; i++) {
                    if (abs(orientation[i]) >= 1e-8) {
                        return i;
                    }
                }
                return 2;
            }
            inline double NonZeroize(double d) {
                return d == 0.0 ? 1e-6 : d;
            }
        }

        int VertexDof(const PIGraph & mg, int vert) {
            auto & v = mg.verts[vert];
            if (v.isSeg) {
                return mg.seg2control[v.id].dof();
            } else {
                return mg.lines[v.id].claz == -1 ? 2 : 1;
            }
        }

        std::vector<double> InverseDepthCoefficientsOfSegAtDirection(const PIGraph & mg,
            int seg, const Vec3 & direction) {
            auto & control = mg.seg2control[seg];
            if (control.dof() == 1) {
                assert(IsFuzzyZero(norm(mg.seg2center[seg]) - 1.0, 1e-2));
                Plane3 plane(mg.seg2center[seg], mg.vps[control.orientationClaz]);
                // variable is 1.0/centerDepth
                // corresponding coeff is 1.0/depthRatio
                // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                double depthRatio = norm(IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), direction), plane).position);
                return std::vector<double>{1.0 / depthRatio};
            } else if (control.dof() == 2) {
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                // -> 1.0/depth = ax + by + (- o1a - o2b)/o3 * z = a(x-o1*z/o3) + b(y-o2*z/o3)
                auto orientation = normalize(mg.vps[control.orientationNotClaz]);
                int sc = SwappedComponent(orientation); // find a non zero component
                Vec3 forientation = orientation;
                std::swap(forientation[sc], forientation[2]);
                Vec3 fdirection = direction;
                std::swap(fdirection[sc], fdirection[2]);
                return std::vector<double>{
                    fdirection[0] - forientation[0] * fdirection[2] / forientation[2],
                        fdirection[1] - forientation[1] * fdirection[2] / forientation[2]
                };
            } else {
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                return std::vector<double>{direction[0], direction[1], direction[2]};
            }
        }

        std::vector<double> InverseDepthCoefficientsOfLineAtDirection(const PIGraph & mg,
            int line, const Vec3 & direction) {
            int claz = mg.lines[line].claz;
            auto & l = mg.lines[line].component;
            if (claz >= 0 && claz < mg.vps.size()) {
                Ray3 infLine(normalize(l.center()), mg.vps[claz]);
                // variable is 1.0/centerDepth
                // corresponding coeff is 1.0/depthRatio
                // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                double depthRatio = norm(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), direction), infLine).second.first);
                return std::vector<double>{1.0 / depthRatio};
            } else {
                double theta = AngleBetweenDirections(normalize(l.first), normalize(l.second));
                double phi = AngleBetweenDirections(normalize(l.first), direction);
                /*           | sin(theta) | | p | | q |
                len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
                | p sin(phi) - q sin(phi - theta) |
                */
                // variables[0] -> 1/p
                // variables[1] -> 1/q
                double coeffFor1_p = -sin(phi - theta) / sin(theta);
                double coeffFor1_q = sin(phi) / sin(theta);
                assert(!IsInfOrNaN(coeffFor1_p) && !IsInfOrNaN(coeffFor1_q));
                return std::vector<double>{coeffFor1_p, coeffFor1_q};
            }
        }


        std::vector<double> InverseDepthCoefficientsAtDirection(const PIGraph & mg, 
            const Vertex & v, const Vec3 & direction) {
            assert(IsFuzzyZero(norm(direction) - 1.0, 1e-2));
            if (v.isSeg) {
                int seg = v.id;
                return InverseDepthCoefficientsOfSegAtDirection(mg, seg, direction);
            } else {
                int line = v.id;
                return InverseDepthCoefficientsOfSegAtDirection(mg, line, direction);
            }
        }

        DenseMatd SegPlaneEquationCoefficients(const PIGraph & mg,
            int seg) {
            auto & control = mg.seg2control[seg];
            if (control.dof() == 1) {
                auto & center = mg.seg2center[seg];
                auto & vp = mg.vps[control.orientationClaz];
                assert(IsFuzzyZero(norm(center) - 1.0, 1e-2));
                assert(IsFuzzyZero(norm(vp) - 1.0, 1e-2));
                // variable is 1.0/centerDepth
                //  let center = [cx, cy, cz], vp = [vx, vy, vz]
                //  then the equation should be k vx x + k vy y + k vz z = 1
                //  if the depth of center is Dc, then
                //      k vx Dc cx + k vy Dc cy + k vz Dc cz = 1
                //  and then k = 1.0 / (Dc*dot(v, c))
                //  the equation coefficients thus are
                //      k v* = v* / (Dc*dot(v, c))
                //  therefore the equation coefficients corresponding to the variable (inverse center depth, 1.0/Dc) are
                //      v* / dot(v, c)
                DenseMatd coeffs(3, 1, 0.0);
                for (int i = 0; i < 3; i++) {
                    coeffs(i, 0) = vp[i] / vp.dot(center);
                }
                return coeffs;
            } else if (control.dof() == 2) {
                auto orientation = normalize(mg.vps[control.orientationNotClaz]);
                int sc = SwappedComponent(orientation); // find a non zero component
                Vec3 forientation = orientation;
                std::swap(forientation[sc], forientation[2]);
                // (a, b, c) \perp o => o1a + o2b + o3c = 0 => c =  (- o1a - o2b)/o3
                DenseMatd coeffs(3, 2, 0.0);
                coeffs(0, 0) = coeffs(1, 1) = 1;
                coeffs(2, 0) = -forientation[0] / forientation[2];
                coeffs(2, 1) = -forientation[1] / forientation[2];
                for (int k = 0; k < 2; k++) {
                    std::swap(coeffs(sc, k), coeffs(2, k));
                }
                return coeffs;
            } else if (control.dof() == 3) {
                return DenseMatd::eye(3, 3);
            }
        }


        Plane3 SegInstance(const PIGraph & mg, const double * variables, int seg) {
            auto & c = mg.seg2control[seg];
            if (c.orientationClaz >= 0) {
                return Plane3(mg.seg2center[seg] / NonZeroize(variables[0]), mg.vps[c.orientationClaz]);
            } else if (c.orientationClaz == -1 && c.orientationNotClaz >= 0) {
                double vs[] = { variables[0], variables[1], 0.0 }; // fake vs
                // v1 * o1 + v2 * o2 + v3 * o3 = 0, since region.orientation is orthogonal to normal
                auto orientation = normalize(mg.vps[c.orientationNotClaz]);
                int c = SwappedComponent(orientation);
                std::swap(orientation[c], orientation[2]); // now fake orientation
                vs[2] = (-vs[0] * orientation[0] - vs[1] * orientation[1])
                    / orientation[2];
                std::swap(vs[c], vs[2]); // now real vs
                return Plane3FromEquation(vs[0], vs[1], vs[2]);
            } else /*if (region.type == MGUnary::RegionWithFixedNormal)*/ {
                return Plane3FromEquation(variables[0], variables[1], variables[2]);
            }
        }

        Line3 LineInstance(const PIGraph & mg, const double * variables, int line) {
            auto & l = mg.lines[line];
            if (l.claz >= 0 && l.claz < mg.vps.size()) {
                Ray3 infLine(normalize(l.component.center()) / NonZeroize(variables[0]), mg.vps[l.claz]);
                return Line3(DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), normalize(l.component.first)), infLine).second.second,
                    DistanceBetweenTwoLines(Ray3(Point3(0, 0, 0), normalize(l.component.second)), infLine).second.second);
            } else /*if (line.type == MGUnary::LineFree)*/ {
                /*           | sin(theta) | | p | | q |
                len:---------------------------------   => 1/len: [(1/q)sin(phi) - (1/p)sin(phi-theta)] / sin(theta)
                | p sin(phi) - q sin(phi - theta) |
                */
                // variables[0] -> 1/p
                // variables[1] -> 1/q
                return Line3(normalize(l.component.first) / NonZeroize(variables[0]), normalize(l.component.second) / NonZeroize(variables[1]));
            }
        }
        



        void SolvePIGraph(int ccid, PIGraph & mg, misc::Matlab & matlab, int tryNum) {

            auto & vertsInCC = mg.cc2verts.at(ccid);

            /// build varriables and equations

            // vert start position in variable vector
            std::vector<int> vert2varPosition(mg.verts.size(), -1);
            std::vector<int> vert2nvar(mg.verts.size(), -1);
            int nvars = 0;
            for (int i = 0; i < vertsInCC.size(); i++) {
                int vertid = vertsInCC[i];
                auto & vert = mg.verts[vertid];
                int nvar = 0;
                if (vert.isSeg) {
                    int seg = vert.id;
                    auto & control = mg.seg2control[seg];
                    nvar = control.dof();
                } else {
                    int line = vert.id;
                    int claz = mg.lines[line].claz;
                    nvar = claz == -1 ? 2 : 1;
                }
                vert2nvar[vertid] = nvar;
                vert2varPosition[vertid] = nvars;
                nvars += nvar;
            }

            std::vector<int> edgesInCC;
            for (int i = 0; i < mg.edges.size(); i++) {
                auto & edge = mg.edges[i];
                if (mg.vert2cc[edge.vert1] == ccid && mg.vert2cc[edge.vert2] == ccid) {
                    edgesInCC.push_back(i);
                }
            }

            // connectivity term
            // the elements of sparse matrix A
            // inverse depths of vert1s on each anchor A1 * X
            // inverse depths of vert2s on each anchor A2 * X
            // weight on each anchor W
            std::vector<SparseMatElementd> A1triplets, A2triplets;        
            std::vector<double> WA;
            int eidA = 0;
            for (int i = 0; i < edgesInCC.size(); i++) {
                auto & edge = mg.edges[edgesInCC[i]];
                auto & vert1 = mg.verts[edge.vert1];
                auto & vert2 = mg.verts[edge.vert2];
                int varpos1 = vert2varPosition[edge.vert1];
                int varpos2 = vert2varPosition[edge.vert2];
                
                double eachWeight = sqrt(1.0 - Gaussian(edge.weight, 0.1) / edge.anchors.size());
                for (auto & anchor : edge.anchors) {
                    auto coeffs1 = InverseDepthCoefficientsAtDirection(mg, vert1, anchor);
                    for (int k = 0; k < coeffs1.size(); k++) {
                        A1triplets.emplace_back(eidA, varpos1 + k, coeffs1[k]);
                    }
                    auto coeffs2 = InverseDepthCoefficientsAtDirection(mg, vert2, anchor);
                    for (int k = 0; k < coeffs2.size(); k++) {
                        A2triplets.emplace_back(eidA, varpos2 + k, coeffs2[k]);
                    }
                    WA.push_back(eachWeight * 1e3);
                    eidA++;
                }
            }
            int neqsA = eidA;
            auto A1 = MakeSparseMatFromElements(neqsA, nvars, A1triplets.begin(), A1triplets.end());
            auto A2 = MakeSparseMatFromElements(neqsA, nvars, A2triplets.begin(), A2triplets.end());
     


            // coplanarity term
            std::vector<SparseMatElementd> C1triplets, C2triplets;
            std::vector<double> WC;
            int eidC = 0;
            for (int i = 0; i < edgesInCC.size(); i++) {
                auto & edge = mg.edges[edgesInCC[i]];
                auto & vert1 = mg.verts[edge.vert1];
                auto & vert2 = mg.verts[edge.vert2];
                int varpos1 = vert2varPosition[edge.vert1];
                int varpos2 = vert2varPosition[edge.vert2];
                if (vert1.isSeg && vert2.isSeg) {
                    auto segPlaneCoeffs1 = SegPlaneEquationCoefficients(mg, vert1.id);
                    assert(segPlaneCoeffs1.rows == 3);
                    int dof1 = segPlaneCoeffs1.cols;
                    auto segPlaneCoeffs2 = SegPlaneEquationCoefficients(mg, vert2.id);
                    assert(segPlaneCoeffs2.rows == 3);
                    int dof2 = segPlaneCoeffs2.cols;

                    for (int k = 0; k < 3; k++) {
                        for (int s = 0; s < segPlaneCoeffs1.cols; s++) {
                            C1triplets.emplace_back(eidC, varpos1 + s, segPlaneCoeffs1(k, s));
                        }
                        for (int s = 0; s < segPlaneCoeffs2.cols; s++) {
                            C2triplets.emplace_back(eidC, varpos2 + s, segPlaneCoeffs2(k, s));
                        }
                        WC.push_back(Gaussian(edge.weight, 0.1) * (dof1 - 1) * (dof2 - 1));
                        eidC++;
                    }
                }
            }
            int neqsC = eidC;
            auto C1 = MakeSparseMatFromElements(neqsC, nvars, C1triplets.begin(), C1triplets.end());
            auto C2 = MakeSparseMatFromElements(neqsC, nvars, C2triplets.begin(), C2triplets.end());


            matlab << "clear;";
            matlab.setVar("A1", A1);
            matlab.setVar("A2", A2);
            matlab.setVar("C1", C1);
            matlab.setVar("C2", C2);
            matlab.setVar("WA", cv::Mat(WA));
            matlab.setVar("WC", cv::Mat(WC));

            matlab << "m = size(A1, 1);"; // number of equations
            matlab << "n = size(A1, 2);"; // number of variables

            matlab << "scale = 10;";           

            double minE = std::numeric_limits<double>::infinity();

            std::vector<double> X;
            for (int i = 0; i < tryNum; i++) {
                matlab << "D1D2 = ones(m, 1);"; //  current depths of anchors
                int repNum = 0;
                while(true) {
                    matlab << "K = (A1 - A2) .* repmat(D1D2 .* WA, [1, n]);";
                    matlab << "R = (C1 - C2) .* repmat(WC, [1, n]);";
                    matlab
                        << "cvx_begin"
                        << "variable X(n);"
                        << "minimize sum_square(K * X) + sum_square(R * X)"
                        << "subject to"
                        << "    ones(m, 1) / scale <= A1 * X <= ones(m, 1) * scale;"
                        << "    ones(m, 1) / scale <= A2 * X <= ones(m, 1) * scale;"
                        << "cvx_end";
                    matlab << "D1D2 = 1./ (A1 * X) ./ (A2 * X);";
                    matlab << "D1D2 = abs(D1D2) / norm(D1D2);";
                    //matlab << "P1P2 = "
                    matlab << "K = (A1 - A2) .* repmat(D1D2 .* WA, [1, n]);";

                    matlab << "e = sum_square(K * X) + sum_square((C1 - C2) * 0.1 * X);";
                    double curE = matlab.var("e").scalar();
                    std::cout << "e = " << curE << std::endl;
                    if (IsInfOrNaN(curE)) {
                        break;
                    }
                    if (curE < minE) {
                        minE = curE;
                        matlab << "X = 2 * X ./ median((A1 + A2) * X);";
                        X = matlab.var("X").toCVMat();
                    } else /*if(repNum > 10)*/ {
                        break;
                    }
                    repNum++;
                }
                if (IsInfOrNaN(minE)) {
                    matlab << "scale = scale * 2;";
                } else {
                    break;
                }
            }
         

            // install results
            for (int vert : vertsInCC) {
                int varpos = vert2varPosition[vert];
                auto & v = mg.verts[vert];
                if (v.isSeg) {
                    int seg = v.id;
                    mg.seg2recPlanes[seg] = SegInstance(mg, X.data() + varpos, seg);
                } else {
                    int line = v.id;
                    mg.line2recLines[line] = LineInstance(mg, X.data() + varpos, line);
                }
            }

        }







        std::vector<double> VariableCoefficientsForInverseDepthAtDirection(const std::vector<Vec3> & vps, const SegControl & control, 
            const Vec3 & center, const Vec3 & direction) {
            if (control.dof() == 1) {
                Plane3 plane(center, vps[control.orientationClaz]);
                // variable is 1.0/centerDepth
                // corresponding coeff is 1.0/depthRatio
                // so that 1.0/depth = 1.0/centerDepth * 1.0/depthRatio -> depth = centerDepth * depthRatio
                double depthRatio = norm(IntersectionOfLineAndPlane(Ray3(Point3(0, 0, 0), direction), plane).position);
                return std::vector<double>{1.0 / depthRatio};
            } else if (control.dof() == 2) {
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                // -> 1.0/depth = ax + by + (- o1a - o2b)/o3 * z = a(x-o1*z/o3) + b(y-o2*z/o3)
                auto orientation = normalize(vps[control.orientationNotClaz]);
                int sc = SwappedComponent(orientation);
                Vec3 forientation = orientation;
                std::swap(forientation[sc], forientation[2]);
                Vec3 fdirection = direction;
                std::swap(fdirection[sc], fdirection[2]);
                return std::vector<double>{
                    fdirection[0] - forientation[0] * fdirection[2] / forientation[2],
                        fdirection[1] - forientation[1] * fdirection[2] / forientation[2]
                };
            } else {
                // depth = 1.0 / (ax + by + cz) where (x, y, z) = direction, (a, b, c) = variables
                // -> 1.0/depth = ax + by + cz
                return std::vector<double>{direction[0], direction[1], direction[2]};
            }
        }


        Plane3 SegInstance(const std::vector<Vec3> & vps, const SegControl & control, const Vec3 & center, const double * variables) {
            auto & c = control;
            if (c.orientationClaz >= 0) {
                return Plane3(center / NonZeroize(variables[0]), vps[c.orientationClaz]);
            } else if (c.orientationClaz == -1 && c.orientationNotClaz >= 0) {
                double vs[] = { variables[0], variables[1], 0.0 }; // fake vs
                // v1 * o1 + v2 * o2 + v3 * o3 = 0, since region.orientation is orthogonal to normal
                auto orientation = normalize(vps[c.orientationNotClaz]);
                int c = SwappedComponent(orientation);
                std::swap(orientation[c], orientation[2]); // now fake orientation
                vs[2] = (-vs[0] * orientation[0] - vs[1] * orientation[1])
                    / orientation[2];
                std::swap(vs[c], vs[2]); // now real vs
                return Plane3FromEquation(vs[0], vs[1], vs[2]);
            } else {
                return Plane3FromEquation(variables[0], variables[1], variables[2]);
            }
        }


        void ReconstructLayoutAnnotation(PILayoutAnnotation & anno, misc::Matlab & matlab) {
            
            // group faces as planes
            std::vector<int> face2plane(anno.nfaces());
            std::map<int, std::set<int>> plane2faces;

            std::vector<int> faces(anno.nfaces());
            std::iota(faces.begin(), faces.end(), 0);
            int nplanes = core::ConnectedComponents(faces.begin(), faces.end(), [&anno](int face) {
                std::vector<int> coplanarFaces;
                for (auto & p : anno.coplanarFacePairs) {
                    if (p.first == face) {
                        coplanarFaces.push_back(p.second);
                    } else if (p.second == face) {
                        coplanarFaces.push_back(p.first);
                    }
                }
                return coplanarFaces;
            }, [&face2plane, &plane2faces](int face, int ccid) {
                face2plane[face] = ccid;
                plane2faces[ccid].insert(face);
            });


            // get plane properties
            std::vector<SegControl> plane2control(nplanes);
            std::vector<Vec3> plane2center(nplanes);
            std::vector<int> plane2varPosition(nplanes);
            std::vector<int> plane2nvar(nplanes);

            int nvars = 0;
            for (int i = 0; i < nplanes; i++) {
                SegControl control = { -1, -1, true };
                Vec3 center;
                for (int face : plane2faces[i]) {
                    auto & ctrl = anno.face2control[face];
                    if (ctrl.dof() < control.dof()) {
                        control = ctrl;                        
                    }
                    for (int c : anno.face2corners[face]) {
                        center += normalize(anno.corners[c]);
                    }
                }
                center /= norm(center);
                plane2control[i] = control;
                plane2center[i] = center;
                int nvar = control.dof();
                plane2nvar[i] = nvar;
                plane2varPosition[i] = nvars;
                nvars += nvar;
            }

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

            // the elements of sparse matrix A
            // inverse depths of vert1s on each anchor A1 * X
            // inverse depths of vert2s on each anchor A2 * X
            // weight on each anchor W
            std::vector<SparseMatElementd> A1triplets, A2triplets;
            int eid = 0;
            for (int i = 0; i < anno.nborders(); i++) {
                int plane1 = face2plane[border2face[i].first];
                int plane2 = face2plane[border2face[i].second];
                assert(plane1 != -1 && plane2 != -1);
                if (!anno.border2connected[i]) {
                    continue;
                }
                int varpos1 = plane2varPosition[plane1];
                int varpos2 = plane2varPosition[plane2];

                Vec3 corners[] = { anno.corners[anno.border2corners[i].first], anno.corners[anno.border2corners[i].second] };
                for (auto & anchor : corners) {
                    auto coeffs1 = VariableCoefficientsForInverseDepthAtDirection(anno.vps, 
                        plane2control[plane1], plane2center[plane1], anchor);
                    for (int k = 0; k < coeffs1.size(); k++) {
                        A1triplets.emplace_back(eid, varpos1 + k, coeffs1[k]);
                    }
                    auto coeffs2 = VariableCoefficientsForInverseDepthAtDirection(anno.vps,
                        plane2control[plane2], plane2center[plane2], anchor);
                    for (int k = 0; k < coeffs2.size(); k++) {
                        A2triplets.emplace_back(eid, varpos2 + k, coeffs2[k]);
                    }
                    eid++;
                }
            }

            int neqs = eid;

            auto A1 = MakeSparseMatFromElements(neqs, nvars, A1triplets.begin(), A1triplets.end());
            auto A2 = MakeSparseMatFromElements(neqs, nvars, A2triplets.begin(), A2triplets.end());

            matlab << "clear;";
            matlab.setVar("A1", A1);
            matlab.setVar("A2", A2);

            matlab << "m = size(A1, 1);"; // number of equations
            matlab << "n = size(A1, 2);"; // number of variables

            //matlab << "scale = 100.0;";

            double minE = std::numeric_limits<double>::infinity();

            std::vector<double> X;
            for (int i = 0; i < 100; i++) {
                matlab << "D1D2 = ones(m, 1);"; //  current depths of anchors
                while (true) {
                    matlab << "K = (A1 - A2) .* repmat(D1D2, [1, n]);";
                    matlab
                        << "cvx_begin"
                        << "variable X(n);"
                        << "minimize norm(K * X)"
                        << "subject to"
                        //<< "    ones(m, 1) <= A1 * X <= ones(m, 1) * scale;"
                        //<< "    ones(m, 1) <= A2 * X <= ones(m, 1) * scale;"
                        << "    ones(m, 1) <= A1 * X;"
                        << "    ones(m, 1) <= A2 * X;"
                        << "cvx_end";
                    matlab << "e = norm(K * X);";
                    double curE = matlab.var("e").scalar();
                    std::cout << "e = " << curE << std::endl;
                    if (IsInfOrNaN(curE)) {
                        break;
                    }
                    if (curE < minE) {
                        minE = curE;
                        matlab << "X = 2 * X ./ median((A1 + A2) * X);";
                        X = matlab.var("X").toCVMat();
                    } else {
                        break;
                    }
                    matlab << "D1D2 = 1./ (A1 * X) ./ (A2 * X);";
                    matlab << "D1D2 = abs(D1D2) / norm(D1D2);";
                }   
                if (IsInfOrNaN(minE)) {
                    matlab << "scale = scale * 2;";
                } else {
                    break;
                }
            }

            // install back as planes
            for (int i = 0; i < nplanes; i++) {
                Plane3 inst = SegInstance(anno.vps, plane2control[i], plane2center[i], X.data() + plane2varPosition[i]);
                for (int face : plane2faces[i]) {
                    anno.face2plane[face] = inst;
                }
            }          
        }











        void SolvePIGraphIteratively(PIGraph & mg, misc::Matlab & matlab) {

            // set all 

        }


#endif
    }
}