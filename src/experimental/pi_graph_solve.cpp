#include "../gui/scene.hpp"
#include "../core/algorithms.hpp"

#include "pi_graph_solve.hpp"

namespace pano {
    namespace experimental {

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
                auto & anchors = edge.anchors;
                double & weight = edge.weight;
                for (int bp : mg.bnd2bndPieces[bnd]) {
                    if (mg.bndPiece2occlusion[bp] != OcclusionRelation::Connected) {
                        continue;
                    }
                    anchors.push_back(mg.bndPiece2dirs[bp].front());
                    anchors.push_back(mg.bndPiece2dirs[bp].back());
                    weight += mg.bndPiece2length[bp];
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
                OcclusionRelation occ = mg.bndPiece2occlusion[bp];
                if ((occ == OcclusionRelation::Connected || occ == OcclusionRelation::LeftIsFront) && seg2vert[seg1] != -1) {
                    // connect seg1 and line
                    auto & edge = edgesMap[MakeOrderedPair(seg2vert[seg1], line2vert[line])];
                    edge.anchors.push_back(mg.bndPiece2dirs[bp].front());
                    edge.anchors.push_back(mg.bndPiece2dirs[bp].back());
                    edge.weight += mg.bndPiece2length[bp];
                }
                if ((occ == OcclusionRelation::Connected || occ == OcclusionRelation::RightIsFront) && seg2vert[seg2] != -1) {
                    // connect seg2 and line
                    auto & edge = edgesMap[MakeOrderedPair(seg2vert[seg2], line2vert[line])];
                    edge.anchors.push_back(mg.bndPiece2dirs[bp].front());
                    edge.anchors.push_back(mg.bndPiece2dirs[bp].back());
                    edge.weight += mg.bndPiece2length[bp];
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
                AttachmentRelation att = mg.linePiece2attachment[lp];
                if (att == AttachmentRelation::Attached) {
                    // connect seg and line
                    auto & edge = edgesMap[MakeOrderedPair(seg2vert[seg], line2vert[line])];
                    edge.anchors.push_back(mg.linePiece2samples[lp].front());
                    edge.anchors.push_back(mg.linePiece2samples[lp].back());
                    edge.weight += mg.linePiece2length[lp];
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
                edge.weight += mg.lineRelation2weight[lr];
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

        std::vector<double> VariableCoefficientsForInverseDepthAtDirection(const PIGraph & mg, 
            const Vertex & v, const Vec3 & direction) {
            if (v.isSeg) {
                int seg = v.id;
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
            } else {
                int line = v.id;
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

            // the elements of sparse matrix A
            // inverse depths of vert1s on each anchor A1 * X
            // inverse depths of vert2s on each anchor A2 * X
            // weight on each anchor W
            std::vector<SparseMatElementd> A1triplets, A2triplets;
            std::vector<double> W;
            int eid = 0;
            for (int i = 0; i < edgesInCC.size(); i++) {
                auto & edge = mg.edges[edgesInCC[i]];
                auto & vert1 = mg.verts[edge.vert1];
                auto & vert2 = mg.verts[edge.vert2];
                int varpos1 = vert2varPosition[edge.vert1];
                int varpos2 = vert2varPosition[edge.vert2];
                
                for (auto & anchor : edge.anchors) {
                    auto coeffs1 = VariableCoefficientsForInverseDepthAtDirection(mg, vert1, anchor);
                    for (int k = 0; k < coeffs1.size(); k++) {
                        A1triplets.emplace_back(eid, varpos1 + k, coeffs1[k]);
                    }
                    auto coeffs2 = VariableCoefficientsForInverseDepthAtDirection(mg, vert2, anchor);
                    for (int k = 0; k < coeffs2.size(); k++) {
                        A2triplets.emplace_back(eid, varpos2 + k, coeffs2[k]);
                    }
                    W.push_back(edge.weight);
                    eid++;
                }
            }

            int neqs = eid;

            auto A1 = MakeSparseMatFromElements(neqs, nvars, A1triplets.begin(), A1triplets.end());
            auto A2 = MakeSparseMatFromElements(neqs, nvars, A2triplets.begin(), A2triplets.end());

            matlab << "clear;";
            matlab.setVar("A1", A1);
            matlab.setVar("A2", A2);
            matlab.setVar("W", cv::Mat(W));

            matlab << "m = size(A1, 1);"; // number of equations
            matlab << "n = size(A1, 2);"; // number of variables

            matlab << "scale = 1.0;";
            matlab << "D1D2 = ones(m, 1);"; //  current depths of anchors

            std::vector<double> X;
            for (int i = 0; i < tryNum; i++) {
                matlab << "K = (A1 - A2) .* repmat(D1D2 .* W, [1, n]);";
                matlab
                    << "cvx_begin"
                    << "variable X(n);"
                    << "minimize norm(K * X)"
                    << "subject to"
                    << "    ones(m, 1) <= A1 * X <= ones(m, 1) * scale;"
                    << "    ones(m, 1) <= A2 * X <= ones(m, 1) * scale;"
                    << "cvx_end";

                matlab << "scale = scale * 1.5;";
                //matlab << "D1D2 = 1./ (A1 * X) ./ (A2 * X);";
                //matlab << "K = (A1 - A2) .* repmat(D1D2 .* W, [1, n]);";
                matlab << "e = norm(K * X);";
                double e = matlab.var("e").scalar();
                if (!IsInfOrNaN(e)) {
                    std::cout << "e = " << e << std::endl;
                    X = matlab.var("X").toCVMat();
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

    }
}