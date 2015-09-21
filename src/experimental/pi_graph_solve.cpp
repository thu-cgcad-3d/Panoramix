#include "../core/algorithms.hpp"
#include "pi_graph_solve.hpp"

namespace pano {
    namespace experimental {

        void SolvePIGraph(PIGraph & mg, misc::Matlab & matlab) {
            
            // cc decompose, according to 1) occlusion status of bndPieces and 2) weight of lineRelations
            
            mg.seg2ccid.resize(mg.nsegs, -1);
            mg.line2ccid.resize(mg.nlines(), -1);

            
            // build the constraint graph
            struct Vertex {
                bool isSeg; int id; 
                bool operator < (const Vertex & u) const { return std::tie(isSeg, id) < std::tie(u.isSeg, u.id); }
                bool operator == (const Vertex & u) const { return std::tie(isSeg, id) == std::tie(u.isSeg, u.id); }
            };

            // verts
            std::vector<Vertex> verts;
            verts.reserve(mg.nsegs + mg.nlines()); 
            std::vector<int> seg2vertid(mg.nsegs, -1);
            for (int seg = 0; seg < mg.nsegs; seg++) { 
                if (mg.seg2control[seg].used) {
                    verts.push_back({ true, seg });
                    seg2vertid[seg] = verts.size() - 1;
                }
            }
            std::vector<int> line2vertid(mg.nlines(), -1);
            for (int line = 0; line < mg.nlines(); line++) { 
                verts.push_back({ false, line }); 
                line2vertid[line] = verts.size() - 1;
            }

            // edges
            std::map<std::pair<int, int>, std::vector<Vec3>> edges;
            // seg - bnd - seg
            for (int bnd = 0; bnd < mg.nbnds(); bnd++) {

            }

            // seg - bndpiece - linepiece - line


            // seg - linepiece - line


            // line - linerelation - line




            /*mg.nccs = core::ConnectedComponents(verts.begin(), verts.end(), [&mg](const Vertex & vert) {
                std::vector<Vertex> neighbors;
                if (vert.isSeg) {
                    int seg = vert.id;
                    auto & linePieces = mg.seg2linePieces[seg];
                    for (int lp : linePieces) {
                        if (!mg.linePiece2used[lp]) continue;
                        int line = mg.linePiece2line[lp];
                        neighbors.push_back({ false, line });
                    }
                    auto & bnds = mg.seg2bnds[seg];
                    for (int bnd : bnds) {
                        int anotherSeg = mg.bnd2segs[bnd].first;
                        if (anotherSeg == seg) {
                            anotherSeg = mg.bnd2segs[bnd].second;
                        }
                        auto & bps = mg.bnd2bndPieces[bnd];
                        if (std::any_of(bps.begin(), bps.end(), [&mg](int bp) {
                            return mg.bndPiece2occlusion[bp] == OcclusionRelation::Connected; 
                        })) {
                            neighbors.push_back({ true, anotherSeg });
                        }
                    }
                } else {

                }
            }, [&mg](const Vertex & vert, int ccid) {
            
            });*/

        }

    }
}