#pragma once

#include "../core/basic_types.hpp"
#include "../core/cameras.hpp"

#include "../gui/basic_types.hpp"

namespace pano {

    namespace experimental {

        using namespace pano::core;

        struct SegControl {
            int orientationClaz, orientationNotClaz; bool used;
            int dof() const {
                if (orientationClaz != -1) return 1;
                if (orientationNotClaz != -1) return 2;
                return 3;
            }
            bool operator == (const SegControl & c) const {
                return std::tie(orientationClaz, orientationNotClaz, used) == 
                    std::tie(c.orientationClaz, c.orientationNotClaz, used); 
            }
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(orientationClaz, orientationNotClaz, used);
            }
        };

        // whether a bnd is an occlusion
        enum class OcclusionRelation {
            Connected,
            LeftIsFront,
            RightIsFront,
            Unknown
        };

        // whether a line attaches a seg
        enum class AttachmentRelation {
            Attached,
            Detached,
            Unknown
        };

        // vertex in a constraint graph
        struct Vertex {
            bool isSeg; int id;
            Vertex() : isSeg(false), id(-1) {}
            Vertex(bool s, int i) : isSeg(s), id(i) {}
            bool operator < (const Vertex & u) const { return std::tie(isSeg, id) < std::tie(u.isSeg, u.id); }
            bool operator == (const Vertex & u) const { return std::tie(isSeg, id) == std::tie(u.isSeg, u.id); }
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(isSeg, id);
            }
        };

        // edge in a constraint graph
        struct Edge {
            int vert1, vert2;
            std::vector<Vec3> anchors;
            double weight;
            Edge() : vert1(-1), vert2(-1), weight(0.0) {}
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(vert1, vert2, anchors, weight);
            }
        };


        struct PIGraph {
            PanoramicView view;
            std::vector<Vec3> vps;
            int verticalVPId;

            // seg
            Imagei segs;
            int nsegs;
            std::vector<std::vector<int>> seg2bnds;
            std::vector<std::vector<int>> seg2linePieces;            
            std::vector<SegControl> seg2control;
            std::vector<double> seg2area;
            double fullArea;
            std::vector<Vec3> seg2center;
            std::vector<std::vector<std::vector<Vec3>>> seg2contours;

            // linePiece
            std::vector<std::vector<Vec3>> linePiece2samples;
            std::vector<double> linePiece2length;
            std::vector<int> linePiece2line;
            std::vector<int> linePiece2seg; // could be -1
            std::vector<AttachmentRelation> linePiece2attachment; // for linePice2seg only, as for linePIece2bndPiece, see the bndPiece2occlusion
            std::vector<int> linePiece2bndPiece; // could be -1 (either 2seg or 2bndPiece)
            std::vector<bool> linePiece2bndPieceInSameDirection;
            int nlinePieces() const { return linePiece2samples.size(); }

            // line
            std::vector<Classified<Line3>> lines;
            std::vector<std::vector<int>> line2linePieces;
            std::vector<std::vector<int>> line2lineRelations;
            int nlines() const { return lines.size(); }

            // lineRelation
            std::vector<Vec3> lineRelation2anchor;
            std::vector<std::pair<int, int>> lineRelation2lines;
            std::vector<double> lineRelation2weight;
            std::vector<bool> lineRelation2IsIncidence;
            int nlineRelations() const { return lineRelation2anchor.size(); }

            // bndPiece (a STRAIGHT boundary piece in a bnd)
            std::vector<std::vector<Vec3>> bndPiece2dirs; // continuous
            std::vector<double> bndPiece2length;
            std::vector<int> bndPiece2classes;
            std::vector<int> bndPiece2bnd;
            std::vector<std::vector<int>> bndPiece2linePieces;
            std::vector<OcclusionRelation> bndPiece2occlusion;
            int nbndPieces() const { return bndPiece2dirs.size(); }

            // bnd (a CONTINUOUS boundary between TWO segs)
            std::vector<std::vector<int>> bnd2bndPieces; // continuously connected
            std::vector<std::pair<int, int>> bnd2segs; // left , right
            std::vector<std::pair<int, int>> bnd2juncs; // from, to
            int nbnds() const { return bnd2bndPieces.size(); }

            // junc (junction of >=3 bnds)
            std::vector<Vec3> junc2positions;
            std::vector<std::vector<int>> junc2bnds;
            int njuncs() const { return junc2positions.size(); }

            // constraint graph
            std::vector<int> seg2vert;
            std::vector<int> line2vert;
            std::vector<Vertex> verts;
            std::vector<Edge> edges;
            int nccs;
            std::vector<int> vert2cc;
            std::map<int, std::vector<int>> cc2verts;
            std::vector<int> ccidsBigToSmall;

            // reconstructed
            std::vector<Plane3> seg2recPlanes;
            std::vector<Line3> line2recLines;


            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(view, vps, verticalVPId);
                ar(segs, nsegs, seg2bnds, seg2linePieces, seg2control, seg2area, fullArea, seg2center, seg2contours);
                ar(linePiece2samples, linePiece2length, linePiece2line, linePiece2seg, linePiece2attachment, linePiece2bndPiece, linePiece2bndPieceInSameDirection);
                ar(lines, line2linePieces, line2lineRelations);
                ar(lineRelation2anchor, lineRelation2lines, lineRelation2weight, lineRelation2IsIncidence);
                ar(bndPiece2dirs, bndPiece2length, bndPiece2classes, bndPiece2bnd, bndPiece2linePieces, bndPiece2occlusion);
                ar(bnd2bndPieces, bnd2segs, bnd2juncs);
                ar(junc2positions, junc2bnds);
                ar(seg2vert, line2vert, verts, edges, nccs, vert2cc, cc2verts, ccidsBigToSmall);
                ar(seg2recPlanes, line2recLines);
            }
        };

        PIGraph BuildPIGraph(const PanoramicView & view, const std::vector<Vec3> & vps, int verticalVPId,
            const Imagei & segs, const std::vector<Classified<Line3>> & lines,
            double bndPieceSplitAngleThres,
            double bndPieceClassifyAngleThres,
            double bndPieceBoundToLineAngleThres,
            double intersectionAngleThreshold,
            double incidenceAngleAlongDirectionThreshold,
            double incidenceAngleVerticalDirectionThreshold);








        // PerfectSegMaskView
        View<PartialPanoramicCamera, Imageub> PerfectSegMaskView(const PIGraph & mg, int seg, double focal = 100.0);

        // CollectFeatureMeanOnSegs
        template <class T, int N, class CameraT>
        std::vector<Vec<T, N>> CollectFeatureMeanOnSegs(const PIGraph & mg,
            const CameraT & pcam, const ImageOf<Vec<T, N>> & feature) {
            std::vector<Vec<T, N>> featureMeanTable(mg.nsegs);
            for (int i = 0; i < mg.nsegs; i++) {
                auto regionMaskView = PerfectSegMaskView(mg, i, 100.0);
                if (regionMaskView.image.empty()) {
                    continue;
                }
                auto sampler = MakeCameraSampler(regionMaskView.camera, pcam);
                auto featureOnRegion = sampler(feature);
                int votes = 0;
                Vec<T, N> featureSum;
                for (auto it = regionMaskView.image.begin(); it != regionMaskView.image.end(); ++it) {
                    if (!*it) {
                        continue;
                    }
                    featureSum += featureOnRegion(it.pos());
                    votes += 1;
                }
                auto featureMean = featureSum / std::max(votes, 1);
                featureMeanTable[i] = featureMean;
            }
            return featureMeanTable;
        }
    }

}