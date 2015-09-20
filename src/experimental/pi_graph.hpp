#pragma once

#include "../core/basic_types.hpp"
#include "../core/cameras.hpp"

namespace pano {

    namespace experimental {

        using namespace pano::core;

        struct SegControl {
            int orientationClaz, orientationNotClaz; bool used;
            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(orientationClaz, orientationNotClaz, used);
            }
        };

        enum class OcclusionRelation {
            FirstIsFront,
            SecondIsFront,
            Connected,
            Unknown
        };

        struct PIGraph {
            PanoramicView view;
            std::vector<Vec3> vps;

            // seg
            Imagei segs;
            int nsegs;
            std::vector<std::vector<int>> seg2bnds;
            std::vector<std::vector<int>> seg2linePieces;            
            std::vector<SegControl> seg2control;
            std::vector<double> seg2area;
            std::vector<Vec3> seg2center;
            std::vector<Plane3> seg2plane;
            std::vector<std::vector<std::vector<Vec3>>> seg2contours;

            // linePiece
            std::vector<std::vector<Vec3>> linePiece2samples;
            std::vector<int> linePiece2line;
            std::vector<int> linePiece2seg; // could be -1
            std::vector<bool> linePiece2used; // for linePice2seg only, as for linePIece2bndPiece, see the bndPiece2occlusion
            std::vector<int> linePiece2bndPiece; // could be -1 (either 2seg or 2bndPiece)
            std::vector<bool> linePiece2bndPieceInSameDirection;

            // line
            std::vector<Classified<Line3>> lines;
            std::vector<std::vector<int>> line2linePieces;
            std::vector<std::vector<int>> line2lineRelations;
            std::vector<Line3> line2reconstructed;

            // lineRelation
            std::vector<Vec3> lineRelation2anchor;
            std::vector<std::pair<int, int>> lineRelation2lines;
            std::vector<double> lineRelation2weight;
            std::vector<bool> lineRelation2IsIncidence;
            std::vector<bool> lineRelation2used;

            // bndPiece (a STRAIGHT boundary piece in a bnd)
            std::vector<std::vector<Vec3>> bndPiece2dirs; // continuous
            std::vector<double> bndPiece2length;
            std::vector<int> bndPiece2classes;
            std::vector<int> bndPiece2bnd;
            std::vector<std::vector<int>> bndPiece2linePieces;
            std::vector<std::vector<Vec3>> bndPiece2anchors;
            std::vector<OcclusionRelation> bndPiece2occlusion;

            // bnd (a CONTINUOUS boundary between TWO segs)
            std::vector<std::vector<int>> bnd2bndPieces; // continuously connected
            std::vector<std::pair<int, int>> bnd2segs; // left , right
            std::vector<std::pair<int, int>> bnd2juncs; // from, to

            // junc (junction of >=3 bnds)
            std::vector<Vec3> junc2positions;
            std::vector<std::vector<int>> junc2bnds;

            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(view, vps);
                ar(segs, nsegs, seg2bnds, seg2linePieces, seg2control, seg2area, seg2center, seg2plane, seg2contours);
                ar(linePiece2samples, linePiece2line, linePiece2seg, linePiece2used, linePiece2bndPiece, linePiece2bndPieceInSameDirection);
                ar(lines, line2linePieces, line2lineRelations, line2reconstructed);
                ar(lineRelation2anchor, lineRelation2lines, lineRelation2weight, lineRelation2IsIncidence, lineRelation2used);
                ar(bndPiece2dirs, bndPiece2length, bndPiece2classes, bndPiece2bnd, bndPiece2linePieces, bndPiece2anchors, bndPiece2occlusion);
                ar(bnd2bndPieces, bnd2segs, bnd2juncs);
                ar(junc2positions, junc2bnds);
            }

        };

        PIGraph BuildPIGraph(const PanoramicView & view, const std::vector<Vec3> & vps,
            const Imagei & segs, const std::vector<Classified<Line3>> & lines,
            double bndPieceSplitAngleThres,
            double bndPieceClassifyAngleThres,
            double bndPieceBoundToLineAngleThres,
            double intersectionAngleThreshold,
            double incidenceAngleAlongDirectionThreshold,
            double incidenceAngleVerticalDirectionThreshold);


    }

}