#pragma once

#include "../core/basic_types.hpp"
#include "../core/cameras.hpp"

#include "../gui/basic_types.hpp"

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







        // Print PIGraph
        template <
            class SegColorerT = core::ConstantFunctor<gui::ColorTag>,
            class LinePieceColorerT = core::ConstantFunctor<gui::ColorTag>,
            class BndPieceColorerT = core::ConstantFunctor<gui::ColorTag>
        >
        inline Image3f Print(const PIGraph & mg,
            SegColorerT && segColor = gui::Transparent,
            LinePieceColorerT && lpColor = gui::Transparent,
            BndPieceColorerT && bpColor = gui::Transparent,
            int boundaryWidth = 1, int lineWidth = 2) {
            Image3f rendered = Image3f::zeros(mg.segs.size());
            // segs
            for (auto it = rendered.begin(); it != rendered.end(); ++it) {
                int seg = mg.segs(it.pos());
                gui::Color color = segColor(seg);
                *it = Vec3f(color.bluef(), color.greenf(), color.redf());
            }
            // lines
            if (lineWidth > 0) {
                for (int lp = 0; lp < mg.linePiece2line.size(); lp++) {
                    gui::Color color = lpColor(lp);
                    if (color.isTransparent())
                        continue;
                    auto & ps = mg.linePiece2samples[lp];
                    for (int i = 1; i < ps.size(); i++) {
                        auto p1 = ToPixel(mg.view.camera.toScreen(ps[i - 1]));
                        auto p2 = ToPixel(mg.view.camera.toScreen(ps[i]));
                        if (Distance(p1, p2) >= rendered.cols / 2) {
                            continue;
                        }
                        cv::clipLine(cv::Rect(0, 0, rendered.cols, rendered.rows), p1, p2);
                        cv::line(rendered, p1, p2, (cv::Scalar)color / 255.0, lineWidth);
                    }
                }
            }
            // region boundary
            if (boundaryWidth > 0) {
                for (int bp = 0; bp < mg.bndPiece2bnd.size(); bp ++) {
                    gui::Color color = bpColor(bp);
                    if (color.isTransparent())
                        continue;
                    auto & e = mg.bndPiece2dirs[bp];
                    for (int i = 1; i < e.size(); i++) {
                        auto p1 = core::ToPixel(mg.view.camera.toScreen(e[i - 1]));
                        auto p2 = core::ToPixel(mg.view.camera.toScreen(e[i]));
                        if (Distance(p1, p2) >= rendered.cols / 2) {
                            continue;
                        }
                        cv::clipLine(cv::Rect(0, 0, rendered.cols, rendered.rows), p1, p2);
                        cv::line(rendered, p1, p2, (cv::Scalar)color / 255.0, boundaryWidth);
                    }
                }
            }
            return rendered;
        }


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