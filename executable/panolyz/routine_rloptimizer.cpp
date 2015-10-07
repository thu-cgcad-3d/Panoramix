#include "../../src/core/meta.hpp"
#include "../../src/experimental/rloptimizer.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/canvas.hpp"
#include "../../src/gui/utility.hpp"
#include "../../src/core/cameras.hpp"
#include "../../src/core/clock.hpp"

#include "routines.hpp"

namespace panolyz {

    using namespace pano;
    using namespace core;
    using namespace experimental;

    // output a 2d image of this rl graph
    template <
        class CameraT,
        class RhColorerT = core::ConstantFunctor<gui::ColorTag>,
        class LhColorerT = core::ConstantFunctor<gui::ColorTag>,
        class RRhColorerT = core::ConstantFunctor<gui::ColorTag>
    >
    Image3f Print(const RLGraph & mg, const Imagei & segmentedRegions, const CameraT & camera, const std::vector<RHandle> & rid2rhs,
        RhColorerT && rhColor = gui::Transparent,
        LhColorerT && lhColor = gui::Transparent,
        RRhColorerT && rbColor = gui::Transparent,
        int boundaryWidth = 1, int lineWidth = 2) {
        Image3f rendered = Image3f::zeros(segmentedRegions.size());
        // regions
        for (auto it = rendered.begin(); it != rendered.end(); ++it) {
            RHandle rh = rid2rhs[segmentedRegions(it.pos())];
            if (rh.invalid())
                continue;
            gui::Color color = rhColor(rh);
            *it = Vec3f(color.bluef(), color.greenf(), color.redf());
        }
        // lines
        if (lineWidth > 0) {
            for (auto & ld : mg.components<LData>()) {
                gui::Color color = lhColor(ld.topo.hd);
                if (color.isTransparent())
                    continue;
                static const double sampleAngle = M_PI / 100.0;
                auto & line = ld.data.normalizedLine;
                double spanAngle = AngleBetweenDirections(line.first, line.second);
                std::vector<Pixel> ps; ps.reserve(spanAngle / sampleAngle);
                for (double angle = 0.0; angle <= spanAngle; angle += sampleAngle) {
                    Vec3 dir = RotateDirection(line.first, line.second, angle);
                    ps.push_back(ToPixel(camera.toScreen(dir)));
                }
                for (int i = 1; i < ps.size(); i++) {
                    auto & p1 = ps[i - 1];
                    auto & p2 = ps[i];
                    cv::clipLine(cv::Rect(0, 0, rendered.cols, rendered.rows), p1, p2);
                    cv::line(rendered, p1, p2, (cv::Scalar)color, lineWidth);
                }
            }
        }
        // region boundary
        if (boundaryWidth > 0) {
            for (auto & b : mg.constraints<RRData>()) {
                gui::Color color = rbColor(b.topo.hd);
                if (color.isTransparent())
                    continue;
                for (auto & e : b.data.normalizedEdges) {
                    if (e.size() <= 1) continue;
                    for (int i = 1; i < e.size(); i++) {
                        auto p1 = core::ToPixel(camera.toScreen(e[i - 1]));
                        auto p2 = core::ToPixel(camera.toScreen(e[i]));
                        cv::clipLine(cv::Rect(0, 0, rendered.cols, rendered.rows), p1, p2);
                        cv::line(rendered, p1, p2, (cv::Scalar)color, boundaryWidth);
                    }
                }
            }
        }
        return rendered;
    }



    namespace RLOpt {


        void Run() {

            bool refresh = false;

            std::string path;
            Image3ub image = gui::PickAnImage("H:\\GitHub\\Panoramix\\data\\panorama\\indoor", &path);
            ResizeToMakeHeightUnder(image, 600);

            misc::Matlab matlab;

            View<PanoramicCamera, Image3ub> view = CreatePanoramicView(image);

            std::vector<PerspectiveCamera> cams;
            std::vector<std::vector<Classified<Line2>>> lines;
            std::vector<Vec3> vps;
            std::vector<DenseMatd> lineVPScores;
            if (refresh || !misc::LoadCache(path, "vplines", cams, lines, vps, lineVPScores)) {
                cams = CreateCubicFacedCameras(view.camera);
                lines.resize(cams.size());
                for (int i = 0; i < cams.size(); i++) {
                    auto pim = view.sampled(cams[i]).image;
                    LineSegmentExtractor lineExtractor;
                    lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
                    auto ls = lineExtractor(pim, 2, 300); // use pyramid
                    lines[i].reserve(ls.size());
                    for (auto & l : ls) {
                        lines[i].push_back(ClassifyAs(l, -1));
                    }
                }
                // estimate vp            
                vps = EstimateVanishingPointsAndClassifyLines(cams, lines, &lineVPScores);
                OrderVanishingPoints(vps, view.camera.up());
                if (0) {
                    auto ctable = gui::CreateRandomColorTableWithSize(vps.size());
                    for (int i = 0; i < cams.size(); i++) {
                        auto pim = view.sampled(cams[i]).image;
                        gui::AsCanvas(pim).thickness(3).colorTable(ctable).add(lines[i]).show();
                    }
                }
                misc::SaveCache(path, "vplines", cams, lines, vps, lineVPScores);
            }

            // segs
            Imagei segs;
            int nsegs;
            if (refresh || !misc::LoadCache(path, "segs", segs, nsegs)) {
                SegmentationExtractor seger;
                seger.params().algorithm = SegmentationExtractor::GraphCut;
                std::tie(segs, nsegs) = seger(image, true);
                misc::SaveCache(path, "segs", segs, nsegs);
            }       
            

            std::vector<PerspectiveCamera> hcams;
            if (refresh || !misc::LoadCache(path, "hcams", hcams)) {
                hcams = CreatePanoContextCameras(view.camera, 500, 400, 300);
                misc::SaveCache(path, "hcams", hcams);
            }

            // occ
            std::vector<Scored<Chain3>> occbnds;
            if (refresh || !misc::LoadCache(path, "occ", occbnds)) {
                for (int i = 0; i < hcams.size(); i++) {
                    auto pim = view.sampled(hcams[i]).image;
                    auto occ = AsDimensionConvertor(hcams[i]).toSpace(DetectOcclusionBoundary(matlab, pim));
                    for (auto & soc : occ) {
                        occbnds.push_back(soc);
                    }
                }
                misc::SaveCache(path, "occ", occbnds);
            }

            // extract gcs
            std::vector<Weighted<View<PerspectiveCamera, Image5d>>> gcs;
            if (refresh || !misc::LoadCache(path, "gcs", gcs)) {
                gcs.resize(hcams.size());
                for (int i = 0; i < hcams.size(); i++) {
                    auto pim = view.sampled(hcams[i]);
                    auto pgc = ComputeGeometricContext(matlab, pim.image, false, true);
                    gcs[i].component.camera = hcams[i];
                    gcs[i].component.image = pgc;
                    gcs[i].score = sin(AngleBetweenUndirectedVectors(hcams[i].forward(), view.camera.up()));
                }
                misc::SaveCache(path, "gcs", gcs);
            }


            // gc
            Image5d gc;
            if (refresh || !misc::LoadCache(path, "gc", gc)) {
                gc = Combine(view.camera, gcs).image;
                misc::SaveCache(path, "gc", gc);
            }

            // graph
            OrderVanishingPoints(vps, view.camera.up());

            RLGraph g;
            std::vector<RHandle> segId2Rhs;
            if (refresh || !misc::LoadCache(path, "g", g, segId2Rhs)) {
                RLGraphBuilder gb;
                std::vector<Chain3> occs;
                for (auto & soc : occbnds) {
                    if (soc.score > 0) {
                        occs.push_back(soc.component);
                    }
                }
                g = gb(vps, lines, lineVPScores, cams, segs, gc, occs, view.camera, &segId2Rhs);
                misc::SaveCache(path, "g", g, segId2Rhs);
            }
            {
                auto pim = Print(g, segs, view.camera, segId2Rhs, 
                    core::ConstantFunctor<gui::ColorTag>(gui::Transparent),
                    core::ConstantFunctor<gui::ColorTag>(gui::Transparent), 
                    [&g](RRHandle rrh) {
                    return g.data(rrh).occDetectionResult == RRData::NotConnected ? gui::Yellow : gui::Transparent;
                }, 2);
                cv::imshow("occ", pim);
                cv::waitKey();
            }

            // first time reconstruct
            auto pp = MakeASinglePatch(g);
          
            RLGraphPatchDict<RLGraphLabel> initialLabels;
            FillSingleLabels(g, vps, pp, initialLabels);
            std::unordered_map<RHandle, Plane3> recPlanes;
            std::unordered_map<LHandle, Line3> recLines;
            Reconstruct(g, vps, pp, initialLabels, recPlanes, recLines);
            Visualize(g, view, recPlanes, recLines);


            // make patches
            std::vector<RLGraphPatch> patches;
            MakeIcosahedron([&patches, &g, &view, &segId2Rhs, &segs](double x, double y, double z) {
                Vec3 dir(x, y, z);
                auto p2 = ToPixel(view.camera.toScreen(dir));
                if (Contains(segs, p2)) {
                    RHandle rh = segId2Rhs[segs(p2)];
                    patches.push_back(SamplePatch(g, rh, 10, DegreesToRadians(20)));
                }
                return 1;
            });


            RLOptimizer opt;
            opt.setup(g, vps, patches);
            opt.preprocess();
            misc::SaveCache(path, "opt_fea", opt.graphFeatureDictTable(), opt.patchFeatureDictTable());


            //opt.inference(opt.suggestedWeights(), recPlanes, recLines);


        }


    }
}