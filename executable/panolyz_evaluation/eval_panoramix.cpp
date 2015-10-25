
#include <thread>
#include "eval.hpp"

#include "../../src/experimental/rl_graph_solver.hpp"
#include "../../src/experimental/rl_graph_occlusion.hpp"
#include "../../src/experimental/rl_graph_modeling.hpp"
#include "../../src/experimental/rl_graph_vis.hpp"

#include "../../src/experimental/pi_graph_control.hpp"
#include "../../src/experimental/pi_graph_occlusion.hpp"
#include "../../src/experimental/pi_graph_cg.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"

namespace panolyz {

    template <class FunT>
    void ParallelRun(int n, FunT && fun) {
        std::vector<std::thread> threads;
        for (int i = 0; i < n; i++) {
            threads.emplace_back(std::forward<FunT>(fun), i);
        }
        for (auto & t : threads) {
            t.join();
        }
    }



    std::string PanoramixResultFilePath(const std::string & impath, int version) {
        return impath + ".panoramix." + std::to_string(version) + ".cereal";
    }

    // pi graph

    struct PIGraphModel : ReconstructedModel {
        PIGraph mg;
        explicit PIGraphModel(PIGraph && g) : mg(std::move(g)) {}

        virtual double depthAt(const Vec3 & direction) const {
            auto p = ToPixel(mg.view.camera.toScreen(direction));
            p.x = WrapBetween(p.x, 0, mg.segs.cols);
            p.y = BoundBetween(p.y, 0, mg.segs.rows - 1);
            int seg = mg.segs(p);
            NOT_IMPLEMENTED_YET();
            //return norm(IntersectionOfLineAndPlane(Ray3(Origin(), direction), mg.seg2recPlanes[seg]).position);
        }
        virtual void visualize(const std::vector<Vec3> & directions) const {
            //VisualizeReconstruction({ mg.ccidsBigToSmall.front() }, mg);
        }
    };


    std::unique_ptr<ReconstructedModel> PredictionOfPanoramix(const std::string & impath,
        const PredictOptions & options, misc::Matlab & matlab, const PILayoutAnnotation & anno) {

        std::cout << "folder is " << misc::FolderOfFile(impath) << std::endl;

        auto resultPath = PanoramixResultFilePath(impath, 0);

        PIGraph mg;
        if (1 || !LoadFromDisk(resultPath, mg)) {
            // rebuild pigraph
            misc::Matlab matlab;

            Image3ub image;
            bool extendedOnTop, extendedOnBottom;
            static const bool automaticallyRectifyIfNeeded = false;
            if (0 || !misc::LoadCache(impath, "rectified", image, extendedOnTop, extendedOnBottom)) {
                image = cv::imread(impath);
                if (!automaticallyRectifyIfNeeded) {
                    if (!gui::MakePanoramaByHand(image, &extendedOnTop, &extendedOnBottom)) {
                        WARNNING("failed making panorama");
                    }
                } else {
                    if (!MakePanorama(image, -1, &extendedOnTop, &extendedOnBottom)) {
                        WARNNING("failed making panorama");
                    }
                }
                misc::SaveCache(impath, "rectified", image, extendedOnTop, extendedOnBottom);
            }



            //Image3ub originalImage = anno.rectifiedImage;

            //auto image = originalImage.clone();
            ResizeToHeight(image, 700);

            /// prepare things!
            View<PanoramicCamera, Image3ub> view;
            std::vector<PerspectiveCamera> cams;
            std::vector<std::vector<Classified<Line2>>> rawLine2s;
            std::vector<Classified<Line3>> line3s;
            std::vector<Vec3> vps;
            int vertVPId;
            Imagei segs;
            int nsegs;


            bool refresh_preparation = false;
            if (refresh_preparation || !misc::LoadCache(impath, "preparation", view, cams, rawLine2s, line3s, vps, vertVPId, segs, nsegs)) {
                view = CreatePanoramicView(image);

                // collect lines in each view
                cams = CreateCubicFacedCameras(view.camera, image.rows, image.rows, image.rows * 0.4);
                std::vector<Line3> rawLine3s;
                rawLine2s.resize(cams.size());
                for (int i = 0; i < cams.size(); i++) {
                    auto pim = view.sampled(cams[i]).image;
                    LineSegmentExtractor lineExtractor;
                    lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
                    auto ls = lineExtractor(pim); // use pyramid
                    rawLine2s[i] = ClassifyEachAs(ls, -1);
                    for (auto & l : ls) {
                        rawLine3s.emplace_back(normalize(cams[i].toSpace(l.first)),
                            normalize(cams[i].toSpace(l.second)));
                    }
                }
                rawLine3s = MergeLines(rawLine3s, DegreesToRadians(3), DegreesToRadians(5));

                // estimate vp
                line3s = ClassifyEachAs(rawLine3s, -1);
                vps = EstimateVanishingPointsAndClassifyLines(line3s, nullptr, true);
                vertVPId = NearestDirectionId(vps, Vec3(0, 0, 1));


                if (0) {
                    gui::ColorTable ctable = gui::RGBGreys;
                    for (int i = 0; i < cams.size(); i++) {
                        auto & cam = cams[i];
                        std::vector<Classified<Line2>> lines;
                        for (auto & l3 : line3s) {
                            if (!cam.isVisibleOnScreen(l3.component.first) || !cam.isVisibleOnScreen(l3.component.second)) {
                                continue;
                            }
                            auto p1 = cam.toScreen(l3.component.first);
                            auto p2 = cam.toScreen(l3.component.second);
                            lines.push_back(ClassifyAs(Line2(p1, p2), l3.claz));
                        }
                        auto pim = view.sampled(cams[i]).image;
                        gui::AsCanvas(pim).thickness(3).colorTable(ctable).add(lines).show();
                    }
                }


                // estimate segs
                nsegs = SegmentationForPIGraph(view, line3s, segs, DegreesToRadians(1));
                RemoveThinRegionInSegmentation(segs, 1, true);
                nsegs = DensifySegmentation(segs, true);
                assert(IsDenseSegmentation(segs));

                if (0) {
                    auto ctable = gui::CreateGreyColorTableWithSize(nsegs);
                    ctable.randomize();
                    gui::ColorTable rgb = gui::RGBGreys;
                    auto canvas = gui::MakeCanvas(view.image).alpha(0.9).add(ctable(segs));
                    for (auto & l : line3s) {
                        static const double sampleAngle = M_PI / 100.0;
                        auto & line = l.component;
                        double spanAngle = AngleBetweenDirections(line.first, line.second);
                        std::vector<Point2> ps; ps.reserve(spanAngle / sampleAngle);
                        for (double angle = 0.0; angle <= spanAngle; angle += sampleAngle) {
                            Vec3 dir = RotateDirection(line.first, line.second, angle);
                            ps.push_back(view.camera.toScreen(dir));
                        }
                        for (int i = 1; i < ps.size(); i++) {
                            auto & p1 = ps[i - 1];
                            auto & p2 = ps[i];
                            if (Distance(p1, p2) >= view.image.cols / 2) {
                                continue;
                            }
                            canvas.thickness(2);
                            canvas.colorTable(rgb).add(gui::ClassifyAs(Line2(p1, p2), l.claz));
                        }
                    }
                    canvas.show();
                }

                // save
                misc::SaveCache(impath, "preparation", view, cams, rawLine2s, line3s, vps, vertVPId, segs, nsegs);
            }



            // om
            Image3ub omap;
            bool refresh_omap = false;
            if (false/*refresh_omap || !misc::LoadCache(impath, "omap", omap)*/) {
                // estimate omaps
                std::vector<View<PerspectiveCamera, Image3ub>> omaps(cams.size());
                ParallelRun(cams.size(), [&](int i) {
                    std::vector<HPoint2> hvps(vps.size());
                    for (int k = 0; k < vps.size(); k++) {
                        hvps[k] = cams[i].toScreenInHPoint(vps[k]);
                    }
                    // collect line 2ds in this view
                    ClassifyLines(rawLine2s[i], hvps);

                    Imagei omap = ComputeOrientationMaps(rawLine2s[i], hvps, cams[i].screenSize());
                    std::vector<Imageub> omapChannels = { omap == 0, omap == 1, omap == 2 };
                    Image3ub merged;
                    cv::merge(omapChannels, merged);
                    omaps[i] = View<PerspectiveCamera, Image3ub>(merged, cams[i]);
                });
                auto panoOmap = Combine(view.camera, omaps);
                for (auto & p : panoOmap.image) {
                    if ((p.val[0] > 0) + (p.val[1] > 0) + (p.val[2] > 0) > 1) {
                        std::fill(p.val, std::end(p.val), 0);
                    }
                }
                omap = panoOmap.image;
                misc::SaveCache(impath, "omap", omap);
            }



            // gc !!!!
            std::vector<PerspectiveCamera> hcams;
            std::vector<Weighted<View<PerspectiveCamera, Image5d>>> gcs;
            Image5d gc;
            static const int hcamNum = 16;
            static const Sizei hcamScreenSize(500, 500);
            //static const Sizei hcamScreenSize(500, 700);
            static const int hcamFocal = 200;
            std::string hcamsgcsFileName;
            {
                std::stringstream ss;
                ss << "hcamsgcs_" << hcamNum << "_" << hcamScreenSize.width << "_" << hcamScreenSize.height << "_" << hcamFocal;
                hcamsgcsFileName = ss.str();
            }
            if (0 || !misc::LoadCache(impath, hcamsgcsFileName, hcams, gcs)) {
                // extract gcs
                hcams = CreateHorizontalPerspectiveCameras(view.camera, hcamNum, hcamScreenSize.width, hcamScreenSize.height, hcamFocal);
                gcs.resize(hcams.size());
                for (int i = 0; i < hcams.size(); i++) {
                    auto pim = view.sampled(hcams[i]);
                    auto pgc = ComputeGeometricContext(matlab, pim.image, false, true);
                    gcs[i].component.camera = hcams[i];
                    gcs[i].component.image = pgc;
                    gcs[i].score = abs(1.0 - normalize(hcams[i].forward()).dot(normalize(view.camera.up())));
                }
                misc::SaveCache(impath, hcamsgcsFileName, hcams, gcs);
            }
            std::string gcmergedFileName;
            {
                std::stringstream ss;
                ss << "gc_" << hcamNum << "_" << hcamScreenSize.width << "_" << hcamScreenSize.height << "_" << hcamFocal;
                gcmergedFileName = ss.str();
            }
            if (0 || !misc::LoadCache(impath, gcmergedFileName, gc)) {
                gc = Combine(view.camera, gcs).image;
                misc::SaveCache(impath, gcmergedFileName, gc);
            }


            if (0) {
                std::vector<Imaged> gcChannels;
                cv::split(gc, gcChannels);
                gui::AsCanvas(ConvertToImage3d(gc)).show(1, "gc");
            }



            // build pigraph!
            bool refresh_mg_init = refresh_preparation || false;
            if (refresh_mg_init || !misc::LoadCache(impath, "mg_init", mg)) {
                std::cout << "########## refreshing mg init ###########" << std::endl;
                mg = BuildPIGraph(view, vps, vertVPId, segs, line3s,
                    DegreesToRadians(1), DegreesToRadians(1), DegreesToRadians(1),
                    0.04, DegreesToRadians(15), DegreesToRadians(2));
                misc::SaveCache(impath, "mg_init", mg);
            }

            bool refresh_line2leftRightSegs = refresh_mg_init || false;
            std::vector<std::array<std::set<int>, 2>> line2leftRightSegs;
            static const double angleDistForSegLineNeighborhood = DegreesToRadians(5);
            if (refresh_line2leftRightSegs || !misc::LoadCache(impath, "line2leftRightSegs", line2leftRightSegs)) {
                std::cout << "########## refreshing line2leftRightSegs ###########" << std::endl;
                line2leftRightSegs = CollectSegsNearLines(mg, angleDistForSegLineNeighborhood);
                misc::SaveCache(impath, "line2leftRightSegs", line2leftRightSegs);
            }


            const auto printLines = [&mg, &impath](int delay, const std::string & saveAs) {
                static gui::ColorTable rgbColors = gui::RGBGreys;
                rgbColors.exceptionalColor() = gui::Gray;
                Image3ub image = mg.view.image.clone();
                auto canvas = gui::AsCanvas(image);
                for (auto & l : mg.lines) {
                    static const double sampleAngle = M_PI / 100.0;
                    auto & line = l.component;
                    int claz = l.claz;
                    if (claz >= mg.vps.size()) {
                        claz = -1;
                    }
                    double spanAngle = AngleBetweenDirections(line.first, line.second);
                    std::vector<Point2> ps; ps.reserve(spanAngle / sampleAngle);
                    for (double angle = 0.0; angle <= spanAngle; angle += sampleAngle) {
                        Vec3 dir = RotateDirection(line.first, line.second, angle);
                        ps.push_back(mg.view.camera.toScreen(dir));
                    }
                    for (int i = 1; i < ps.size(); i++) {
                        auto p1 = (ps[i - 1]);
                        auto p2 = (ps[i]);
                        if (Distance(p1, p2) >= mg.view.image.cols / 2) {
                            continue;
                        }
                        canvas.thickness(3);
                        canvas.colorTable(rgbColors).add(gui::ClassifyAs(Line2(p1, p2), claz));
                    }
                }
                canvas.show(delay, "lines");
                if (saveAs != "") {
                    cv::imwrite(misc::FolderOfFile(impath) + "/" + saveAs + ".png", canvas.image());
                }
            };

            const auto printOrientedSegs = [&mg, &impath](int delay, const std::string & saveAs) {
                static const gui::ColorTable randColors = gui::CreateRandomColorTableWithSize(mg.nsegs);
                static gui::ColorTable rgbColors = gui::RGBGreys;
                rgbColors.exceptionalColor() = gui::Gray;
                auto pim = Print2(mg,
                    [&mg](int seg, Pixel pos) -> gui::Color {
                    static const gui::ColorTable ctable = gui::ColorTableDescriptor::RGBGreys;
                    auto & c = mg.seg2control[seg];
                    if (!c.used) {
                        return gui::Black;
                    }
                    if (c.orientationClaz != -1) {
                        return ctable[c.orientationClaz].blendWith(gui::White, 0.3);
                    }
                    if (c.orientationNotClaz != -1) {
                        static const int w = 10;
                        if (IsBetween((pos.x + pos.y) % w, 0, w / 2 - 1)) {
                            return ctable[c.orientationNotClaz].blendWith(gui::White, 0.3);
                        } else {
                            return gui::White;
                        }
                    }
                    return gui::White;
                },
                    [&mg](int lp) {
                    return gui::Transparent;
                },
                    [&mg](int bp) -> gui::Color {
                    return gui::Gray;
                }, 1, 0);
                auto canvas = gui::AsCanvas(pim);
                canvas.show(delay, "oriented segs");
                if (saveAs != "") {
                    cv::imwrite(misc::FolderOfFile(impath) + "/" + saveAs + ".png", Image3ub(canvas.image() * 255));
                }
                return canvas.image();
            };
            
            const auto printPIGraph = [&mg, &impath](int delay, const std::string & saveAs) {
                static const gui::ColorTable randColors = gui::CreateRandomColorTableWithSize(mg.nsegs);
                static gui::ColorTable rgbColors = gui::RGBGreys;
                rgbColors.exceptionalColor() = gui::Gray;
                auto pim = Print2(mg,
                    [&mg](int seg, Pixel pos) -> gui::Color {
                    static const gui::ColorTable ctable = gui::ColorTableDescriptor::RGBGreys;
                    auto & c = mg.seg2control[seg];
                    if (!c.used) {
                        return gui::Black;
                    }
                    if (c.orientationClaz != -1) {
                        return ctable[c.orientationClaz].blendWith(gui::White, 0.3);
                    }
                    if (c.orientationNotClaz != -1) {
                        static const int w = 10;
                        if (IsBetween((pos.x + pos.y) % w, 0, w/2-1)) {
                            return ctable[c.orientationNotClaz].blendWith(gui::White, 0.3);
                        } else {
                            return gui::White;
                        }
                    }
                    return gui::White;
                },
                    [&mg](int lp) {
                    return gui::Transparent;
                },
                    [&mg](int bp) -> gui::Color {
                    return gui::Gray;
                }, 1, 0);
                auto canvas = gui::AsCanvas(pim);
                for (auto & l : mg.lines) {
                    static const double sampleAngle = M_PI / 100.0;
                    auto & line = l.component;
                    int claz = l.claz;
                    if (claz >= mg.vps.size()) {
                        claz = -1;
                    }
                    double spanAngle = AngleBetweenDirections(line.first, line.second);
                    std::vector<Point2> ps; ps.reserve(spanAngle / sampleAngle);
                    for (double angle = 0.0; angle <= spanAngle; angle += sampleAngle) {
                        Vec3 dir = RotateDirection(line.first, line.second, angle);
                        ps.push_back(mg.view.camera.toScreen(dir));
                    }
                    for (int i = 1; i < ps.size(); i++) {
                        auto p1 = ToPixel(ps[i - 1]);
                        auto p2 = ToPixel(ps[i]);
                        if (Distance(p1, p2) >= mg.view.image.cols / 2) {
                            continue;
                        }
                        gui::Color color = rgbColors[claz];
                        cv::clipLine(cv::Rect(0, 0, canvas.image().cols, canvas.image().rows), p1, p2);
                        cv::line(canvas.image(), p1, p2, (cv::Scalar)color / 255.0, 2);
                    }
                }
                canvas.show(delay, "pi graph");
                if (saveAs != "") {
                    cv::imwrite(misc::FolderOfFile(impath) + "/" + saveAs + ".png", Image3ub(canvas.image() * 255));
                }
                return canvas.image();
            };

            auto linesSegsIm = printPIGraph(0, "initial_lines_segs");

            if(false){
                Image3ub lsim = linesSegsIm * 255;
                ReverseRows(lsim);
                gui::SceneBuilder sb;
                gui::ResourceStore::set("tex", lsim);
                Sphere3 sphere;
                sphere.center = Origin();
                sphere.radius = 1.0;
                sb.begin(sphere).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("tex").end();
                sb.show(true, true);
            }



            // attach orientation constraints
            bool refresh_mg_oriented = refresh_mg_init || false;
            if (refresh_mg_oriented || !misc::LoadCache(impath, "mg_oriented", mg)) {
                std::cout << "########## refreshing mg oriented ###########" << std::endl;
                AttachPrincipleDirectionConstraints(mg);
                AttachWallConstraints(mg, M_PI / 60.0);
                AttachGCConstraints(mg, gc, 0.7, 0.7, true);
                misc::SaveCache(impath, "mg_oriented", mg);
            }

            //printLines(0, "oriented_lines");
            //printOrientedSegs(0, "oriented_segs");
            printPIGraph(0, "oriented_lines_segs");

            // detect occlusions
            bool refresh_lsw = refresh_mg_oriented || true;
            std::vector<LineSidingWeight> lsw;
            if (refresh_lsw || !misc::LoadCache(impath, "lsw", lsw)) {
                std::cout << "########## refreshing lsw ###########" << std::endl;
                lsw = ComputeLinesSidingWeights(mg, DegreesToRadians(3), 0.2, 0.1, angleDistForSegLineNeighborhood);
                misc::SaveCache(impath, "lsw", lsw);
            }



            if (true) {
                auto pim = Print(mg, ConstantFunctor<gui::Color>(gui::White), 
                    ConstantFunctor<gui::Color>(gui::Transparent), 
                    ConstantFunctor<gui::Color>(gui::Gray), 2, 0);
                auto drawLine = [&mg, &pim](const Line3 & line, const std::string & text, const gui::Color & color, bool withTeeth) {
                    double angle = AngleBetweenDirections(line.first, line.second);
                    std::vector<Pixel> ps;
                    for (double a = 0.0; a <= angle; a += 0.01) {
                        ps.push_back(ToPixel(mg.view.camera.toScreen(RotateDirection(line.first, line.second, a))));
                    }
                    for (int i = 1; i < ps.size(); i++) {
                        auto p1 = ps[i - 1];
                        auto p2 = ps[i];
                        if (Distance(p1, p2) >= pim.cols / 2) {
                            continue;
                        }
                        cv::clipLine(cv::Rect(0, 0, pim.cols, pim.rows), p1, p2);
                        cv::line(pim, p1, p2, (cv::Scalar)color / 255.0, 2);
                        if (withTeeth && i % 3 == 0) {
                            auto teethp = ToPixel(RightPerpendicularDirectiion(ecast<double>(p2 - p1))) * 2 + p1;
                            auto tp1 = p1, tp2 = teethp;
                            cv::clipLine(cv::Rect(0, 0, pim.cols, pim.rows), tp1, tp2);
                            cv::line(pim, tp1, tp2, (cv::Scalar)color / 255.0, 1);
                        }
                    }
                    //cv::circle(pim, ps.back(), 2.0, (cv::Scalar)color / 255.0, 2);
                    if (!text.empty()) {
                        cv::putText(pim, text, ps.back() + Pixel(5, 0), 1, 0.7, color);
                    }
                };
                for (int line = 0; line < mg.nlines(); line++) {
                    auto & ws = lsw[line];
                    auto l = mg.lines[line].component;
                    if (!ws.isOcclusion()) {
                        drawLine(l, "", gui::Black, false);
                    } else if (ws.onlyConnectLeft()) {
                        drawLine(l, std::to_string(line), gui::Blue, true);
                    } else if (ws.onlyConnectRight()) {
                        drawLine(l.reversed(), std::to_string(line), gui::Blue, true);
                    } else {
                        drawLine(l, std::to_string(line), gui::Red, false);
                    }
                }
                gui::AsCanvas(pim).show(0, "line labels");
            }


            bool refresh_mg_occdetected = refresh_lsw || refresh_line2leftRightSegs || false;
            if (refresh_mg_occdetected || !misc::LoadCache(impath, "mg_occdetected", mg)) {
                std::cout << "########## refreshing mg occdetected ###########" << std::endl;
                ApplyLinesSidingWeights(mg, lsw, line2leftRightSegs, true);
                if (anno.extendedOnTop) {
                    DisableTopSeg(mg);
                }
                //if (anno.extendedOnBottom) {
                //    DisableBottomSeg(mg);
                //}
                misc::SaveCache(impath, "mg_occdetected", mg);
            }


            PIConstraintGraph cg;
            PICGDeterminablePart dp;
            bool refresh_mg_reconstructed = refresh_mg_occdetected || true;
            if (refresh_mg_reconstructed || !misc::LoadCache(impath, "mg_reconstructed", mg, cg, dp)) {
                std::cout << "########## refreshing mg reconstructed ###########" << std::endl;
                cg = BuildPIConstraintGraph(mg, /*lsw, line2leftRightSegs, */DegreesToRadians(1));
                
                for(int i = 0; i < 2; i++) {
                    dp = LocateDeterminablePart(cg, DegreesToRadians(3));
                    double energy = Solve(dp, cg, matlab, 5);
                    if (IsInfOrNaN(energy)) {
                        std::cout << "solve failed" << std::endl;
                        return nullptr;
                    }
                    auto depthMap = DepthMap(dp, cg, mg);
                    gui::AsCanvas(depthMap).show();

                    VisualizeReconstruction(dp, cg, mg, true);

                   /* int disabledNum = DisableUnsatisfiedConstraints(dp, cg, [](double distRankRatio, double avgDist, double maxDist) {
                        return maxDist > 0.05;
                    });
                    std::cout << "disabled constraint num = " << disabledNum << std::endl;
                    if (disabledNum == 0) {
                        break;
                    }*/

                    DisorientDanglingLines2(dp, cg, mg, 0.5);
                    DisorientDanglingSegs(dp, cg, mg, 0.1);
                }

                misc::SaveCache(impath, "mg_reconstructed", mg, cg);
            }
         
            VisualizeReconstruction(dp, cg, mg, false);


            // save to disk
            SaveToDisk(resultPath, mg);
        }

        return std::make_unique<PIGraphModel>(std::move(mg));
    }

}
