#include <QApplication>
#include <QFileDialog>
#include <QtCore>

#include "../../src/core/parallel.hpp"
#include "../../src/misc/cache.hpp"

#include "../../src/gui/singleton.hpp"
#include "../../src/gui/utility.hpp"
#include "../../src/gui/canvas.hpp"

#include "../../src/experimental/pi_graph_control.hpp"
#include "../../src/experimental/pi_graph_occlusion.hpp"
#include "../../src/experimental/pi_graph_cg.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"


using namespace pano;
using namespace pano::core;
using namespace pano::experimental;



std::vector<Image3d> GT(const PILayoutAnnotation & anno, const std::vector<PerspectiveCamera> & testCams) {
    std::vector<Polygon3> polygons(anno.nfaces());
    for (int i = 0; i < anno.nfaces(); i++) {
        auto & plane = anno.face2plane[i];
        auto & poly = polygons[i];
        poly.normal = plane.normal;
        for (int c : anno.face2corners[i]) {
            Ray3 ray(Origin(), anno.corners[c]);
            poly.corners.push_back(IntersectionOfLineAndPlane(ray, plane).position);
        }
    }

    std::vector<Image3d> gts(testCams.size());
    ParallelRun(testCams.size(), std::thread::hardware_concurrency(), [&](int i) {
        std::cout << "computing gt surface normals - " << i << std::endl;
        auto & cam = testCams[i];
        Image3d gt(cam.screenSize());
        for (auto it = gt.begin(); it != gt.end(); ++it) {
            Vec3 direction = normalize(cam.toSpace(it.pos()));
            Ray3 ray(Origin(), direction);
            double depth = std::numeric_limits<double>::infinity();
            for (int i = 0; i < anno.nfaces(); i++) {
                auto & poly = polygons[i];
                auto inter = IntersectionOfLineAndPolygon(ray, poly);
                if (inter.failed()) {
                    continue;
                }
                *it = normalize(poly.normal);
            }
        }
        gts[i] = gt;

        for (auto it = gt.begin(); it != gt.end(); ++it) {
            if (*it == Origin()) {
                for (int dx = -1; dx <= 1; dx++) {
                    if (*it != Origin()) {
                        break;
                    }
                    for (int dy = -1; dy <= 1; dy++) {
                        auto p = it.pos() + Pixel(dx, dy);
                        if (Contains(gt, p) && gt(p) != Origin()) {
                            *it = gt(p);
                            break;
                        }
                    }
                }
            }
        }
    });

    return gts;
}




// compute normal maps for comparison
std::vector<Imagei> OrientationMaps(const PILayoutAnnotation & anno, const std::vector<PerspectiveCamera> & testCams) {
    std::vector<Imagei> oms(testCams.size());
    ParallelRun(testCams.size(), std::thread::hardware_concurrency(), [&](int i) {
        std::vector<HPoint2> hvps(anno.vps.size());
        for (int k = 0; k < anno.vps.size(); k++) {
            hvps[k] = testCams[i].toScreenInHPoint(anno.vps[k]);
        }
        // collect line 2ds in this view
        auto pim = anno.view.sampled(testCams[i]).image;
        LineSegmentExtractor lineExtractor;
        lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
        auto ls = lineExtractor(pim); // use pyramid
        
        auto rawLines2 = ClassifyEachAs(ls, -1);
        ClassifyLines(rawLines2, hvps);

        Imagei omap = ComputeOrientationMaps(rawLines2, hvps, testCams[i].screenSize());
        oms[i] = omap;
    });
    return oms;
}

std::vector<Image7d> GeometricContext(const PILayoutAnnotation & anno, const std::vector<PerspectiveCamera> & testCams, misc::Matlab & matlab) {
    auto up = normalize(anno.vps[anno.vertVPId]);
    if (up.dot(-anno.view.camera.up()) < 0) {
        up = -up;
    }
    std::vector<Image7d> rawgcs(testCams.size());
    for (int i = 0; i < testCams.size(); i++) {
        auto pim = anno.view.sampled(testCams[i]).image;
        rawgcs[i] = ComputeRawGeometricContext(matlab, pim, false, true);
    }
    return rawgcs;
}


// proposed method
std::vector<Image3d> Panoramix(const PILayoutAnnotation & anno, const std::vector<PerspectiveCamera> & testCams, misc::Matlab & matlab, 
    bool useGTOcclusions) {
    
    const auto & impath = anno.impath;

    auto image = anno.rectifiedImage.clone();
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
    PIGraph mg;
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
        line2leftRightSegs = CollectSegsNearLines(mg, angleDistForSegLineNeighborhood * 2);
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


    if (false) {
        Image3ub lsim = printPIGraph(0, "initial_lines_segs") * 255;
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
    //printPIGraph(0, "oriented_lines_segs");

    // detect occlusions
    bool refresh_lsw = refresh_mg_oriented || true;
    std::vector<LineSidingWeight> lsw;
    //std::vector<std::map<int, double>> line2leftSegsWithWeight, line2rightSegsWithWeight;
    if (!useGTOcclusions) {
        if (refresh_lsw || !misc::LoadCache(impath, "lsw", lsw)) {
            std::cout << "########## refreshing lsw ###########" << std::endl;
            lsw = ComputeLinesSidingWeights(mg, DegreesToRadians(3), 0.2, 0.1,
                angleDistForSegLineNeighborhood);
            misc::SaveCache(impath, "lsw", lsw);
        }
    } else {
        if (refresh_lsw || !misc::LoadCache(impath, "lsw_from_gt", lsw)) {
            std::cout << "########## refreshing lsw from gt ###########" << std::endl;
            lsw = ComputeLinesSidingWeightsFromAnnotation(mg, anno, DegreesToRadians(0.5), DegreesToRadians(8), 0.6);
            misc::SaveCache(impath, "lsw_from_gt", lsw);
        }
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


    bool refresh_mg_occdetected = refresh_lsw || refresh_line2leftRightSegs || true;
    if (refresh_mg_occdetected || !misc::LoadCache(impath, "mg_occdetected", mg)) {
        std::cout << "########## refreshing mg occdetected ###########" << std::endl;
        ApplyLinesSidingWeights(mg, lsw, line2leftRightSegs, true);
        if (anno.extendedOnTop && !anno.topIsPlane) {
            DisableTopSeg(mg);
        }
        if (anno.extendedOnBottom && !anno.bottomIsPlane) {
            DisableBottomSeg(mg);
        }
        misc::SaveCache(impath, "mg_occdetected", mg);
    }


    PIConstraintGraph cg;
    PICGDeterminablePart dp;
    bool refresh_mg_reconstructed = refresh_mg_occdetected || true;
    if (refresh_mg_reconstructed || !misc::LoadCache(impath, "mg_reconstructed", mg, cg, dp)) {
        std::cout << "########## refreshing mg reconstructed ###########" << std::endl;
        cg = BuildPIConstraintGraph(mg, DegreesToRadians(1), 0.2);

        for (int i = 0; i < 2; i++) {
            dp = LocateDeterminablePart(cg, DegreesToRadians(3));
            double energy = Solve(dp, cg, matlab, 5, 1e6);
            if (IsInfOrNaN(energy)) {
                std::cout << "solve failed" << std::endl;
                return {};
            }
            
            std::pair<double, double> validRange;
            auto depthMap = DepthMap(dp, cg, mg, &validRange);
            Imagei depthMapDisc = (depthMap - validRange.first) / (validRange.second - validRange.first) * 255;
            auto jetctable = gui::CreateJetColorTableWithSize(MinMaxValOfImage(depthMapDisc).second + 1);
            gui::AsCanvas(jetctable(depthMapDisc)).show(1, "depth map");

            auto surfaceNormalMap = SurfaceNormalMap(dp, cg, mg);
            auto surfaceNormalMapForShow = surfaceNormalMap.clone();
            for (auto & n : surfaceNormalMapForShow) {
                auto nn = n;
                for (int i = 0; i < mg.vps.size(); i++) {
                    n[i] = abs(nn.dot(vps[i]));
                }
            }
            gui::AsCanvas(surfaceNormalMapForShow).show(0, "surface normal map");

            VisualizeReconstruction(dp, cg, mg, false);

            DisorientDanglingLines3(dp, cg, mg, 0.05, 0.1);
            DisorientDanglingSegs3(dp, cg, mg, 0.01, 0.1);
            OverorientSkewSegs(dp, cg, mg, DegreesToRadians(3), DegreesToRadians(45), 0.5);
        }

        misc::SaveCache(impath, "mg_reconstructed", mg, cg, dp);
    }  


    std::vector<Image3d> surfaceNormalMaps(testCams.size());
    ParallelRun(testCams.size(), std::thread::hardware_concurrency(), [&](int i){
        std::cout << "computing surface normal map for panoramix on testCamera - " << i << std::endl;
        auto & map = surfaceNormalMaps[i];
        auto & cam = testCams[i];
        map = Image3d(cam.screenSize());
        for (auto it = map.begin(); it != map.end(); ++it) {
            auto p = it.pos();
            auto dir = normalize(cam.toSpace(p));
            auto pp = ToPixel(anno.view.camera.toScreen(dir));
            pp.x = WrapBetween(pp.x, 0, mg.segs.cols);
            if (!Contains(mg.segs, pp)) {
                continue;
            }
            int seg = mg.segs(pp);
            int ent = cg.seg2ent[seg];
            if (ent == -1) {
                continue;
            }
            if (!Contains(dp.determinableEnts, ent)) {
                continue;
            }
            auto & plane = cg.entities[ent].supportingPlane.reconstructed;
            auto n = normalize(plane.normal);
            *it = n;
        }
    });


    return surfaceNormalMaps;

}










// evaluate helpers
double ErrorOfSurfaceNormal(const Vec3 & gt, const Vec3 & cand) {
    if (norm(cand) == 0) {
        return M_PI_2;
    }
    return AngleBetweenUndirectedVectors(gt, cand);
}


int GeometricContextToLabel(const Vec7 & label) {
    // 0: front, 1: left, 2: right, 3: floor, 4: ceiling, 5: clutter, 6: unknown
    int labels[5];
    std::iota(labels, labels + 5, 0);
    std::sort(labels, labels + 5, [&label](int a, int b) {return label[a] > label[b]; });
    int maxid = labels[0];
    if (label[maxid] == 0.0) {
        return -1;
    }
    return maxid;
}

int SurfaceNormalToLabel(const Vec3 & normal, const Vec3 & up, const Vec3 & forward, const Vec3 & left) {
    if (norm(normal) == 0) {
        return -1;
    }
    // 0: front, 1: left, 2: right, 3: floor, 4: ceiling, 5: clutter, 6: unknown
    Vec3 dirs[] = { forward, left, -left, -up, up };
    double angleToDirs[5];
    for (int i = 0; i < 5; i++) {
        angleToDirs[i] = AngleBetweenDirections(dirs[i], normal);
    }
    return std::min_element(std::begin(angleToDirs), std::end(angleToDirs)) - std::begin(angleToDirs);
}



int main(int argc, char ** argv) {
    
    gui::Singleton::InitGui();
    misc::Matlab matlab("", true);
    
    std::vector<std::string> impaths;
    gui::PickImages("H:\\DataSet\\pi\\dataset\\selected\\", &impaths);    
    
    bool recache_gc = false;
    bool recache_om = false;
    bool recache_pn = false;
    bool recache_gt = false;

    bool pn_vs_om = false;
    bool pn_vs_gc = true;

    bool consider_horizontal_cams_only = true;

    // cache mode
    // cache gt
    for (auto & impath : impaths) {
        std::vector<PerspectiveCamera> testCams1;
        std::vector<Image3d> gt1;
        if (recache_gt || !misc::LoadCache(impath, "testCams1_gt1", testCams1, gt1)) {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            if (anno.nfaces() == 0) {
                while (true) {
                    pano::experimental::EditLayoutAnnotation(impath, anno);
                    pano::experimental::ReconstructLayoutAnnotation(anno, matlab);
                    pano::experimental::VisualizeLayoutAnnotation(anno);
                    int selected = pano::gui::SelectFrom({ "Accept", "Edit Again", "Abandon" },
                        "Your decision?",
                        "Accept the edit, or edit it again, or just abandon the edit this time?", 0, 2);
                    if (selected == 0) {
                        pano::experimental::SaveLayoutAnnotation(impath, anno);
                        break;
                    } else if (selected == 2) {
                        break;
                    }
                }
            }
            testCams1 = CreatePanoContextCameras(anno.view.camera);
            gt1 = GT(anno, testCams1);
            misc::SaveCache(impath, "testCams1_gt1", testCams1, gt1);
        }
    }

    // cache gc
    for (auto & impath : impaths) {
        std::vector<PerspectiveCamera> testCams1;
        std::vector<Image7d> gc1;
        if (recache_gc || !misc::LoadCache(impath, "testCams1_gc1", testCams1, gc1)) {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            testCams1 = CreatePanoContextCameras(anno.view.camera);
            gc1 = GeometricContext(anno, testCams1, matlab);
            misc::SaveCache(impath, "testCams1_gc1", testCams1, gc1);
        }
    }

    // cache om
    for (auto & impath : impaths) {
        std::vector<PerspectiveCamera> testCams1;
        std::vector<Imagei> om1;
        if (recache_om || !misc::LoadCache(impath, "testCams1_om1", testCams1, om1)) {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            testCams1 = CreatePanoContextCameras(anno.view.camera);
            om1 = OrientationMaps(anno, testCams1);
            misc::SaveCache(impath, "testCams1_om1", testCams1, om1);
        }
    }

    // cache pn
    for (auto & impath : impaths) {
        std::vector<PerspectiveCamera> testCams1;
        std::vector<Image3d> pn1;
        if (recache_pn || !misc::LoadCache(impath, "testCams1_pn1", testCams1, pn1)) {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            testCams1 = CreatePanoContextCameras(anno.view.camera);
            pn1 = Panoramix(anno, testCams1, matlab, false);
            misc::SaveCache(impath, "testCams1_pn1", testCams1, pn1);
        }
    }

    // cache pn_gtocc
    for (auto & impath : impaths) {
        std::vector<PerspectiveCamera> testCams1;
        std::vector<Image3d> pngtocc1;
        if (recache_pn || !misc::LoadCache(impath, "testCams1_pngtocc1", testCams1, pngtocc1)) {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            testCams1 = CreatePanoContextCameras(anno.view.camera);
            pngtocc1 = Panoramix(anno, testCams1, matlab, true);
            misc::SaveCache(impath, "testCams1_pngtocc1", testCams1, pngtocc1);
        }
    }


    // evaluate mode
    if (pn_vs_om) {
        for (auto & impath : impaths) {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            std::vector<PerspectiveCamera> testCams1 = CreatePanoContextCameras(anno.view.camera);
            
            std::vector<PerspectiveCamera> testCamsA, testCamsB, testCamsC;
            std::vector<Image3d> pn1;
            misc::LoadCache(impath, "testCams1_pn1", testCamsA, pn1);
            std::vector<Imagei> om1;
            misc::LoadCache(impath, "testCams1_om1", testCamsB, om1);
            std::vector<Image3d> gt1;
            misc::LoadCache(impath, "testCams1_gt1", testCamsC, gt1);

            assert(testCams1 == testCamsA && testCams1 == testCamsB && testCams1 == testCamsC);

            double error_pn = 0.0;
            double error_om = 0.0;

            int npixels = 0;
            for (int i = 0; i < testCams1.size(); i++) {
                auto & pn = pn1[i];
                auto & om = om1[i];
                auto & gt = gt1[i];

                for (auto it = gt.begin(); it != gt.end(); ++it) {
                    Vec3 gtNormal = *it;
                    if (norm(gtNormal) == 0) {
                        continue;
                    }
                    auto p = it.pos();

                    Vec3 pnNormal = pn(p);
                    //if (norm(pnNormal) == 0) {
                    //    pnNormal = gtNormal;
                    //}
                    if (om(p) == -1 || norm(pnNormal) == 0) {
                        continue;
                    }
                    Vec3 omNormal = om(p) == -1 ? Origin() : anno.vps[om(p)];

                    error_pn += ErrorOfSurfaceNormal(gtNormal, pnNormal);
                    error_om += ErrorOfSurfaceNormal(gtNormal, omNormal);

                    npixels++;
                }
            }

            error_pn /= npixels;
            error_om /= npixels;

            std::cout << "{{{PN[" << error_pn << "] vs OM[" << error_om << "]}}} on image\"" << impath << "\"" << std::endl;
        }      
    }


    if (pn_vs_gc) {

        auto ctable = gui::CreateGreyColorTableWithSize(4);
        std::map<std::string, double> impath2errorpn, impath2errorgc;

        for (auto & impath : impaths) {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            std::vector<PerspectiveCamera> testCams1 = CreatePanoContextCameras(anno.view.camera);

            std::vector<PerspectiveCamera> testCamsA, testCamsB, testCamsC;
            std::vector<Image3d> pn1;
            misc::LoadCache(impath, "testCams1_pn1", testCamsA, pn1);
            std::vector<Image7d> gc1;
            misc::LoadCache(impath, "testCams1_gc1", testCamsB, gc1);
            std::vector<Image3d> gt1;
            misc::LoadCache(impath, "testCams1_gt1", testCamsC, gt1);

            assert(testCams1 == testCamsA && testCams1 == testCamsB && testCams1 == testCamsC);

            double error_pn = 0.0;
            double error_gc = 0.0;

            int npixels = 0;
            std::vector<Imagei> gcLabelsTable(testCams1.size()), pnLabelsTable(testCams1.size()), gtLabelsTable(testCams1.size());
            for (int i = 0; i < testCams1.size(); i++) {

                auto & pn = pn1[i];
                auto & gc = gc1[i];
                auto & gt = gt1[i];

                auto & cam = testCams1[i];
                if (consider_horizontal_cams_only &&
                    AngleBetweenUndirectedVectors(anno.view.camera.up(), cam.forward()) < DegreesToRadians(60)) {
                    continue;
                }

                Imagei gcLabels(cam.screenSize(), -1);
                Imagei pnLabels(cam.screenSize(), -1);
                Imagei gtLabels(cam.screenSize(), -1);

                for (auto it = gt.begin(); it != gt.end(); ++it) {
                    auto p = it.pos();
                    Vec3 gtNormal = *it;
                    if (norm(gtNormal) == 0) {
                        continue;
                    }

                    int gtLabel = SurfaceNormalToLabel(gtNormal, cam.up(), cam.forward(), cam.leftward());
                    gtLabels(p) = gtLabel;

                    Vec3 pnNormal = pn(p); 
                    int pnLabel = SurfaceNormalToLabel(pnNormal, cam.up(), cam.forward(), cam.leftward());
                    Vec7 gcScores = gc(p);
                    int gcLabel = GeometricContextToLabel(gcScores);

                    pnLabels(p) = pnLabel;
                    gcLabels(p) = gcLabel;

                    error_pn += pnLabel != gtLabel;
                    error_gc += gcLabel != gtLabel;

                    npixels++;
                }

                gcLabelsTable[i] = gcLabels;
                pnLabelsTable[i] = pnLabels;
                gtLabelsTable[i] = gtLabels;

                gui::MakeCanvas(ctable(gcLabels)).show(1, "gcLabels");
                gui::MakeCanvas(ctable(pnLabels)).show(1, "pnLabels");
                gui::MakeCanvas(ctable(gtLabels)).show(0, "gtLabels");
            }

            error_pn /= npixels;
            error_gc /= npixels;

            std::cout << std::endl << std::endl;
            std::cout << "{{{{{{{{{  PN[" << error_pn << "] vs GC[" << error_gc << "]  }}}}}}}}}}}} on image\"" << impath << "\"" << std::endl;
            std::cout << std::endl << std::endl;

            impath2errorpn[impath] = error_pn;
            impath2errorgc[impath] = error_gc;
        }

        double mean_error_pn = 0.0;
        double mean_error_gc = 0.0;
        for (auto & impath : impaths) {
            double error_pn = impath2errorpn[impath];
            mean_error_pn += error_pn;
            double error_gc = impath2errorgc[impath];
            mean_error_gc += error_gc;
            std::cout << "{{{{{{{{{  PN[" << error_pn << "] vs GC[" << error_gc << "]  }}}}}}}}}}}} on image\"" << impath << "\"" << std::endl;
        }

        mean_error_pn /= impaths.size();
        mean_error_gc /= impaths.size();

        std::cout << "{{{{{{{{{  PN[" << mean_error_pn << "] vs GC[" << mean_error_gc << "]  }}}}}}}}}}}} on average" << std::endl;

    }


    return 0;

}