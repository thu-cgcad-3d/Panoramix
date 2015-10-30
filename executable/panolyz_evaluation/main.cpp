
#include "panoramix.hpp"


std::vector<Imagei> GTFaceLabels(const PILayoutAnnotation & anno, const std::vector<PerspectiveCamera> & testCams) {
    auto impath = anno.impath;
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

    std::vector<Imagei> faceLabels(testCams.size());
    ParallelRun(testCams.size(), std::thread::hardware_concurrency(), [&](int i) {
        std::cout << "computing gt face label map - " << i << std::endl;
        auto & cam = testCams[i];
        Imagei faceLabelMap(cam.screenSize(), -1);
        for (auto it = faceLabelMap.begin(); it != faceLabelMap.end(); ++it) {
            Vec3 direction = normalize(cam.toSpace(it.pos()));
            Ray3 ray(Origin(), direction);
            double depth = std::numeric_limits<double>::infinity();
            for (int i = 0; i < anno.nfaces(); i++) {
                auto & poly = polygons[i];
                auto inter = IntersectionOfLineAndPolygon(ray, poly);
                if (!inter.failed()) {
                    *it = i;
                    break;
                }
            }
        }
        for (auto it = faceLabelMap.begin(); it != faceLabelMap.end(); ++it) {
            if (*it == -1) {
                for (int dx = -1; dx <= 1; dx++) {
                    if (*it != -1) {
                        break;
                    }
                    for (int dy = -1; dy <= 1; dy++) {
                        auto p = it.pos() + Pixel(dx, dy);
                        if (Contains(faceLabelMap, p) && faceLabelMap(p) != -1) {
                            *it = faceLabelMap(p);
                            break;
                        }
                    }
                }
            }
        }
        faceLabels[i] = faceLabelMap;
    });

    return faceLabels;
}


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










// evaluation helpers
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
    misc::Matlab matlab("", false);

    using TaskQueue = std::vector<std::function<Any(const std::string & impath)>>;
    TaskQueue activeQ;

    {
        TaskQueue Q1;

        // set tasks
        // default
        Q1.push_back([&matlab](const std::string & impath) -> Any {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            anno.impath = impath;

            PanoramixOptions options;
            options.useWallPrior = true;
            options.usePrincipleDirectionPrior = true;
            options.useGeometricContextPrior = true;

            options.useGTOcclusions = false;
            options.looseLinesSecondTime = false;
            options.looseSegsSecondTime = false;
            options.restrictSegsSecondTime = false;

            options.refresh_preparation = false;
            options.refresh_mg_init = options.refresh_preparation || false;
            options.refresh_mg_oriented = options.refresh_mg_init || false;
            options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
            options.refresh_lsw = options.refresh_mg_oriented || false;
            options.refresh_mg_occdetected = options.refresh_lsw || options.refresh_line2leftRightSegs || false;
            options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

            RunPanoramix(anno, options, matlab, false);
            return options;
        });
        // use gt occ
        Q1.push_back([&matlab](const std::string & impath) -> Any {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            anno.impath = impath;

            PanoramixOptions options;
            options.useWallPrior = true;
            options.usePrincipleDirectionPrior = true;
            options.useGeometricContextPrior = true;

            options.useGTOcclusions = true;
            options.looseLinesSecondTime = false;
            options.looseSegsSecondTime = false;
            options.restrictSegsSecondTime = false;

            options.refresh_preparation = false;
            options.refresh_mg_init = options.refresh_preparation || false;
            options.refresh_mg_oriented = options.refresh_mg_init || false;
            options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
            options.refresh_lsw = options.refresh_mg_oriented || false;
            options.refresh_mg_occdetected = options.refresh_lsw || options.refresh_line2leftRightSegs || false;
            options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

            RunPanoramix(anno, options, matlab, false);
            return options;
        });
        // without wall
        Q1.push_back([&matlab](const std::string & impath) -> Any {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            anno.impath = impath;

            PanoramixOptions options;
            options.useWallPrior = false;
            options.usePrincipleDirectionPrior = true;
            options.useGeometricContextPrior = true;

            options.useGTOcclusions = false;
            options.looseLinesSecondTime = false;
            options.looseSegsSecondTime = false;
            options.restrictSegsSecondTime = false;

            options.refresh_preparation = false;
            options.refresh_mg_init = options.refresh_preparation || false;
            options.refresh_mg_oriented = options.refresh_mg_init || false;
            options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
            options.refresh_lsw = options.refresh_mg_oriented || false;
            options.refresh_mg_occdetected = options.refresh_lsw || options.refresh_line2leftRightSegs || false;
            options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

            RunPanoramix(anno, options, matlab, false);
            return options;
        });
        // without principle direction
        Q1.push_back([&matlab](const std::string & impath) -> Any {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            anno.impath = impath;

            PanoramixOptions options;
            options.useWallPrior = true;
            options.usePrincipleDirectionPrior = false;
            options.useGeometricContextPrior = true;

            options.useGTOcclusions = false;
            options.looseLinesSecondTime = false;
            options.looseSegsSecondTime = false;
            options.restrictSegsSecondTime = false;

            options.refresh_preparation = false;
            options.refresh_mg_init = options.refresh_preparation || false;
            options.refresh_mg_oriented = options.refresh_mg_init || false;
            options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
            options.refresh_lsw = options.refresh_mg_oriented || false;
            options.refresh_mg_occdetected = options.refresh_lsw || options.refresh_line2leftRightSegs || false;
            options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

            RunPanoramix(anno, options, matlab, false);
            return options;
        });
        // without gc!
        Q1.push_back([&matlab](const std::string & impath) -> Any {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            anno.impath = impath;

            PanoramixOptions options;
            options.useWallPrior = true;
            options.usePrincipleDirectionPrior = true;
            options.useGeometricContextPrior = false;

            options.useGTOcclusions = false;
            options.looseLinesSecondTime = false;
            options.looseSegsSecondTime = false;
            options.restrictSegsSecondTime = false;

            options.refresh_preparation = false;
            options.refresh_mg_init = options.refresh_preparation || false;
            options.refresh_mg_oriented = options.refresh_mg_init || false;
            options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
            options.refresh_lsw = options.refresh_mg_oriented || false;
            options.refresh_mg_occdetected = options.refresh_lsw || options.refresh_line2leftRightSegs || false;
            options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

            RunPanoramix(anno, options, matlab, false);
            return options;
        });
        // without gc + gt occ !
        Q1.push_back([&matlab](const std::string & impath) -> Any {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            anno.impath = impath;

            PanoramixOptions options;
            options.useWallPrior = true;
            options.usePrincipleDirectionPrior = true;
            options.useGeometricContextPrior = false;

            options.useGTOcclusions = true;
            options.looseLinesSecondTime = false;
            options.looseSegsSecondTime = false;
            options.restrictSegsSecondTime = false;

            options.refresh_preparation = false;
            options.refresh_mg_init = options.refresh_preparation || false;
            options.refresh_mg_oriented = options.refresh_mg_init || false;
            options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
            options.refresh_lsw = options.refresh_mg_oriented || false;
            options.refresh_mg_occdetected = options.refresh_lsw || options.refresh_line2leftRightSegs || false;
            options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

            RunPanoramix(anno, options, matlab, false);
            return options;
        });
        // only gc!
        Q1.push_back([&matlab](const std::string & impath) -> Any {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            anno.impath = impath;

            PanoramixOptions options;
            options.useWallPrior = false;
            options.usePrincipleDirectionPrior = false;
            options.useGeometricContextPrior = true;

            options.useGTOcclusions = false;
            options.looseLinesSecondTime = false;
            options.looseSegsSecondTime = false;
            options.restrictSegsSecondTime = false;

            options.refresh_preparation = false;
            options.refresh_mg_init = options.refresh_preparation || false;
            options.refresh_mg_oriented = options.refresh_mg_init || false;
            options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
            options.refresh_lsw = options.refresh_mg_oriented || false;
            options.refresh_mg_occdetected = options.refresh_lsw || options.refresh_line2leftRightSegs || false;
            options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

            RunPanoramix(anno, options, matlab, false);
            return options;
        });
        // only gc + gt occ!
        Q1.push_back([&matlab](const std::string & impath) -> Any {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            anno.impath = impath;

            PanoramixOptions options;
            options.useWallPrior = false;
            options.usePrincipleDirectionPrior = false;
            options.useGeometricContextPrior = true;

            options.useGTOcclusions = true;
            options.looseLinesSecondTime = false;
            options.looseSegsSecondTime = false;
            options.restrictSegsSecondTime = false;

            options.refresh_preparation = false;
            options.refresh_mg_init = options.refresh_preparation || false;
            options.refresh_mg_oriented = options.refresh_mg_init || false;
            options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
            options.refresh_lsw = options.refresh_mg_oriented || false;
            options.refresh_mg_occdetected = options.refresh_lsw || options.refresh_line2leftRightSegs || false;
            options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

            RunPanoramix(anno, options, matlab, false);
            return options;
        });






        TaskQueue Q2;
        Q2.push_back([&matlab](const std::string & impath) -> Any {
            auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
            anno.impath = impath;

            PanoramixOptions options;
            options.useWallPrior = true;
            options.usePrincipleDirectionPrior = true;
            options.useGeometricContextPrior = true;

            options.useGTOcclusions = false;
            options.looseLinesSecondTime = false;
            options.looseSegsSecondTime = false;
            options.restrictSegsSecondTime = false;

            options.refresh_preparation = false;
            options.refresh_mg_init = options.refresh_preparation || false;
            options.refresh_mg_oriented = options.refresh_mg_init || false;
            options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
            options.refresh_lsw = options.refresh_mg_oriented || false;
            options.refresh_mg_occdetected = options.refresh_lsw || options.refresh_line2leftRightSegs || false;
            options.refresh_mg_reconstructed = options.refresh_mg_occdetected || false;

            RunPanoramix(anno, options, matlab, true);

            return options;
        });



        activeQ = Q2;
    }


    std::vector<std::string> impaths;
    gui::PickImages("H:\\DataSet\\pi\\dataset\\selected\\", &impaths);
    for (int i = 0; i < activeQ.size(); i++) {
        auto & task = activeQ[i];
        std::cout << "[[[[[[[[[ TASK " << i << "]]]]]]]]" << std::endl;
        std::vector<double> results(impaths.size(), 0.0);
        for (auto & impath : impaths) {
            Clock clock = "Task " + std::to_string(i) + " on \"" + impath + "\"";
            try {
                task(impath);
            } catch (...) {
                std::cout << "############### ERROR #############" << std::endl;
            }
        }
    }







    //bool show_gui = false;
    //
    //bool recache_gc = false;
    //bool recache_om = false;
    //bool recache_pn = false;

    //bool recache_gt = false;

    //bool consider_horizontal_cams_only = true;


    //// cache mode
    //// cache gt
    //auto getCachedGT = [&recache_gt, &matlab](const std::string & impath, const std::vector<PerspectiveCamera> & testCams) {
    //    std::vector<PerspectiveCamera> testCams1;
    //    std::vector<Imagei> gt1;
    //    if (recache_gt || !misc::LoadCache(impath, "testCams1_gt1", testCams1, gt1)) {
    //        auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    //        if (anno.nfaces() == 0) {
    //            while (true) {
    //                pano::experimental::EditLayoutAnnotation(impath, anno);
    //                pano::experimental::ReconstructLayoutAnnotation(anno, matlab);
    //                pano::experimental::VisualizeLayoutAnnotation(anno);
    //                int selected = pano::gui::SelectFrom({ "Accept", "Edit Again", "Abandon" },
    //                    "Your decision?",
    //                    "Accept the edit, or edit it again, or just abandon the edit this time?", 0, 2);
    //                if (selected == 0) {
    //                    pano::experimental::SaveLayoutAnnotation(impath, anno);
    //                    break;
    //                } else if (selected == 2) {
    //                    break;
    //                }
    //            }
    //        }
    //        testCams1 = CreatePanoContextCameras(anno.view.camera);
    //        gt1 = GTFaceLabels(anno, testCams1);
    //        misc::SaveCache(impath, "testCams1_gt1", testCams1, gt1);
    //    } 
    //    assert(testCams == testCams1);
    //    return gt1;
    //};

    //// cache gc
    //auto getCachedGC = [&recache_gc, &matlab](const std::string & impath, const std::vector<PerspectiveCamera> & testCams) {
    //    std::vector<PerspectiveCamera> testCams1;
    //    std::vector<Image7d> gc1;
    //    if (recache_gc || !misc::LoadCache(impath, "testCams1_gc1", testCams1, gc1)) {
    //        auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    //        testCams1 = CreatePanoContextCameras(anno.view.camera);
    //        gc1 = GeometricContext(anno, testCams1, matlab);
    //        misc::SaveCache(impath, "testCams1_gc1", testCams1, gc1);
    //    }
    //    assert(testCams == testCams1);
    //    return gc1;
    //};

    //// cache om
    //auto getCachedOM = [&recache_om](const std::string & impath, const std::vector<PerspectiveCamera> & testCams) {
    //    std::vector<PerspectiveCamera> testCams1;
    //    std::vector<Imagei> om1;
    //    if (recache_om || !misc::LoadCache(impath, "testCams1_om1", testCams1, om1)) {
    //        auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    //        testCams1 = CreatePanoContextCameras(anno.view.camera);
    //        om1 = OrientationMaps(anno, testCams1);
    //        misc::SaveCache(impath, "testCams1_om1", testCams1, om1);
    //    }
    //    assert(testCams == testCams1);
    //    return om1;
    //};

    //// cache pn
    //auto getCachedPN = [&recache_pn, &matlab, &show_gui](const std::string & impath, const std::vector<PerspectiveCamera> & testCams) {
    //    std::vector<PerspectiveCamera> testCams1;
    //    std::vector<Image3d> pn1;
    //    if (recache_pn || !misc::LoadCache(impath, "testCams1_pn1", testCams1, pn1)) {
    //        auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    //        testCams1 = CreatePanoContextCameras(anno.view.camera);
    //        
    //        PanoramixOptions options;
    //        options.useWallPrior = true;
    //        options.usePrincipleDirectionPrior = true;
    //        options.useGeometricContextPrior = true;

    //        options.useGTOcclusions = false;
    //        options.looseLinesSecondTime = false;
    //        options.looseSegsSecondTime = false;
    //        options.restrictSegsSecondTime = false;

    //        options.refresh_preparation = true;
    //        options.refresh_mg_init = options.refresh_preparation || false;
    //        options.refresh_mg_oriented = options.refresh_mg_init || false;

    //        options.refresh_line2leftRightSegs = options.refresh_mg_init || false;
    //        options.refresh_lsw = options.refresh_mg_init || false;
    //        options.refresh_mg_occdetected = options.refresh_lsw || options.refresh_line2leftRightSegs || true;
    //        options.refresh_mg_reconstructed = options.refresh_mg_occdetected || true;

    //        pn1 = Panoramix(anno, testCams1, matlab, options, show_gui);
    //        misc::SaveCache(impath, "testCams1_pn1", testCams1, pn1);
    //    }
    //    assert(testCams == testCams1);
    //    return pn1;
    //};







    //bool pn_vs_om = false;
    //bool pn_vs_gc = true;


    //// evaluate mode
    //if (pn_vs_om) {
    //    for (auto & impath : impaths) {
    //        auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    //        std::vector<PerspectiveCamera> testCams1 = CreatePanoContextCameras(anno.view.camera);
    //        
    //        std::vector<PerspectiveCamera> testCamsA, testCamsB, testCamsC;
    //        std::vector<Image3d> pn1;
    //        misc::LoadCache(impath, "testCams1_pn1", testCamsA, pn1);
    //        std::vector<Imagei> om1;
    //        misc::LoadCache(impath, "testCams1_om1", testCamsB, om1);
    //        std::vector<Image3d> gt1;
    //        misc::LoadCache(impath, "testCams1_gt1", testCamsC, gt1);

    //        assert(testCams1 == testCamsA && testCams1 == testCamsB && testCams1 == testCamsC);

    //        double error_pn = 0.0;
    //        double error_om = 0.0;

    //        int npixels = 0;
    //        for (int i = 0; i < testCams1.size(); i++) {
    //            auto & pn = pn1[i];
    //            auto & om = om1[i];
    //            auto & gt = gt1[i];

    //            for (auto it = gt.begin(); it != gt.end(); ++it) {
    //                Vec3 gtNormal = *it;
    //                if (norm(gtNormal) == 0) {
    //                    continue;
    //                }
    //                auto p = it.pos();

    //                Vec3 pnNormal = pn(p);
    //                //if (norm(pnNormal) == 0) {
    //                //    pnNormal = gtNormal;
    //                //}
    //                if (om(p) == -1 || norm(pnNormal) == 0) {
    //                    continue;
    //                }
    //                Vec3 omNormal = om(p) == -1 ? Origin() : anno.vps[om(p)];

    //                error_pn += ErrorOfSurfaceNormal(gtNormal, pnNormal);
    //                error_om += ErrorOfSurfaceNormal(gtNormal, omNormal);

    //                npixels++;
    //            }
    //        }

    //        error_pn /= npixels;
    //        error_om /= npixels;

    //        std::cout << "{{{PN[" << error_pn << "] vs OM[" << error_om << "]}}} on image\"" << impath << "\"" << std::endl;
    //    }      
    //}


    //if (pn_vs_gc) {

    //    auto ctable = gui::CreateGreyColorTableWithSize(4);
    //    std::map<std::string, double> impath2errorpn, impath2errorgc;

    //    for (auto & impath : impaths) {
    //        auto anno = LoadOrInitializeNewLayoutAnnotation(impath);
    //        std::vector<PerspectiveCamera> testCams1 = CreatePanoContextCameras(anno.view.camera);

    //        std::vector<PerspectiveCamera> testCamsA, testCamsB, testCamsC;
    //        std::vector<Image3d> pn1 = getCachedPN(impath, testCams1);
    //        std::vector<Image7d> gc1 = getCachedGC(impath, testCams1);
    //        std::vector<Imagei> gt1 = getCachedGT(impath, testCams1);

    //        double error_pn = 0.0;
    //        double error_gc = 0.0;

    //        int npixels = 0;
    //        std::vector<Imagei> gcLabelsTable(testCams1.size()), pnLabelsTable(testCams1.size()), gtLabelsTable(testCams1.size());
    //        for (int i = 0; i < testCams1.size(); i++) {

    //            auto & pn = pn1[i];
    //            auto & gc = gc1[i];
    //            auto & gt = gt1[i];

    //            auto & cam = testCams1[i];
    //            if (consider_horizontal_cams_only &&
    //                AngleBetweenUndirectedVectors(anno.view.camera.up(), cam.forward()) < DegreesToRadians(60)) {
    //                continue;
    //            }

    //            Imagei gcLabels(cam.screenSize(), -1);
    //            Imagei pnLabels(cam.screenSize(), -1);
    //            Imagei gtLabels(cam.screenSize(), -1);

    //            for (auto it = gt.begin(); it != gt.end(); ++it) {
    //                auto p = it.pos();
    //                int gtFaceId = *it;
    //                if (gtFaceId == -1) {
    //                    std::cout << "!!!gt face id is -1!!!!" << std::endl;
    //                    continue;
    //                }                    

    //                Vec3 gtNormal = anno.face2plane[gtFaceId].normal;
    //                if (norm(gtNormal) == 0) {
    //                    continue;
    //                }

    //                int gtLabel = SurfaceNormalToLabel(gtNormal, cam.up(), cam.forward(), cam.leftward());
    //                gtLabels(p) = gtLabel;

    //                Vec3 pnNormal = pn(p); 
    //                int pnLabel = SurfaceNormalToLabel(pnNormal, cam.up(), cam.forward(), cam.leftward());
    //                Vec7 gcScores = gc(p);
    //                int gcLabel = GeometricContextToLabel(gcScores);

    //                pnLabels(p) = pnLabel;
    //                gcLabels(p) = gcLabel;

    //                error_pn += pnLabel != gtLabel;
    //                error_gc += gcLabel != gtLabel;

    //                npixels++;
    //            }

    //            gcLabelsTable[i] = gcLabels;
    //            pnLabelsTable[i] = pnLabels;
    //            gtLabelsTable[i] = gtLabels;

    //            if (show_gui) {
    //                gui::MakeCanvas(ctable(gcLabels)).show(1, "gcLabels");
    //                gui::MakeCanvas(ctable(pnLabels)).show(1, "pnLabels");
    //                gui::MakeCanvas(ctable(gtLabels)).show(0, "gtLabels");
    //            }
    //        }

    //        error_pn /= npixels;
    //        error_gc /= npixels;

    //        std::cout << std::endl << std::endl;
    //        std::cout << "{{{{{{{{{  PN[" << error_pn << "] vs GC[" << error_gc << "]  }}}}}}}}}}}} on image\"" << impath << "\"" << std::endl;
    //        std::cout << std::endl << std::endl;

    //        impath2errorpn[impath] = error_pn;
    //        impath2errorgc[impath] = error_gc;
    //    }

    //    double mean_error_pn = 0.0;
    //    double mean_error_gc = 0.0;
    //    for (auto & impath : impaths) {
    //        double error_pn = impath2errorpn[impath];
    //        mean_error_pn += error_pn;
    //        double error_gc = impath2errorgc[impath];
    //        mean_error_gc += error_gc;
    //        std::cout << "{{{{{{{{{  PN[" << error_pn << "] vs GC[" << error_gc << "]  }}}}}}}}}}}} on image\"" << impath << "\"" << std::endl;
    //    }

    //    mean_error_pn /= impaths.size();
    //    mean_error_gc /= impaths.size();

    //    std::cout << "{{{{{{{{{  PN[" << mean_error_pn << "] vs GC[" << mean_error_gc << "]  }}}}}}}}}}}} on average" << std::endl;

    //}


    return 0;

}