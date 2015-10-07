
#include "eval.hpp"

#include "../../src/experimental/rl_graph_solver.hpp"
#include "../../src/experimental/rl_graph_occlusion.hpp"
#include "../../src/experimental/rl_graph_modeling.hpp"
#include "../../src/experimental/tools.hpp"

#include "../../src/experimental/pi_graph_control.hpp"
#include "../../src/experimental/pi_graph_occlusion.hpp"
//#include "../../src/experimental/pi_graph_solve.hpp"
#include "../../src/experimental/pi_graph_vis.hpp"

namespace panolyz {

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
            return norm(IntersectionOfLineAndPlane(Ray3(Origin(), direction), mg.seg2recPlanes[seg]).position);
        }
        virtual void visualize(const std::vector<Vec3> & directions) const {
            //VisualizeReconstruction({ mg.ccidsBigToSmall.front() }, mg);
        }
    };


    std::unique_ptr<ReconstructedModel> PredictionOfPanoramix(const std::string & impath,
        const PredictOptions & options, misc::Matlab & matlab, const PILayoutAnnotation & anno) {

        auto resultPath = PanoramixResultFilePath(impath, 0);

        PIGraph mg;
        if (1 || !LoadFromDisk(resultPath, mg)) {
            // rebuild pigraph
            misc::Matlab matlab;

            Image3ub originalImage = anno.rectifiedImage;

            auto image = originalImage.clone();
            ResizeToHeight(image, 700);

            /// prepare things!
            View<PanoramicCamera, Image3ub> view;
            std::vector<PerspectiveCamera> cams;
            std::vector<Classified<Line3>> line3s;
            std::vector<Vec3> vps;
            int vertVPId;
            Imagei segs;
            int nsegs;


            if (1 || !misc::LoadCache(impath, "preparation", view, cams, line3s, vps, vertVPId, segs, nsegs)) {
                view = CreatePanoramicView(image);

                // collect lines in each view
                cams = CreateCubicFacedCameras(view.camera, image.rows, image.rows, image.rows * 0.4);
                std::vector<Line3> rawLine3s;
                std::vector<Line2> rawLine2s;
                for (int i = 0; i < cams.size(); i++) {
                    auto pim = view.sampled(cams[i]).image;
                    LineSegmentExtractor lineExtractor;
                    lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
                    auto ls = lineExtractor(pim, 2, 300); // use pyramid
                    for (auto & l : ls) {
                        rawLine2s.push_back(l);
                        rawLine3s.emplace_back(normalize(cams[i].toSpace(l.first)),
                            normalize(cams[i].toSpace(l.second)));
                    }
                }
                rawLine3s = MergeLines(rawLine3s, DegreesToRadians(3), DegreesToRadians(5));

                // estimate vp
                line3s = ClassifyEachAs(rawLine3s, -1);
                vps = EstimateVanishingPointsAndClassifyLines(line3s);
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
                SegmentationExtractor segmenter;
                segmenter.params().algorithm = SegmentationExtractor::GraphCut;
                segmenter.params().sigma = 10.0;
                segmenter.params().c = 1.0;
                segmenter.params().superpixelSizeSuggestion = 2000;
                segmenter.params().useYUVColorSpace = true;
                std::tie(segs, nsegs) = segmenter(view.image, rawLine3s, view.camera, DegreesToRadians(1));
                RemoveThinRegionInSegmentation(segs, true);
                nsegs = DensifySegmentation(segs, true);
                assert(IsDenseSegmentation(segs));

                if (1) {
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
                misc::SaveCache(impath, "preparation", view, cams, line3s, vps, vertVPId, segs, nsegs);
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
                if (1) {
                    std::vector<Imaged> gcChannels;
                    cv::split(gc, gcChannels);
                    gui::AsCanvas(ConvertToImage3d(gc)).show();
                }
                misc::SaveCache(impath, gcmergedFileName, gc);
            }



            // build pigraph!
            mg = BuildPIGraph(view, vps, vertVPId, segs, line3s,
                DegreesToRadians(1), DegreesToRadians(1), DegreesToRadians(2),
                DegreesToRadians(5), DegreesToRadians(60), DegreesToRadians(5));

            const auto printPIGraph = [&mg](int delay){
                static const gui::ColorTable randColors = gui::CreateRandomColorTableWithSize(mg.nsegs);
                static const gui::ColorTable rgbGrayColors = gui::RGBGreys;
                auto pim = Print(mg,
                    [&mg](int seg) -> gui::Color {
                    return gui::White;
                },
                    [&mg](int lp) {
                    return mg.linePiece2bndPiece[lp] == -1 ? gui::Red : gui::Black;
                },
                    [&mg](int bp) -> gui::Color {
                    return gui::LightGray;
                }, 1, 3);
                cv::imshow("pigraph", Image3f(pim * 0.99f) + Image3f(mg.view.image / 255.0f * 0.01f));
                cv::waitKey(delay);
            };

            printPIGraph(0);

            // attach orientation constraints
            AttachPrincipleDirectionConstraints(mg);
            AttachWallConstraints(mg, M_PI / 60.0);
            AttachGCConstraints(mg, gc);


            const auto printPIGraphControls = [&mg](int delay) {
                const auto bpColor = [&mg](int bp) -> gui::Color {
                    auto occ = mg.bndPiece2segRelation[bp];
                    if (occ == SegRelation::Coplanar) {
                        return gui::LightGray;
                    } else if (occ == SegRelation::Connected) {
                        return gui::Gray;
                    } else if (occ == SegRelation::LeftIsFront) {
                        return gui::Color(gui::Red) * 0.6;
                    } else {
                        return gui::Color(gui::Blue) * 0.6;
                    }
                };
                const auto lpColor = [&mg, bpColor](int lp) -> gui::Color {
                    if (mg.linePiece2bndPiece[lp] != -1) {
                        return bpColor(mg.linePiece2bndPiece[lp]) * 0.9;
                    } else {
                        return mg.linePiece2segLineRelation[lp] == SegLineRelation::Detached ? gui::Black : gui::LightGray;
                    }
                };
                const auto segColor = [&mg](int seg) -> gui::Color {
                    static const gui::ColorTable ctable = gui::ColorTableDescriptor::RGBGreys;
                    auto & c = mg.seg2control[seg];
                    if (!c.used) {
                        return gui::Black;
                    }
                    if (c.orientationClaz != -1) {
                        return ctable[c.orientationClaz];
                    }
                    if (c.orientationNotClaz != -1) {
                        return ctable[c.orientationNotClaz] * 0.7;
                    }
                    return gui::White;
                };
                auto pim = Print(mg, segColor, ConstantFunctor<gui::Color>(gui::Transparent), bpColor, 2, 3);
                cv::imshow("pigraph_controled", Image3f(pim * 0.99f) + Image3f(mg.view.image / 255.0f * 0.01f));
                cv::waitKey(delay);
            };

            printPIGraphControls(0);
            DetectOcclusions2(mg);
           

            //// initialize connectivities
            //for (int i = 0; i < mg.nbndPieces(); i++) {
            //    mg.bndPiece2segRelation[i] = SegRelation::Coplanar;
            //}
            //for (int i = 0; i < mg.nlinePieces(); i++) {
            //    mg.linePiece2segLineRelation[i] = SegLineRelation::Attached;
            //}
            //struct Vertex {
            //    enum Type {
            //        IsSeg,
            //        IsLine
            //    } type;
            //    bool isSeg() const { return type == IsSeg; }
            //    bool isLine() const { return type == IsLine; }
            //    int id;
            //};
            //struct Edge {
            //    int vert1, vert2;
            //    enum Type {
            //        SegSeg,
            //        SegLine,
            //        LineLine
            //    } type;
            //    int id; // SegSeg->bndPiece, SegLine->linePiece, LineLine->lineRelation
            //    double weight;
            //    double width;
            //    std::vector<Vec3> anchors;
            //};
            //std::vector<Vertex> verts;
            //std::vector<int> seg2vert(mg.nsegs, -1), line2vert(mg.nlines(), -1);
            //std::vector<Edge> edges;
            //std::vector<int> bndPiece2edge(mg.nbndPieces(), -1), linePiece2edge(mg.nlinePieces(), -1), lineRelation2edge(mg.nlineRelations(), -1);
            //
            //// build constraint graph
            //{
            //    // add verts
            //    // add segs
            //    for (int i = 0; i < mg.nsegs; i++) {
            //        Vertex v;
            //        v.type = Vertex::IsSeg;
            //        v.id = i;
            //        verts.push_back(v);
            //        int vert = verts.size() - 1;
            //        seg2vert[i] = vert;
            //    }
            //    // add lines
            //    for (int i = 0; i < mg.nlines(); i++) {
            //        Vertex v;
            //        v.type = Vertex::IsLine;
            //        v.id = i;
            //        verts.push_back(v);
            //        int vert = verts.size() - 1;
            //        line2vert[i] = vert;
            //    }
            //    // add edges
            //    // add bndPieces
            //    static const double minAngleThresholdForAWideEdge = DegreesToRadians(5);
            //    for (int i = 0; i < mg.nbndPieces(); i++) {
            //        Edge e;
            //        e.type = Edge::SegSeg;
            //        e.id = i;
            //        auto & segPair = mg.bnd2segs[mg.bndPiece2bnd[i]];
            //        e.vert1 = seg2vert[segPair.first];
            //        e.vert2 = seg2vert[segPair.second];
            //        double len = mg.bndPiece2length[i];
            //        e.width = len;
            //        if (len >= minAngleThresholdForAWideEdge) {
            //            e.anchors = { mg.bndPiece2dirs[i].front(), mg.bndPiece2dirs[i].back() };
            //        } else {
            //            e.anchors = { normalize(mg.bndPiece2dirs[i].front() + mg.bndPiece2dirs[i].back()) };
            //        }
            //        e.weight = 1.0;
            //        edges.push_back(e);
            //        int edge = edges.size() - 1;
            //        bndPiece2edge[i] = edge;
            //    }
            //    // add linePieces
            //    for (int i = 0; i < mg.nlinePieces(); i++) {
            //        int bndPiece = mg.linePiece2bndPiece[i];
            //        int line = mg.linePiece2line[i];
            //        if (bndPiece == -1) {
            //            int seg = mg.linePiece2seg[i];
            //            assert(seg != -1);
            //            Edge e;
            //            e.type = Edge::SegLine;
            //            e.id = i;
            //            e.vert1 = seg2vert[seg];
            //            e.vert2 = line2vert[line];
            //            double len = mg.linePiece2length[i];
            //            e.width = len;
            //            if (len >= minAngleThresholdForAWideEdge) {
            //                e.anchors = { mg.linePiece2samples[i].front(), mg.linePiece2samples[i].back() };
            //            } else {
            //                e.anchors = { normalize(mg.linePiece2samples[i].front() + mg.linePiece2samples[i].back()) };
            //            }
            //            e.weight = 1.0;
            //            edges.push_back(e);
            //            int edge = edges.size() - 1;
            //            linePiece2edge[i] = edge;
            //        } else {
            //            int segPair[] = { -1, -1 };
            //            std::tie(segPair[0], segPair[1]) = mg.bnd2segs[mg.bndPiece2bnd[bndPiece]];
            //            for (int seg : segPair) {
            //                Edge e;
            //                e.type = Edge::SegLine;
            //                e.id = i;
            //                e.vert1 = seg2vert[seg];
            //                e.vert2 = line2vert[line];
            //                double len = mg.linePiece2length[i];
            //                e.width = len;
            //                if (len >= minAngleThresholdForAWideEdge) {
            //                    e.anchors = { mg.linePiece2samples[i].front(), mg.linePiece2samples[i].back() };
            //                } else {
            //                    e.anchors = { normalize(mg.linePiece2samples[i].front() + mg.linePiece2samples[i].back()) };
            //                }
            //                e.weight = 1.0;
            //                edges.push_back(e);
            //                int edge = edges.size() - 1;
            //                linePiece2edge[i] = edge;
            //            }
            //        }
            //    }
            //    // add lineRelations
            //    for (int i = 0; i < mg.nlineRelations(); i++) {
            //        Edge e;
            //        e.type = Edge::LineLine;
            //        e.id = i;
            //        auto & linePair = mg.lineRelation2lines[i];
            //        e.vert1 = line2vert[linePair.first];
            //        e.vert2 = line2vert[linePair.second];
            //        e.width = 1.0;
            //        e.weight = 1.0;
            //        edges.push_back(e);
            //        int edge = edges.size() - 1;
            //        lineRelation2edge[i] = edge;
            //    }
            //}
            //// adjust weights 
            //{

            //}




            while (true) {
               

                // fetch the largest cc(?)
                // and optimze it
                {


                }

            }





            // save to disk
            SaveToDisk(resultPath, mg);
        }

        return std::make_unique<PIGraphModel>(std::move(mg));
    }

}

 //{
 //    auto cctable = gui::CreateRandomColorTableWithSize(mg.nccs);
 //    auto pim = Print(mg,
 //        [&mg, &cctable](int seg) -> gui::Color {
 //        int vert = mg.seg2vert[seg];
 //        if (vert == -1) return gui::Transparent;
 //        return cctable[mg.vert2cc[vert]];
 //    },
 //        [&mg, &cctable](int lp) -> gui::Color {
 //        int vert = mg.line2vert[mg.linePiece2line[lp]];
 //        if (vert == -1) return gui::Transparent;
 //        return cctable[mg.vert2cc[vert]];
 //    },
 //        [](int bp) -> gui::Color {
 //        return gui::LightGray;
 //    }, 1, 2);
 //    cv::imshow("pigraph_ccs", Image3f(pim * 0.99f) + Image3f(mg.view.image / 255.0f * 0.01f));
 //    cv::waitKey();
 //}