#include <boost/filesystem.hpp>

#include "../../src/core/basic_types.hpp"
#include "../../src/core/single_view.hpp"
#include "../../src/experimental/rl_graph_solver.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/experimental/projective_solver.hpp"
#include "../../src/misc/matlab_api.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/gui/canvas.hpp"
#include "../../src/gui/utility.hpp"

#include "routines.hpp"


namespace panolyz {

    namespace PanoContext {


        using namespace pano;
        using namespace pano::core;
        using namespace pano::experimental;

        static const std::string root = "H:\\DataSet\\PanoContext\\bedroom";

        struct Task;
        struct GlobalConfig {
            misc::MXA anno;
            struct AnnoRule {
                int type;
                std::vector<int> linefroms, linetos;
            };
            std::map<int, AnnoRule> annorules;

            explicit GlobalConfig(const std::string & root);
            Failable<Task> newTask(int id) const;
        };

    
        struct Task {
            std::string name;
            std::string path;
            Image3ub originalImage;
            View<PanoramicCamera, Image3ub> view;
            std::vector<Vec3> vps;
            int vertVPId;
            struct Object {
                std::vector<Point3> originalPoints;
                int type;
                std::string name;

                std::vector<Point3> worldPoints;
                TransformedIn3D<Box3> worldBox;
                View<PerspectiveCamera, Imagef> maskView;
                std::vector<Point3> biggerMaskCorners;
                View<PerspectiveCamera, Imagef> biggerMaskView;
            };            
            std::vector<Object> objects;


            std::vector<gui::PenConfig> penConfigs;

            // features
            Chain3 floorChain, ceilingChain;
            Image5d gc;
            std::vector<Classified<Line3>> lines;
            
            Imagei segs; 
            int nsegs; 
            std::vector<int> segIdsOld2New;
            std::vector<std::vector<std::vector<Vec3>>> segContours;
            std::vector<Vec3> segCenters;
            std::vector<double> segAreas;
            Image3ub segimWithObjectsMask;
            enum SegPosition {
                InFloor, InCeiling, VerticalPlane, FreePlane, InObject
            };
            std::vector<SegPosition> segPositions;


            void updateChains();
            void calcGC(bool refresh = false);
            void calcLines(bool refresh = false);
            //void calcSegs(bool refresh = false);

            void visualize();

            void solve();
        };



        namespace util {


            void Visualize(const PanoramicView & view, const Imagei & segs, int nsegs,
                const std::vector<Vec3> & vps,
                const std::vector<Task::Object> & objects,
                const  std::map<int, GlobalConfig::AnnoRule> & annorules) {
                // visualize objects
                gui::SceneBuilder sb;
                // draw image
                gui::ResourceStore::set("tex", view.image);
                auto segim = gui::CreateRandomColorTableWithSize(nsegs)(segs);
                gui::ResourceStore::set("segim", segim);
                sb.begin(MakeTransformableIn3D(UnitSphere()).scale(100))
                    .shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("segim").end();

                // draw objects (as lines)
                for (auto & o : objects) {
                    if (o.originalPoints.empty() || o.name == "room")
                        continue;

                    // obj visible lines
                    auto & pp = o.worldPoints;
                    auto & ar = annorules.at(o.type);
                    std::vector<Line3> objLines;
                    objLines.reserve(ar.linefroms.size());
                    for (int i = 0; i < ar.linefroms.size(); i++) {
                        objLines.emplace_back(pp.at(ar.linefroms[i]), pp.at(ar.linetos[i]));
                        objLines.back() = normalize(objLines.back()) / 2.0;
                    }

                    // mask contour lines
                    std::vector<Line3> chlines;
                    for (int i = 0; i < o.biggerMaskCorners.size(); i++) {
                        chlines.emplace_back(o.biggerMaskCorners[i] / 2,
                            o.biggerMaskCorners[(i + 1) % o.biggerMaskCorners.size()] / 2);
                    }

                    sb.begin(o.worldBox, [o](gui::InteractionID iid, const TransformedIn3D<Box3> & b) {
                        std::cout << "name: " << o.name << "\t type: " << o.type << std::endl;
                    }).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("tex").end();
                    sb.installingOptions().discretizeOptions.color = gui::Yellow;
                    sb.begin(objLines).lineWidth(5).shaderSource(gui::OpenGLShaderSourceDescriptor::XLines).end();
                    sb.begin(chlines).lineWidth(10).shaderSource(gui::OpenGLShaderSourceDescriptor::XLines).end();
                    sb.begin(pp).shaderSource(gui::OpenGLShaderSourceDescriptor::XPoints).end();
                }

                sb.installingOptions().pointSize = 50.0;
                sb.installingOptions().discretizeOptions.color = gui::Red;
                sb.add(Point3(0.5, 0, 0));
                sb.installingOptions().discretizeOptions.color = gui::Green;
                sb.add(Point3(0, 0.5, 0));
                sb.installingOptions().discretizeOptions.color = gui::Blue;
                sb.add(Point3(0, 0, 0.5));

                sb.installingOptions().pointSize = 30.0;
                sb.installingOptions().discretizeOptions.color = gui::Red;
                sb.begin(vps[0] / 2.0).shaderSource(gui::OpenGLShaderSourceDescriptor::XPoints).end();
                sb.installingOptions().discretizeOptions.color = gui::Green;
                sb.begin(vps[1] / 2.0).shaderSource(gui::OpenGLShaderSourceDescriptor::XPoints).end();
                sb.installingOptions().discretizeOptions.color = gui::Blue;
                sb.begin(vps[2] / 2.0).shaderSource(gui::OpenGLShaderSourceDescriptor::XPoints).end();

                for (int i = 0; i < 3; i++) {
                    std::cout << "vps[" << i << "] = " << vps[i] << std::endl;
                }


                gui::VisualizeWithPanoramicOperation(sb.scene(),
                    gui::RenderOptions().bwColor(0).bwTexColor(1.0)
                    .backgroundColor(gui::White).camera(PerspectiveCamera(800, 800, Point2(400, 400), 500)));
            }

            std::vector<bool> ComputeObjectMaskedSegs(const View<PanoramicCamera, Image3ub> & view,
                const std::vector<Task::Object> & objects,
                const Imagei & segs, int nsegs,
                Image3ub & segimWithObjectsMask) {

                auto segim = gui::CreateRandomColorTableWithSize(nsegs)(segs);
                std::vector<bool> segsMaskedByObjects(nsegs, true);
                std::vector<View<PerspectiveCamera, Imagef>> masks;
                for (auto & o : objects) {
                    if (o.biggerMaskCorners.empty())
                        continue;
                    if (o.name == "room" || o.name == "doorway")
                        continue;
                    masks.push_back(o.biggerMaskView);
                }
                auto objectsMask = Combine(view.camera, masks).image;
                assert(segs.size() == objectsMask.size());
                for (auto it = segs.begin(); it != segs.end(); ++it) {
                    if (objectsMask(it.pos()) == 0) { // not masked
                        segsMaskedByObjects[*it] = false;
                    }
                }

                segimWithObjectsMask = segim.clone();
                for (auto it = segimWithObjectsMask.begin(); it != segimWithObjectsMask.end(); ++it) {
                    if (segsMaskedByObjects[segs(it.pos())]) {
                        *it = view.image(it.pos());
                    } else {
                        *it = (Vec3(*it) * 1.0 + Vec3(view.image(it.pos()))) / 2.0;
                    }
                }

                return segsMaskedByObjects;
            }

        }




        GlobalConfig::GlobalConfig(const std::string & root) {
            misc::MAT annotations(root + "\\ANNO_ALL.mat", misc::MAT::Read);
            anno = annotations.var("ANNO_ALL");
            int annoLength = anno.length();
            // load anno rules
            misc::MAT annorulesFile(root + "\\ANNORULE.mat", misc::MAT::Read);
            auto annorulesData = annorulesFile.var("annorule");
            int nannorules = annorulesData.length();
            for (int i = 0; i < nannorules; i++) {
                int type = annorulesData.field("type", i);
                auto & ar = annorules[type];
                ar.type = type;
                auto lineData = annorulesData.field("line", i);
                if (lineData.null())
                    continue;
                for (int j = 0; j < lineData.m(); j++) {
                    ar.linefroms.push_back(lineData.at<double>(j, 0) - 1);
                    ar.linetos.push_back(lineData.at<double>(j, 1) - 1);
                }
            }
        }

        Failable<Task> GlobalConfig::newTask(int id) const {
            Task task;
            task.name = anno.field("name", id);

            // image
            task.path = root + "\\" + task.name + "\\" + task.name + ".jpg";
            task.originalImage = gui::ImageRead(task.path);
            auto image = task.originalImage.clone();
            ResizeToHeight(image, 1600);

            task.view = CreatePanoramicView(image);

            /// load groundtruth data
            task.vps.resize(3);
            DenseMatd vpsData = anno.field("vp", id);
            assert(vpsData.rows == 6 && vpsData.cols == 3);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    task.vps[i][j] = vpsData(i, j);
                }
            }
            task.vps = {
                { task.vps[2][1], task.vps[2][0], task.vps[2][2] },
                { task.vps[1][1], task.vps[1][0], task.vps[1][2] },
                { -task.vps[0][1], -task.vps[0][0], -task.vps[0][2] }
            };
            for (auto & vp : task.vps) {
                vp = normalize(vp);
            }

            // objects
            auto anno3d = anno.field("ANNO3D", id);
            if (!anno3d) {
                return nullptr;
            }
            auto objects3DData = anno3d.field("objects3D");
            if (!objects3DData) {
                return nullptr;
            }

            int nobjects = objects3DData.length();


            task.objects.resize(nobjects);
            for (int i = 0; i < nobjects; i++) {
                DenseMatd points = objects3DData.field("out_points_w", i);
                task.objects[i].originalPoints.resize(points.rows);
                for (int j = 0; j < points.rows; j++) {
                    for (int k = 0; k < 3; k++) {
                        task.objects[i].originalPoints[j][k] = points(j, k);
                    }
                }
                task.objects[i].type = objects3DData.field("type", i);
                task.objects[i].name = objects3DData.field("name", i);
                task.objects[i].worldBox =
                    AsInLocalCoordinates(BoundingBoxOfContainer(task.objects[i].originalPoints), 
                    task.vps[0], task.vps[1], task.vps[2], Origin()).scale(0.001);
                task.objects[i].worldPoints.reserve(task.objects[i].originalPoints.size());
                for (auto & p : task.objects[i].originalPoints) {
                    task.objects[i].worldPoints.push_back(task.objects[i].worldBox.toWorld(normalize(p)) / 2);
                }

                if (task.objects[i].name == "room" || task.objects[i].worldPoints.empty())
                    continue;

                // compute mask corners
                // project partially
                std::vector<int> ch;
                {
                    Vec3 center;
                    for (auto & p : task.objects[i].worldPoints) {
                        center += p;
                    }
                    center /= norm(center);
                    double angle = 0.0;
                    for (auto & p : task.objects[i].worldPoints) {
                        angle = std::max(angle, AngleBetweenDirections(center, p));
                    }

                    // project on to a persp cam
                    static const double focal = 500;
                    double sz = focal * tan(angle) * 2;
                    task.objects[i].maskView.camera = PerspectiveCamera(sz, sz,
                        Point2(sz, sz) / 2.0, focal,
                        Origin(), center);
                    std::vector<Point2f> pp2s(task.objects[i].worldPoints.size());
                    for (int j = 0; j < task.objects[i].worldPoints.size(); j++) {
                        pp2s[j] = task.objects[i].maskView.camera.toScreen(task.objects[i].worldPoints[j]);
                    }
                    // get convex hull
                    cv::convexHull(pp2s, ch);
                    task.objects[i].maskView.image = Imagef(task.objects[i].maskView.camera.screenSize(), 0.0f);
                    std::vector<Point2i> pp2is;
                    for (int id : ch) {
                        pp2is.push_back(pp2s[id]);
                    }
                    cv::fillConvexPoly(task.objects[i].maskView.image, pp2is, 1.0);
                }

                {
                    // get space mask
                    Vec3 center;
                    for (int id : ch) {
                        center += task.objects[i].worldPoints[id];
                    }
                    center /= norm(center);

                    // expanded hull points
                    auto & chpoints = task.objects[i].biggerMaskCorners;
                    chpoints.resize(ch.size());
                    for (int j = 0; j < ch.size(); j++) {
                        chpoints[j] = normalize(task.objects[i].worldPoints[ch[j]]);
                        auto v = normalize(chpoints[j] - center);
                        chpoints[j] += (v * tan(0.1));
                        chpoints[j] /= norm(chpoints[j]);
                    }

                    double angle = 0.0;
                    for (auto & p : chpoints) {
                        angle = std::max(angle, AngleBetweenDirections(center, p));
                    }
                    // project on to a persp cam
                    static const double focal = 500;
                    double sz = focal * tan(angle) * 2;
                    task.objects[i].biggerMaskView.camera = PerspectiveCamera(sz, sz,
                        Point2(sz, sz) / 2.0, focal,
                        Origin(), center);
                    std::vector<Point2i> pp2s(task.objects[i].biggerMaskCorners.size());
                    for (int j = 0; j < task.objects[i].biggerMaskCorners.size(); j++) {
                        pp2s[j] = task.objects[i].biggerMaskView.camera.toScreen(task.objects[i].biggerMaskCorners[j]);
                    }
                    task.objects[i].biggerMaskView.image = Imagef(task.objects[i].biggerMaskView.camera.screenSize(), 0.0f);
                    cv::fillConvexPoly(task.objects[i].biggerMaskView.image, pp2s, 1.0);
                }
            }


            task.penConfigs = {
                { "ceiling boundary", "the boundary between ceiling and wall", 5.0, gui::Blue, gui::DashLine },
                { "floor boundary", "the boundary between floor and wall", 5.0, gui::Yellow, gui::DashLine }
            };

            return std::move(task);
        }

        void Task::updateChains() {
            // assign chains
            // draw chains
            Load(path, "chains", floorChain, ceilingChain);
            std::vector<Chain3> chains = { std::move(ceilingChain), std::move(floorChain) };
            gui::DrawChainsInPanorama(view, penConfigs, chains);
            ceilingChain = std::move(chains[0]);
            floorChain = std::move(chains[1]);
            Save(path, "chains", floorChain, ceilingChain);
        }

        void Task::calcGC(bool refresh) {
            misc::Matlab matlab;

            // gc
            std::vector<PerspectiveCamera> hcams;
            if (refresh || !Load(path, "hcams", hcams)) {
                hcams = CreatePanoContextCameras(view.camera, 500, 400, 300);
                Save(path, "hcams", hcams);
            }

            // extract gcs
            std::vector<Weighted<View<PerspectiveCamera, Image5d>>> gcs;
            if (refresh || !Load(path, "gcs", gcs)) {
                gcs.resize(hcams.size());
                for (int i = 0; i < hcams.size(); i++) {
                    auto pim = view.sampled(hcams[i]);
                    auto pgc = ComputeGeometricContext(matlab, pim.image, false, true);
                    gcs[i].component.camera = hcams[i];
                    gcs[i].component.image = pgc;
                    gcs[i].score = sin(AngleBetweenUndirectedVectors(hcams[i].forward(), view.camera.up()));
                }
                Save(path, "gcs", gcs);
            }

            if (refresh || !Load(path, "gc", gc)) {
                gc = Combine(view.camera, gcs).image;
                Save(path, "gc", gc);
            }
            gui::MakeCanvas(gc).maxHeight(500).show();
        }

        void Task::calcLines(bool refresh) {
            // collect lines in each view
            auto cams = CreateCubicFacedCameras(view.camera, 
                view.image.rows, view.image.rows, view.image.rows * 0.4);
            std::vector<std::vector<Classified<Line2>>> clines(cams.size());
            for (int i = 0; i < cams.size(); i++) {
                auto pim = view.sampled(cams[i]).image;
                LineSegmentExtractor lineExtractor;
                lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
                auto ls = lineExtractor(pim, 2, 300); // use pyramid
                clines[i].reserve(ls.size());
                for (auto & l : ls) {
                    clines[i].push_back(ClassifyAs(l, -1));
                }
            }
            vertVPId = NearestDirectionId(vps, Z());

            // extract lines from segmentated region boundaries and classify them using estimated vps
            // make 3d lines
            for (int i = 0; i < cams.size(); i++) {
                for (auto & l : clines[i]) {
                    lines.push_back(ClassifyAs(Line3(
                        normalize(cams[i].toSpace(l.component.first)),
                        normalize(cams[i].toSpace(l.component.second))),
                        -1));
                }
            }

            // classify lines
            ClassifyLines(lines, vps, M_PI / 3.0, 0.1, 0.8, M_PI / 18.0);
        }

        //void Task::calcSegs(bool refresh) {
        //    // extract segments
        //    if (refresh || 
        //        !Load(path, "segs", segs, nsegs, validSegIds, segContours, segCenters, segAreas, segimWithObjectsMask, segPositions)) {
        //        SegmentationExtractor segmenter;
        //        segmenter.params().algorithm = SegmentationExtractor::GraphCut;
        //        std::vector<Line3> chainLines;
        //        for (auto & l : lines) {
        //            chainLines.push_back(normalize(l.component));
        //        }
        //        for (auto & c : {ceilingChain, floorChain}) {
        //            for (int i = 0; i < c.size(); i++) {
        //                chainLines.push_back(c.edge(i));
        //            }
        //        }
        //        std::tie(segs, nsegs) = segmenter(view.image, chainLines, view.camera);

        //        validSegIds = ComputeSpatialRegionProperties(segs, 
        //            view.camera, &segContours, &segCenters, &segAreas);
        //        std::vector<bool> segsMaskedByObjects =
        //            util::ComputeObjectMaskedSegs(view, objects, segs, nsegs, segimWithObjectsMask);

        //        // classify segments
        //        segPositions.resize(validSegIds.size());
        //        auto up = vps[vertVPId], down = -vps[vertVPId];
        //        /*  auto floorView = PerfectRegionMaskView(floorChain.points, down, 100);
        //          auto ceilingView = PerfectRegionMaskView(ceilingChain.points, up, 100);
        //          for (int i = 0; i < validSegIds.size(); i++) {
        //          if (segsMaskedByObjects[validSegIds[i]]) {
        //          segPositions[i] = InObject;
        //          continue;
        //          }
        //          auto & contour = segContours[i];
        //          bool allInFloor = true, allInCeiling = true;
        //          IterateOver(contour, [&allInFloor, &floorView](const Vec3 & c) {
        //          auto & in = allInFloor;
        //          auto & v = floorView;
        //          in = in &&
        //          v.camera.isVisibleOnScreen(c);
        //          if (in) {
        //          auto proj = ToPixel(v.camera.toScreen(c));
        //          in = in && Contains(v.image, proj) && v.image(proj);
        //          }
        //          });
        //          if (allInFloor) {
        //          segPositions[i] = InFloor;
        //          continue;
        //          }
        //          IterateOver(contour, [&allInCeiling, &ceilingView](const Vec3 & c) {
        //          auto & in = allInCeiling;
        //          auto & v = ceilingView;
        //          in = in &&
        //          v.camera.isVisibleOnScreen(c);
        //          if (in) {
        //          auto proj = ToPixel(v.camera.toScreen(c));
        //          in = in && Contains(v.image, proj) && v.image(proj);
        //          }
        //          });
        //          if (allInCeiling) {
        //          segPositions[i] = InCeiling;
        //          continue;
        //          }

        //          bool hasUpperPart = false, hasLowerPart = false;
        //          IterateOver(contour, [&up, &hasLowerPart, &hasUpperPart](const Vec3 & c) {
        //          double d = c.dot(up);
        //          if (d < 0) {
        //          hasLowerPart = true;
        //          } else if (d > 0) {
        //          hasUpperPart = true;
        //          }
        //          });
        //          if (hasUpperPart && hasLowerPart) {
        //          segPositions[i] = VerticalPlane;
        //          continue;
        //          }
        //          segPositions[i] = FreePlane;
        //          }*/

        //        Save(path, "segs", segs, nsegs, validSegIds, segContours, segCenters, segAreas, segimWithObjectsMask, segPositions);
        //    }
        //}

        void Task::visualize() {
            gui::VisualizeAll(view, lines, segs, nsegs, gc);
        }

        void Task::solve() {

            //// create projective solver
            //ProjectiveSolver solver;

            //// collect planes
            //Plane3 floor(-vps[vertVPId], vps[vertVPId]), ceiling(vps[vertVPId], -vps[vertVPId]);
            //auto cfloor = solver.bindPlaneDoF3(floor);
            //auto cceiling = solver.bindPlaneDoF3(ceiling);

            //std::vector<Plane3> segPlanes(nsegs);
            //std::vector<int> segCompIds(nsegs, -1);
            //for (int i = 0; i < validSegIds.size(); i++) {
            //    auto position = segPositions[i];
            //    if (position == FreePlane) {
            //        segCompIds[validSegIds[i]] = solver.bindPlaneDoF3(segPlanes[validSegIds[i]]);
            //    } else if (position == VerticalPlane) {
            //        segCompIds[validSegIds[i]] = solver.bindPlaneDoF2(segPlanes[validSegIds[i]], vps[vertVPId]);
            //    }
            //}

            //// anchors
            //auto boundaries = FindRegionBoundaries(segs, 2);

            //std::set<std::pair<int, int>> disconnectedSegIdPairs;
        }







        void Run(){

            GlobalConfig cc(root);

            while (true){

                std::cout << "INPUT ID [0~441 except 12]: " << std::endl;
                int id = 0;
                std::cin >> id;
                std::cout << "ID = " << id << std::endl;
                
                auto t = cc.newTask(id);
                if (t.failed())
                    continue;

                Task task = t;
                task.calcLines(true);
                //task.calcSegs(false);
                task.calcGC(false);
                             
                task.visualize();
               


                

                continue;
            }

        }
    }
}