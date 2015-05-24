#include <QtCore>

#include "../../src/core/basic_types.hpp"
#include "../../src/experimental/rl_graph.hpp"
#include "../../src/experimental/tools.hpp"
#include "../../src/experimental/rl_graph_occlusion.hpp"
#include "../../src/gui/canvas.hpp"
#include "../../src/gui/scene.hpp"
#include "../../src/core/cameras.hpp"
#include "../../src/core/clock.hpp"

#include "routines.hpp"

namespace panolyz {

    namespace Prepare {

        using namespace panoramix;
        using namespace core;
        using namespace experimental;

        using Image7d = ImageOf<Vec<double, 7>>;

        void RunNormal(std::string folder, bool indoor){
            QString qfolder = QString::fromStdString(folder);
            for (auto & finfo : QDir(qfolder + "\\image").entryInfoList(QStringList() << "*.jpg" << "*.png")){
                qDebug() << finfo.absoluteFilePath();

                std::string path = finfo.absoluteFilePath().toStdString();
                std::string filename = finfo.baseName().toStdString();

                Image3ub image = cv::imread(path);
                ResizeToHeight(image, 600);

                View<PerspectiveCamera> view;
                std::vector<Classified<Line2>> lines;
                std::vector<Vec3> vps;
                double focal;
                std::vector<Classified<Line3>> line3s;

                qDebug() << "computing view";
                VanishingPointsDetector::Params vpdParams(VanishingPointsDetector::TardifSimplified);
                view = CreatePerspectiveView(image, Point3(0, 0, 0), Point3(1, 0, 0), Point3(0, 0, -1),
                    LineSegmentExtractor(), VanishingPointsDetector(vpdParams), &line3s, &lines, &vps, &focal).unwrap();
                std::vector<HPoint2> hvps(vps.size());
                for (int i = 0; i < vps.size(); i++){
                    hvps[i] = view.camera.toScreenInHPoint(vps[i]);
                }

                SaveToDisk(folder + "\\feature\\" + filename + "_view_lines_vps_focal",
                    view, lines, vps, focal);
                {
                    gui::AsCanvas(image).colorTable(gui::ColorTableDescriptor::RGBGreys).thickness(2.0).add(lines).show();
                }


                qDebug() << "computing omap";
                Imagei om = ComputeOrientationMaps(lines, hvps, view.image.size());
                SaveToDisk(folder + "\\feature\\" + filename + "_omap", om);
                {
                    gui::ColorTable rgb = gui::ColorTableDescriptor::RGBGreys;
                    rgb.exceptoinalColor() = gui::Transparent;
                    cv::imwrite(folder + "\\snapshot\\" + filename + "_omap.png", rgb(om));
                }

                qDebug() << "computing segmentations";
                SegmentationExtractor segmenter;
                segmenter.params().algorithm = SegmentationExtractor::GraphCut;
                segmenter.params().sigma = 10.0;
                segmenter.params().c = 1.0;
                segmenter.params().superpixelSizeSuggestion = 2000;
                int segmentsNum = 0;
                Imagei segmentedImage;
                std::tie(segmentedImage, segmentsNum) = segmenter(view.image, false);
                SaveToDisk(folder + "\\feature\\" + filename + "_segs_segsnum", segmentedImage, segmentsNum);
                {
                    gui::ColorTable rand = gui::CreateRandomColorTableWithSize(segmentsNum);
                    cv::imwrite(folder + "\\snapshot\\" + filename + "_segs.png", rand(segmentedImage));
                }
                {
                    gui::ColorTable rand = gui::CreateGreyColorTableWithSize(segmentsNum);
                    rand.randomize();
                    auto segim = rand(segmentedImage);
                    gui::AsCanvas(segim).colorTable(gui::ColorTableDescriptor::RGBGreys).thickness(2.0).add(lines).show();
                }

                qDebug() << "computing geometric contexts";
                auto gc = ComputeGeometricContext(image, indoor, false);
                SaveToDisk(folder + "\\feature\\" + filename + "_gc7", gc);
                {
                    Image3f gcim(gc.size(), 0.0);
                    for (auto it = gcim.begin(); it != gcim.end(); ++it){
                        auto & v = gc(it.pos());
                        *it =
                            (v[0]) * Vec3f(1, 0, 0) +
                            (v[1] + v[3]) * Vec3f(0, 1, 0) +
                            (v[2]) * Vec3f(0, 0, 1) +
                            v[4] * normalize(Vec3f(1, 1, 0)) +
                            v[5] * normalize(Vec3f(1, 0, 1)) +
                            v[6] * normalize(Vec3f(0, 1, 1));
                        *it /= norm(*it);
                    }
                    cv::imwrite(folder + "\\snapshot\\" + filename + "_gc7.png", gcim * 255);
                }

            }
        }

        void RunPano(std::string folder){

            QString qfolder = QString::fromStdString(folder);

            for (auto & finfo : QDir(qfolder + "\\image").entryInfoList(QStringList() << "*.jpg" << "*.png")){

                qDebug() << finfo.absoluteFilePath();

                std::string path = finfo.absoluteFilePath().toStdString();
                std::string filename = finfo.baseName().toStdString();

                Image3ub image = cv::imread(path);
                MakePanorama(image);
                ResizeToHeight(image, 700);
                
                View<PanoramicCamera, Image3ub> view;
                std::vector<PerspectiveCamera> cams;
                std::vector<std::vector<Classified<Line2>>> lines;
                std::vector<Vec3> vps;

                view = CreatePanoramicView(image);

                // collect lines in each view
                qDebug() << "computing lines and vps";
                cams = CreatePanoContextCameras(view.camera, image.rows, image.rows, image.rows * 0.4);
                lines.resize(cams.size());
                for (int i = 0; i < cams.size(); i++){
                    auto pim = view.sampled(cams[i]).image;
                    LineSegmentExtractor lineExtractor;
                    lineExtractor.params().algorithm = LineSegmentExtractor::LSD;
                    auto ls = lineExtractor(pim, 1, 300); // no pyramid
                    lines[i].reserve(ls.size());
                    for (auto & l : ls){
                        lines[i].push_back(ClassifyAs(l, -1));
                    }
                }
                vps = EstimateVanishingPointsAndClassifyLines(cams, lines);
                SaveToDisk(folder + "\\feature\\" + filename + "_view_cams_lines_vps",
                    view, cams, lines, vps);
                {
                    gui::ColorTable ctable = gui::ColorTableDescriptor::RGBGreys;
                    std::vector<View<PerspectiveCamera, Image3f>> lineviews(cams.size());
                    for (int i = 0; i < cams.size(); i++){
                        auto pim = view.sampled(cams[i]).image;
                        auto vis = gui::MakeCanvas(pim);
                        vis.colorTable(ctable).thickness(3.0).add(lines[i]);
                        lineviews[i].camera = cams[i];
                        lineviews[i].image = vis.image() / 255.0f;
                    }
                    auto pv = Combine(view.camera, lineviews);
                    cv::imwrite(folder + "\\snapshot\\" + filename + "_lines.png", pv.image * 255);
                }


                // estimate omap
                qDebug() << "computing omap";
                std::vector<Classified<Line3>> line3s;
                for (int i = 0; i < cams.size(); i++){
                    for (auto & line : lines[i]){
                        line3s.push_back(ClassifyAs(
                            Line3(cams[i].toSpace(line.component.first), 
                            cams[i].toSpace(line.component.second)), line.claz));
                    }
                }
                std::vector<View<PerspectiveCamera, Image3ub>> omaps;
                for (int i = 0; i < cams.size(); i++){
                    std::vector<HPoint2> hvps(vps.size());
                    for (int k = 0; k < vps.size(); k++){
                        hvps[k] = cams[i].toScreenInHPoint(vps[k]);
                    }
                    Imagei omap = ComputeOrientationMaps(lines[i], hvps, cams[i].screenSize());
                    std::vector<Imageub> omapChannels = { omap == 0, omap == 1, omap == 2 };
                    Image3ub merged;
                    cv::merge(omapChannels, merged);
                    omaps.emplace_back(merged, cams[i]);
                }
                auto panoOmap = Combine(view.camera, omaps);
                for (auto & p : panoOmap.image){
                    if ((p.val[0] > 0) + (p.val[1] > 0) + (p.val[2] > 0) > 1){
                        std::fill(p.val, std::end(p.val), 0);
                    }
                }             
                SaveToDisk(folder + "\\feature\\" + filename + "_omap", panoOmap.image);
                {
                    gui::ColorTable rgb = gui::ColorTableDescriptor::RGBGreys;
                    rgb.exceptoinalColor() = gui::Transparent;
                    cv::imwrite(folder + "\\snapshot\\" + filename + "_omap.png", panoOmap.image * 255);
                }

                // computing segmentations
                qDebug() << "computing segmentations";
                SegmentationExtractor segmenter;
                segmenter.params().algorithm = SegmentationExtractor::GraphCut;
                segmenter.params().sigma = 10.0;
                segmenter.params().c = 1.0;
                segmenter.params().superpixelSizeSuggestion = 2000;
                int segmentsNum = 0;
                Imagei segmentedImage;
                std::tie(segmentedImage, segmentsNum) = segmenter(view.image, true);
                SaveToDisk(folder + "\\feature\\" + filename + "_segs_segsnum", segmentedImage, segmentsNum);
                {
                    gui::ColorTable rand = gui::CreateRandomColorTableWithSize(segmentsNum);
                    cv::imwrite(folder + "\\snapshot\\" + filename + "_segs.png", rand(segmentedImage));
                    gui::ColorTable greys = gui::CreateGreyColorTableWithSize(segmentsNum);
                    greys.randomize();
                    gui::ColorTable ctable = gui::ColorTableDescriptor::RGBGreys;
                    std::vector<View<PerspectiveCamera, Image3f>> lineviews(cams.size());
                    for (int i = 0; i < cams.size(); i++){
                        auto pim = MakeCameraSampler(cams[i], view.camera)(greys(segmentedImage));
                        auto vis = gui::MakeCanvas(pim);
                        vis.thickness(3).colorTable(ctable).add(lines[i]);
                        lineviews[i].camera = cams[i];
                        lineviews[i].image = vis.image() / 255.0f;
                    }
                    auto pv = Combine(view.camera, lineviews);
                    cv::imwrite(folder + "\\snapshot\\" + filename + "_segs_lines.png", pv.image * 255);
                }

                // compute gc
                //qDebug() << "computing geometric contexts";
                //std::vector<View<PerspectiveCamera, Image7d>> omaps;
                //for (int i = 0; i < cams.size(); i++){
                //    std::vector<HPoint2> hvps(vps.size());
                //    for (int k = 0; k < vps.size(); k++){
                //        hvps[k] = cams[i].toScreenInHPoint(vps[k]);
                //    }
                //    auto pview = view.sampled(cams[i]);
                //    // 0: front, 1: left, 2: right, 3: floor, 4: ceiling, 5: clutter, 6: unknown
                //    Image7d gc = ComputeGeometricContext(pview.image, SceneClass::Indoor, false);
                //    
                //}
                //auto panoOmap = Combine(view.camera, omaps);

               
                //SaveToDisk(folder + "\\feature\\" + filename + "_gc7", gc);
                //{
                //    Image3f gcim(gc.size(), 0.0);
                //    for (auto it = gcim.begin(); it != gcim.end(); ++it){
                //        auto & v = gc(it.pos());
                //        *it =
                //            (v[0]) * Vec3f(1, 0, 0) +
                //            (v[1] + v[3]) * Vec3f(0, 1, 0) +
                //            (v[2]) * Vec3f(0, 0, 1) +
                //            v[4] * normalize(Vec3f(1, 1, 0)) +
                //            v[5] * normalize(Vec3f(1, 0, 1)) +
                //            v[6] * normalize(Vec3f(0, 1, 1));
                //        *it /= norm(*it);
                //    }
                //    cv::imwrite(folder + "\\snapshot\\" + filename + "_gc7.png", gcim * 255);
                //}
            }
        }


        void Run() {

            //RunNormal("F:\\DataSets\\yanghao.panoramix\\normal\\outdoor", false);
            RunPano("F:\\DataSets\\yanghao.panoramix\\pano\\indoor");

        }



    }
}