#include "../../src/core/basic_types.hpp"
#include "../../src/core/mixed_graph.hpp"
#include "../../src/vis/visualize2d.hpp"

#include "routines.hpp"

namespace panolyz {

    using namespace panoramix;

    ROUTINE_FOR_ALGORITHM(Panorama_v1){

        core::View<core::PanoramicCamera> view;

        std::vector<core::PerspectiveCamera> cams;
        std::vector<std::vector<core::Classified<core::Line2>>> lines;
        std::vector<core::Vec3> vps;

        std::vector<core::PerspectiveCamera> hCams; // horizontal cams
        std::vector<core::GeometricContextEstimator::Feature> gcs;

        core::Imagei segmentedImage;

        core::MixedGraph mg;
        core::MixedGraphPropertyTable props;

        if (1){
            view = core::CreatePanoramicView(image);

            // collect lines in each view
            cams = core::CreateCubicFacedCameras(view.camera, image.rows * 0.7, image.rows * 0.7, image.rows * 0.3);
            lines.resize(cams.size());
            for (int i = 0; i < cams.size(); i++){
                auto pim = view.sampled(cams[i]).image;
                core::LineSegmentExtractor lineExtractor;
                lineExtractor.params().algorithm = core::LineSegmentExtractor::LSD;
                auto ls = lineExtractor(pim, 2, 300); // use pyramid
                lines[i].reserve(ls.size());
                for (auto & l : ls){
                    lines[i].push_back(core::ClassifyAs(l, -1));
                }
                vis::Visualizer2D(pim) << ls << vis::manip2d::Show();
            }

            // estimate vp
            vps = core::EstimateVanishingPointsAndClassifyLines(cams, lines);
            // extract lines from segmentated region boundaries and classify them using estimated vps
            // make 3d lines
            std::vector<core::Line3> line3ds;
            for (int i = 0; i < cams.size(); i++){
                for (auto & l : lines[i]){
                    line3ds.emplace_back(normalize(cams[i].spatialDirection(l.component.first)),
                        normalize(cams[i].spatialDirection(l.component.second)));
                }
            }
            core::SegmentationExtractor segmenter;
            segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
            segmenter.params().c = 100.0;
            int segmentsNum = 0;
            std::tie(segmentedImage, segmentsNum) = segmenter(view.image, line3ds, view.camera, M_PI / 36.0);
            //auto colorTable = vis::CreateRandomColorTableWithSize(segmentsNum);
            ////vis::Visualizer2D(segmentedImage) << vis::manip2d::Show();
            ////int segmentsNum = core::MinMaxValOfImage(segmentedImage).second + 1;
            //for (int i = 0; i < cams.size(); i++){
            //    auto pim = core::MakeCameraSampler(cams[i], view.camera)(segmentedImage);
            //    core::LineSegmentExtractor lineExtractor;
            //    lineExtractor.params().algorithm = core::LineSegmentExtractor::LSD;
            //    lineExtractor.params().minLength = 50;
            //    lineExtractor.params().xBorderWidth = lineExtractor.params().yBorderWidth = 10;
            //    auto ls = lineExtractor(colorTable(pim));
            //    
            //    vis::Visualizer2D(pim) << ls << vis::manip2d::Show();
            //}


            // collect gcs in each horizontal view
            /*hCams = core::CreateHorizontalPerspectiveCameras(view.camera, vps, 500, 500, 200);
            gcs.resize(hCams.size());
            for (int i = 0; i < hCams.size(); i++){
            auto pim = view.sampled(hCams[i]).image;
            core::GeometricContextEstimator gce;
            gcs[i] = gce(pim, isOutdoor ? core::SceneClass::Outdoor : core::SceneClass::Indoor);
            std::vector<core::Image> channels = {
            cv::max(gcs[i][core::GeometricContextLabel::Left], gcs[i][core::GeometricContextLabel::Right]),
            gcs[i][core::GeometricContextLabel::Front],
            cv::max(gcs[i][core::GeometricContextLabel::Ceiling],  gcs[i][core::GeometricContextLabel::Floor])
            };
            core::Image gc = core::Imaged3::zeros(channels.front().size());
            cv::mixChannels(channels, gc, {0, 0, 1, 1, 2, 2});
            vis::Visualizer2D(gc) << vis::manip2d::Show();
            }*/


            // append lines
            for (int i = 0; i < cams.size(); i++){
                core::AppendLines(mg, lines[i], cams[i], vps);
            }

            // append regions
            core::AppendRegions(mg, segmentedImage, view.camera, 0.01, 0.02, 3, 2);

            core::SaveToDisk("./cache/all", view, cams, lines, vps, hCams, gcs, segmentedImage, mg, props);
        }
        else{
            core::LoadFromDisk("./cache/all", view, cams, lines, vps, hCams, gcs, segmentedImage, mg, props);
        }

        // optimize
        if (1){
            props = core::MakeMixedGraphPropertyTable(mg, vps);
            core::AttachPrincipleDirectionConstraints(mg, props, M_PI / 15.0);
            core::AttachWallConstriants(mg, props, M_PI / 30.0);
            //core::AttachGeometricContextConstraints(mg, props, gcs, hCams, 2);

            for (int i = 0; i < 10; i++){

                // visualize using segmented image
                {
                    std::vector<vis::Color> colors(mg.internalComponents<core::RegionData>().size());
                    std::vector<vis::Color> oColors = { vis::ColorTag::Red, vis::ColorTag::Green, vis::ColorTag::Blue };
                    std::vector<vis::Color> noColors = { vis::ColorTag::DimGray, vis::ColorTag::Gray, vis::ColorTag::DarkGray };
                    for (auto & r : mg.components<core::RegionData>()){
                        auto & p = props[r.topo.hd];
                        auto & color = colors[r.topo.hd.id];
                        if (!p.used){
                            color = vis::ColorTag::Black;
                        }
                        else{
                            if (p.orientationClaz != -1){
                                color = oColors[p.orientationClaz];
                            }
                            else if (p.orientationNotClaz != -1){
                                color = noColors[p.orientationNotClaz];
                            }
                            else{
                                color = vis::ColorTag::White;
                            }
                        }
                    }
                    vis::Visualizer2D::Params params;
                    params.colorTable = vis::ColorTable(colors, vis::ColorTag::White);
                    vis::Visualizer2D(segmentedImage)
                        << vis::manip2d::Show();
                    vis::Visualizer2D(segmentedImage, params)
                        << vis::manip2d::Show();
                }

                core::SolveVariablesUsingInversedDepths(mg, props);
                core::NormalizeVariables(mg, props);
                std::cout << "score = " << core::ComputeScore(mg, props) << std::endl;
                core::Visualize(view, mg, props);
                core::LooseOrientationConstraintsOnComponents(mg, props, 0.2, 0.02, 0.1);

                core::SolveVariablesUsingInversedDepths(mg, props);
                core::NormalizeVariables(mg, props);
                std::cout << "score = " << core::ComputeScore(mg, props) << std::endl;
                core::Visualize(view, mg, props);
                core::AttachFloorAndCeilingConstraints(mg, props, 0.1, 0.6);
            }

            core::SaveToDisk("./cache/mgp", mg, props);
        }
        else{
            core::LoadFromDisk("./cache/mgp", mg, props);
        }
    }


    ROUTINE_FOR_ALGORITHM(Normal_v1){

        core::View<core::PerspectiveCamera> view;
        std::vector<core::Classified<core::Line2>> lines;
        std::vector<core::Vec3> vps;

        double focal;
        core::MixedGraph mg;
        
        core::Imagei segmentedImage;
        core::MixedGraphPropertyTable props;

        core::VanishingPointsDetector::Params vpdParams(core::VanishingPointsDetector::TardifSimplified, image.size());
        view = core::CreatePerspectiveView(image, core::Point3(0, 0, 0), core::Point3(1, 0, 0), core::Point3(0, 0, -1),
            core::LineSegmentExtractor(), core::VanishingPointsDetector(vpdParams), nullptr, &lines, &vps, &focal);

        {
            std::vector<core::Classified<core::InfiniteLine2>> vpRays;
            for (int i = 0; i < 3; i++){
                std::cout << "vp[" << i << "] = " << view.camera.screenProjection(vps[i]) << std::endl;
                for (double a = 0; a <= M_PI * 2.0; a += 0.1){
                    core::Point2 p = core::Point2(image.cols / 2, image.rows / 2) + core::Vec2(cos(a), sin(a)) * 1000.0;
                    vpRays.push_back(core::ClassifyAs(core::InfiniteLine2(p, (view.camera.screenProjectionInHPoint(vps[i]) - core::HPoint2(p, 1.0)).numerator), i));
                }
            }
            vis::Visualizer2D(image)
                << vis::manip2d::SetColorTable(vis::ColorTable(vis::ColorTableDescriptor::RGB).appendRandomizedGreyColors(vps.size() - 3))
                << vis::manip2d::SetThickness(1)
                << vpRays
                << vis::manip2d::SetThickness(2)
                << lines
                << vis::manip2d::Show(0);
        }


        // append lines
        core::AppendLines(mg, lines, view.camera, vps);

        // append regions
        core::SegmentationExtractor segmenter;
        segmenter.params().algorithm = core::SegmentationExtractor::GraphCut;
        segmenter.params().c = 100.0;
        std::vector<core::Line2> pureLines(lines.size());
        for (int i = 0; i < lines.size(); i++)
            pureLines[i] = lines[i].component;
        int segmentsNum = 0;
        std::tie(segmentedImage, segmentsNum) = segmenter(image, pureLines, image.cols / 100.0);

        vis::Visualizer2D::Params vParams;
        vParams.colorTable = vis::CreateRandomColorTableWithSize(segmentsNum);
        vis::Visualizer2D(segmentedImage, vParams)
            << vis::manip2d::Show();

        core::AppendRegions(mg, segmentedImage, view.camera, 0.001, 0.001, 3, 1);

        // optimize
        props = core::MakeMixedGraphPropertyTable(mg, vps);
        core::AttachPrincipleDirectionConstraints(mg, props, M_PI / 15.0);
        core::AttachWallConstriants(mg, props, M_PI / 30.0);

        for (int i = 0; i < 10; i++){

            // visualize using segmented image
            {
                std::vector<vis::Color> colors(mg.internalComponents<core::RegionData>().size());
                std::vector<vis::Color> oColors = { vis::ColorTag::Red, vis::ColorTag::Green, vis::ColorTag::Blue };
                std::vector<vis::Color> noColors = { vis::ColorTag::DimGray, vis::ColorTag::Gray, vis::ColorTag::DarkGray };
                for (auto & r : mg.components<core::RegionData>()){
                    auto & p = props[r.topo.hd];
                    auto & color = colors[r.topo.hd.id];
                    if (!p.used){
                        color = vis::ColorTag::Black;
                    }
                    else{
                        if (p.orientationClaz != -1){
                            color = oColors[p.orientationClaz];
                        }
                        else if (p.orientationNotClaz != -1){
                            color = noColors[p.orientationNotClaz];
                        }
                        else{
                            color = vis::ColorTag::White;
                        }
                    }
                }
                vis::Visualizer2D::Params params;
                params.colorTable = vis::ColorTable(colors, vis::ColorTag::White);
                vis::Visualizer2D(segmentedImage)
                    << vis::manip2d::Show();
                vis::Visualizer2D(segmentedImage, params)
                    << vis::manip2d::Show();
            }

            core::SolveVariablesUsingInversedDepths(mg, props);
            core::NormalizeVariables(mg, props);
            std::cout << "score = " << core::ComputeScore(mg, props) << std::endl;
            core::Visualize(view, mg, props);
            core::LooseOrientationConstraintsOnComponents(mg, props, 0.2, 0, 0.05);

            core::SolveVariablesUsingInversedDepths(mg, props);
            core::NormalizeVariables(mg, props);
            std::cout << "score = " << core::ComputeScore(mg, props) << std::endl;
            core::Visualize(view, mg, props);
            core::AttachFloorAndCeilingConstraints(mg, props, 0.1, 0.6);
        }

    }

}