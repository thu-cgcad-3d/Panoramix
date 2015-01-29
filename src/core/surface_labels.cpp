#include "surface_labels.hpp"

namespace panoramix {
    namespace core { 


        SurfaceLabelNames MostLikelyLabelAtPosition(const SurfaceLabelDistribution & distribution, const PixelLoc & p) {
            double maxProb = 1e-8;
            int n = (int)SurfaceLabelNames::None;
            for (int i = 0; i < distribution.size(); i++){
                double prob = distribution[i](p);
                if (prob > maxProb){
                    n = i;
                    maxProb = prob;
                }
            }
            return (SurfaceLabelNames)(n);
        }

        SurfaceLabelNames MostLikelySurfaceLabelInRegion(const SurfaceLabelDistribution & distribution, const std::vector<PixelLoc> & contour) {
            auto size = distribution.front().size();
            cv::Rect roi(cv::boundingRect(contour));
            if (roi.x < 0) roi.x = 0;
            if (roi.y < 0) roi.y = 0;
            if (roi.x + roi.width >= size.width) roi.width = size.width - 1 - roi.x;
            if (roi.y + roi.height >= size.height) roi.height = size.height - 1 - roi.y;

            double maxProb = 1e-8;
            int n = (int)SurfaceLabelNames::None;
            for (int o = 0; o < distribution.size(); o++){
                cv::Mat cropProb(distribution[o], roi);
                cv::Mat mask(cv::Mat::zeros(cropProb.rows, cropProb.cols, CV_8UC1));

                typedef cv::vector<cv::Point> TContour;
                typedef cv::vector<TContour> TContours;

                cv::drawContours(mask, TContours(1, contour), -1, cv::Scalar(255), CV_FILLED, CV_AA, cv::noArray(), 1, -roi.tl());
                double prob = cv::mean(cropProb, mask).val[0];
                if (prob > maxProb){
                    n = o;
                    maxProb = prob;
                }
            }
            return (SurfaceLabelNames)(n);
        }



        namespace {

            inline SurfaceLabelNames ToSurfaceLabels(GeometricContextLabel label){
                switch (label){
                case GeometricContextLabel::Ceiling: return SurfaceLabelNames::Horizontal;
                case GeometricContextLabel::Floor: return SurfaceLabelNames::Horizontal;
                case GeometricContextLabel::Front:
                case GeometricContextLabel::Left:
                case GeometricContextLabel::Right:
                    return SurfaceLabelNames::Vertical;
                case GeometricContextLabel::Furniture: return SurfaceLabelNames::Planar;
                case GeometricContextLabel::Ground: return SurfaceLabelNames::Horizontal;
                case GeometricContextLabel::Sky: return SurfaceLabelNames::Void;
                case GeometricContextLabel::Vertical: return SurfaceLabelNames::Vertical;
                case GeometricContextLabel::NotPlanar: return SurfaceLabelNames::NonPlanar;
                default:
                    return SurfaceLabelNames::None;
                }
            }

            inline std::array<Imaged, (unsigned)SurfaceLabelNames::Count> ComputeSurfaceLabelDistributions(SizeI size, 
                const std::vector<GeometricContextEstimator::Feature> & panoGCs){
               
                std::array<Imaged, (unsigned)SurfaceLabelNames::Count> distributions;
                for (auto & d : distributions){
                    d = Imaged::zeros(size);
                }

                // suppress non max in panoGCs and get masks for each panoGC
                for (int x = 0; x < size.width; x++){
                    for (int y = 0; y < size.height; y++){
                        std::vector<GeometricContextLabel> bestLabels;
                        bestLabels.reserve(panoGCs.size());
                        for (int i = 0; i < panoGCs.size(); i++){
                            GeometricContextLabel bestLabel = GeometricContextLabel::None;
                            double maxScore = 1e-5;
                            for (auto & gcChannel : panoGCs[i]){
                                double score = gcChannel.second(y, x);
                                if (score > maxScore){
                                    bestLabel = gcChannel.first;
                                    maxScore = score;
                                }
                            }
                            if (bestLabel != GeometricContextLabel::None){
                                bestLabels.push_back(bestLabel);
                            }
                        }

                        if (bestLabels.empty()){
                            continue;
                        }
                        for (auto label : bestLabels){
                            distributions[(unsigned)ToSurfaceLabels(label)](y, x) += 1.0 / bestLabels.size();
                        }

                    }
                }

                /// FIXME: Fill the area of unknowns !!!!!!!!!!!!!!!!!!!!!!!!!!!
                // maybe we can use GraphCut!!!!!!!! MRF???
                NOT_IMPLEMENTED_YET();

                return distributions;
            }

        }



        SurfaceLabels<PanoramicCamera> ComputeSurfaceLabels(const View<PanoramicCamera> & view,
            const std::vector<View<PerspectiveCamera, GeometricContextEstimator::Feature>> & perspectiveGCs){
            std::vector<GeometricContextEstimator::Feature> gcs(perspectiveGCs.size());
            for (int i = 0; i < perspectiveGCs.size(); i++){
                auto gc = perspectiveGCs[i].image;
                auto sampleBack = MakeCameraSampler(view.camera, perspectiveGCs[i].camera);
                for (auto & gcc : gc){
                    gcc.second = sampleBack(gcc.second, cv::BORDER_CONSTANT, 0);
                }
                gcs[i] = gc;
            }
            return SurfaceLabels<PanoramicCamera>(view.camera, ComputeSurfaceLabelDistributions(view.image.size(), gcs));
        }

        SurfaceLabels<PanoramicCamera> ComputeSurfaceLabels(const View<PanoramicCamera> & view, SceneClass sceneClaz,
            const GeometricContextEstimator & gcEstimator,
            int num, int width, int height, double focal){
            auto cams = core::CreateHorizontalPerspectiveCameras(view.camera, num, width, height, focal);
            std::vector<View<PerspectiveCamera, GeometricContextEstimator::Feature>> perspectiveGCs(cams.size());
            for (int i = 0; i < cams.size(); i++){
                perspectiveGCs[i].image = gcEstimator(view.sampled(cams[i]).image, sceneClaz);
                perspectiveGCs[i].camera = cams[i];
            }
            return ComputeSurfaceLabels(view, perspectiveGCs);
        }

        SurfaceLabels<PerspectiveCamera> ComputeSurfaceLabels(const View<PerspectiveCamera> & view, SceneClass sceneClaz,
            const GeometricContextEstimator & gcEstimator){
            auto gc = gcEstimator(view.image, sceneClaz);
            return SurfaceLabels<PerspectiveCamera>(view.camera, ComputeSurfaceLabelDistributions(view.image.size(), {gc}));
        }
    	
    }
}
 