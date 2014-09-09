#include "../core/debug.hpp"
#include "../vis/visualize2d.hpp"

#include "lines_net.hpp"

namespace panoramix {
    namespace rec {

        LinesNet::Params::Params() :
            intersectionDistanceThreshold(8),
            incidenceDistanceAlongDirectionThreshold(15),
            incidenceDistanceVerticalDirectionThreshold(3) {
        }

        LinesNet::LinesNet(const Image & image, const Params & params)
            : _image(image), _params(params), _lineSegments(params.lineSegmentExtractor(image)) {
            _lineSegmentIntersections = ComputeLineIntersections(_lineSegments, &_lineSegmentIntersectionIds, true);
        }

        namespace {

            template <class T, int N>
            inline Line<T, N> LineRelativeToPoint(const Line<T, N> & line, const Point<T, N> & p){
                return Line < T, N > {line.first - p, line.second - p};
            }
            inline Vec3 ComputeLineCoeffsRelativeToPoint(const Line2 & line, const Point2 & imCenter, double fakeFocal){
                return normalize(VectorFromHPoint(HPoint2(line.first - imCenter, fakeFocal))
                    .cross(VectorFromHPoint(HPoint2(line.second - imCenter, fakeFocal))));
            }

        }

        void LinesNet::buildNetAndComputeFeaturesUsingVanishingPoints (
            const std::array<HPoint2, 3> & vps, const std::vector<int> & lineSegmentClasses) {

            const auto & lines = _lineSegments;
            if (!lineSegmentClasses.empty())
                assert(lineSegmentClasses.size() == lines.size());

            // build lines net and classify lines
            std::vector<LineHandle> handles;
            handles.reserve(lines.size());
            _lines.clear();
            _lines.internalElements<0>().reserve(lines.size());

            static const double angleThreshold = M_PI / 32;
            static const double sigma = 0.1;

            auto lineSegmentClassesBegin = lineSegmentClasses.begin();
            for (auto & line : lines){

                // classify lines
                LineData ld;
                ld.line.component = line;
                ld.line.claz = -1;

                if (!lineSegmentClasses.empty()){
                    ld.line.claz = *lineSegmentClassesBegin;
                    ++ lineSegmentClassesBegin;
                }
                else{

                    // classify
                    std::vector<double> lineangles(vps.size());
                    std::vector<double> linescores(vps.size());

                    for (int j = 0; j < vps.size(); j++){
                        auto & point = vps[j];
                        double angle = std::min(
                            AngleBetweenDirections(line.direction(), (point - HPoint2(line.center())).numerator),
                            AngleBetweenDirections(-line.direction(), (point - HPoint2(line.center())).numerator));
                        lineangles[j] = angle;
                    }

                    // get score based on angle
                    for (int j = 0; j < vps.size(); j++){
                        double angle = lineangles[j];
                        double score = exp(-(angle / angleThreshold) * (angle / angleThreshold) / sigma / sigma / 2);
                        linescores[j] = (angle > angleThreshold) ? 0 : score;
                    }


                    double curscore = 0.8;
                    for (int j = 0; j < vps.size(); j++){
                        if (linescores[j] > curscore){
                            ld.line.claz = j;
                            curscore = linescores[j];
                        }
                    }
                }

                if (ld.line.claz == -1) {// ignore cluttered lines
                    handles.push_back(LineHandle());
                }
                else{
                    handles.push_back(_lines.add(ld));
                }
            }


            // construct incidence/intersection relations
            _lines.internalElements<1>().reserve(lines.size() * (lines.size()-1) / 2);
            for (int i = 0; i < lines.size(); i++){
                auto & linei = lines[i];
                if (handles[i].isInvalid())
                    continue;
                int clazi = _lines.data(handles[i]).line.claz;
                assert(clazi != -1);

                for (int j = i + 1; j < lines.size(); j++){                    
                    auto & linej = lines[j];
                    if (handles[j].isInvalid())
                        continue;
                    int clazj = _lines.data(handles[j]).line.claz;
                    assert(clazj != -1);
                    
                    auto nearest = DistanceBetweenTwoLines(linei, linej);
                    double d = nearest.first;                    
           
                    if (clazi == clazj){     
                        auto conCenter = (nearest.second.first.position + nearest.second.second.position) / 2.0;
                        auto conDir = (nearest.second.first.position - nearest.second.second.position);
                        auto & vp = vps[clazi];

                        if (Distance(vp.value(), conCenter) < _params.intersectionDistanceThreshold)
                            continue;

                        auto dir = normalize((vp - HPoint2(conCenter)).numerator);
                        double dAlong = abs(conDir.dot(dir));
                        double dVert = sqrt(Square(norm(conDir)) - dAlong*dAlong);

                        if (dAlong < _params.incidenceDistanceAlongDirectionThreshold &&
                            dVert < _params.incidenceDistanceVerticalDirectionThreshold){
                            LineRelationData lrd;
                            lrd.type = LineRelationData::Type::Incidence;
                            lrd.relationCenter = conCenter;
                            _lines.add<1>({ handles[i], handles[j] }, lrd);
                        }
                    }
                    else {
                        if (d < _params.intersectionDistanceThreshold){
                            auto conCenter = HPointFromVector(GetLine2Coeffs(linei.infinieLine())
                                .cross(GetLine2Coeffs(linej.infinieLine()))).value();
                            LineRelationData lrd;
                            lrd.type = LineRelationData::Type::Intersection;
                            lrd.relationCenter = conCenter;
                            _lines.add<1>({ handles[i], handles[j] }, lrd);
                        }
                    }
                }
            }

            // compute voting distributions
            _lineVotingDistribution = ImageWithType<Mat<float, 3, 2>>::zeros(_image.size());
            
            for (auto it = _lineVotingDistribution.begin();
                it != _lineVotingDistribution.end(); 
                ++it){
                Point2 pos(it.pos().x, it.pos().y);
                for (auto & ld : _lines.elements<0>()){
                    auto & line = ld.data.line.component;
                    int claz = ld.data.line.claz;
                    auto & vp = vps[claz];
                    Point2 center = line.center();

                    Vec2 center2vp = vp.value() - center;
                    Vec2 center2pos = pos - center;

                    if (norm(center2pos) <= 1)
                        continue;
                    
                    double angle = AngleBetweenDirections(center2pos, center2vp);
                    double angleSmall = angle > M_PI_2 ? (M_PI - angle) : angle;
                    assert(angleSmall >= 0 && angleSmall <= M_PI_2);

                    double angleScore = 
                        exp(-(angleSmall / angleThreshold) * (angleSmall / angleThreshold) / sigma / sigma / 2);
                    
                    auto proj = ProjectionOfPointOnLine(pos, line);
                    double projRatio = BoundBetween(proj.ratio, 0.0, 1.0);

                    auto & votingData = *it;
                    if (AngleBetweenDirections(center2vp, line.direction()) < M_PI_2){ // first-second-vp
                        votingData(claz, TowardsVanishingPoint) += angleScore * line.length() * (1 - projRatio);
                        votingData(claz, TowardsOppositeOfVanishingPoint) += angleScore * line.length() * projRatio;
                    }
                    else{ // vp-first-second
                        votingData(claz, TowardsOppositeOfVanishingPoint) += angleScore * line.length() * (1 - projRatio);
                        votingData(claz, TowardsVanishingPoint) += angleScore * line.length() * projRatio;
                    }
                }
            }


            IF_DEBUG_USING_VISUALIZERS {

                float distributionMaxVal = std::numeric_limits<float>::lowest();
                cv::Mat_<Vec<float, 3>> distributionImages[2];
                for (int i = 0; i < 2; i++){
                    auto & dim = distributionImages[i];
                    dim.create(_image.size());
                    for (int x = 0; x < _image.cols; x++){
                        for (int y = 0; y < _image.rows; y++){
                            auto & pixel = _lineVotingDistribution(y, x);
                            dim(y, x) = Vec<float, 3>(pixel(0, i), pixel(1, i), pixel(2, i));
                            distributionMaxVal = 
                                std::max({ distributionMaxVal, pixel(0, i), pixel(1, i), pixel(2, i) });
                        }
                    }
                }

                for (auto & dim : distributionImages){
                    dim /= distributionMaxVal;
                }

                //vis::Visualizer2D(distributionImages[0]) << vis::manip2d::Show();
                //vis::Visualizer2D(distributionImages[1]) << vis::manip2d::Show();
                //cv::imshow("distribution 0", distributionImages[0]);
                //cv::imshow("distribution 1", distributionImages[1]);

                vis::Visualizer2D viz(_image);
                viz.params.thickness = 2;
                viz.params.colorTableDescriptor = vis::ColorTableDescriptor::RGB;
                for (auto & ld : _lines.elements<0>()){
                    viz = viz << ld.data.line;
                }
                viz.params.thickness = 1;
                viz = viz << vis::manip2d::SetColor(vis::ColorTag::Red);
                for (auto & rd : _lines.elements<1>()){
                    auto & l1 = _lines.data(rd.topo.lowers[0]).line;
                    auto & l2 = _lines.data(rd.topo.lowers[1]).line;
                    if (rd.data.type == LineRelationData::Type::Incidence){
                        auto nearest = DistanceBetweenTwoLines(l1.component, l2.component);
                        Line2 conLine = { nearest.second.first.position, nearest.second.second.position };
                        viz = viz << conLine;
                    }
                    else if (rd.data.type == LineRelationData::Type::Intersection){
                        auto nearest = DistanceBetweenTwoLines(l1.component, l2.component);
                        auto & relationCenter = rd.data.relationCenter;
                        viz = viz
                            << Line2(nearest.second.first.position, relationCenter)
                            << Line2(relationCenter, nearest.second.second.position);
                    }
                }

                viz << vis::manip2d::Show();
                /*cv::imshow("lines", viz.image());
                cv::waitKey();*/
            }
            

        }
    
    }
}