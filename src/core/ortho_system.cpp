
#include "utilities.hpp"
#include "ortho_system.hpp"

namespace panoramix {
    namespace core {
      
        std::pair<SemiOrthogonalSystem::DirectionIndex, double> SemiOrthogonalSystem::nearestDirectionId(const Vec3 & dir) const{
            double angle = std::numeric_limits<double>::max();
            DirectionIndex di = -1;
            for (DirectionIndex i = 0; i < _directions.size(); i++){
                auto a = core::AngleBetweenUndirectedVectors(_directions[i], dir);
                if (a < angle){
                    di = i;
                    angle = a;
                }
            }
            return std::make_pair(di, angle);
        }

        std::pair<SemiOrthogonalSystem::DirectionIndex, double> SemiOrthogonalSystem::nearestViceDirectionId(const Vec3 & dir) const{
            double angle = std::numeric_limits<double>::max();
            DirectionIndex di = -1;
            for (DirectionIndex i = 1; i < _directions.size(); i++){
                auto a = core::AngleBetweenUndirectedVectors(_directions[i], dir);
                if (a < angle){
                    di = i;
                    angle = a;
                }
            }
            return std::make_pair(di, angle);
        }

        std::pair<SemiOrthogonalSystem::DirectionIndex, SemiOrthogonalSystem::DirectionIndex> 
            SemiOrthogonalSystem::insertOrthoPair(const Vec3 & a, const Vec3 b, double angleErrorThreshold){
            assert(abs(core::AngleBetweenUndirectedVectors(a, b) - M_PI_2) <= angleErrorThreshold);
            assert(abs(core::AngleBetweenUndirectedVectors(a, archDirection()) - M_PI_2) <= angleErrorThreshold);
            assert(abs(core::AngleBetweenUndirectedVectors(b, archDirection()) - M_PI_2) <= angleErrorThreshold);
            
            DirectionIndex dia;
            double anglea;
            std::tie(dia, anglea) = nearestViceDirectionId(a);
            if (anglea <= angleErrorThreshold){
                assert(dia != 0);
                if (Contains(_viceDirectionIdPairs, dia)){
                    DirectionIndex diaOppo = _viceDirectionIdPairs.at(dia);
                    if (AngleBetweenDirections(b, _directions[diaOppo]) <= angleErrorThreshold){
                        return std::make_pair(dia, diaOppo);
                    }
                    else{
                        _directions.push_back(a);
                        _directions.push_back(b);
                        DirectionIndex ia = _directions.size() - 2;
                        DirectionIndex ib = _directions.size() - 1;
                        _viceDirectionIdPairs.emplace(ia, ib);
                        _viceDirectionIdPairs.emplace(ib, ia);
                        return std::make_pair(ia, ib);
                    }
                }
                else{
                    _directions.push_back(b);
                    DirectionIndex ib = _directions.size() - 1;
                    _viceDirectionIdPairs.emplace(dia, ib);
                    _viceDirectionIdPairs.emplace(ib, dia);
                    return std::make_pair(dia, ib);
                }
            }
            
            DirectionIndex dib;
            double angleb;
            std::tie(dib, angleb) = nearestViceDirectionId(b);
            if (angleb <= angleErrorThreshold){
                assert(dib != 0);
                if (Contains(_viceDirectionIdPairs, dib)){
                    DirectionIndex dibOppo = _viceDirectionIdPairs.at(dib);
                    if (AngleBetweenDirections(a, _directions[dibOppo]) <= angleErrorThreshold){
                        return std::make_pair(dibOppo, dib);
                    }
                    else{
                        _directions.push_back(a);
                        _directions.push_back(b);
                        DirectionIndex ia = _directions.size() - 2;
                        DirectionIndex ib = _directions.size() - 1;
                        _viceDirectionIdPairs.emplace(ia, ib);
                        _viceDirectionIdPairs.emplace(ib, ia);
                        return std::make_pair(ia, ib);
                    }
                }
                else{
                    _directions.push_back(a);
                    DirectionIndex ia = _directions.size() - 1;
                    _viceDirectionIdPairs.emplace(ia, dib);
                    _viceDirectionIdPairs.emplace(dib, ia);
                    return std::make_pair(ia, dib);
                }
            }

            _directions.push_back(a);
            _directions.push_back(b);
            DirectionIndex ia = _directions.size() - 2;
            DirectionIndex ib = _directions.size() - 1;
            _viceDirectionIdPairs.emplace(ia, ib);
            _viceDirectionIdPairs.emplace(ib, ia);

        }

        
        SemiOrthogonalSystem ConstructSemiOrthogonalSystemAndClassifyLines(std::vector<Classified<Line3>> & lines,
            const Vec3 & archDirectionRangeCenter,
            double archDirectionRangeRadius) {

            static const double dotThresholdForArchDir = 0.006;
            static const double angleThresholdForNeighborLinePairJudge = M_PI / 300.0;

            static const int longitudeDivideNum = 1000, latitudeDivideNum = 500;

            std::vector<Vec3> lineNormals(lines.size());
            for (int i = 0; i < lines.size(); i++){
                lineNormals[i] = lines[i].component.first.cross(lines[i].component.second);
            }
            std::vector<Vec3> intersections;
            intersections.reserve(lines.size() * lines.size() / 2);
            for (int i = 0; i < lines.size(); i++){
                auto & lineNormalI = lineNormals[i];
                for (int j = i + 1; j < lines.size(); j++){
                    auto & lineNormalJ = lineNormals[j];
                    auto intersection = normalize(lineNormalI.cross(lineNormalJ));
                    if (AngleBetweenDirections(archDirectionRangeCenter, intersection) <= archDirectionRangeRadius)
                        intersections.push_back(intersection);
                }
            }
            if (intersections.empty()){
                return SemiOrthogonalSystem();
            }

            // collect votes of intersection directions
            Imagef votePanel = Imagef::zeros(longitudeDivideNum, latitudeDivideNum);
            size_t pn = intersections.size();
            for (const Vec3 & p : intersections){
                PixelLoc pixel = PixelLocFromGeoCoord(GeoCoord(p), longitudeDivideNum, latitudeDivideNum);
                votePanel(pixel.x, pixel.y) += 1.0;
            }
            cv::GaussianBlur(votePanel, votePanel, cv::Size((longitudeDivideNum / 50) * 2 + 1, (latitudeDivideNum / 50) * 2 + 1),
                4, 4, cv::BORDER_REPLICATE);

            // set the direction with the max votes as the first vanishing point
            double minVal = 0, maxVal = 0;
            int maxIndex[] = { -1, -1 };
            cv::minMaxIdx(votePanel, &minVal, &maxVal, 0, maxIndex);
            cv::Point maxPixel(maxIndex[0], maxIndex[1]);

            SemiOrthogonalSystem sos(normalize(GeoCoordFromPixelLoc(maxPixel, longitudeDivideNum, latitudeDivideNum).toVector()));
            const Vec3 & vec0 = sos.archDirection();

            for (int i = 0; i < lines.size(); i++){
                if (abs(lineNormals[i].dot(vec0)) < dotThresholdForArchDir)
                    lines[i].claz = 0;
                else
                    lines[i].claz = -1;
            }


            std::vector<std::pair<int, int>> pairs;
            for (int i = 0; i < lines.size(); i++){
                if (lines[i].claz == 0)
                    continue;
                if (abs(lineNormals[i].dot(vec0)) < 0.01)
                    continue;
                for (int j = i + 1; j < lines.size(); j++){
                    if (lines[j].claz == 0)
                        continue;
                    if (abs(lineNormals[j].dot(vec0)) < 0.01)
                        continue;
                    double angle = DistanceBetweenTwoLines(normalize(lines[i].component), normalize(lines[j].component)).first;
                    auto & n1 = lineNormals[i];
                    auto & n2 = lineNormals[j];
                    auto inter = n1.cross(n2);
                    if (angle < angleThresholdForNeighborLinePairJudge &&
                        DistanceFromPointToLine(inter, lines[i].component).first < angleThresholdForNeighborLinePairJudge &&
                        DistanceFromPointToLine(inter, lines[j].component).first < angleThresholdForNeighborLinePairJudge){
                        pairs.emplace_back(i, j);
                    }
                }
            }

            for (auto & p : pairs){
                auto & n1 = lineNormals[p.first];
                auto & n2 = lineNormals[p.second];
                auto p1 = normalize(n1.cross(vec0));
                auto p2 = normalize(n2.cross(vec0));
                if (abs(p1.dot(p2)) < 0.02){
                    sos.insertOrthoPair(p1, p2, 0.01);
                }
            }

            return sos;

        }



        void SemiOrthogonalSystem::classify(std::vector<Classified<Line3>> & lines) const {
            NOT_IMPLEMENTED_YET();
        }


    }
}