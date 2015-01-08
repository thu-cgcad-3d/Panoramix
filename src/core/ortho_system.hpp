#ifndef PANORAMIX_CORE_ORTHO_SYSTEM_HPP
#define PANORAMIX_CORE_ORTHO_SYSTEM_HPP

#include "feature.hpp"

namespace panoramix {
    namespace core {


        // semi-orthogonal system
        class SemiOrthogonalSystem {
        public:
            using DirectionIndex = int;
            explicit SemiOrthogonalSystem(const Vec3 & archD = Vec3(0, 0, 1))
                : _directions(1, archD), _archDirectionId(0) {}

            inline const Vec3 & archDirection() const { return _directions[_archDirectionId]; }
            inline Vec3 & archDirection() { return _directions[_archDirectionId]; }

            inline const Vec3 & operator[](DirectionIndex di) const { return _directions[di]; }
            inline Vec3 & operator[](DirectionIndex di) { return _directions[di]; }

            std::pair<DirectionIndex, double> nearestDirectionId(const Vec3 & dir) const;
            std::pair<DirectionIndex, double> nearestViceDirectionId(const Vec3 & dir) const;
            std::pair<DirectionIndex, DirectionIndex> insertOrthoPair(const Vec3 & a, const Vec3 b, 
                double angleErrorThreshold = M_PI_4 / 100.0);
            void classify(std::vector<Classified<Line3>> & lines) const;

            template <class Archiver>
            void serialize(Archiver & ar) {
                ar(_directions, _archDirectionId, _viceDirectionIdPairs);
            }
        private:
            std::vector<Vec3> _directions;
            int _archDirectionId;
            std::map<DirectionIndex, DirectionIndex> _viceDirectionIdPairs;
        };


        SemiOrthogonalSystem ConstructSemiOrthogonalSystemAndClassifyLines(std::vector<Classified<Line3>> & lines,
            const Vec3 & archDirectionRangeCenter = Vec3(0, 0, 1),
            double archDirectionRangeRadius = M_PI_4 / 2.0);


    }
}
 
#endif

