#include "score.hpp"


namespace panolyz {

    auto costFun = [](const RLGraph& mg, const RLGraphControls & controls) -> double {
        auto & vps = controls.vanishingPoints;
        auto vars = SolveVariables(mg, controls, false);
        NormalizeVariables(mg, controls, vars);
        auto planes = Instances<RegionData>(mg, controls, vars);
        auto lines = Instances<LineData>(mg, controls, vars);

        double infPunishment = 0.0;
        double originPunishment = 0.0;
        double nonManhattanPunishment = 0.0;
        double boundaryDistancePunishment = 0.0;
        double occlusionPunishment = 0.0;

        for (auto it = planes.begin(); it != planes.end(); ++it){
            if (!controls[it.hd()].used){
                continue;
            }
            Plane3 & plane = *it;

            // inf punishment
            if (HasValue(plane, IsInfOrNaN<double>)){
                infPunishment += 10.0;
            }

            // origin punishment
            if (!(plane.distanceTo(Point3(0, 0, 0)) > 0.05)){
                originPunishment += 10.0;
            }

            // non manhattan punishment
            double minAngleToVPs = std::min({
                AngleBetweenUndirectedVectors(plane.normal, vps[0]),
                AngleBetweenUndirectedVectors(plane.normal, vps[1]),
                AngleBetweenUndirectedVectors(plane.normal, vps[2])
            });
            if (minAngleToVPs > DegreesToRadians(10)){
                nonManhattanPunishment += mg.data(it.hd()).area * DepthAt(mg.data(it.hd()).normalizedCenter, plane) /
                    std::max(1e-4, cos(AngleBetweenUndirectedVectors(mg.data(it.hd()).normalizedCenter, plane.normal)));
            }
        }

        for (auto it = lines.begin(); it != lines.end(); ++it){
            if (!controls[it.hd()].used){
                continue;
            }
            Line3 & line = *it;

            // inf punishment
            if (HasValue(line, IsInfOrNaN<double>)){
                infPunishment += 10.0;
            }

            // origin punishment
            if (!(DistanceFromPointToLine(Point3(0, 0, 0), line).first > 0.05)){
                originPunishment += 10.0;
            }

            // non manhattan punishment
            double minAngleToVPs = std::min({
                AngleBetweenUndirectedVectors(line.direction(), vps[0]),
                AngleBetweenUndirectedVectors(line.direction(), vps[1]),
                AngleBetweenUndirectedVectors(line.direction(), vps[2])
            });
            if (minAngleToVPs > DegreesToRadians(10)){
                nonManhattanPunishment += line.length();
            }
        }

        // connection consistency
        // occlusion punishment
        for (auto & b : mg.constraints<RegionBoundaryData>()){
            if (!controls[b.topo.hd].used){
                occlusionPunishment += 1;
            }
            double distanceSum = 0.0;
            auto & plane1 = planes[b.topo.component<0>()];
            auto & plane2 = planes[b.topo.component<1>()];
            for (auto & ss : b.data.normalizedSampledPoints){
                for (auto & s : ss){
                    double d1 = DepthAt(s, plane1);
                    double d2 = DepthAt(s, plane2);
                    if (IsInfOrNaN(d1) || IsInfOrNaN(d2)){
                        boundaryDistancePunishment += 100.0;
                        continue;
                    }
                    assert(d1 >= 0 && d2 >= 0);
                    boundaryDistancePunishment += abs(d1 - d2);
                }
            }
        }
        for (auto & b : mg.constraints<RegionLineConnectionData>()){
            if (!controls[b.topo.hd].used){
                occlusionPunishment += 1;
            }
            double distanceSum = 0.0;
            auto & inst1 = planes[b.topo.component<0>()];
            auto & inst2 = lines[b.topo.component<1>()];
            for (auto & ss : b.data.normalizedAnchors){
                double d1 = DepthAt(ss, inst1);
                double d2 = DepthAt(ss, inst2);
                if (IsInfOrNaN(d1) || IsInfOrNaN(d2)){
                    boundaryDistancePunishment += 100.0;
                    continue;
                }
                assert(d1 >= 0 && d2 >= 0);
                boundaryDistancePunishment += abs(d1 - d2);
            }
        }
        for (auto & b : mg.constraints<LineRelationData>()){
            if (!controls[b.topo.hd].used){
                occlusionPunishment += 1;
            }
            double distanceSum = 0.0;
            auto & inst1 = lines[b.topo.component<0>()];
            auto & inst2 = lines[b.topo.component<1>()];
            auto & ss = b.data.normalizedRelationCenter;
            double d1 = DepthAt(ss, inst1);
            double d2 = DepthAt(ss, inst2);
            if (IsInfOrNaN(d1) || IsInfOrNaN(d2)){
                boundaryDistancePunishment += 100.0;
                continue;
            }
            assert(d1 >= 0 && d2 >= 0);
            boundaryDistancePunishment += abs(d1 - d2) * 2;
        }

        infPunishment *= 10.0;
        originPunishment *= 0.01;
        nonManhattanPunishment *= 30.0;
        boundaryDistancePunishment *= 0.05;
        occlusionPunishment *= 1.0;

        std::cout << "infPunishment: " << infPunishment << std::endl;
        std::cout << "originPunishment: " << originPunishment << std::endl;
        std::cout << "nonManhattanPunishment: " << nonManhattanPunishment << std::endl;
        std::cout << "boundaryDistancePunishment: " << boundaryDistancePunishment << std::endl;
        std::cout << "occlusionPunishment: " << occlusionPunishment << std::endl;

        return infPunishment + originPunishment + nonManhattanPunishment + boundaryDistancePunishment + occlusionPunishment;
    };


    double RegionManhattanCost(const RLGraph & mg, const RLGraphControls & controls, const Cache & cache, const Features & features){
        double manhattanCost = 0.0;
        auto & vps = controls.vanishingPoints;
        auto & planes = cache.planes;
        for (auto it = planes.begin(); it != planes.end(); ++it){
            if (!controls[it.hd()].used){
                continue;
            }
            const Plane3 & plane = *it;

            // non manhattan punishment
            double minAngleToVPs = std::min({
                AngleBetweenUndirectedVectors(plane.normal, vps[0]),
                AngleBetweenUndirectedVectors(plane.normal, vps[1]),
                AngleBetweenUndirectedVectors(plane.normal, vps[2])
            });

            auto & gc = features.gc[it.hd()];
            if (minAngleToVPs > DegreesToRadians(10)){
                double gcNonManhattan = gc[2] + gc[3] + gc[4];
                manhattanCost += mg.data(it.hd()).area * DepthAt(mg.data(it.hd()).normalizedCenter, plane) /
                    std::max(1e-4, cos(AngleBetweenUndirectedVectors(mg.data(it.hd()).normalizedCenter, plane.normal))) *
                    gcNonManhattan;
            }
        }
        return manhattanCost;
    }

    double LineManhattanCost(const RLGraph & mg, const RLGraphControls & controls, const Cache & cache, const Features & features){
        double lineCost = 0.0;
        auto & vps = controls.vanishingPoints;
        auto & lines = cache.lines;
        for (auto it = lines.begin(); it != lines.end(); ++it){
            if (!controls[it.hd()].used){
                continue;
            }
            const Line3 & line = *it;

            // non manhattan punishment
            double minAngleToVPs = std::min({
                AngleBetweenUndirectedVectors(line.direction(), vps[0]),
                AngleBetweenUndirectedVectors(line.direction(), vps[1]),
                AngleBetweenUndirectedVectors(line.direction(), vps[2])
            });
            if (minAngleToVPs > DegreesToRadians(10)){
                lineCost += line.length();
            }
        }
        return lineCost;
    }


    double RegionBoundaryConnectionCost(const RLGraph & mg, const RLGraphControls & controls,
        const Cache & cache, const Features & features){
        double bCost = 0.0;
        auto & planes = cache.planes;
        for (auto & b : mg.constraints<RegionBoundaryData>()){
            if (!controls[b.topo.hd].used){
                continue;
            }
            double distanceSum = 0.0;
            auto & plane1 = planes[b.topo.component<0>()];
            auto & plane2 = planes[b.topo.component<1>()];
            for (auto & ss : b.data.normalizedSampledPoints){
                for (auto & s : ss){
                    double d1 = DepthAt(s, plane1);
                    double d2 = DepthAt(s, plane2);
                    if (IsInfOrNaN(d1) || IsInfOrNaN(d2)){
                        bCost += 10.0;
                        continue;
                    }
                    assert(d1 >= 0 && d2 >= 0);
                    bCost += abs(d1 - d2);
                }
            }
        }
        return bCost;
    }


    double RegionLineConnectionCost(const RLGraph & mg, const RLGraphControls & controls,
        const Cache & cache, const Features & features){
        double bCost = 0.0;
        auto & planes = cache.planes;
        auto & lines = cache.lines;
        for (auto & b : mg.constraints<RegionLineConnectionData>()){
            if (!controls[b.topo.hd].used){
                continue;
            }
            double distanceSum = 0.0;
            auto & inst1 = planes[b.topo.component<0>()];
            auto & inst2 = lines[b.topo.component<1>()];
            for (auto & ss : b.data.normalizedAnchors){
                double d1 = DepthAt(ss, inst1);
                double d2 = DepthAt(ss, inst2);
                if (IsInfOrNaN(d1) || IsInfOrNaN(d2)){
                    bCost += 10.0;
                    continue;
                }
                assert(d1 >= 0 && d2 >= 0);
                bCost += abs(d1 - d2);
            }
        }
        return bCost;

    }


    double LineConnectionCost(const RLGraph & mg, const RLGraphControls & controls,
        const Cache & cache, const Features & features){
        double bCost = 0.0;
        auto & planes = cache.planes;
        auto & lines = cache.lines;
        for (auto & b : mg.constraints<LineRelationData>()){
            if (!controls[b.topo.hd].used){
                continue;
            }
            double distanceSum = 0.0;
            auto & inst1 = lines[b.topo.component<0>()];
            auto & inst2 = lines[b.topo.component<1>()];
            auto & ss = b.data.normalizedRelationCenter;
            double d1 = DepthAt(ss, inst1);
            double d2 = DepthAt(ss, inst2);
            if (IsInfOrNaN(d1) || IsInfOrNaN(d2)){
                bCost += 100.0;
                continue;
            }
            assert(d1 >= 0 && d2 >= 0);
            bCost += abs(d1 - d2) * 2;
        }
        return bCost;
    }



}