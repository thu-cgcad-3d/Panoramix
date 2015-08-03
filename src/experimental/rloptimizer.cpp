#include "../core/utility.hpp"
#include "../core/containers.hpp"
#include "../core/clock.hpp"
#include "../gui/scene.hpp"
#include "projective_solver.hpp"
#include "rloptimizer.hpp"

namespace pano {
    namespace experimental {


        RLGraphPatch SamplePatch(const RLGraph & g, RHandle centerRh, int maxRNum, double angleRadius) {
            RLGraphPatch patch;
            MaxHeap<RHandle, double> heapr;
            MaxHeap<LHandle, double> heapl;
            heapr.push(centerRh, 1.0);

            auto & center = g.data(centerRh).normalizedCenter;

            while ((!heapr.empty() || !heapl.empty()) && patch.size() < maxRNum) {
                RHandle curRh;
                double curRhScore = -1.0;
                if (!heapr.empty()) {
                    curRh = heapr.top();
                    curRhScore = heapr.topScore();
                }
                LHandle curLh;
                double curLhScore = -1.0;
                if (!heapl.empty()) {
                    curLh = heapl.top();
                    curLhScore = heapl.topScore();
                }

                std::vector<RHandle> rhcands;
                std::vector<LHandle> lhcands;
                if (curLhScore < curRhScore) { // insert a r
                    assert(curRh.valid());
                    patch.insert(curRh);
                    // insert related rr rl rrl
                    for (auto & h : g.topo(curRh).constraints<RRData>()) {
                        if (patch.contains(g.topo(h).component<0>()) && patch.contains(g.topo(h).component<1>())) {
                            patch.insert(h);
                        }
                    }
                    for (auto & h : g.topo(curRh).constraints<RLData>()) {
                        if (patch.contains(g.topo(h).component<0>()) && patch.contains(g.topo(h).component<1>())) {
                            patch.insert(h);
                        }
                    }
                    for (auto & h : g.topo(curRh).constraints<RRLData>()) {
                        if (patch.contains(g.topo(h).component<0>()) &&
                            patch.contains(g.topo(h).component<1>()) &&
                            patch.contains(g.topo(h).component<2>())) {
                            patch.insert(h);
                        }
                    }

                    heapr.pop();
                    for (auto rrh : g.topo(curRh).constraints<RRData>()) {
                        auto anotherRh = g.topo(rrh).component<0>();
                        if (anotherRh == curRh) {
                            anotherRh = g.topo(rrh).component<1>();
                        }
                        rhcands.push_back(anotherRh);
                    }
                    for (auto rlh : g.topo(curRh).constraints<RLData>()) {
                        LHandle lh = g.topo(rlh).component<1>();
                        lhcands.push_back(lh);
                    }
                    for (auto rrlh : g.topo(curRh).constraints<RRLData>()) {
                        RHandle anotherRh = g.topo(rrlh).component<0>();
                        if (anotherRh == curRh) {
                            anotherRh = g.topo(rrlh).component<1>();
                        }
                        rhcands.push_back(anotherRh);
                        LHandle lh = g.topo(rrlh).component<2>();
                        lhcands.push_back(lh);
                    }
                } else { // insert a l
                    assert(curLh.valid());
                    patch.insert(curLh);
                    // insert related rl rrl ll
                    for (auto & h : g.topo(curLh).constraints<RLData>()) {
                        if (patch.contains(g.topo(h).component<0>()) && patch.contains(g.topo(h).component<1>())) {
                            patch.insert(h);
                        }
                    }
                    for (auto & h : g.topo(curLh).constraints<RRLData>()) {
                        if (patch.contains(g.topo(h).component<0>()) &&
                            patch.contains(g.topo(h).component<1>())) {
                            patch.insert(g.topo(h).component<2>());
                            patch.insert(h);
                        }
                        if (patch.contains(g.topo(h).component<1>()) &&
                            patch.contains(g.topo(h).component<2>())) {
                            patch.insert(g.topo(h).component<0>());
                            patch.insert(h);
                        }
                        if (patch.contains(g.topo(h).component<2>()) &&
                            patch.contains(g.topo(h).component<0>())) {
                            patch.insert(g.topo(h).component<1>());
                            patch.insert(h);
                        }
                    }
                    for (auto & h : g.topo(curLh).constraints<LLData>()) {
                        if (patch.contains(g.topo(h).component<0>()) && patch.contains(g.topo(h).component<1>())) {
                            patch.insert(h);
                        }
                    }

                    heapl.pop();

                    for (auto llh : g.topo(curLh).constraints<LLData>()) {
                        auto anotherLh = g.topo(llh).component<0>();
                        if (anotherLh == curLh) {
                            anotherLh = g.topo(llh).component<1>();
                        }
                        lhcands.push_back(anotherLh);
                    }
                }

                for (auto rhcand : rhcands) {
                    if (patch.contains(rhcand)) {
                        continue;
                    }
                    double angleDist = AngleBetweenDirections(center, g.data(rhcand).normalizedCenter);
                    if (angleDist > angleRadius) {
                        continue;
                    }
                    double score = 1.0 / angleDist;
                    heapr.push(rhcand, score);
                }
                for (auto lhcand : lhcands) {
                    if (patch.contains(lhcand)) {
                        continue;
                    }
                    auto nearest = DistanceFromPointToLine(center, g.data(lhcand).normalizedLine).second.position;
                    double angleDist = AngleBetweenDirections(center, nearest);
                    if (angleDist > angleRadius) {
                        continue;
                    }
                    double score = 1.0 / angleDist;
                    heapl.push(lhcand, score);
                }
            }

            for (auto & c : g.constraints<RRLData>()) {
                auto h = c.topo.hd;
                if (patch.contains(g.topo(h).component<0>()) &&
                    patch.contains(g.topo(h).component<1>())) {
                    patch.insert(g.topo(h).component<2>());
                    patch.insert(h);
                }
                if (patch.contains(g.topo(h).component<1>()) &&
                    patch.contains(g.topo(h).component<2>())) {
                    patch.insert(g.topo(h).component<0>());
                    patch.insert(h);
                }
                if (patch.contains(g.topo(h).component<2>()) &&
                    patch.contains(g.topo(h).component<0>())) {
                    patch.insert(g.topo(h).component<1>());
                    patch.insert(h);
                }

                if (!patch.contains(h))
                    continue;

                auto rh1 = g.topo(h).component<0>();
                auto rh2 = g.topo(h).component<1>();
                auto lh = g.topo(h).component<2>();

                // find boundary between rh1 and rh2
                int morerradded = 0;
                for (const auto & rrd : g.constraints<RRData>()) {
                    auto hh = rrd.topo.hd;
                    if (g.topo(hh).component<0>() == rh1 && g.topo(hh).component<1>() == rh2) {
                        if (!patch.contains(hh)) {
                            patch.insert(hh);
                            morerradded++;
                        }
                    }
                }    
                int morerldadded = 0;
                for (const auto & rld : g.constraints<RLData>()) {
                    auto hh = rld.topo.hd;
                    if ((g.topo(hh).component<0>() == rh1 || g.topo(hh).component<0>() == rh2) && 
                        g.topo(hh).component<1>() == lh) {
                        if (!patch.contains(hh)) {
                            patch.insert(hh);
                            morerldadded++;
                        }
                    }
                }

                std::cout << morerldadded << std::endl;

            }

            patch.centerRh = centerRh;
            return patch;
        }


        std::pair<std::vector<RLGraphPatch::HandleType>, std::vector<int>> RLGraphPatch::handleInfo() const {
            std::vector<RLGraphPatch::HandleType> htypes;
            std::vector<int> hidx;
            htypes.reserve(size());
            hidx.reserve(size());
            for (auto & h : container<RHandle>()) {
                htypes.push_back(RH);
                hidx.push_back(h.id);
            }
            for (auto & h : container<LHandle>()) {
                htypes.push_back(LH);
                hidx.push_back(h.id);
            }
            for (auto & h : container<RRHandle>()) {
                htypes.push_back(RRH);
                hidx.push_back(h.id);
            }
            for (auto & h : container<LLHandle>()) {
                htypes.push_back(LLH);
                hidx.push_back(h.id);
            }
            return std::make_pair(std::move(htypes), std::move(hidx));
        }

        RLGraphPatchDict<int> RLGraphPatch::handlePositions() const {
            RLGraphPatchDict<int> pos;
            int p = 0;
            for (auto & h : container<RHandle>()) {
                pos[h] = p++;
            }
            for (auto & h : container<LHandle>()) {
                pos[h] = p++;
            }
            for (auto & h : container<RRHandle>()) {
                pos[h] = p++;
            }
            for (auto & h : container<LLHandle>()) {
                pos[h] = p++;
            }
            return pos;
        }



        void Visualize(const RLGraph & g, const PanoramicView & view,
            const std::unordered_map<RHandle, Plane3> & recplanes,
            const std::unordered_map<LHandle, Line3> & reclines) {
            gui::ResourceStore::set("texture", view.image);

            std::vector<double> centerDepths;

            std::vector<gui::SpatialProjectedPolygon> polygons;
            std::vector<Line3> lines;
            for (auto & p : recplanes) {
                centerDepths.push_back(p.second.distanceTo(view.camera.center()));
                auto & rd = g.data(p.first);
                gui::SpatialProjectedPolygon spp;
                for (auto & cs : rd.normalizedContours) {
                    for (auto & c : cs) {
                        spp.corners.push_back(c);
                    }
                }
                spp.plane = p.second;
                if (HasValue(spp.plane, IsInfOrNaN<double>))
                    continue;
                spp.projectionCenter = view.camera.center();
                polygons.push_back(spp);
            }
            for (auto & l : reclines) {
                if (HasValue(l.second, IsInfOrNaN<double>))
                    continue;
                centerDepths.push_back(DistanceFromPointToLine(view.camera.center(), l.second).first);
                lines.push_back(l.second);
            }

            std::nth_element(centerDepths.begin(), centerDepths.begin() + centerDepths.size() / 2, centerDepths.end());
            double d = centerDepths[centerDepths.size() / 2];
            for (auto & spp : polygons) {
                spp.plane.anchor /= d;
                assert(!HasValue(spp.plane.anchor, IsInfOrNaN<double>));
            }
            for (auto & l : lines) {
                l /= d;
                assert(!HasValue(l, IsInfOrNaN<double>));
            }

            gui::SceneBuilder viz;
            viz.begin(polygons).shaderSource(gui::OpenGLShaderSourceDescriptor::XPanorama).resource("texture").end();
            viz.begin(lines).shaderSource(gui::OpenGLShaderSourceDescriptor::XLines).end();
            viz.show(true, false, gui::RenderOptions()
                .renderMode(gui::RenderModeFlag::Triangles | gui::RenderModeFlag::Lines)
                .backgroundColor(gui::White)
                .bwColor(0.0)
                .bwTexColor(1.0)
                .cullBackFace(false)
                .cullFrontFace(false)
                .camera(PerspectiveCamera(1000, 800, Point2(500, 400), 800, Point3(-1, 1, 1), Point3(0, 0, 0))));
        }




        void AppliedRLCons(const RLGraph & g, const RLGraphLabelsProvider & labelsProv, const RLGraphPatch & patch,
            const RLGraphPatchDict<RLGraphLabel> & cl,
            std::vector<std::pair<RHandle, LHandle>> & rlcons, std::vector<const std::vector<Vec3>*> & nanchorsPtrTable) {

            rlcons.clear();
            rlcons.reserve(patch.container<RLHandle>().size() + patch.container<RRLHandle>().size());

            nanchorsPtrTable.clear();
            nanchorsPtrTable.reserve(rlcons.capacity());

            for (auto & c : patch.container<RLHandle>()) {
                auto rlh = c;
                auto rh = g.topo(rlh).component<0>();
                if (cl.at(rh) == labelsProv.rNotPlanar()) { // rh is nonplanar
                    continue;
                }
                auto lh = g.topo(rlh).component<1>();
                rlcons.emplace_back(rh, lh);
                nanchorsPtrTable.push_back(&g.data(rlh).normalizedAnchors);
            }

            for (auto & c : patch.container<RRLHandle>()) {
                auto rrlh = c;
                auto rh1 = g.topo(rrlh).component<0>();
                auto rh2 = g.topo(rrlh).component<1>();
                auto lh = g.topo(rrlh).component<2>();
                // find boundary between rh1 and rh2
                RRHandle rrh;
                for (const RRHandle & h : patch.container<RRHandle>()) {
                    if (g.topo(h).component<0>() == rh1 && g.topo(h).component<1>() == rh2) {
                        rrh = h;
                        break;
                    }
                }
                assert(rrh.valid());

                int rlabel1 = cl.at(rh1);
                int rlabel2 = cl.at(rh2);
                bool r1isnonplanar = rlabel1 == labelsProv.rNotPlanar();
                bool r2isnonplanar = rlabel2 == labelsProv.rNotPlanar();
                if (r1isnonplanar && r2isnonplanar) { // all of them are nonplanar
                    continue;
                }

                int rrlabel = cl.at(rrh);
                bool connect_rh1_lh = (rrlabel == labelsProv.rrConnected() || rrlabel == labelsProv.rrFirstIsCloser()) && !r1isnonplanar;
                bool connect_rh2_lh = (rrlabel == labelsProv.rrConnected() || rrlabel == labelsProv.rrSecondIsCloser()) && !r2isnonplanar;

                if (connect_rh1_lh) {
                    rlcons.emplace_back(rh1, lh);
                    nanchorsPtrTable.push_back(&g.data(rrlh).normalizedAnchors);
                }
                if (connect_rh2_lh) {
                    rlcons.emplace_back(rh2, lh);
                    nanchorsPtrTable.push_back(&g.data(rrlh).normalizedAnchors);
                }
            }

        }

        RLGraphLabel OrientationToRLabel(const std::vector<Vec3> & vps, 
            const RLGraphLabelsProvider & provideLabels, const Vec3 & rcenter) {
            for (int i = 0; i < vps.size(); i++) {
                if (AngleBetweenDirections(vps[i], rcenter) < DegreesToRadians(5)) {
                    return provideLabels.rTowardVP(i);
                }
            }
           /* if (AngleBetweenUndirectedVectors(vps[0], rcenter) > M_PI_2 - DegreesToRadians(3)) {
                return provideLabels.rVertical();
            }*/
            return provideLabels.rFree();
        }
        RLGraphLabel GCMeanToRLabel(const RLGraphLabelsProvider & provideLabels, const Vec5 & gcMean) {
            auto maxInd = std::max_element(std::begin(gcMean), std::end(gcMean)) - std::begin(gcMean);
            if (gcMean[maxInd] > 0.7) {
                switch (GeometricContextIndex(maxInd)) {
                case GeometricContextIndex::CeilingOrSky:
                case GeometricContextIndex::FloorOrGround:
                    return provideLabels.rTowardVP(0);
                case GeometricContextIndex::Vertical:
                    return provideLabels.rVertical();
                default:
                    break;
                }
            }
            return provideLabels.lFree();
        }


        RLGraphPatch MakeASinglePatch(const RLGraph & g) {
            RLGraphPatch p;
            for (auto & c : g.components<RData>()) {
                p.insert(c.topo.hd);
            }
            for (auto & c : g.components<LData>()) {
                p.insert(c.topo.hd);
            }
            for (auto & c : g.constraints<RRData>()) {
                p.insert(c.topo.hd);
            }
            for (auto & c : g.constraints<LLData>()) {
                p.insert(c.topo.hd);
            }
            for (auto & c : g.constraints<RLData>()) {
                p.insert(c.topo.hd);
            }
            for (auto & c : g.constraints<RRLData>()) {
                p.insert(c.topo.hd);
            }
            // center anchor
            
            // or try largest region
            RHandle largestRh;
            double maxArea = 0.0;
            for (auto & rh : p.container<RHandle>()) {
                double a = g.data(rh).area;
                if (a > maxArea) {
                    largestRh = rh;
                    maxArea = a;
                }
            }
            p.centerRh = largestRh;

            return p;
        }


        void FillSingleLabels(const RLGraph & g, const std::vector<Vec3> & vps, const RLGraphPatch & patch, 
            RLGraphPatchDict<RLGraphLabel> & labels) {
            RLGraphLabelsProvider labelsProv{ vps.size() };
            for (auto h : patch.container<RHandle>()) {
                auto & c = g.data(h);
                auto rlabel = OrientationToRLabel(vps, labelsProv, c.normalizedCenter);
                if (rlabel == labelsProv.rFree()) {
                    rlabel = GCMeanToRLabel(labelsProv, c.gcMean);
                }
                labels[h] = rlabel;
            }
            for (auto h : patch.container<LHandle>()) {
                auto & c = g.data(h);
                if (c.initialClaz == -1) {
                    labels[h] = labelsProv.lFree();
                } else {
                    auto llabel = labelsProv.lTowardVP(c.initialClaz);
                    labels[h] = llabel;
                }
            }
            for (auto h : patch.container<RRHandle>()) {
                auto rrlabel = //g.data(h).occDetectionResult == RRData::Connected ? labelsProv.rrConnected() : labelsProv.rrDisconnected();
                    labelsProv.rrConnected();
                labels[h] = rrlabel;
            }
            for (auto h : patch.container<LLHandle>()) {
                auto & c = g.data(h);
                /*labels[h] = (c.mjType == LLData::I || c.mjType == LLData::IAcrossView ||
                    c.mjType == LLData::K || c.mjType == LLData::Y || c.mjType == LLData::Hex) ? 
                    labelsProv.llConnected() : labelsProv.llDisconnected();*/
                labels[h] = labelsProv.llDisconnected();
            }
        }




        //double MedianCenterDepth(const RLGraph & mg, const RLGraphControls & controls,
        //    const RLGraphVars & vars) {

        //    std::vector<double> centerDepths;
        //    centerDepths.reserve(mg.internalComponents<RegionData>().size() + mg.internalComponents<LineData>().size());
        //    for (auto & c : mg.components<RegionData>()) {
        //        if (!controls[c.topo.hd].used)
        //            continue;
        //        double d = DepthAt(c.data.normalizedCenter, Instance(mg, controls, vars, c.topo.hd));
        //        if (!IsInfOrNaN(d)) {
        //            centerDepths.push_back(d);
        //        }
        //        //assert(!IsInfOrNaN(centerDepths.back()));
        //    }
        //    for (auto & c : mg.components<LineData>()) {
        //        if (!controls[c.topo.hd].used)
        //            continue;
        //        double d = DepthAt(normalize(c.data.line.center()), Instance(mg, controls, vars, c.topo.hd));
        //        if (!IsInfOrNaN(d)) {
        //            centerDepths.push_back(d);
        //        }
        //        //assert(!IsInfOrNaN(centerDepths.back()));
        //    }
        //    std::nth_element(centerDepths.begin(), centerDepths.begin() + centerDepths.size() / 2, centerDepths.end());
        //    double medianCenterDepth = centerDepths[centerDepths.size() / 2];
        //    return medianCenterDepth;
        //}




        bool Reconstruct(const RLGraph & g, const std::vector<Vec3> & vps, 
            const RLGraphPatch & patch, 
            const RLGraphPatchDict<RLGraphLabel> & labels,
            std::unordered_map<RHandle, Plane3> & planes, 
            std::unordered_map<LHandle, Line3> & lines) {

            RLGraphLabelsProvider labelsProv{ vps.size() };

            if (labels.at(patch.centerRh) == labelsProv.rNotPlanar())
                return false;
            
            auto up = vps[0];

            // setup solver
            // entities
            ProjectiveSolver solver;
            std::unordered_map<RHandle, int> rh2solverEntId;
            for (auto & h : patch.container<RHandle>()) {
                int label = labels.at(h);
                auto & plane = planes[h];
                if (label == labelsProv.rFree()) { // free
                    rh2solverEntId[h] = solver.bindPlaneDoF3(plane);
                } else if (labelsProv.rIsToVP(label)) { // to vp
                    plane = Plane3(vps[labelsProv.rToVPId(label)], vps[labelsProv.rToVPId(label)]);
                    rh2solverEntId[h] = solver.bindPlaneDoF1(plane);
                } else if (label == labelsProv.rVertical()) { // vertical
                    rh2solverEntId[h] = solver.bindPlaneDoF2(plane, vps[0]);
                } else if (label == labelsProv.rNotPlanar()) { // nonplanar                    
                } else {
                    assert(false);
                }
            }
            std::unordered_map<LHandle, int> lh2solverEntId;
            for (auto & h : patch.container<LHandle>()) {
                int label = labels.at(h);
                auto & line = lines[h];
                line = g.data(h).normalizedLine;
                if (label == labelsProv.lFree()) { // free
                    lh2solverEntId[h] = solver.bindLineDoF2(line);
                } else if (labelsProv.lIsToVP(label)) { // to vp
                    Ray3 ray(line.center(), vps[labelsProv.lToVPId(label)]);
                    auto p1 = DistanceBetweenTwoLines(ray, Ray3(Origin(), line.first)).second.first;
                    auto p2 = DistanceBetweenTwoLines(ray, Ray3(Origin(), line.second)).second.second;
                    line = Line3(p1, p2);
                    lh2solverEntId[h] = solver.bindLineDoF1(line);
                } else {
                    assert(false);
                }
            }

            // anchors
            static const double depthLB = 1e-3;
            static const double depthUB = 1e3;
            static const double depthAnchor = 1.0;

            // rr
            for (auto & h : patch.container<RRHandle>()) {
                int label = labels.at(h);
                auto rrh = h;
                auto rh1 = g.topo(rrh).component<0>();
                auto rh2 = g.topo(rrh).component<1>();
                if (labels.at(rh1) == labelsProv.rNotPlanar() || labels.at(rh2) == labelsProv.rNotPlanar()) { // one of them is nonplanar
                    continue;
                }
                int ent1 = rh2solverEntId.at(rh1);
                int ent2 = rh2solverEntId.at(rh2);
                auto & points = g.data(rrh).normalizedEdges;
                // register anchors
                std::vector<int> anchorIds;
                for (auto & ps : points) {
                    for (auto & p : ps) {
                        anchorIds.push_back(solver.makeNormalAnchor(p));
                    }
                }
                if (label == labelsProv.rrConnected()) { // connected
                    for (int aid : anchorIds) {
                        solver.makeAEqualToBAt(ent1, ent2, aid);
                    }
                } else if (label == labelsProv.rrFirstIsCloser()) { // 1 closer
                    for (int aid : anchorIds) {
                        solver.makeACloserThanBAt(ent1, ent2, aid);
                    }
                } else if (label == labelsProv.rrSecondIsCloser()) { // 2 closer
                    for (int aid : anchorIds) {
                        solver.makeAFartherThanDepthAt(ent1, ent2, aid);
                    }
                } else if(label == labelsProv.rrDisconnected()){
                    // do nothing
                } else{
                    assert(false);
                }

                // depth lb/ub
                for (int aid : anchorIds) {
                    solver.makeACloserThanDepthAt(ent1, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent1, depthLB, aid);
                    solver.makeACloserThanDepthAt(ent2, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent2, depthLB, aid);
                }
            }

            // ll
            for (LLHandle h : patch.container<LLHandle>()) {
                int label = labels.at(h);
                auto llh = h;
                auto lh1 = g.topo(llh).component<0>();
                auto lh2 = g.topo(llh).component<1>();
                int ent1 = lh2solverEntId.at(lh1);
                int ent2 = lh2solverEntId.at(lh2);
                auto & point = g.data(llh).normalizedRelationCenter;
                // register anchor
                int anchorId = solver.makeNormalAnchor(point);
                if (label == labelsProv.llConnected()) { // connected
                    solver.makeAEqualToBAt(ent1, ent2, anchorId);
                } else if (label == labelsProv.llDisconnected()) { // disconnected
                } else {
                    assert(false);
                }

                // depth lb/ub
                {
                    int aid = anchorId;
                    solver.makeACloserThanDepthAt(ent1, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent1, depthLB, aid);
                    solver.makeACloserThanDepthAt(ent2, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent2, depthLB, aid);
                }
            }

            // rl, rrl
            std::vector<std::pair<RHandle, LHandle>> rlcons;
            std::vector<const std::vector<Vec3>*> nanchorsPtrTable;
            AppliedRLCons(g, labelsProv, patch, labels, rlcons, nanchorsPtrTable);

            assert(rlcons.size() == nanchorsPtrTable.size());

            for (int i = 0; i < rlcons.size(); i++) {
                auto rh = rlcons[i].first;
                auto lh = rlcons[i].second;
                int ent1 = rh2solverEntId.at(rh);
                int ent2 = lh2solverEntId.at(lh);
                auto & as = *nanchorsPtrTable[i];
                for (auto & a : as) {
                    int aid = solver.makeNormalAnchor(a);
                    solver.makeAEqualToBAt(ent1, ent2, aid);

                    solver.makeACloserThanDepthAt(ent1, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent1, depthLB, aid);
                    solver.makeACloserThanDepthAt(ent2, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent2, depthLB, aid);
                }
            }

            // how to put arbitrary depth anchors!!!?
            {
                auto center = normalize(g.data(patch.centerRh).normalizedCenter);
                int aid = solver.makeNormalAnchor(center);
                solver.makeAEqualToDepthAt(rh2solverEntId.at(patch.centerRh), depthAnchor, aid);
            }

            // solve
            double nanOrInfRatio = 0.0;
            bool solvable = solver.solve(&nanOrInfRatio, false);
            assert(solvable);
            return solvable;

        }








        namespace {
            bool NextSub(int * subs, const int * dims, int len) {
                if (len == 0)
                    return false;
                subs[0]++;
                int k = 0;
                while (k < len && subs[k] >= dims[k]) {
                    subs[k] = 0;
                    if (k + 1 < len) subs[k + 1] ++;
                    k++;
                }
                if (k == len) {
                    return false; 
                }
                return true;
            }
            bool NextSub(int * subs, int dim, int len) {
                if (len == 0)
                    return false;
                subs[0]++;
                int k = 0;
                while (k < len && subs[k] >= dim) {
                    subs[k] = 0;
                    if (k + 1 < len) subs[k + 1] ++;
                    k++;
                }
                if (k == len) {
                    return false;
                }
                return true;
            }
            int Sub2Ind(const int * subs, const int * dims, int len) {
                int ind = 0;
                for (int k = len - 1; k >= 0; k--) {
                    ind = ind * dims[k] + subs[k];
                }
                return ind;
            }
            int Sub2Ind(const int * subs, int dim, int len) {
                int ind = 0;
                for (int k = len - 1; k >= 0; k--) {
                    ind = ind * dim + subs[k];
                }
                return ind;
            }


        }


#define INFINITY std::numeric_limits<double>::infinity();


        // feature layout
#pragma pack(push)
#pragma pack(1)
#define Part(name, nbytes) struct {uint8_t content[(nbytes)];} name
        struct FeatureLayout {
            Part(R, 8);
            Part(L, 1);
            Part(RR, 6);
            Part(LL, LLData::ManhattanJunctionTypeNum * 2);
            Part(RRRorRRRR, 5);
            Part(Patch, 21);
        };
#undef Part
#pragma pack(pop)
        static const int FullFeatureLength = sizeof(FeatureLayout);

#define FeatureBeginPos(name) (offsetof(FeatureLayout, name))
#define FeatureEndPos(name) (offsetof(FeatureLayout, name) + sizeof(std::declval<FeatureLayout>().name))
#define FeatureLength(name) (sizeof(std::declval<FeatureLayout>().name))


        static const double LineVPScoreThreshold = 0.8;

        void RLOptimizer::determinFactorGraphVars() {

            SetClock();

            // determine allowed labels
            _allowedLabels = RLGraphAllowedLabelsTable({
                _g.internalComponents<RData>().size(),
                _g.internalComponents<LData>().size(),
                _g.internalConstraints<RRData>().size(),
                _g.internalConstraints<LLData>().size()
            });

            // r
            findPeakyRhs(DegreesToRadians(6));
            findHorizonRhs(DegreesToRadians(5));
            for (auto && c : _g.components<RData>()) {
                auto rh = c.topo.hd;
                assert(_peakyRhs.size() == _vps.size());
                bool determined = false;
                for (int i = 0; i < _vps.size(); i++) {
                    if (Contains(_peakyRhs[i], rh)) {
                        _allowedLabels[rh] = { provideLabels().rTowardVP(i), provideLabels().rNotPlanar() };
                        determined = true;
                        break;
                    }
                }
                if (determined) {
                    continue;
                }
                if (Contains(_horizonRhs, rh)) {
                    _allowedLabels[rh] = { provideLabels().rVertical(), provideLabels().rNotPlanar() };
                    continue;
                }
                auto & gcMean = _g.data(rh).gcMean;
                auto maxInd = std::max_element(std::begin(gcMean), std::end(gcMean)) - std::begin(gcMean);
                if (gcMean[maxInd] > 0.5) {
                    switch (GeometricContextIndex(maxInd)) {
                    case GeometricContextIndex::CeilingOrSky:
                    case GeometricContextIndex::FloorOrGround:
                        _allowedLabels[rh] = { provideLabels().rTowardVP(0) };
                        continue;
                    case GeometricContextIndex::Vertical:
                        _allowedLabels[rh] = { provideLabels().rVertical() };
                        continue;
                    default:
                        break;
                    }
                }
                _allowedLabels[rh].push_back(provideLabels().rFree());
                for (int i = 0; i < _vps.size(); i++) {
                    _allowedLabels[rh].push_back(provideLabels().rTowardVP(i));
                }
                _allowedLabels[rh].push_back(provideLabels().rNotPlanar());
            }
            // l
            for (auto && c : _g.components<LData>()) {
                auto lh = c.topo.hd;
                auto & vpScores = _g.data(lh).vpScores;
                assert(vpScores.size() == _vps.size());
                std::vector<int> ids(vpScores.size());
                std::iota(ids.begin(), ids.end(), 0);
                for (int i = 0; i < ids.size(); i++) {
                    if (vpScores[i] > LineVPScoreThreshold) {
                        _allowedLabels[lh].push_back(provideLabels().lTowardVP(i));
                    }
                }
                if (_allowedLabels[lh].empty()) {
                    _allowedLabels[lh].push_back(provideLabels().lFree());
                }
            }
            // rr
            for (auto && c : _g.constraints<RRData>()) {
                auto rrh = c.topo.hd;
                _allowedLabels[rrh] = { provideLabels().rrConnected(), provideLabels().rrFirstIsCloser(), provideLabels().rrSecondIsCloser() };
            }
            // ll
            for (auto && c : _g.constraints<LLData>()) {
                auto llh = c.topo.hd;
                _allowedLabels[llh] = { provideLabels().llConnected(), provideLabels().llDisconnected() };
            }


            _varhs = RLGraphVarHandleTable({
                _g.internalComponents<RData>().size(),
                _g.internalComponents<LData>().size(),
                _g.internalConstraints<RRData>().size(),
                _g.internalConstraints<LLData>().size()
            }, ml::FactorGraph::VarHandle());


            // register (discrete) var cats
            for (auto && c : _g.components<RData>()) {
                auto vc = _fg.addVarCategory(_allowedLabels[c.topo.hd].size(), 1.0);
                _varhs[c.topo.hd] = _fg.addVar(vc);
            }
            for (auto && c : _g.components<LData>()) {
                auto vc = _fg.addVarCategory(_allowedLabels[c.topo.hd].size(), 1.0);
                _varhs[c.topo.hd] = _fg.addVar(vc);
            }
            for (auto && c : _g.constraints<RRData>()) {
                auto vc = _fg.addVarCategory(_allowedLabels[c.topo.hd].size(), 1.0);
                _varhs[c.topo.hd] = _fg.addVar(vc);
            }
            for (auto && c : _g.constraints<LLData>()) {
                auto vc = _fg.addVarCategory(_allowedLabels[c.topo.hd].size(), 1.0);
                _varhs[c.topo.hd] = _fg.addVar(vc);
            }

        }


        void RLOptimizer::determinFactorGraphFactors() {

            SetClock();

            _graphFeatureDict = RLFeatureDictTable({
            	_g.internalComponents<RData>().size(),
				_g.internalComponents<LData>().size(),
				_g.internalConstraints<RRData>().size(),
				_g.internalConstraints<LLData>().size(),
				_g.internalConstraints<RRRData>().size(),
				_g.internalConstraints<RRRRData>().size()
            });
            _patchFeatureDict.resize(_patches.size());

            // r
            for (auto & r : _g.components<RData>()) {
                // compute feature
                auto h = r.topo.hd;
                auto & feaDict = _graphFeatureDict[r.topo.hd];
                auto & allowedLabels = _allowedLabels[r.topo.hd];
                feaDict = Dictionary<std::vector<double>>({ provideLabels().rNum() });
                // add fea vecs
                double areaRatio = r.data.area / _rAreaSum;
                for (RLGraphLabel label : allowedLabels) {
                    static_assert(FeatureLength(R) == 8, "");
                    RLFeatureVec fea(FeatureLength(R), 0.0);
                    if (Contains(_peakyRhs[0], r.topo.hd)) { // tend to be horizontal
                        if (label == provideLabels().rTowardVP(0)) { // judging as horizontal
                            fea[0] = areaRatio;
                        } else if (label == provideLabels().rVertical()) { // judging as vertical
                            fea[1] = areaRatio;
                        } else if (label == provideLabels().rNotPlanar()) { // ... as not planar
                            fea[2] = areaRatio;
                        } else { // judging as free
                            fea[3] = areaRatio;
                        }
                    } else if (Contains(_horizonRhs, r.topo.hd)) { // tend to be vertical
                        if (label == provideLabels().rVertical()) { // juding as vertical
                            fea[4] = areaRatio;
                        } else if (label == provideLabels().rTowardVP(0)) { // judging as horizontal
                            fea[5] = areaRatio;
                        } else if (label == provideLabels().rNotPlanar()) { // .. as not planar
                            fea[6] = areaRatio;
                        } else { // judging as free
                            fea[7] = areaRatio;
                        }
                    }
                    feaDict.insert({ label }, std::move(fea));
                }
                // add factor
                ml::FactorGraph::VarHandle vh = _varhs[r.topo.hd];
                auto cost = [this, h](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == 1);
                    int labelid = varlabels[0];
                    RLGraphLabel label = _allowedLabels[h][labelid];
                    assert(label >= 0 && label < provideLabels().rNum());
                    auto weights = static_cast<const double*>(givenData);
                    if (!_graphFeatureDict[h].contains({ label })) {
                        return INFINITY;
                    }
                    auto & fea = _graphFeatureDict[h].at({ label });
                    assert(fea.size() == FeatureLength(R));
                    return std::inner_product(weights + FeatureBeginPos(R), weights + FeatureEndPos(R), fea.begin(), 0.0);
                };
                _fg.addFactor({ vh }, _fg.addFactorCategory(ml::FactorGraph::FactorCategory{ std::move(cost), 1.0 }));
            }

            // l
            for (auto & l : _g.components<LData>()) {
                // compute feature
                auto h = l.topo.hd;
                auto & feaDict = _graphFeatureDict[l.topo.hd];
                auto & allowedLabels = _allowedLabels[l.topo.hd];
                feaDict = Dictionary<std::vector<double>>({ provideLabels().lNum() });
                // add fea vecs
                double lenRatio = l.data.normalizedLine.length() / _lLengthSum;
                auto & vpScores = l.data.vpScores;
                for (RLGraphLabel label : allowedLabels) {
                    static_assert(FeatureLength(L) == 1, "");
                    RLFeatureVec fea(FeatureLength(L), 0.0);
                    if (provideLabels().lIsToVP(label)) {
                        fea[0] = vpScores[provideLabels().lToVPId(label)] * lenRatio;
                    } else {
                        fea[0] = lenRatio;
                    }
                    feaDict.insert({ label }, std::move(fea));
                }
                // add factor
                ml::FactorGraph::VarHandle vh = _varhs[l.topo.hd];
                auto cost = [this, h](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == 1);
                    int labelid = varlabels[0];
                    RLGraphLabel label = _allowedLabels[h][labelid];
                    assert(label >= 0 && label < provideLabels().lNum());
                    auto weights = static_cast<const double*>(givenData);
                    if (!_graphFeatureDict[h].contains({ label })) {
                        return INFINITY;
                    }
                    auto & fea = _graphFeatureDict[h].at({ label });
                    assert(fea.size() == FeatureLength(L));
                    return std::inner_product(weights + FeatureBeginPos(L), weights + FeatureEndPos(L), fea.begin(), 0.0);
                };
                _fg.addFactor({ vh }, _fg.addFactorCategory(ml::FactorGraph::FactorCategory{ std::move(cost), 1.0 }));
            }

            // rr
            for (auto & rr : _g.constraints<RRData>()) {
                // compute feature
                auto h = rr.topo.hd;
                auto & feaDict = _graphFeatureDict[rr.topo.hd];
                auto & allowedLabels = _allowedLabels[rr.topo.hd];
                feaDict = Dictionary<std::vector<double>>({ provideLabels().rrNum() });
                // add fea vecs
                double lenRatio = rr.data.length / _rrLengthSum;
                int resp = ToUnderlying(rr.data.occDetectionResult);
                for (RLGraphLabel label : allowedLabels) {
                    static_assert(FeatureLength(RR) == 6, "");
                    RLFeatureVec fea(FeatureLength(RR), 0.0);
                    assert(label >= 0 && label < 3 && resp >= 0 && resp < 2);                    
                    fea[resp * 3 + label] = lenRatio;
                    feaDict.insert({ label }, std::move(fea));
                }
                // add factor
                ml::FactorGraph::VarHandle vh = _varhs[rr.topo.hd];
                auto cost = [this, h](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == 1);
                    int labelid = varlabels[0];
                    RLGraphLabel label = _allowedLabels[h][labelid];
                    assert(label >= 0 && label < provideLabels().rrNum());
                    auto weights = static_cast<const double*>(givenData);
                    if (!_graphFeatureDict[h].contains({ label })) {
                        return INFINITY;
                    }
                    auto & fea = _graphFeatureDict[h].at({ label });
                    assert(fea.size() == FeatureLength(RR));
                    return std::inner_product(weights + FeatureBeginPos(RR), weights + FeatureEndPos(RR), fea.begin(), 0.0);
                };
                _fg.addFactor({ vh }, _fg.addFactorCategory(ml::FactorGraph::FactorCategory{ std::move(cost), 1.0 }));
            }

            // ll 
            std::vector<int> mjCounts(LLData::ManhattanJunctionTypeNum, 0);
            for (auto & ll : _g.constraints<LLData>()) {
                mjCounts[ll.data.mjType] ++;
            }
            for (auto & ll : _g.constraints<LLData>()) {
                auto h = ll.topo.hd;
                auto lh1 = ll.topo.component<0>();
                auto lh2 = ll.topo.component<1>();
                auto & line1 = _g.data(lh1).normalizedLine;
                auto & line2 = _g.data(lh2).normalizedLine;

                // compute feature
                auto & feaDict = _graphFeatureDict[ll.topo.hd];
                auto & allowedLabels = _allowedLabels[ll.topo.hd];
                feaDict = Dictionary<std::vector<double>>({ provideLabels().llNum() });
                // add fea vecs
                for (RLGraphLabel label : allowedLabels) {
                    static_assert(FeatureLength(LL) == LLData::ManhattanJunctionTypeNum * 2, "");
                    RLFeatureVec fea(FeatureLength(LL), 0.0);
                    fea[ToUnderlying(ll.data.mjType) * 2 + label] = 1.0 / std::max(mjCounts[ll.data.mjType], 1);
                    feaDict.insert({ label }, std::move(fea));
                }
                // add factor
                ml::FactorGraph::VarHandle vh = _varhs[ll.topo.hd];
                auto cost = [this, h](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == 1);
                    int labelid = varlabels[0];
                    RLGraphLabel label = _allowedLabels[h][labelid];
                    assert(label >= 0 && label < provideLabels().llNum());
                    auto weights = static_cast<const double*>(givenData);
                    if (!_graphFeatureDict[h].contains({ label })) {
                        return INFINITY;
                    }
                    auto & fea = _graphFeatureDict[h].at({ label });
                    assert(fea.size() == FeatureLength(LL));
                    return std::inner_product(weights + FeatureBeginPos(LL), weights + FeatureEndPos(LL), fea.begin(), 0.0);
                };
                _fg.addFactor({ vh }, _fg.addFactorCategory(ml::FactorGraph::FactorCategory{ std::move(cost), 1.0 }));
            }

            // rrr and rrrr
            std::vector<std::vector<RHandle>> rhsInRRRorRRRRs;
            std::vector<Dictionary<RLFeatureVec> *> feaDictPtrInRRRorRRRRs;
            for (auto & rrr : _g.constraints<RRRData>()) {
                rhsInRRRorRRRRs.push_back({
                    rrr.topo.component<0>(),
                    rrr.topo.component<1>(),
                    rrr.topo.component<2>()
                });
                feaDictPtrInRRRorRRRRs.push_back(&(_graphFeatureDict[rrr.topo.hd]));
            }
            for (auto & rrrr : _g.constraints<RRRRData>()) {
                rhsInRRRorRRRRs.push_back({
                    rrrr.topo.component<0>(),
                    rrrr.topo.component<1>(),
                    rrrr.topo.component<2>(),
                    rrrr.topo.component<3>()
                });
                feaDictPtrInRRRorRRRRs.push_back(&(_graphFeatureDict[rrrr.topo.hd]));
            }

            for (int k = 0; k < rhsInRRRorRRRRs.size(); k++) {
                auto & rhs = rhsInRRRorRRRRs[k];
                auto & feaDict = *feaDictPtrInRRRorRRRRs[k];

                // get rrhs and rrh directions
                std::vector<RRHandle> rrhs(rhs.size());
                std::vector<int> rrhascend(rhs.size(), true);
                for (int i = 0; i < rhs.size(); i++) {
                    auto & nextrh = rhs[(i + 1) % 3];
                    for (const RRHandle & rrh : _g.topo(rhs[i]).constraints<RRData>()) {
                        if (_g.topo(rrh).component<0>() == nextrh || _g.topo(rrh).component<1>() == nextrh) {
                            rrhs[i] = rrh;
                            rrhascend[i] = _g.topo(rrh).component<1>() == nextrh;
                            break;
                        }
                    }
                }
                for (auto & rrh : rrhs) {
                    assert(rrh.valid());
                }             

                // count allowed labels num in each rr
                std::vector<int> rrlabelNums(rrhs.size());
                for (int i = 0; i < rrhs.size(); i++) {
                    rrlabelNums[i] = _allowedLabels[rrhs[i]].size();
                }

                // store rrlabelids, rrlabels and std rrlabels
                std::vector<int> rrlabelids(rrhs.size(), 0);
                std::vector<RLGraphLabel> rrlabels(rrhs.size());
                std::vector<RLGraphLabel> rrlabelsStd(rrhs.size());

                // feature dict for this rrr or rrrr
                feaDict = Dictionary<std::vector<double>>(std::vector<size_t>(provideLabels().rrNum(), rrhs.size()));
                // for each allowed combination of rrlabel
                do {
                    // fill rrlabels and std rrlabels
                    for (int i = 0; i < rrhs.size(); i++) {
                        rrlabels[i] = _allowedLabels[rrhs[i]][rrlabelids[i]];
                        rrlabelsStd[i] = rrhascend[i] ? rrlabels[i] : provideLabels().rrReverse(rrlabels[i]);
                    }

                    // fea vec
                    static_assert(FeatureLength(RRRorRRRR) == 5, "");
                    RLFeatureVec fea(FeatureLength(RRRorRRRR), 0.0);
                    int nRRRorRRRR = _g.internalConstraints<RRRData>().size() +
                        _g.internalConstraints<RRRRData>().size();
                    if (rrlabelsStd.size() == 0) {
                        fea[0] = 1.0;
                    } else if (rrlabelsStd.size() == 1) {
                        fea[1] = 1.0;
                    } else if (rrlabelsStd.size() == 2) {
                        if (rrlabelsStd[0] == rrlabelsStd[1]) {
                            fea[2] = 1.0;
                        } else {
                            fea[3] = 1.0;
                        }
                    } else if (rrlabelsStd.size() >= 3) {
                        fea[4] = 1.0;
                    }
                    for (auto & f : fea) {
                        f /= nRRRorRRRR;
                    }
                    feaDict.insert(rrlabels.data(), std::move(fea));
                } while (NextSub(rrlabelids.data(), rrlabelNums.data(), rrlabelids.size()));

                auto cost = [this, rrhs, &feaDict](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    assert(nvar == rrhs.size());
                    auto weights = static_cast<const double*>(givenData);
                    std::vector<RLGraphLabel> rrlabels(nvar);
                    for (int i = 0; i < nvar; i++) {
                        rrlabels[i] = _allowedLabels[rrhs[i]][varlabels[i]];
                    }
                    if (!feaDict.contains(rrlabels.data())) {
                        return INFINITY;
                    }
                    auto & fea = feaDict.at(rrlabels.data());
                    assert(fea.size() == FeatureLength(RRRorRRRR));
                    return std::inner_product(weights + FeatureBeginPos(RRRorRRRR), weights + FeatureEndPos(RRRorRRRR), fea.begin(), 0.0);
                };
                ml::FactorGraph::FactorCategory fc;
                fc.costs = std::move(cost);
                fc.c_alpha = 1.0;

                // get vhs
                std::vector<ml::FactorGraph::VarHandle> vhs(rrhs.size());
                for (int i = 0; i < vhs.size(); i++) {
                    vhs[i] = _varhs[rrhs[i]];
                }
                _fg.addFactor(vhs.begin(), vhs.end(), _fg.addFactorCategory(std::move(fc)));
            }

            // reconstructed patches            
            HandledTable<RRHandle, std::map<std::pair<int, int>, double>> 
                rrViolationOnVPPairs(_g.internalConstraints<RRData>().size());
            for (auto & rr : _g.constraints<RRData>()) {
                std::map<std::pair<int, int>, double> violates;
                for (int i = 0; i < _vps.size(); i++) {
                    for (int j = i + 1; j < _vps.size(); j++) {
                        Vec3 dir = normalize(_vps[i].cross(_vps[j]));
                        double maxVio = 0.0;
                        auto & ps = rr.data.normalizedSampledPoints;
                        for (int a1 = 0; a1 < ps.size(); a1++) {
                            for (int a2 = 0; a2 < ps[a1].size(); a2++) {
                                auto & pa = ps[a1][a2];
                                for (int b1 = a1; b1 < ps.size(); b1++) {
                                    for (int b2 = 0; b2 < ps[b1].size(); b2++) {
                                        auto & pb = ps[b1][b2];
                                        auto center = normalize(pa + pb);
                                        auto vertDir = normalize(dir.cross(center));
                                        double vio = abs(pa.dot(vertDir) - pb.dot(vertDir));
                                        if (vio > maxVio) {
                                            maxVio = vio;
                                        }
                                    }
                                }
                            }
                        }
                        rrViolationOnVPPairs[rr.topo.hd][std::make_pair(i, j)] = maxVio;
                    }
                }
            }
            for (int i = 0; i < _patches.size(); i++) {
                auto & patch = _patches[i];
               
                std::vector<RLGraphPatch::HandleType> handleTypes;
                std::vector<int> handleIndices;
                std::tie(handleTypes, handleIndices) = patch.handleInfo();
                auto handlePositions = patch.handlePositions();

                std::vector<size_t> labelNums;
                for (auto & h : patch.container<RHandle>()) {
                    labelNums.push_back(provideLabels().rNum());
                }
                for (auto & h : patch.container<LHandle>()) {
                    labelNums.push_back(provideLabels().lNum());
                }
                for (auto & h : patch.container<RRHandle>()) {
                    labelNums.push_back(provideLabels().rrNum());
                }
                for (auto & h : patch.container<LLHandle>()) {
                    labelNums.push_back(provideLabels().llNum());
                }
                auto & feaDict = _patchFeatureDict[i];
                feaDict = Dictionary<RLPatchFeatureData>(labelNums);
                
                // iterate over allowed rlabels and llabels
                std::vector<int> allowedRLLabelNums;
                allowedRLLabelNums.reserve(
                    patch.container<RHandle>().size() + 
                    patch.container<LHandle>().size());
                for (auto & h : patch.container<RHandle>()) {
                    allowedRLLabelNums.push_back(_allowedLabels[h].size());
                }
                for (auto & h : patch.container<LHandle>()) {
                    allowedRLLabelNums.push_back(_allowedLabels[h].size());
                }
                std::vector<int> rllabelids(allowedRLLabelNums.size(), 0);
                int solvenum = 0;
                do {
                    int id = 0;
                    // rl label ids -> rl labels
                    std::vector<RLGraphLabel> labels;
                    labels.reserve(patch.size());
                    for (auto & h : patch.container<RHandle>()) {
                        labels.push_back(_allowedLabels[h][rllabelids[id++]]);
                    }
                    for (auto & h : patch.container<LHandle>()) {
                        labels.push_back(_allowedLabels[h][rllabelids[id++]]);
                    }
                    // bredth first search
                    std::queue<std::vector<RLGraphLabel>> Q;
                    Q.push(labels);
                    while (!Q.empty()) {
                        // get current labels
                        std::vector<RLGraphLabel> curlabels = std::move(Q.front());
                        Q.pop();
                        assert(curlabels.size() <= labelNums.size());
                        
                        if (curlabels.size() == labelNums.size()) { // labels completed
                            // use projective solver to reconstruct 3d entities
                            // and extract 3d features
                            auto result = featureInReconstructedPatch(patch, curlabels);
                            solvenum++;
                            if (result.failed())
                                continue;
                            auto patchFeaData = result.unwrap();
                            bool isValid = true;
                            for (auto & p : Get<1>(patchFeaData)) {
                                if (HasValue(p.second, IsInfOrNaN<double>)) {
                                    isValid = false;
                                    break;
                                }
                            }
                            for (auto & l : Get<2>(patchFeaData)) {
                                if (HasValue(l.second, IsInfOrNaN<double>)) {
                                    isValid = false;
                                    break;
                                }
                            }
                            if (!isValid)
                                continue;
                            // record the feature vector if not failed
                            feaDict.insert(curlabels.data(), std::move(patchFeaData));                            
                            std::cout << "[" << i << "]" << "feanum/solvenum = " << feaDict.size() << "/" << solvenum << std::endl;

                            if (StaticStorage::has("view")) {
                                View<PanoramicCamera, Image3ub> view = StaticStorage::get("view");
                                auto & fea = feaDict.at(curlabels.data());
                                Visualize(_g, view, Get<1>(fea), Get<2>(fea));
                            }

                            continue;
                        } else { // labels not completed yet, 
                            int k = curlabels.size();
                            auto htype = handleTypes[k];
                            int hid = handleIndices[k];
                            if (htype == RLGraphPatch::RRH) {
                                RRHandle rrh(hid);
                                auto & allowedLabels = _allowedLabels[rrh];
                                auto rh1 = _g.topo(rrh).component<0>();
                                auto rh2 = _g.topo(rrh).component<1>();
                                RLGraphLabel label1 = curlabels[handlePositions.at(rh1)];
                                RLGraphLabel label2 = curlabels[handlePositions.at(rh2)];

                                bool connectable = true;
                                if (label1 == provideLabels().rNotPlanar() || label2 == provideLabels().rNotPlanar()) {
                                    connectable = false;
                                } else if (provideLabels().rIsToVP(label1) &&
                                    provideLabels().rIsToVP(label2) && 
                                    provideLabels().rToVPId(label1) != provideLabels().rToVPId(label2)) {
                                    int vpids[] = { provideLabels().rToVPId(label1), provideLabels().rToVPId(label2) };
                                    std::sort(vpids, vpids + 2);
                                    double vio = rrViolationOnVPPairs[rrh].at(std::make_pair(vpids[0], vpids[1]));
                                    if (vio > tan(DegreesToRadians(1))) {
                                        connectable = false;
                                    }
                                }
                                curlabels.push_back(-1);
                                for (RLGraphLabel allowedLabel : allowedLabels) {
                                    if (allowedLabel == provideLabels().rrConnected() && !connectable) {
                                        continue;
                                    }
                                    curlabels.back() = allowedLabel;
                                    Q.push(curlabels);
                                }

                            } else if (htype == RLGraphPatch::LLH) {
                                LLHandle llh(hid);
                                auto & allowedLabels = _allowedLabels[llh];
                                curlabels.push_back(-1);
                                for (RLGraphLabel allowedLabel : allowedLabels) {
                                    curlabels.back() = allowedLabel;
                                    Q.push(curlabels);
                                }
                            } else {
                                SHOULD_NEVER_BE_CALLED();
                            }
                        }
                    }

                } while (NextSub(rllabelids.data(), allowedRLLabelNums.data(), rllabelids.size()));

                auto cost = [this, i, &labelNums, &patch, &feaDict](const int * varlabels, size_t nvar,
                    ml::FactorGraph::FactorCategoryId fcid, void * givenData) -> double {
                    auto weights = static_cast<const double*>(givenData);
                    std::vector<RLGraphLabel> labels;
                    labels.reserve(labelNums.size());
                    int id = 0;
                    for (auto & h : patch.container<RHandle>()) {
                        labels.push_back(_allowedLabels[h][varlabels[id++]]);
                    }
                    for (auto & h : patch.container<LHandle>()) {
                        labels.push_back(_allowedLabels[h][varlabels[id++]]);
                    }
                    for (auto & h : patch.container<RRHandle>()) {
                        labels.push_back(_allowedLabels[h][varlabels[id++]]);
                    }
                    for (auto & h : patch.container<LLHandle>()) {
                        labels.push_back(_allowedLabels[h][varlabels[id++]]);
                    }

                    if (!feaDict.contains(labels.data())) {
                        return INFINITY;
                    }
                    auto & fea = Get<0>(feaDict.at(labels.data()));
                    assert(fea.size() == FeatureLength(Patch));
                    return std::inner_product(weights + FeatureBeginPos(Patch), weights + FeatureEndPos(Patch), fea.begin(), 0.0);
                };
                ml::FactorGraph::FactorCategory fc;
                fc.costs = std::move(cost);
                fc.c_alpha = 1.0;

                std::vector<ml::FactorGraph::VarHandle> vhs;
                vhs.reserve(
                    patch.container<RHandle>().size() +
                    patch.container<LHandle>().size() +
                    patch.container<RRHandle>().size() +
                    patch.container<LLHandle>().size()
                    );
                for (auto & h : patch.container<RHandle>()) {
                    vhs.push_back(_varhs[h]);
                }
                for (auto & h : patch.container<LHandle>()) {
                    vhs.push_back(_varhs[h]);
                }
                for (auto & h : patch.container<RRHandle>()) {
                    vhs.push_back(_varhs[h]);
                }
                for (auto & h : patch.container<LLHandle>()) {
                    vhs.push_back(_varhs[h]);
                }

                _fg.addFactor(vhs.begin(), vhs.end(), _fg.addFactorCategory(std::move(fc)));
            }

        }




        void RLOptimizer::preprocess() {
            _fg.clear();
            {
                _rAreaSum = 0.0;
                for (auto & r : _g.components<RData>()) {
                    _rAreaSum += r.data.area;
                }
            }
            {
                _lLengthSum = 0;
                for (auto & l : _g.components<LData>()) {
                    _lLengthSum += l.data.normalizedLine.length();
                }
            }
            {
                _rrLengthSum = 0.0;
                for (auto & rr : _g.constraints<RRData>()) {
                    _rrLengthSum += rr.data.length;
                }
            }
            determinFactorGraphVars();
            determinFactorGraphFactors();
        }


        int RLOptimizer::featureLength() const {
            return FullFeatureLength;
        }


        Failable<RLPatchFeatureData> RLOptimizer::featureInReconstructedPatch(
            const RLGraphPatch & patch, const std::vector<RLGraphLabel> & allLabelsInPatch) const {
            
            // fill labels
            RLGraphPatchDict<RLGraphLabel> labels;
            int id = 0;
            for (auto & h : patch.container<RHandle>()) {
                labels[h] = allLabelsInPatch[id++];
            }
            for (auto & h : patch.container<LHandle>()) {
                labels[h] = allLabelsInPatch[id++];
            }
            for (auto & h : patch.container<RRHandle>()) {
                labels[h] = allLabelsInPatch[id++];
            }
            for (auto & h : patch.container<LLHandle>()) {
                labels[h] = allLabelsInPatch[id++];
            }

            if (labels[patch.centerRh] == provideLabels().rNotPlanar())
                return nullptr;

            int labelsNum = id;

            RLPatchFeatureData feaData;
            auto & fea = Get<0>(feaData);
            auto & planes = Get<1>(feaData);
            auto & lines = Get<2>(feaData);

            // setup solver
            // entities
            ProjectiveSolver solver;
            std::unordered_map<RHandle, int> rh2solverEntId;
            for (auto & c : labels.container<RHandle>()) {
                int label = c.second;
                auto & plane = planes[c.first];
                if (label == provideLabels().rFree()) { // free
                    rh2solverEntId[c.first] = solver.bindPlaneDoF3(plane);
                } else if (provideLabels().rIsToVP(label)) { // to vp
                    plane = Plane3(_vps[provideLabels().rToVPId(label)], _vps[provideLabels().rToVPId(label)]);
                    rh2solverEntId[c.first] = solver.bindPlaneDoF1(plane);
                } else if (label == provideLabels().rVertical()) { // vertical
                    rh2solverEntId[c.first] = solver.bindPlaneDoF2(plane, up());
                } else if (label == provideLabels().rNotPlanar()) { // nonplanar                    
                } else {
                    assert(false);
                }
            }
            std::unordered_map<LHandle, int> lh2solverEntId;
            for (auto & c : labels.container<LHandle>()) {
                int label = c.second;
                auto & line = lines[c.first];
                line = _g.data(c.first).normalizedLine;
                if (label == provideLabels().lFree()) { // free
                    lh2solverEntId[c.first] = solver.bindLineDoF2(line);
                } else if (provideLabels().lIsToVP(label)) { // to vp
                    Ray3 ray(line.center(), _vps[provideLabels().lToVPId(label)]);
                    auto p1 = DistanceBetweenTwoLines(ray, Ray3(Origin(), line.first)).second.first;
                    auto p2 = DistanceBetweenTwoLines(ray, Ray3(Origin(), line.second)).second.second;
                    line = Line3(p1, p2);
                    lh2solverEntId[c.first] = solver.bindLineDoF1(line);
                } else {
                    assert(false);
                }
            }

            // anchors
            static const double depthLB = 1e-2;
            static const double depthUB = 1e2;
            static const double depthAnchor = 1.0;

            // rr
            for (auto & c : labels.container<RRHandle>()) {
                int label = c.second;
                auto rrh = c.first;
                auto rh1 = _g.topo(rrh).component<0>();
                auto rh2 = _g.topo(rrh).component<1>();
                if (labels.at(rh1) == provideLabels().rNotPlanar() || labels.at(rh2) == provideLabels().rNotPlanar()) { // one of them is nonplanar
                    continue;
                }
                int ent1 = rh2solverEntId.at(rh1);
                int ent2 = rh2solverEntId.at(rh2);
                auto & points = _g.data(rrh).normalizedSampledPoints;
                // register anchors
                std::vector<int> anchorIds;
                for (auto & ps : points) {
                    for (auto & p : ps) {
                        anchorIds.push_back(solver.makeNormalAnchor(p));
                    }
                }
                if (label == provideLabels().rrConnected()) { // connected
                    for (int aid : anchorIds) {
                        solver.makeAEqualToBAt(ent1, ent2, aid);
                    }
                } else if (label == provideLabels().rrFirstIsCloser()) { // 1 closer
                    for (int aid : anchorIds) {
                        solver.makeACloserThanBAt(ent1, ent2, aid);
                    }
                } else if (label == provideLabels().rrSecondIsCloser()) { // 2 closer
                    for (int aid : anchorIds) {
                        solver.makeAFartherThanDepthAt(ent1, ent2, aid);
                    }
                } else {
                    assert(false);
                }

                // depth lb/ub
                for (int aid : anchorIds) {
                    solver.makeACloserThanDepthAt(ent1, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent1, depthLB, aid);
                    solver.makeACloserThanDepthAt(ent2, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent2, depthLB, aid);
                }
            }

            // ll
            for (auto & c : labels.container<LLHandle>()) {
                int label = c.second;
                auto llh = c.first;
                auto lh1 = _g.topo(llh).component<0>();
                auto lh2 = _g.topo(llh).component<1>();
                int ent1 = lh2solverEntId.at(lh1);
                int ent2 = lh2solverEntId.at(lh2);
                auto & point = _g.data(llh).normalizedRelationCenter;
                // register anchor
                int anchorId = solver.makeNormalAnchor(point);
                if (label == provideLabels().llConnected()) { // connected
                    solver.makeAEqualToBAt(ent1, ent2, anchorId);
                } else if (label == provideLabels().llDisconnected()) { // disconnected
                } else {
                    assert(false);
                }

                // depth lb/ub
                {
                    int aid = anchorId;
                    solver.makeACloserThanDepthAt(ent1, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent1, depthLB, aid);
                    solver.makeACloserThanDepthAt(ent2, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent2, depthLB, aid);
                }
            }

            // rl, rrl
            std::vector<std::pair<RHandle, LHandle>> rlcons;
            std::vector<const std::vector<Vec3>*> nanchorsPtrTable;
            appliedRLCons(patch, labels, rlcons, nanchorsPtrTable);

            assert(rlcons.size() == nanchorsPtrTable.size());

            for (int i = 0; i < rlcons.size(); i++) {
                auto rh = rlcons[i].first;
                auto lh = rlcons[i].second;
                int ent1 = rh2solverEntId.at(rh);
                int ent2 = lh2solverEntId.at(lh);
                auto & as = *nanchorsPtrTable[i];
                for (auto & a : as) {
                    int aid = solver.makeNormalAnchor(a);
                    solver.makeAEqualToBAt(ent1, ent2, aid);

                    solver.makeACloserThanDepthAt(ent1, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent1, depthLB, aid);
                    solver.makeACloserThanDepthAt(ent2, depthUB, aid);
                    solver.makeAFartherThanDepthAt(ent2, depthLB, aid);
                }
            }

            // how to put arbitrary depth anchors!!!?
            {
                //// try longest line
                //LHandle longestLh;
                //double maxLen = 0.0;
                //for (auto & lh : patch.container<LHandle>()) {
                //    double len = _g.data(lh).normalizedLine.length();
                //    if (len > maxLen) {
                //        maxLen = len;
                //        longestLh = lh;
                //    }
                //}
                //if (longestLh.valid()) {
                //    // put a fixed depth anchor to its center
                //    auto center = normalize(_g.data(longestLh).normalizedLine.center());
                //    int aid = solver.makeNormalAnchor(center);
                //    solver.makeAEqualToDepthAt(lh2solverEntId.at(longestLh), depthAnchor, aid);
                //} else {
                //    // or try largest region
                //    RHandle largestRh;
                //    double maxArea = 0.0;
                //    for (auto & rh : patch.container<RHandle>()) {
                //        if (labels[rh] == provideLabels().rNotPlanar()) {
                //            continue;
                //        }
                //        double a = _g.data(rh).area;
                //        if (a > maxArea) {
                //            largestRh = rh;
                //            maxArea = a;
                //        }
                //    }
                //    // put a fixed depth anchor to its center
                //    auto center = normalize(_g.data(largestRh).normalizedCenter);
                //    int aid = solver.makeNormalAnchor(center);
                //    solver.makeAEqualToDepthAt(rh2solverEntId.at(largestRh), depthAnchor, aid);
                //}
                auto center = normalize(_g.data(patch.centerRh).normalizedCenter);
                int aid = solver.makeNormalAnchor(center);
                solver.makeAEqualToDepthAt(rh2solverEntId.at(patch.centerRh), depthAnchor, aid);
            }

            // solve
            double nanOrInfRatio = 0.0;
            bool solvable = solver.solve(&nanOrInfRatio, true);

            fea = RLFeatureVec(FeatureLength(Patch), 0.0);
            if (!solvable) {
                return nullptr;
            } else {
                // normalize resulted planes and lines
                double centerDepthsSum = 0.0;
                int centerDepthsNum = 0;
                for (auto & rp : planes) {
                    int label = labels[rp.first];
                    if (label == provideLabels().rNotPlanar()) { // non planar
                        // todo
                        continue;
                    }

                    auto & center = _g.data(rp.first).normalizedCenter;
                    const Plane3 & plane = rp.second;
                    double d = norm(IntersectionOfLineAndPlane(Ray3(Origin(), center), plane).position);
                    centerDepthsSum += d;
                    centerDepthsNum++;
                }
                double centerDepthMean = centerDepthsSum / centerDepthsNum;

                // manhattan region area ratio
                double manhattanRegionArea = 0.0;
                double allRegionArea = 0.0;
                for (auto & p : planes) {
                    int label = labels[p.first];
                    if (label == provideLabels().rNotPlanar()) { // non planar
                        // todo
                        continue;
                    }

                    double area = _g.data(p.first).area;
                    allRegionArea += area;

                    auto & plane = p.second;
                    static const double angleThreshold = DegreesToRadians(5);
                    double angle = M_PI;
                    int nearestVPId = -1;
                    for (int k = 0; k < _vps.size(); k++) {
                        double a = AngleBetweenUndirectedVectors(plane.normal, _vps[k]);
                        if (a < angle) {
                            angle = a;
                            nearestVPId = k;
                        }
                    }
                    if (angle < angleThreshold) {
                        manhattanRegionArea += area;
                    }
                }

                // gc confusion mat
                Mat<double, 2, 3> gcWeightedMeanConfusionMat;
                for (auto & p : planes) {
                    int label = labels[p.first];
                    if (label == provideLabels().rNotPlanar()) { // non planar
                        // todo
                        continue;
                    }

                    double area = _g.data(p.first).area;
                    auto & gcMean = _g.data(p.first).gcMean;

                    auto & plane = p.second;
                    double horizscore = abs(normalize(plane.normal).dot(normalize(up())));
                    double vertscore = 1.0 - horizscore;

                    double gcHoriz = gcMean[ToUnderlying(GeometricContextIndex::FloorOrGround)] +
                        gcMean[ToUnderlying(GeometricContextIndex::CeilingOrSky)];
                    double gcVert = gcMean[ToUnderlying(GeometricContextIndex::Vertical)];
                    double gcClutterOrPorous = gcMean[ToUnderlying(GeometricContextIndex::ClutterOrPorous)];

                    Mat<double, 2, 3> gcConfusionMat;

                    gcConfusionMat(0, 0) = horizscore * gcHoriz;
                    gcConfusionMat(0, 1) = horizscore * gcVert;
                    gcConfusionMat(0, 2) = horizscore * gcClutterOrPorous;

                    gcConfusionMat(1, 0) = vertscore * gcHoriz;
                    gcConfusionMat(1, 1) = vertscore * gcVert;
                    gcConfusionMat(1, 2) = vertscore * gcClutterOrPorous;

                    gcWeightedMeanConfusionMat += (area * gcConfusionMat);
                }
                gcWeightedMeanConfusionMat *= (1.0 / allRegionArea);

                // plane orientation/position confusion mat
                // corresponding to the (Note this!)THREE principle vanishing points
                Mat<double, 3, 3> planeWeightedMeanOPConfusionMat;
                assert(_vps.size() >= 3);
                for (auto & p : planes) {
                    int label = labels[p.first];
                    if (label == provideLabels().rNotPlanar()) { // non planar
                        // todo
                        continue;
                    }

                    double area = _g.data(p.first).area;
                    auto & center = _g.data(p.first).normalizedCenter;

                    auto & plane = p.second;
                    double normalDotsToVPs[3];
                    for (int i = 0; i < 3; i++)
                        normalDotsToVPs[i] = abs(normalize(plane.normal).dot(normalize(_vps[i])));
                    double centerDotsToVPs[3];
                    for (int i = 0; i < 3; i++)
                        centerDotsToVPs[i] = abs(center.dot(normalize(_vps[i])));
                    Mat<double, 3, 3> confusionMat;
                    for (int i = 0; i < 3; i++) {
                        for (int j = 0; j < 3; j++) {
                            confusionMat(i, j) = normalDotsToVPs[i] * centerDotsToVPs[j];
                        }
                    }
                    planeWeightedMeanOPConfusionMat += (area * confusionMat);
                }
                planeWeightedMeanOPConfusionMat *= (1.0 / allRegionArea);


                // region connectivity
                double rrDistanceSum = 0;
                int appliedRRAnchorsNum = 0;
                for (auto & c : labels.container<RRHandle>()) {
                    int label = c.second;
                    auto rrh = c.first;
                    auto rh1 = _g.topo(rrh).component<0>();
                    auto rh2 = _g.topo(rrh).component<1>();
                    if (labels.at(rh1) == provideLabels().rNotPlanar() || labels.at(rh2) == provideLabels().rNotPlanar()) { // one of them is nonplanar
                        continue;
                    }

                    auto & plane1 = planes.at(rh1);
                    auto & plane2 = planes.at(rh2);
                    auto & points = _g.data(rrh).normalizedSampledPoints;
                    for (auto & ps : points) {
                        for (auto & p : ps) {
                            double d1 = norm(IntersectionOfLineAndPlane(Ray3(Origin(), p), plane1).position);
                            double d2 = norm(IntersectionOfLineAndPlane(Ray3(Origin(), p), plane2).position);
                            rrDistanceSum += Gaussian(abs(d1 - d2), 1.0);
                            appliedRRAnchorsNum++;
                        }
                    }
                }
                double rrDistanceMean = rrDistanceSum / std::max(appliedRRAnchorsNum, 1);

                // region line connectivity
                double rlDistanceSum = 0.0;
                int appliedRLAnchorsNum = 0;
                for (int i = 0; i < rlcons.size(); i++) {
                    auto rh = rlcons[i].first;
                    auto lh = rlcons[i].second;
                    const Plane3 & plane = planes.at(rh);
                    const Line3 & line = lines.at(lh);
                    auto & anchors = *nanchorsPtrTable[i];
                    for (auto & a : anchors) {
                        double d1 = norm(IntersectionOfLineAndPlane(Ray3(Origin(), a), plane).position);
                        double d2 = norm(DistanceBetweenTwoLines(Ray3(Origin(), a), line.ray()).second.first);
                        rlDistanceSum += Gaussian(abs(d1 - d2), 1.0);
                        appliedRLAnchorsNum++;
                    }
                }
                double rlDistanceMean = rlDistanceSum / std::max(appliedRLAnchorsNum, 1);

                int allAnchorsNum = 0;
                for (auto & c : patch.container<RRHandle>()) {
                    for (auto & ps : _g.data(c).normalizedSampledPoints) {
                        allAnchorsNum += ps.size();
                    }
                }
                for (auto & c : patch.container<RLHandle>()) {
                    allAnchorsNum += _g.data(c).normalizedAnchors.size();
                }
                for (auto & c : patch.container<RRLHandle>()) {
                    allAnchorsNum += _g.data(c).normalizedAnchors.size() * 2;
                }

                int appliedAnchorsNum = appliedRRAnchorsNum + appliedRLAnchorsNum;
                int notAppliedAnchors = allAnchorsNum - appliedAnchorsNum;
                assert(notAppliedAnchors >= 0);

                static_assert(FeatureLength(Patch) == 21, "");
                fea = RLFeatureVec(FeatureLength(Patch), 0.0);
                fea[0] = 0.0; // 0
                fea[1] = manhattanRegionArea / allRegionArea; // 1
                std::copy(gcWeightedMeanConfusionMat.val, gcWeightedMeanConfusionMat.val + 6, fea.begin() + 2); // 2:7
                std::copy(planeWeightedMeanOPConfusionMat.val, planeWeightedMeanOPConfusionMat.val + 9, fea.begin() + 8); // 8:16
                fea[17] = double(appliedAnchorsNum) / std::max(allAnchorsNum, 1); // 17
                fea[18] = double(notAppliedAnchors) / std::max(allAnchorsNum, 1); // 18
                fea[19] = rrDistanceMean / centerDepthMean; // 19
                fea[20] = rlDistanceMean / centerDepthMean; // 20

                for (auto & f : fea) {
                    f /= _patches.size();
                }
            }

            return std::move(feaData);
        }

        void RLOptimizer::findPeakyRhs(double rangeAngle) {
            _peakyRhs.clear();
            _peakyRhs.resize(_vps.size());

            // find peaky regions
            for (auto & r : _g.components<RData>()) {
                auto h = r.topo.hd;

                auto & contours = r.data.normalizedContours;
                double radiusAngle = 0.0;
                for (auto & cs : r.data.normalizedContours) {
                    for (auto & c : cs) {
                        double angle = AngleBetweenDirections(r.data.normalizedCenter, c);
                        if (angle > radiusAngle) {
                            radiusAngle = angle;
                        }
                    }
                }

                bool mayCrossAnyVP = false;
                for (auto & vp : _vps) {
                    double angle = AngleBetweenUndirectedVectors(vp, r.data.normalizedCenter);
                    if (angle < radiusAngle) {
                        mayCrossAnyVP = true;
                        break;
                    }
                }

                if (!mayCrossAnyVP) {
                    continue;
                }

                static const double focal = 100.0;
                auto maskView = PerfectRegionMaskView(r.data.normalizedContours, r.data.normalizedCenter, focal);

                // intersection test
                for (int i = 0; i < _vps.size(); i++) {
                    auto p1 = ToPixel(maskView.camera.toScreen(_vps[i]));
                    auto p2 = ToPixel(maskView.camera.toScreen(-_vps[i]));

                    int dilateSize = focal * rangeAngle;
                    bool intersected = false;
                    for (int x = -dilateSize; x <= dilateSize; x++) {
                        if (intersected)
                            break;
                        for (int y = -dilateSize; y <= dilateSize; y++) {
                            if (intersected)
                                break;
                            auto pp1 = Pixel(p1.x + x, p1.y + y);
                            auto pp2 = Pixel(p2.x + x, p2.y + y);
                            if (Contains(maskView.image, pp1) && maskView.image(pp1)) {
                                _peakyRhs[i].insert(r.topo.hd);
                                intersected = true;
                            } else if (Contains(maskView.image, pp2) && maskView.image(pp2)) {
                                _peakyRhs[i].insert(r.topo.hd);
                                intersected = true;
                            }
                        }
                    }
                }
            }
        }

        void RLOptimizer::findHorizonRhs(double rangeAngle) {
            _horizonRhs.clear();
            for (auto & r : _g.components<RData>()) {
                auto h = r.topo.hd;
                auto & contours = r.data.normalizedContours;
                bool intersected = false;
                for (auto & cs : r.data.normalizedContours) {
                    if (intersected)
                        break;
                    for (auto & c : cs) {
                        double angle = M_PI_2 - AngleBetweenUndirectedVectors(c, up());
                        if (angle <= rangeAngle) {
                            intersected = true;
                            break;
                        }
                    }
                }
                if (intersected) {
                    _horizonRhs.insert(h);
                }
            }
        }



        void RLOptimizer::appliedRLCons(const RLGraphPatch & patch, 
            const RLGraphPatchDict<RLGraphLabel> & cl,
            std::vector<std::pair<RHandle, LHandle>> & rlcons, std::vector<const std::vector<Vec3>*> & nanchorsPtrTable) const {

            rlcons.clear();
            rlcons.reserve(patch.container<RLHandle>().size() + patch.container<RRLHandle>().size());

            nanchorsPtrTable.clear();
            nanchorsPtrTable.reserve(rlcons.capacity());

            for (auto & c : patch.container<RLHandle>()) {
                auto rlh = c;
                auto rh = _g.topo(rlh).component<0>();
                if (cl.at(rh) == provideLabels().rNotPlanar()) { // rh is nonplanar
                    continue;
                }
                auto lh = _g.topo(rlh).component<1>();
                rlcons.emplace_back(rh, lh);
                nanchorsPtrTable.push_back(&_g.data(rlh).normalizedAnchors);
            }

            for (auto & c : patch.container<RRLHandle>()) {
                auto rrlh = c;
                auto rh1 = _g.topo(rrlh).component<0>();
                auto rh2 = _g.topo(rrlh).component<1>();
                auto lh = _g.topo(rrlh).component<2>();
                // find boundary between rh1 and rh2
                RRHandle rrh;
                for (const RRHandle & h : patch.container<RRHandle>()) {
                    if (_g.topo(h).component<0>() == rh1 && _g.topo(h).component<1>() == rh2) {
                        rrh = h;
                        break;
                    }
                }
                assert(rrh.valid());

                int rlabel1 = cl.at(rh1);
                int rlabel2 = cl.at(rh2);
                bool r1isnonplanar = rlabel1 == provideLabels().rNotPlanar();
                bool r2isnonplanar = rlabel2 == provideLabels().rNotPlanar();
                if (r1isnonplanar && r2isnonplanar) { // all of them is nonplanar
                    continue;
                }

                int rrlabel = cl.at(rrh);
                bool connect_rh1_lh = (rrlabel == provideLabels().rrConnected() || rrlabel == provideLabels().rrFirstIsCloser()) && !r1isnonplanar;
                bool connect_rh2_lh = (rrlabel == provideLabels().rrConnected() || rrlabel == provideLabels().rrSecondIsCloser()) && !r2isnonplanar;

                if (connect_rh1_lh) {
                    rlcons.emplace_back(rh1, lh);
                    nanchorsPtrTable.push_back(&_g.data(rrlh).normalizedAnchors);
                }
                if (connect_rh2_lh) {
                    rlcons.emplace_back(rh2, lh);
                    nanchorsPtrTable.push_back(&_g.data(rrlh).normalizedAnchors);
                }
            }

        }


        std::vector<double> RLOptimizer::suggestedWeights() const {
            NOT_IMPLEMENTED_YET();
        }



        void RLOptimizer::inference(const std::vector<double> & weights,
            HandledTable<RHandle, Plane3> & planes, HandledTable<LHandle, Line3> & lines) const {
            assert(weights.size() == featureLength());
            auto predictedLabels = _fg.solveWithSimpleCallback(1e2, 10, nullptr, (void*)weights.data());

            
            // whole graph as a cluster

            // todo

            // todo

        }

        void RLOptimizer::inferenceWithLossFunction(const std::vector<double> & weights, 
            std::function<double(const HandledTable<RHandle, Plane3> & planes, const HandledTable<LHandle, Line3> & lines)> & lossFun) const {


        }

    }
}