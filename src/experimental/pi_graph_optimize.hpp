#pragma once

#include "pi_graph_solve.hpp"

namespace pano {
    namespace experimental {



        // get the CompactModel
        std::vector<Polygon3> CompactModel(const PICGDeterminablePart & dp, const PIConstraintGraph & cg, const PIGraph & mg, double distThres);





        std::vector<Vec3> ComputeSegNormals(const PICGDeterminablePart & dp, const PIConstraintGraph & cg, 
            const PIGraph & mg, bool smoothed);


        template <class CameraT>
        Image3d SurfaceNormalMap(const CameraT & cam, const PICGDeterminablePart & dp, const PIConstraintGraph & cg,
            const PIGraph & mg, bool smoothed) {
            auto seg2normal = ComputeSegNormals(dp, cg, mg, smoothed);
            Image3d snm(cam.screenSize());
            for (auto it = snm.begin(); it != snm.end(); ++it) {
                auto p = it.pos();
                auto dir = normalize(cam.toSpace(p));
                auto pp = ToPixel(mg.view.camera.toScreen(dir));
                pp.x = WrapBetween(pp.x, 0, mg.segs.cols);
                pp.y = BoundBetween(pp.y, 0, mg.segs.rows - 1);
                int seg = mg.segs(pp);
                *it = seg2normal[seg];
            }
            return snm;
        }


        template <class CameraT>
        Imaged SurfaceDepthMap(const CameraT & cam, const PICGDeterminablePart & dp, const PIConstraintGraph & cg,
            const PIGraph & mg, std::pair<double, double> * validDepthRange = nullptr) {
            Imaged depths(cam.screenSize(), 0.0);
            double minv = std::numeric_limits<double>::max();
            double maxv = 0.0;
            for (auto it = depths.begin(); it != depths.end(); ++it) {
                auto pos = it.pos();
                int seg = mg.segs(pos);
                int ent = cg.seg2ent[seg];
                if (ent == -1) {
                    continue;
                }
                if (!Contains(dp.determinableEnts, ent)) {
                    continue;
                }
                auto & plane = cg.entities[ent].supportingPlane.reconstructed;
                Vec3 dir = normalize(cam.toSpace(pos));
                double depth = norm(IntersectionOfLineAndPlane(Ray3(Origin(), dir), plane).position);
                if (depth < minv) {
                    minv = depth;
                }
                if (depth > maxv) {
                    maxv = depth;
                }
                *it = depth;
            }
            if (validDepthRange) {
                validDepthRange->first = minv;
                validDepthRange->second = maxv;
            }
            return depths;
        }

    }
}