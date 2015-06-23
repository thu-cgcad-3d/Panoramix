#ifndef PANORAMIX_GUI_UTILITY_HPP
#define PANORAMIX_GUI_UTILITY_HPP

#include "../core/cameras.hpp"
#include "basic_types.hpp"
#include "scene.hpp"

class QWidget;

namespace panoramix {
    namespace gui {

        int SelectFrom(const std::vector<std::string> & strs,
            const std::string & title = std::string(),
            const std::string & text = std::string(),
            int acceptId = -1, int rejectId = -1);

        core::Image PickAnImage(const std::string & dir = std::string(), std::string * picked = nullptr);

        std::vector<core::Image> PickImages(const std::string & dir = std::string());

        std::vector<core::Image> PickAllImagesFromAFolder(const std::string & dir = std::string());

        void PaintWith(const std::function<core::Image()> & updater,
            const std::vector<PenConfig> & penConfigs,
            const std::function<bool(const std::vector<core::Point2> & polyline, int penId)> & callback);


        void VisualizeWithPanoramicOperation(const Scene & scene, const RenderOptions & options);


        void VisualizeAll(const core::View<core::PanoramicCamera, core::Image3ub> & view,
            const std::vector<core::Classified<core::Line3>> & lines,
            const core::Imagei & segs, int nsegs,
            const core::Image5d & gc);

        void DrawChainsInPanorama(const core::PanoramicView & view,
            const std::vector<PenConfig> & penConfigs,
            std::vector<core::Chain3> & chains);
        

    }
}

#endif