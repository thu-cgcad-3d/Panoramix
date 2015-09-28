#pragma once

#include "pi_graph_solve.hpp"
#include "pi_graph_annotation.hpp"

namespace pano {
    namespace experimental {

        struct PIEvaluation {
            double avgNormalViolation;
            double avgDepthViolation;
        };

        //PIEvaluation Evaluate(const PILayoutAnnotation & anno, const PanoramicCamera & cam, const std::function<>)
        PIEvaluation Evaluate(const PILayoutAnnotation & anno, const PerspectiveCamera & cam, const PIGraph & mg);

        PIEvaluation Evaluate(const PILayoutAnnotation & anno, const View<PerspectiveCamera, Imagei> & omaps);
        PIEvaluation Evaluate(const PILayoutAnnotation & anno, const View<PerspectiveCamera, Image7d> & gc);

    }
}
