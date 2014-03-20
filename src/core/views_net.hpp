#ifndef PANORAMIX_CORE_VIEWS_NET_HPP
#define PANORAMIX_CORE_VIEWS_NET_HPP
 
#include <memory>

#include "basic_types.hpp"
#include "feature.hpp"
#include "mesh.hpp"

namespace panoramix {
	namespace core { 

		// perspective view
		class PerspectiveView {
		public:
			using Camera = PerspectiveCamera;
			struct FeatureConfigures {
				FeatureConfigure<LineSegmentExtractor> lineSegmentConfig;
				FeatureConfigure<CVFeatureExtractor<cv::SIFT>> siftConfig;
				FeatureConfigure<CVFeatureExtractor<cv::SURF>> surfConfig;
				FeatureConfigure<ManhattanJunctionExtractor> manhattanJunctionConfig;
				FeatureConfigure<GeometricContextExtractor> geometricContextConfig;
			};

		public:
			inline explicit PerspectiveView(const Image & image, const PerspectiveCamera & camera)
				: _image(image), _camera(camera), _correctedCamera(camera) {}

			inline const Image & image() const { return _image; }
			inline const Camera & camera() const { return _camera; }
			inline const Camera & correctedCamera() const { return _correctedCamera; }
			inline Camera & correctedCamera() { return _correctedCamera; }

			inline const FeatureConfigures & featureConfigures() const { return _featureConfigs; }
			inline FeatureConfigures & featureConfigures() { return _featureConfigs; }

		private:
			const Image _image;
			const Camera _camera;

			Camera _correctedCamera;
			FeatureConfigures _featureConfigs;
		};

		// panoramic view
		class PanoramicView {
		public:
			using Camera = PanoramicCamera;


			inline const Image & image() const { return _image; }
			inline const Camera & camera() const { return _camera; }

		public:
			const Image _image;
			const Camera _camera;
		};


		// relation between a pair of views
		template <class View1T, class View2T = View1T>
		class ViewRelation {

		};


		// views net
		class ViewsNet {
		public:
			using ViewMesh = Mesh<std::shared_ptr<PerspectiveView>, std::shared_ptr<ViewRelation<PerspectiveView>>>;
			using ViewHandle = ViewMesh::VertHandle;
			using ViewRelationHandle = ViewMesh::HalfHandle;

		public:
			ViewsNet(){}
			
			

		private:
			ViewMesh _views;

		};

 
	}
}
 
#endif