#ifndef PANORAMIX_CORE_PREDICTOR_HPP
#define PANORAMIX_CORE_PREDICTOR_HPP

#include <memory>

#include "basic_types.hpp"
#include "feature.hpp"
#include "mesh.hpp"

namespace panoramix {
	namespace core {

		class PhotoPairAnalyzer;
 

		class SinglePhotoAnalyzer {
		public:
			inline explicit SinglePhotoAnalyzer(const Image & image, const PerspectiveCamera & camera)
				: _image(image), _camera(camera) {}

			void analyze();

		private:
			Image _image;
			PerspectiveCamera _camera;

			LineSegmentExtractor::Feature _lineSegments;
			std::shared_ptr<LineSegmentExtractor> _lineSegmentExtractor;
		
		private:
			friend class PhotoPairAnalyzer;
		};

		class PhotoPairAnalyzer {

		};

		class PhotoDatabase {
		public:
			using PhotoNet = Mesh<std::shared_ptr<SinglePhotoAnalyzer>, std::shared_ptr<PhotoPairAnalyzer>>;
			using SinglePhotoHandle = PhotoNet::VertHandle;
			using PhotoPairHandle = PhotoNet::HalfHandle;

		public:
			PhotoDatabase(){}
			
			

		private:
			PhotoNet _photoConnections;
		};

		
		class Calibrator {

		};

		class Stitcher {

		};

		class Predictor {

		};

 
	}
}
 
#endif