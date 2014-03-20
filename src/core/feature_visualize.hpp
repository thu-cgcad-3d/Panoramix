#ifndef PANORAMIX_CORE_FEATURE_VISUALIZE_HPP
#define PANORAMIX_CORE_FEATURE_VISUALIZE_HPP

#include "basic_types.hpp"

namespace panoramix {
	namespace core {
 
		class ImageFeatureVisualizer {
		public:
			using Image = cv::Mat;
			struct Params {
				inline Params() 
					: winName("Image Feature Visualizer"), 
					  color(255, 255, 255), thickness(1), 
					  lineType(8), shift(0), alphaForNewImage(0.5) {}

				std::string winName;
				Color color;
				int thickness;
				int lineType;
				int shift;
				float alphaForNewImage;
			};
		public:
			inline explicit ImageFeatureVisualizer(const Image & im, 
				const Params & params = Params()) 
				: _params(params) {
				im.copyTo(_image);
			}
			inline ~ ImageFeatureVisualizer() { 
				cv::imshow(_params.winName, _image); 
				cv::waitKey(); 
			}
			inline Image & image() { return _image; }
			inline Params & params() { return _params; }
			inline const Params & params() const { return _params; }
		private:
			Image _image;
			Params _params;
		};

		template <class T>
		inline ImageFeatureVisualizer & operator << (ImageFeatureVisualizer & viz, const Line<T, 2> & line) {
			cv::line(viz.image(), cv::Point(line.first(0), line.first(1)),
				cv::Point(line.second(0), line.second(1)),
				viz.params().color, viz.params().thickness, viz.params().lineType, viz.params().shift);
			return viz;
		}

		inline ImageFeatureVisualizer & operator << (ImageFeatureVisualizer & viz, const KeyPoint & p) {
			cv::circle(viz.image(), p.pt, p.size, viz.params().color, viz.params().thickness, viz.params().lineType, viz.params().shift);
			return viz;
		}

		inline ImageFeatureVisualizer & operator << (ImageFeatureVisualizer & viz, const Image & im) {
			viz.image() = viz.image() * (1 - viz.params().alphaForNewImage) + im * viz.params().alphaForNewImage;
			return viz;
		}

		template <class T, int N>
		inline ImageFeatureVisualizer & operator << (ImageFeatureVisualizer & viz, const std::array<T, N> & a) {
			for (auto & e : c)
				viz << e;
			return viz;
		}

		template <class ContainerT>
		inline ImageFeatureVisualizer & VisualizeAllInContainer(ImageFeatureVisualizer & viz, const ContainerT & c) {
			for (auto & e : c)
				viz << e;
			return viz;
		}

		#define VISUALIZE_AS_CONTAINER(claz) \
			template <class T> \
			inline ImageFeatureVisualizer & operator << (ImageFeatureVisualizer & viz, const claz<T> & c) { \
				return VisualizeAllInContainer(viz, c); \
			}

		VISUALIZE_AS_CONTAINER(std::list)
		VISUALIZE_AS_CONTAINER(std::vector)
		VISUALIZE_AS_CONTAINER(std::deque)
		VISUALIZE_AS_CONTAINER(std::set)
		VISUALIZE_AS_CONTAINER(std::unordered_set)
		VISUALIZE_AS_CONTAINER(std::forward_list)
 
	}
}
 
#endif