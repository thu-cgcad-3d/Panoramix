#ifndef PANORAMIX_CORE_FEATURE_VISUALIZE_HPP
#define PANORAMIX_CORE_FEATURE_VISUALIZE_HPP

#include "basic_types.hpp"

namespace panoramix {
	namespace core {
 
		class ImageFeatureVisualizer {
		public:
			using Image = cv::Mat;
			inline explicit ImageFeatureVisualizer(const Image & im, 
				const std::string & name = "Image Feature Visualizer") 
				: _name(name) {
				im.copyTo(_image);
			}
			inline ~ ImageFeatureVisualizer() { 
				cv::imshow(_name, _image); 
				cv::waitKey(); 
			}
			inline Image & image() { return _image; }
		private:
			Image _image;
			const std::string _name;
		};

		inline ImageFeatureVisualizer & operator << (ImageFeatureVisualizer & viz, const Line2 & line) {
			cv::line(viz.image(), cv::Point(line.first(0), line.first(1)),
				cv::Point(line.second(0), line.second(1)),
				cv::Scalar(255, 255, 255), 2);
			return viz;
		}

		inline ImageFeatureVisualizer & operator << (ImageFeatureVisualizer & viz, const KeyPoint & p) {
			cv::circle(viz.image(), p.pt, p.size, cv::Scalar(255, 200, 200));
			return viz;
		}

		inline ImageFeatureVisualizer & operator << (ImageFeatureVisualizer & viz, const Image & im) {
			viz.image() = viz.image() * 0.75 + im * 0.25;
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