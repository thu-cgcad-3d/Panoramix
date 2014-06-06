#ifndef PANORAMIX_CORE_SERIALIZATION_HPP
#define PANORAMIX_CORE_SERIALIZATION_HPP

#include <fstream>

#include <cereal/types/array.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/bitset.hpp>
#include <cereal/types/chrono.hpp>
#include <cereal/types/common.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/types/deque.hpp>
#include <cereal/types/forward_list.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/queue.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/stack.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/unordered_set.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>


#include <cereal/archives/binary.hpp>
#include <cereal/archives/portable_binary.hpp>

#include <opencv2/opencv.hpp>

namespace cv {

    // MUST be defined in the namespace of the underlying type (cv::XXX), 
    //    definition of alias names in namespace panoramix::core won't work!
    // see http://stackoverflow.com/questions/13192947/argument-dependent-name-lookup-and-typedef

    // Serialization for cv::Mat
    template <class Archive>
    void save(Archive & ar, Mat const & im) {
        ar(im.elemSize(), im.type(), im.cols, im.rows);
        ar(cereal::binary_data(im.data, im.cols * im.rows * im.elemSize()));
    }

    // Serialization for cv::Mat
    template <class Archive>
    void load(Archive & ar, Mat & im) {
        size_t elemSize;
        int type, cols, rows;
        ar(elemSize, type, cols, rows);
        im.create(rows, cols, type);
        ar(cereal::binary_data(im.data, cols * rows * elemSize));
    }

    // Serialization for cv::Matx<T, M, N>
    template <class Archive, class T, int M, int N>
    inline void serialize(Archive & ar, Matx<T, M, N> & m) {
        ar(m.val);
    }

    // Serialization for cv::Point_<T>
    template <class Archive, class T>
    inline void serialize(Archive & ar, Point_<T> & p) {
        ar(p.x, p.y);
    }

    // Serialization for cv::Size_<T>
    template <class Archive, class T>
    inline void serialize(Archive & ar, Size_<T> & s) {
        ar(s.width, s.height);
    }

    // Serialization for cv::KeyPoint
    template <class Archive>
    inline void serialize(Archive & ar, KeyPoint & p) {
        ar(p.pt, p.size, p.angle, p.response, p.octave, p.class_id);
    }

    // Serialization for cv::Moments
    template <class Archive>
    inline void serialize(Archive & ar, Moments & m) {
        ar(m.m00, m.m10, m.m01, m.m20, m.m11, m.m02, m.m30, m.m21, m.m12, m.m03);
        ar(m.mu20, m.mu11, m.mu02, m.mu30, m.mu21, m.mu12, m.mu03);
        ar(m.nu20, m.nu11, m.nu02, m.nu30, m.nu21, m.nu12, m.nu03);
    }

}

#endif