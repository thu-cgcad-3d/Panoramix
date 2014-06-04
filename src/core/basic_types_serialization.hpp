#ifndef PANORAMIX_CORE_BASIC_TYPES_SERIALIZATION_HPP
#define PANORAMIX_CORE_BASIC_TYPES_SERIALIZATION_HPP

#include <fstream>

#include <cereal/types/common.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/string.hpp>

#include <cereal/archives/binary.hpp>

#include "basic_types.hpp"

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

}

namespace panoramix {
    namespace core {

        template <class Archive, class T, int N>
        inline void serialize(Archive & ar, HPoint<T, N> & p) {
            ar(p.coord, p.scalar);
        }

        template <class Archive>
        inline void serialize(Archive & ar, GeoCoord & gc) {
            ar(gc.longitude, gc.latitude);
        }

        template <class Archive, class T, int N>
        inline void serialize(Archive & ar, Line<T, N> & l) {
            ar(l.first, l.second);
        }

        template <class Archive, class T, int N>
        inline void serialize(Archive & ar, PositionOnLine<T, N> & p) {
            ar(p.ratio, p.position);
        }

        template <class Archive, class T, int N>
        inline void serialize(Archive & ar, HLine<T, N> & l) {
            ar(l.first, l.second);
        }

        template <class Archive, class T>
        inline void serialize(Archive & ar, Classified<T> & c) {
            ar(c.claz, c.component);
        }

        template <class Archive, class T, int N>
        inline void serialize(Archive & ar, Box<T, N> & b) {
            ar(b.isNull, b.minCorner, b.maxCorner);
        }

    }
}


#endif