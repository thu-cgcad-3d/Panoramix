#ifndef PANORAMIX_CORE_DATA_TRAITS_HPP
#define PANORAMIX_CORE_DATA_TRAITS_HPP

#include <type_traits>

namespace panoramix {
    namespace deriv {

        // data traits 
        // use this macro to define a catogory tag
#define DEFINE_DATA_CATEGORY(Checker, TagStructName) \
        struct TagStructName {}; \
        template <class T> \
        std::enable_if_t<(Checker), TagStructName> TAG() { return TagStructName(); }

        namespace {
            template <class T>
            struct TagStruct {
                using type = decltype(TAG<T>());
            };
        }
        template <class T>
        using TagType = typename TagStruct<T>::type;

        // use DataTraits<T> to retrieve traits
        template <class T, class Tag = TagType<T>>
        struct DataTraits {};

        namespace {
            template <class T>
            struct DataStorageStruct {
                using type = typename DataTraits<T>::StorageType;
            };
        }
        template <class T>
        using DataStorageType = typename DataStorageStruct<T>::type;

        namespace {
            template <class T>
            struct DataScalarStruct {
                using type = typename DataTraits<T>::ScalarType;
            };
        }
        template <class T>
        using DataScalarType = typename DataScalarStruct<T>::type;

        template <class T>
        struct IsStorageType : public std::is_same<T, DataStorageType<T>> {};

        template <class T>
        struct IsScalarType : public std::is_same<DataStorageType<T>, DataScalarType<T>> {};

        namespace  {
            template<class T> struct RemoveAll { typedef T type; };
            template<class T> struct RemoveAll<T&> : RemoveAll<T>{};
            template<class T> struct RemoveAll<T&&> : RemoveAll<T>{};
            template<class T> struct RemoveAll<T const> : RemoveAll<T>{};
            template<class T> struct RemoveAll<T volatile> : RemoveAll<T>{};
            template<class T> struct RemoveAll<T const volatile> : RemoveAll<T>{};
        };


        template <class T>
        using RemoveAllType = typename RemoveAll<T>::type;






    }
}
 
#endif