#ifndef PANORAMIX_CORE_DATA_TRAITS_HPP
#define PANORAMIX_CORE_DATA_TRAITS_HPP

#include <type_traits>
#include <iterator>

namespace panoramix {
    namespace deriv {

        // data traits 
        // use this macro to define a catogory tag
#define DEFINE_DATA_CATEGORY(Checker, TagStructName) \
        struct TagStructName {}; \
        template <class T> \
        std::enable_if_t<(Checker), TagStructName> TAG() { return TagStructName(); }
#define SATISFIES(TypeName, TagStructName) \
    std::enable_if_t<std::is_same<TagStructName, decltype(TAG<TypeName>())>::value, int> = 0

        namespace {
            template <class T>
            struct TagStruct {
                using type = decltype(TAG<T>());
            };
        }
        template <class T>
        using TagType = typename TagStruct<T>::type;

        // use DataTraits<T> to retrieve traits of data
        template <class T, class Tag = TagType<T>>
        struct DataTraits {};

        namespace {
            template <class T>
            struct DataStorageStruct {
                using type = typename DataTraits<T>::StorageType;
            };
            template <class T>
            struct DataScalarStruct {
                using type = typename DataTraits<T>::ScalarType;
            };
        }
        
        // storage type
        template <class T>
        using DataStorageType = typename DataStorageStruct<T>::type;

        // scalar type
        template <class T>
        using DataScalarType = typename DataScalarStruct<T>::type;

        template <class T>
        struct IsStorageType : public std::is_same<T, DataStorageType<T>> {};

        template <class T>
        struct IsScalarType : public std::is_same<DataStorageType<T>, DataScalarType<T>> {};

        // role in product
        enum class RoleInProduct {
            Scalar,
            Array,
            Matrix
        };

#define SATISFIES_SCALAR(T) std::enable_if_t<DataTraits<T>::roleInProduct == RoleInProduct::Scalar, int> = 0
#define SATISFIES_ARRAY(T) std::enable_if_t<DataTraits<T>::roleInProduct == RoleInProduct::Array, int> = 0
#define SATISFIES_MATRIX(T) std::enable_if_t<DataTraits<T>::roleInProduct == RoleInProduct::Matrix, int> = 0

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


        struct ResultRetrievedByValueTag {};
        struct ResultRetrievedByCacheTag {};
        template <class T>
        struct ResultTag {
            using type =
                std::conditional_t < std::is_lvalue_reference<T>::value, ResultRetrievedByValueTag,
                std::conditional_t < DataTraits<T>::shouldBeCached, ResultRetrievedByCacheTag,
                ResultRetrievedByValueTag >> ;
        };

    }
}
 
#endif