#ifndef PANORAMIX_SANDBOX_FOR_EXPRESSION_HPP
#define PANORAMIX_SANDBOX_FOR_EXPRESSION_HPP

#include <iostream>
 
namespace panoramix {
    namespace sandbox {

        int foo();

#define DEFINE_CATEGORY(Checker, TagStructName) \
        struct TagStructName {}; \
        template <class T> \
        std::enable_if_t<Checker, TagStructName> TAG() { return TagStructName(); }

        template <class T, class Tag = decltype(TAG<T>())>
        struct Traits {};

        template <class T>
        struct DataStorageStruct {
            using type = typename Traits<T>::StorageType;
        };
        template <class T>
        using DataStorageType = typename DataStorageStruct<T>::type;


        


    }
}
 
#endif