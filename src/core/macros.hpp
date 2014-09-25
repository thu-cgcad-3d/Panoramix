#ifndef PANORAMIX_CORE_MACROS_HPP
#define PANORAMIX_CORE_MACROS_HPP
 
namespace panoramix {
    namespace core {

        template <class ...T> 
        struct DelayStaticAssert {
            enum { value = false };
        };

        // not implemented error
#define NOT_IMPLEMENTED_YET() \
    throw std::runtime_error("This feature has not yet been implemented! \n" \
    "in function: "__FUNCSIG__ "\n" \
    "in line: " + std::to_string(__LINE__) + "\n" \
    "in file: "__FILE__)

        // not tested warning
#define NOT_TESTED_YET() \
    std::cout << ("This feature has not yet been tested! \n" \
    "in function: "__FUNCSIG__ "\n" \
    "in line: " + std::to_string(__LINE__) + "\n" \
    "in file: "__FILE__) << std::endl 

        // should never be called error
#define SHOULD_NEVER_BE_CALLED() \
    throw std::runtime_error("This feature should never be called! \n" \
    "in function: "__FUNCSIG__ "\n" \
    "in line: " + std::to_string(__LINE__) + "\n" \
    "in file: "__FILE__)

        // should never be instanciated error
#define SHOULD_NEVER_BE_INSTANCIATED(...) \
    static_assert(panoramix::core::DelayStaticAssert<__VA_ARGS__>::value, "This feature should never be instanciated by compiler! \n" \
    "in function: "__FUNCSIG__ "\n" \
    "in line: " + std::to_string(__LINE__) + "\n" \
    "in file: "__FILE__)


    }
}
 
#endif