#ifndef PANORAMIX_CORE_RESULT_HPP
#define PANORAMIX_CORE_RESULT_HPP

#include <utility>
 
namespace panoramix {
    namespace core {
    
        // failable result
        template <class T>
        class Result {
        public:            
            Result() : _exist(false) {}
            
            Result(nullptr_t) : _exist(false) {}
            Result & operator = (nullptr_t) {
                _exist = false;
                return *this;
            }
            
            Result(T && v) : _exist(true), _temp(std::move(v)) {}
            Result & operator = (T && v) {
                _temp = std::move(v);
                _exist = true;
                return *this;
            }

            Result(Result && opt) { swap(opt); }
            Result & operator = (Result && opt) {
                swap(opt);
                return *this;
            }

            Result(const Result &) = delete;
            Result & operator = (const Result &) = delete;

            inline void swap(Result & opt) {
                std::swap(_exist, opt._exist);
                std::swap(_temp, opt._temp);
            }

            inline bool null() const { return !_exist; }
            inline bool failed() const { return null(); }
            inline const T & ref() const { return _temp; }

            inline T unwrap() {
                if (!_exist){
                    std::cout << "unwrapping a null Result<T>!!!!!!!" << std::endl;
                }
                assert(_exist);
                _exist = false;
                return std::move(_temp);
            }
            template <class TT>
            inline T unwrap(TT && defaultValue) { 
                bool e = _exist;
                _exist = false;
                return e ? std::move(_temp) : std::forward<TT>(defaultValue);
            }

        private:
            bool _exist;
            T _temp;
        };

        template <class T, class = std::enable_if_t<!std::is_reference<T>::value>>
        inline Result<T> AsResult(T && v){
            return Result<T>(std::move(v));
        }


    }
}

namespace std {

    template <class T>
    void swap(panoramix::core::Result<T> & a, panoramix::core::Result<T> & b){
        a.swap(b);
    }

}

 
#endif