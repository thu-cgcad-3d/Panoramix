#ifndef PANORAMIX_CORE_OPTIONAL_HPP
#define PANORAMIX_CORE_OPTIONAL_HPP

#include <utility>
 
namespace panoramix {
    namespace core {
    
        // optional
        template <class T>
        class Optional {
        public:            
            Optional() : _exist(false) {}
            
            Optional(nullptr_t) : _exist(false) {}
            Optional & operator = (nullptr_t) {
                _exist = false;
                return *this;
            }
            
            Optional(T && v) : _exist(true), _temp(std::move(v)) {}
            Optional & operator = (T && v) {
                _temp = std::move(v);
                _exist = true;
                return *this;
            }

            Optional(Optional && opt) { swap(opt); }
            Optional & operator = (Optional && opt) {
                swap(opt);
                return *this;
            }

            Optional(const Optional &) = delete;
            Optional & operator = (const Optional &) = delete;

            inline void swap(Optional & opt) {
                std::swap(_exist, opt._exist);
                std::swap(_temp, opt._temp);
            }

            inline bool null() const { return !_exist; }
            inline const T & ref() const { return _temp; }

            inline T unwrap() {
                if (!_exist){
                    std::cout << "unwrapping a null Optional<T>!!!!!!!" << std::endl;
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
        inline Optional<T> AsOptional(T && v){
            return Optional<T>(std::move(v));
        }


    }
}

namespace std {

    template <class T>
    void swap(panoramix::core::Optional<T> & a, panoramix::core::Optional<T> & b){
        a.swap(b);
    }

}

 
#endif