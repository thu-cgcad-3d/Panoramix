#ifndef PANORAMIX_CORE_ANY_HPP
#define PANORAMIX_CORE_ANY_HPP

#include <typeinfo>
#include <string>
#include <exception>
#include <cassert>

namespace panoramix {
    namespace core {


        // class Any
        struct DataBase {
            virtual DataBase * clone() const = 0;
            virtual const std::type_info & type() const = 0;
        };

        template <class T>
        struct Data : DataBase {
            inline Data() {}
            inline Data(const T & d) : value(d) {}
            inline Data(T && d) : value(std::move(d)) {}

            virtual DataBase * clone() const override { return new Data(value); }
            virtual const std::type_info & type() const override { return typeid(T); }

            template <class Archive> inline void serialize(Archive & ar) { ar(value); }

            T value;
        };


        class Any {
        public:
            inline Any() : _data(nullptr) {}

            // from class Any
            inline Any(const Any & a) : _data(a._data->clone()) {}
            inline Any & operator = (const Any & a) {
                if (this == &a) return *this;
                delete _data;
                _data = a._data->clone();
                return *this;
            }
            //inline Any(Any && a) { swap(a);  ASSERTVALID; }
            inline Any & operator = (Any && a) { swap(a); return *this; }
            inline void swap(Any & a){ std::swap(_data, a._data); }

            // from other types
            template <class T, class = std::enable_if_t<!std::is_same<T, Any>::value>> // accepts const T(&)
            inline Any(const T & v) : _data(new Data<T>(v)) {}
            template <class T, class = std::enable_if_t<!std::is_same<std::decay_t<T>, Any>::value>> // accepts T&& and T&
            inline Any(T && v) : _data(new Data<std::decay_t<T>>(std::forward<T>(v))) {}

            // set to nullptr
            inline Any(nullptr_t) : _data(nullptr) {}
            inline Any & operator = (nullptr_t) {
                if (_data == nullptr) return *this;
                delete _data;
                _data = nullptr;
                return *this;
            }

            ~Any() {
                delete _data;
                _data = nullptr;
            }

            inline bool null() const { return _data == nullptr; }
            inline bool operator == (nullptr_t) const { return _data == nullptr; }
            inline bool operator != (nullptr_t) const { return _data != nullptr; }

            template <class T>
            inline T & ref() const {
                assert(_data->type() == typeid(T));
                return reinterpret_cast<Data<T>*>(_data)->value;
            }

            template <class T>
            inline operator T() const { return ref<T>(); }
            template <class T>
            inline T cast() const { return ref<T>(); }

            template <class T>
            inline bool is() const { return _data->type() == typeid(T); }

        private:
            DataBase * _data;
        };

   	}
}


namespace std {

    inline void swap(panoramix::core::Any & a, panoramix::core::Any & b) { a.swap(b); }

}

 
#endif