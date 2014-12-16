#ifndef PANORAMIX_CORE_MISC_HPP
#define PANORAMIX_CORE_MISC_HPP

#include <stack>
#include <iterator>
#include <cassert>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <iostream>
#include <string>

#include "meta.hpp"
 
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


        class CastError : public std::exception {
        public:
            inline explicit CastError(const std::string & message)
                : std::exception(message.c_str()) {}
            inline explicit CastError(const char * message)
                : std::exception(message) {}
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
                if (_data->type() == typeid(T))
                    return reinterpret_cast<Data<T>*>(_data)->value;
                throw CastError("type not matched! "
                    "intrinsic type is '" + std::string(_data->type().name()) +
                    "', while target type is '" + typeid(T).name() + "'");
            }

            template <class T>
            inline T & uncheckedRef() const {
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

        template <class ...T>
        using AnyOfTypes = Any;



        //// ITERATORS
        template <class IteratorT, class T>
        struct IsIteratorOfType 
            : std::is_same<typename std::iterator_traits<IteratorT>::value_type, T> {
        };

 
        // element of container MUST support PredT(ele) -> bool
        // ConditionalIterator will automatically skip elements which DO NOT satisfy PredT in iteration
        template <class IteratorT, class PredT>
        class ConditionalIterator : public std::iterator<std::forward_iterator_tag, 
            typename std::iterator_traits<IteratorT>::value_type,
            typename std::iterator_traits<IteratorT>::difference_type,
            typename std::iterator_traits<IteratorT>::pointer,
            typename std::iterator_traits<IteratorT>::reference> {

        public:
            using Iterator = IteratorT;

            inline ConditionalIterator(Iterator it_, Iterator end_, PredT pred_ = PredT())
                : _it(it_), _end(end_), _pred(pred_) {
                if (_it != _end && !_pred(*_it))
                    ++ (*this);
            }

            inline ConditionalIterator & operator++() {
                assert(_it != _end);
                ++_it;
                while (_it != _end && !_pred(*_it))
                    ++_it;
                return *this;
            }

            inline reference operator * () const {
                return *_it;
            }

            inline pointer operator -> () const {
                return &(*_it);
            }

            inline bool operator == (const ConditionalIterator & i) const {
                return _it == i._it;
            }

            inline bool operator != (const ConditionalIterator & i) const {
                return !(*this == i);
            }

            inline Iterator internalIterator() const {
                return _it;
            }

        private:
            IteratorT _it;
            IteratorT _end;
            PredT _pred;
        };


        // class ConditionalContainerWrapper
        template <class ContainerT, class ElementPredT>
        class ConditionalContainerWrapper {
        public:
            using OriginalIterator = typename ContainerT::iterator;
            using iterator = ConditionalIterator<OriginalIterator, ElementPredT>;
            using value_type = typename std::iterator_traits<iterator>::value_type;

            inline ConditionalContainerWrapper(ContainerT * cont_, ElementPredT elePred_ = ElementPredT()) 
                : _cont(cont_), _elePred(elePred_){}
            inline iterator begin() { return iterator(std::begin(*_cont), std::end(*_cont), _elePred); }
            inline iterator end() { return iterator(std::end(*_cont), std::end(*_cont), _elePred); }
            inline iterator begin() const { return iterator(std::begin(*_cont), std::end(*_cont), _elePred); }
            inline iterator end() const { return iterator(std::end(*_cont), std::end(*_cont), _elePred); }

        private:
            ContainerT * _cont;
            ElementPredT _elePred;
        };


        // class ConstConditionalContainerWrapper
        template <class ContainerT, class ElementPredT>
        class ConstConditionalContainerWrapper {
        public:
            using OriginalIterator = typename ContainerT::const_iterator;
            using iterator = ConditionalIterator<OriginalIterator, ElementPredT>;
            using value_type = typename std::iterator_traits<iterator>::value_type;

            inline ConstConditionalContainerWrapper(const ContainerT * cont_, ElementPredT elePred_ = ElementPredT()) 
                : _cont(cont_), _elePred(elePred_){}
            inline iterator begin() const { return iterator(std::begin(*_cont), std::end(*_cont), _elePred); }
            inline iterator end() const { return iterator(std::end(*_cont), std::end(*_cont), _elePred); }

        private:
            const ContainerT * _cont;
            ElementPredT _elePred;
        };

 
        // make conditional container
        template <class ContainerT, class ElementPredT>
        inline ConditionalContainerWrapper<ContainerT, ElementPredT> 
            MakeConditionalContainer(ContainerT * cont_, ElementPredT elePred_ = ElementPredT()) {
            return ConditionalContainerWrapper<ContainerT, ElementPredT>(cont_, elePred_);
        }

        template <class ContainerT, class ElementPredT>
        inline ConstConditionalContainerWrapper<ContainerT, ElementPredT> 
            MakeConditionalContainer(const ContainerT * cont_, ElementPredT elePred_ = ElementPredT()) {
            return ConstConditionalContainerWrapper<ContainerT, ElementPredT>(cont_, elePred_);
        }





        // yield
        template <class T, class ProcessorT>
        class YieldIterator : public std::iterator<std::output_iterator_tag, T> {
        public:
            struct Wrapper {
                inline explicit Wrapper(const ProcessorT & p) : processor(p) {}
                inline explicit Wrapper(ProcessorT && p) : processor(std::move(p)) {}
                inline Wrapper & operator = (const T & data){
                    processor(data);
                    return *this;
                }
                ProcessorT processor;
            };
        public:
            inline explicit YieldIterator(const ProcessorT & p) : _w(p) {}
            inline explicit YieldIterator(ProcessorT && p) : _w(std::move(p)) {}
            inline YieldIterator & operator ++() {
                return *this;
            }
            inline Wrapper & operator * () {
                return _w;
            }
        private:
            Wrapper _w;            
        };

        template <class T, class ProcessorT>
        inline YieldIterator<T, std::decay_t<ProcessorT>> MakeYield(ProcessorT && p) {
            return YieldIterator<T, std::decay_t<ProcessorT>>(std::forward<ProcessorT>(p));
        }






        // make an iterator
        template <class CoreDataT, class GetValueT, class SetToNextT, class CompareCoreDataT = std::equal_to<CoreDataT>>
        class EasyForwardIterator : public std::iterator<std::forward_iterator_tag, 
            decltype(std::declval<GetValueT>()(std::declval<CoreDataT>()))> {
        public:
            inline explicit EasyForwardIterator(const CoreDataT & cdata, const GetValueT & getValue, const SetToNextT & setToNext, 
                const CompareCoreDataT & cmpCData = CompareCoreDataT())
                : _coreData(cdata), _getValue(getValue), _setToNext(setToNext), _compareCoreData(cmpCData) {}
            inline EasyForwardIterator & operator++() {
                _setToNext(_coreData);
                return *this;
            }

            inline value_type operator * () const {
                return _getValue(_coreData);
            }

            inline pointer operator -> () const {
                return &(_getValue(_coreData));
            }

            inline bool operator == (const EasyForwardIterator & i) const {
                return _compareCoreData(_coreData, i._coreData);
            }

            inline bool operator != (const EasyForwardIterator & i) const {
                return !(*this == i);
            }

            inline const CoreDataT & coreData() const {
                return _coreData;
            }
        private:
            CoreDataT _coreData;
            GetValueT _getValue;
            SetToNextT _setToNext;
            CompareCoreDataT _compareCoreData;
        };

        // make easy forward iterator
        template <class CoreDataT, class GetValueT, class SetToNextT, class CompareCoreDataT = std::equal_to<CoreDataT>>
        inline EasyForwardIterator<CoreDataT, GetValueT, SetToNextT, CompareCoreDataT> 
            MakeEasyForwardIterator(const CoreDataT & cdata, const GetValueT & getValue, const SetToNextT & setToNext, 
            const CompareCoreDataT & cmpCData = CompareCoreDataT()) {
            return EasyForwardIterator<CoreDataT, GetValueT, SetToNextT, CompareCoreDataT>(cdata, getValue, setToNext, cmpCData);
        }




        // tick tock
        inline std::pair<std::chrono::system_clock::time_point, std::string> 
            Tick(const std::string & taskName){
            return std::make_pair(std::chrono::system_clock::now(), taskName);
        }

        inline std::chrono::system_clock::duration 
            Tock(std::pair<std::chrono::system_clock::time_point, std::string> startInfo) {
            auto duration = std::chrono::system_clock::now() - startInfo.first;
            std::cout << "[" << startInfo.second << "] Time Elapsed: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() 
                << " ms" << std::endl;
            return duration;
        }
        
    }
}


namespace std {

    inline void swap(panoramix::core::Any & a, panoramix::core::Any & b) { a.swap(b); }

}
 
#endif