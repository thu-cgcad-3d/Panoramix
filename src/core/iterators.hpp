#ifndef PANORAMIX_CORE_ITERATORS_HPP
#define PANORAMIX_CORE_ITERATORS_HPP

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


        //// ITERATORS
        template <class IteratorT, class T>
        struct IsIteratorOfType 
            : std::is_same<typename std::iterator_traits<IteratorT>::value_type, T> {
        };

        template <class ContainerT, class T>
        struct IsContainerOfType
            : IsIteratorOfType<decltype(std::begin(std::declval<ContainerT>())), T> {
        };


        // range
        template <class IteratorT>
        struct Range {
            IteratorT b, e;
            Range(IteratorT bb, IteratorT ee) : b(bb), e(ee) {}
            template <class ContainerT>
            explicit Range(ContainerT && cont) : b(std::begin(cont)), e(std::end(cont)) {}
            IteratorT begin() const { return b; }
            IteratorT end() const { return e; }
            
            template <class FunT>
            inline void forEach(FunT && fun) const {
                IteratorT i = b;
                while (i != e){
                    fun(*i);
                    ++i;
                }
            }
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



        //template <class WalkerT>
        //class WalkerContainerWrapper {
        //public:
        //    using reference_type = decltype(std::declval<NexterT>().next());
        //    using value_type = std::decay_t<reference_type>;
        //    struct iterator : public std::iterator<std::forward_iterator_tag, value_type> {
        //        const bool isEnd;
        //        WalkerT nexter;
        //        reference_type * curValuePtr;
        //        explicit iterator(bool e, WalkerT n) : isEnd(e), nexter(n), curValuePtr(nullptr) {}
        //        inline reference_type operator *() const { return *curValuePtr; }
        //        inline reference_type * operator -> () const { return curValuePtr; }
        //        inline iterator & operator ++() { curValuePtr = &(nexter.next()); }
        //        inline bool operator == (iterator it) const { 
        //            assert(isEnd || it.isEnd);
        //            return isEnd ? it.nexter.hasNext() : nexter.hasNext();
        //        }
        //    };
        //    
        //    explicit WalkerContainerWrapper(NexterT nexter) : _nexter(nexter) {}
        //    inline iterator begin() { return iterator(false, _nexter); }
        //    inline iterator end() { return iterator(true, _nexter); }

        //private:
        //    WalkerT _nexter;
        //};

        //template <class WalkerT>
        //inline WalkerContainerWrapper<WalkerT> WrapWalker(const WalkerT & w) {
        //    return WalkerContainerWrapper<WalkerT>(w);
        //}

        
    }
}

#endif