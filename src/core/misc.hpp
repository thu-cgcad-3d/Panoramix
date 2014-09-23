#ifndef PANORAMIX_CORE_MISC_HPP
#define PANORAMIX_CORE_MISC_HPP

#include <stack>
#include <iterator>
#include <cassert>
#include <map>
#include <set>
 
namespace panoramix {
    namespace core {    


        template <class ContainerT, class ValueT>
        inline bool Contains(const ContainerT & c, const ValueT & v) {
            return std::find(c.cbegin(), c.cend(), v) != c.cend();
        }

        template <class KeyT, class ValueT, class PredT, class AllocT>
        inline bool Contains(const std::map<KeyT, ValueT, PredT, AllocT> & m, const KeyT & k) {
            return m.find(k) != m.end();
        }

        template <class KeyT, class PredT, class AllocT>
        inline bool Contains(const std::set<KeyT, PredT, AllocT> & m, const KeyT & k) {
            return m.find(k) != m.end();
        }

 
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


        // make an iterator
        template <class CoreDataT, class GetValueT, class SetToNextT, class CompareCoreDataT = std::equal_to<CoreDataT>>
        class EasyForwardIterator : public std::iterator<std::forward_iterator_tag, 
            decltype(std::declval<GetValueT>()(std::declval<CoreDataT>()))> {
        public:
            inline explicit EasyForwardIterator(const CoreDataT & cdata, GetValueT getValue, SetToNextT setToNext, 
                CompareCoreDataT cmpCData = CompareCoreDataT())
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
            MakeEasyForwardIterator(const CoreDataT & cdata, GetValueT getValue, SetToNextT setToNext, CompareCoreDataT cmpCData = CompareCoreDataT()) {
            return EasyForwardIterator<CoreDataT, GetValueT, SetToNextT, CompareCoreDataT>(cdata, getValue, setToNext, cmpCData);
        }
       
    }
}
 
#endif