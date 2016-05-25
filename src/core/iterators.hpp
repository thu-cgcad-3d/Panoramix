#pragma once

#include <cassert>
#include <chrono>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "meta.hpp"

namespace pano {
namespace core {

//// ITERATORS
template <class IterT, class T>
struct IsIteratorOfType
    : std::is_same<typename std::iterator_traits<IterT>::value_type, T> {};

template <class ContainerT, class T>
struct IsContainerOfType
    : IsIteratorOfType<decltype(std::begin(std::declval<ContainerT>())), T> {};

// IsIterator
template <class T, class = void> struct IsIterator : no {};
template <class T>
struct IsIterator<T, void_t<typename T::iterator_category,
                            typename T::value_type, typename T::difference_type,
                            typename T::pointer, typename T::reference>> : yes {
};
template <class T> struct IsIterator<T *> : yes {};

// Range
template <class IterT> struct Range {
  IterT b, e;
  Range(IterT bb, IterT ee) : b(bb), e(ee) {}
  template <class ContainerT>
  explicit Range(ContainerT &&cont) : b(std::begin(cont)), e(std::end(cont)) {}
  IterT begin() const { return b; }
  IterT end() const { return e; }

  template <class FunT> void forEach(FunT &&fun) const {
    IterT i = b;
    while (i != e) {
      fun(*i);
      ++i;
    }
  }

  bool operator==(const Range &r) const { return b == r.b && e == r.e; }
  bool operator!=(const Range &r) const { return !(*this == r); }
};

// MakeRange
template <class IterT> Range<IterT> MakeRange(IterT b, IterT e) {
  return Range<IterT>(b, e);
}

// TransformIterator
template <class T, class IterT, class FunT>
class TransformIterator
    : public std::iterator<
          typename std::iterator_traits<IterT>::iterator_category, T,
          typename std::iterator_traits<IterT>::difference_type> {
  static_assert(IsIterator<IterT>::value, "IterT must be an iterator type");

public:
  using difference_type = typename std::iterator_traits<IterT>::difference_type;

  constexpr TransformIterator() : current(), fun() {}
  template <class F>
  constexpr TransformIterator(IterT it, F &&f)
      : current(it), fun(std::forward<F>(f)) {}
  constexpr IterT base() const { return current; }

  decltype(auto) operator*() const { return fun(*current); }
  decltype(auto) operator-> () const { return &(**this); }

  TransformIterator &operator++() {
    ++current;
    return *this;
  }
  TransformIterator &operator++(int) {
    auto tmp = *this;
    ++current;
    return tmp;
  }
  TransformIterator &operator--() {
    --current;
    return *this;
  }
  TransformIterator &operator--(int) {
    auto tmp = *this;
    --current;
    return tmp;
  }

  TransformIterator &operator+=(difference_type off) { // increment by integer
    current += off;
    return (*this);
  }
  TransformIterator
  operator+(difference_type off) const { // return this + integer
    return (TransformIterator(current + off, fun));
  }

  TransformIterator &operator-=(difference_type off) { // decrement by integer
    current -= off;
    return (*this);
  }
  TransformIterator
  operator-(difference_type off) const { // return this - integer
    return (TransformIterator(current - off, fun));
  }

protected:
  IterT current;
  FunT fun;
};

template <class T, class IterT, class FunT>
constexpr bool operator==(const TransformIterator<T, IterT, FunT> &i1,
                          const TransformIterator<T, IterT, FunT> &i2) {
  return i1.base() == i2.base();
}
template <class T, class IterT, class FunT>
constexpr bool operator!=(const TransformIterator<T, IterT, FunT> &i1,
                          const TransformIterator<T, IterT, FunT> &i2) {
  return !(i1 == i2);
}
template <class T, class IterT, class FunT>
constexpr bool operator<(const TransformIterator<T, IterT, FunT> &i1,
                         const TransformIterator<T, IterT, FunT> &i2) {
  return i1.base() < i2.base();
}
template <class T, class IterT, class FunT>
constexpr auto operator-(const TransformIterator<T, IterT, FunT> &i1,
                         const TransformIterator<T, IterT, FunT> &i2) {
  return i1.base() - i2.base();
}

// MakeTransformIterator
template <class IterT, class FunT>
constexpr auto MakeTransformIterator(IterT it, FunT f) {
  using value_t = std::decay_t<decltype(f(*it))>;
  return TransformIterator<value_t, IterT, FunT>(it, f);
}

// ConcatedIterator
template <class T, class IterT1, class IterT2>
class ConcatedIterator : public std::iterator<std::forward_iterator_tag, T> {
  static_assert(
      std::is_same<typename std::iterator_traits<IterT1>::value_type,
                   typename std::iterator_traits<IterT2>::value_type>::value,
      "Value types mismatched!");

public:
  constexpr ConcatedIterator(IterT1 i1, IterT1 e1, IterT2 i2, IterT2 e2)
      : iter1(i1), end1(e1), iter2(i2), end2(e2) {}

  ConcatedIterator &operator++() {
    if (iter1 == end1) {
      assert(iter2 != end2);
      ++iter2;
    } else {
      ++iter1;
    }
    return *this;
  }

  T &operator*() const { return iter1 == end1 ? *iter2 : *iter1; }
  T *operator->() const { return iter1 == end1 ? &(*iter2) : &(*iter1); }
  bool operator==(const ConcatedIterator &i) const {
    return std::make_tuple(iter1, end1, iter2, end2) ==
           std::make_tuple(i.iter1, i.end1, i.iter2, i.end2);
  }
  bool operator!=(const ConcatedIterator &i) const { return !(*this == i); }

protected:
  IterT1 iter1, end1;
  IterT2 iter2, end2;
};

// MakeConcatedRange
template <class IterT1, class IterT2>
constexpr auto MakeConcatedRange(IterT1 b1, IterT1 e1, IterT2 b2, IterT2 e2) {
  using ValueT = typename std::iterator_traits<IterT1>::value_type;
  return MakeRange(ConcatedIterator<ValueT, IterT1, IterT2>(b1, e1, b2, e2),
                   ConcatedIterator<ValueT, IterT1, IterT2>(e1, e1, e2, e2));
}

// element of container MUST support PredT(ele) -> bool
// ConditionalIterator will automatically skip elements which DO NOT satisfy
// PredT in iteration
template <class IterT, class PredT>
class ConditionalIterator
    : public std::iterator<
          std::forward_iterator_tag,
          typename std::iterator_traits<IterT>::value_type,
          typename std::iterator_traits<IterT>::difference_type,
          typename std::iterator_traits<IterT>::pointer,
          typename std::iterator_traits<IterT>::reference> {

public:
  using Iterator = IterT;
  ConditionalIterator(Iterator it_, Iterator end_, PredT pred_ = PredT())
      : _it(it_), _end(end_), _pred(pred_) {
    if (_it != _end && !_pred(*_it))
      ++(*this);
  }
  ConditionalIterator &operator++() {
    assert(_it != _end);
    ++_it;
    while (_it != _end && !_pred(*_it))
      ++_it;
    return *this;
  }

  reference operator*() const { return *_it; }
  pointer operator->() const { return &(*_it); }
  bool operator==(const ConditionalIterator &i) const { return _it == i._it; }
  bool operator!=(const ConditionalIterator &i) const { return !(*this == i); }
  Iterator internalIterator() const { return _it; }

private:
  IterT _it;
  IterT _end;
  PredT _pred;
};

// MakeConditionalRange
template <class IterT, class PredT>
constexpr auto MakeConditionalRange(IterT b, IterT e, PredT pred) {
  return MakeRange(ConditionalIterator<IterT, PredT>(b, e, pred),
                   ConditionalIterator<IterT, PredT>(e, e, pred));
}


//// MakeConditionalIterator
//template <class IterT, class PredT>
//constexpr auto MakeConditionalIterator(IterT iter, IterT end, PredT pred) {
//  return ConditionalIterator<IterT, PredT>(iter, end, pred);
//}

//// class ConditionalContainerWrapper
//template <class ContainerT, class ElementPredT>
//class ConditionalContainerWrapper {
//public:
//  using OriginalIterator = typename ContainerT::iterator;
//  using iterator = ConditionalIterator<OriginalIterator, ElementPredT>;
//  using value_type = typename std::iterator_traits<iterator>::value_type;
//
//  ConditionalContainerWrapper(ContainerT *cont_,
//                              ElementPredT elePred_ = ElementPredT())
//      : _cont(cont_), _elePred(elePred_) {}
//  iterator begin() {
//    return iterator(std::begin(*_cont), std::end(*_cont), _elePred);
//  }
//  iterator end() {
//    return iterator(std::end(*_cont), std::end(*_cont), _elePred);
//  }
//  iterator begin() const {
//    return iterator(std::begin(*_cont), std::end(*_cont), _elePred);
//  }
//  iterator end() const {
//    return iterator(std::end(*_cont), std::end(*_cont), _elePred);
//  }
//
//private:
//  ContainerT *_cont;
//  ElementPredT _elePred;
//};
//
//// class ConstConditionalContainerWrapper
//template <class ContainerT, class ElementPredT>
//class ConstConditionalContainerWrapper {
//public:
//  using OriginalIterator = typename ContainerT::const_iterator;
//  using iterator = ConditionalIterator<OriginalIterator, ElementPredT>;
//  using value_type = typename std::iterator_traits<iterator>::value_type;
//
//  ConstConditionalContainerWrapper(const ContainerT *cont_,
//                                   ElementPredT elePred_ = ElementPredT())
//      : _cont(cont_), _elePred(elePred_) {}
//  iterator begin() const {
//    return iterator(std::begin(*_cont), std::end(*_cont), _elePred);
//  }
//  iterator end() const {
//    return iterator(std::end(*_cont), std::end(*_cont), _elePred);
//  }
//
//private:
//  const ContainerT *_cont;
//  ElementPredT _elePred;
//};
//
//// make conditional container
//template <class ContainerT, class ElementPredT>
//ConditionalContainerWrapper<ContainerT, ElementPredT>
//MakeConditionalContainer(ContainerT *cont, ElementPredT elePred) {
//  return ConditionalContainerWrapper<ContainerT, ElementPredT>(cont, elePred);
//}
//
//template <class ContainerT, class ElementPredT>
//ConstConditionalContainerWrapper<ContainerT, ElementPredT>
//MakeConditionalContainer(const ContainerT *cont, ElementPredT elePred) {
//  return ConstConditionalContainerWrapper<ContainerT, ElementPredT>(cont,
//                                                                    elePred);
//}
}
}
