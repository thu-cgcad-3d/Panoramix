#pragma once

#include "handle.hpp"

namespace pano {
namespace core {

// helper classes for HomogeneousGraph
template <class Tag, int L, int ChildN> struct Topo {
  static const int Level = L;
  enum { ChildrenNum = ChildN };

  std::array<HandleOfTypeAtLevel<Tag, Level - 1>, ChildN>
      lowers; // use std::array
  std::set<HandleOfTypeAtLevel<Tag, Level + 1>> uppers;
  HandleOfTypeAtLevel<Tag, Level> hd;

  explicit inline Topo(int id = -1) : hd(id) {}
  explicit inline Topo(
      int id, std::initializer_list<HandleOfTypeAtLevel<Tag, Level - 1>> ls)
      : hd(id) {
    assert(ls.size() == lowers.size());
    std::copy(ls.begin(), ls.end(), lowers.begin());
  }
  template <class Archive> inline void serialize(Archive &ar) {
    ar(lowers, uppers, hd);
  }
};

// dynamic sized
template <class Tag, int L> struct Topo<Tag, L, Dynamic> {
  static const int Level = L;
  enum { ChildrenNum = Dynamic };

  std::vector<HandleOfTypeAtLevel<Tag, Level - 1>> lowers; // use std::array
  std::set<HandleOfTypeAtLevel<Tag, Level + 1>> uppers;
  HandleOfTypeAtLevel<Tag, Level> hd;

  explicit inline Topo(int id = -1) : hd(id) {}
  explicit inline Topo(
      int id, std::initializer_list<HandleOfTypeAtLevel<Tag, Level - 1>> ls)
      : hd(id), lowers(ls) {}

  template <class IteratorT>
  explicit inline Topo(int id, IteratorT lsBegin, IteratorT lsEnd)
      : hd(id), lowers(lsBegin, lsEnd) {}

  template <class Archive> inline void serialize(Archive &ar) {
    ar(lowers, uppers, hd);
  }
};

// zero sized
template <class Tag, int L> struct Topo<Tag, L, 0> {
  static const int Level = L;
  enum { ChildrenNum = 0 };

  std::set<HandleOfTypeAtLevel<Tag, Level + 1>> uppers;
  HandleOfTypeAtLevel<Tag, Level> hd;

  explicit inline Topo(int id = -1) : hd(id) {}
  explicit inline Topo(
      int id, std::initializer_list<HandleOfTypeAtLevel<Tag, Level - 1>> ls)
      : hd(id) {
    assert(ls.size() == 0);
  }
  template <class Archive> inline void serialize(Archive &ar) {
    ar(uppers, hd);
  }
};

// configuration for each layer
template <class DataT, int ChildN> struct LayerConfig {
  using DataType = DataT;
  enum { ChildNum = ChildN };
};

template <class T> struct IsLayerConfig : public std::false_type {};

template <class DataT, int ChildN>
struct IsLayerConfig<LayerConfig<DataT, ChildN>> : public std::true_type {};

namespace {

// content of each layer
template <class Tag, int Level, class T> struct LayerContent {};

template <class Tag, int Level, class DataT, int ChildN>
struct LayerContent<Tag, Level, LayerConfig<DataT, ChildN>> {
  using TopoType = Topo<Tag, Level, ChildN>;
  using DataType = DataT;
  using TripletType = Triplet<TopoType, DataT>;
  using TableType = TripletArray<TopoType, DataT>;
  TableType contentTable;
  template <class Archive> inline void serialize(Archive &ar) {
    ar(contentTable);
  }
};

// make a layer content tuple with levels
template <class Tag, class CConfTupleT, int Level>
struct LayerContentFromConfigTuple {
  using type =
      LayerContent<Tag, Level, typename std::tuple_element<size_t(Level),
                                                           CConfTupleT>::type>;
};

template <class Tag, class CConfTupleT, class SequenceT>
struct LayerContentTuple {};

template <class Tag, class CConfTupleT, int... Levels>
struct LayerContentTuple<Tag, CConfTupleT, Sequence<Levels...>> {
  using type = std::tuple<
      typename LayerContentFromConfigTuple<Tag, CConfTupleT, Levels>::type...>;
};
}

// Homogeneous Graph
template <class BaseDataT, class... CConfs> class HomogeneousGraph {
  using Tag = HomogeneousGraph;

  // number of layers
  static const unsigned LayerNum = sizeof...(CConfs) + 1;
  // tuple of all layers
  using ContentsType = typename LayerContentTuple<
      Tag, std::tuple<LayerConfig<BaseDataT, 0>, CConfs...>,
      typename SequenceGenerator<LayerNum>::type>::type;
  // layer content type at level Level
  template <int Level> struct LayerContentTypeStruct {
    using type = typename std::tuple_element<Level, ContentsType>::type;
  };
  // predicate of element existance at level Level
  template <int Level>
  using ElementExistsPred =
      TripletExistsPred<typename LayerContentTypeStruct<Level>::type::TopoType,
                        typename LayerContentTypeStruct<Level>::type::DataType>;

public:
  // element triplet type at level Level
  template <int Level>
  using TripletType =
      typename LayerContentTypeStruct<Level>::type::TableType::value_type;

  // internal elements
  template <int Level>
  inline const typename LayerContentTypeStruct<Level>::type::TableType &
  internalElements() const {
    return std::get<Level>(_contents).contentTable;
  }

  template <int Level>
  inline typename LayerContentTypeStruct<Level>::type::TableType &
  internalElements() {
    return std::get<Level>(_contents).contentTable;
  }

  // topo
  template <int Level>
  inline const typename LayerContentTypeStruct<Level>::type::TopoType &
  topo(HandleOfTypeAtLevel<Tag, Level> h) const {
    return internalElements<Level>()[h.id].topo;
  }

  template <int Level>
  inline typename LayerContentTypeStruct<Level>::type::TopoType &
  topo(HandleOfTypeAtLevel<Tag, Level> h) {
    return internalElements<Level>()[h.id].topo;
  }

  // data
  template <int Level>
  inline const typename LayerContentTypeStruct<Level>::type::DataType &
  data(HandleOfTypeAtLevel<Tag, Level> h) const {
    return internalElements<Level>()[h.id].data;
  }

  template <int Level>
  inline typename LayerContentTypeStruct<Level>::type::DataType &
  data(HandleOfTypeAtLevel<Tag, Level> h) {
    return internalElements<Level>()[h.id].data;
  }

  // traverse
  template <int Level>
  inline ConditionalContainerWrapper<
      typename LayerContentTypeStruct<Level>::type::TableType,
      ElementExistsPred<Level>>
  elements() {
    return ConditionalContainerWrapper<
        typename LayerContentTypeStruct<Level>::type::TableType,
        ElementExistsPred<Level>>(&(internalElements<Level>()));
  }

  template <int Level>
  inline ConstConditionalContainerWrapper<
      typename LayerContentTypeStruct<Level>::type::TableType,
      ElementExistsPred<Level>>
  elements() const {
    return ConstConditionalContainerWrapper<
        typename LayerContentTypeStruct<Level>::type::TableType,
        ElementExistsPred<Level>>(&(internalElements<Level>()));
  }

  // add element
  template <int Level>
  HandleOfTypeAtLevel<Tag, Level>
      add(std::initializer_list<HandleOfTypeAtLevel<Tag, Level - 1>> depends,
          const typename LayerContentTypeStruct<Level>::type::DataType &d =
              typename LayerContentTypeStruct<Level>::type::DataType()) {
    int id = static_cast<int>(internalElements<Level>().size());
    internalElements<Level>().emplace_back(
        typename LayerContentTypeStruct<Level>::type::TopoType(id, depends), d,
        true);
    for (const HandleOfTypeAtLevel<Tag, Level - 1> &lowh : depends) {
      topo(lowh).uppers.insert(HandleOfTypeAtLevel<Tag, Level>(id));
    }
    return HandleOfTypeAtLevel<Tag, Level>(id);
  }

  template <int Level>
  HandleOfTypeAtLevel<Tag, Level>
      add(std::initializer_list<HandleOfTypeAtLevel<Tag, Level - 1>> depends,
          typename LayerContentTypeStruct<Level>::type::DataType &&d) {
    int id = static_cast<int>(internalElements<Level>().size());
    internalElements<Level>().emplace_back(
        typename LayerContentTypeStruct<Level>::type::TopoType(id, depends),
        std::move(d), true);
    for (const HandleOfTypeAtLevel<Tag, Level - 1> &lowh : depends) {
      topo(lowh).uppers.insert(HandleOfTypeAtLevel<Tag, Level>(id));
    }
    return HandleOfTypeAtLevel<Tag, Level>(id);
  }

  template <int Level, class IteratorT>
  HandleOfTypeAtLevel<Tag, Level>
  add(IteratorT dependsBegin, IteratorT dependsEnd,
      const typename LayerContentTypeStruct<Level>::type::DataType &d) {
    int id = static_cast<int>(internalElements<Level>().size());
    internalElements<Level>().emplace_back(
        typename LayerContentTypeStruct<Level>::type::TopoType(id, dependsBegin,
                                                               dependsEnd),
        d, true);
    while (dependsBegin != dependsEnd) {
      topo(*dependsBegin).uppers.insert(HandleOfTypeAtLevel<Tag, Level>(id));
      ++dependsBegin;
    }
    return HandleOfTypeAtLevel<Tag, Level>(id);
  }

  template <int Level, class IteratorT>
  HandleOfTypeAtLevel<Tag, Level>
  add(IteratorT dependsBegin, IteratorT dependsEnd,
      typename LayerContentTypeStruct<Level>::type::DataType &&d) {
    int id = static_cast<int>(internalElements<Level>().size());
    internalElements<Level>().emplace_back(
        typename LayerContentTypeStruct<Level>::type::TopoType(id, dependsBegin,
                                                               dependsEnd),
        std::move(d), true);
    while (dependsBegin != dependsEnd) {
      topo(*dependsBegin).uppers.insert(HandleOfTypeAtLevel<Tag, Level>(id));
      ++dependsBegin;
    }
    return HandleOfTypeAtLevel<Tag, Level>(id);
  }

  // add element of the lowest level
  HandleOfTypeAtLevel<Tag, 0>
  add(const typename LayerContentTypeStruct<0>::type::DataType &d =
          LayerContentTypeStruct<0>::type::DataType()) {
    int id = static_cast<int>(internalElements<0>().size());
    internalElements<0>().emplace_back(
        typename LayerContentTypeStruct<0>::type::TopoType(id), d, true);
    return HandleOfTypeAtLevel<Tag, 0>(id);
  }

  HandleOfTypeAtLevel<Tag, 0>
      add(typename LayerContentTypeStruct<0>::type::DataType &&d) {
    int id = static_cast<int>(internalElements<0>().size());
    internalElements<0>().emplace_back(
        typename LayerContentTypeStruct<0>::type::TopoType(id),
        std::forward<typename LayerContentTypeStruct<0>::type::DataType>(d),
        true);
    return HandleOfTypeAtLevel<Tag, 0>(id);
  }

  // removed
  template <int Level>
  inline bool removed(HandleOfTypeAtLevel<Tag, Level> h) const {
    return !internalElements<Level>()[h.id].exists;
  }

  // remove
  template <int Level> void remove(HandleOfTypeAtLevel<Tag, Level> h) {
    if (h.invalid() || removed(h))
      return;
    cleanLowers<Level>(h, std::integral_constant<bool, (Level > 0)>());
    internalElements<Level>()[h.id].exists = false; // mark as deleted
    cleanUppers<Level>(h,
                       std::integral_constant<bool, (Level < LayerNum - 1)>());
  }

private:
  template <int Level>
  inline void cleanLowers(HandleOfTypeAtLevel<Tag, Level> h, std::false_type) {}
  template <int Level>
  void cleanLowers(HandleOfTypeAtLevel<Tag, Level> h, std::true_type) {
    auto &c = internalElements<Level>()[h.id];
    auto &clowerTable = internalElements<Level - 1>();

    for (auto &lowh :
         c.topo.lowers) { // remove this h from all lowers' uppers set
      if (lowh.invalid() || removed(lowh))
        continue;
      auto &low = clowerTable[lowh.id];
      low.topo.uppers.erase(h);
    }
  }

  template <int Level>
  inline void cleanUppers(HandleOfTypeAtLevel<Tag, Level> h, std::false_type) {}
  template <int Level>
  void cleanUppers(HandleOfTypeAtLevel<Tag, Level> h, std::true_type) {
    auto &c = internalElements<Level>()[h.id];
    for (auto &uph : c.topo.uppers) {
      remove(uph);
    }
  }

public:
  // garbage collection
  inline void gc() {
    if (hasGarbage())
      gcUsingSequence(typename SequenceGenerator<LayerNum>::type());
  }

  inline void clear() {
    clearUsingSequence(typename SequenceGenerator<LayerNum>::type());
  }

  inline bool hasGarbage() const {
    return hasGarbageUsingSequence(
        typename SequenceGenerator<LayerNum>::type());
  }

  inline bool isDense() const { return !hasGarbage(); }

private:
  template <int... S> void gcUsingSequence(Sequence<S...>) {
    std::tuple<std::vector<HandleOfTypeAtLevel<Tag, S>>...> nlocs;
    int dummy[] = {RemoveAndMap(internalElements<S>(), std::get<S>(nlocs))...};
    int dummy2[] = {updateEachLayerHandles<S>(nlocs)...};
  }

  template <int Level, class NLocTupleT>
  int updateEachLayerHandles(const NLocTupleT &nlocs) {
    auto &eles = internalElements<Level>();
    for (size_t i = 0; i < eles.size(); i++) {
      UpdateOldHandle(std::get<Level>(nlocs), eles[i].topo.hd);
    }
    updateLowers<Level>(nlocs, std::integral_constant<bool, (Level > 0)>());
    updateUppers<Level>(nlocs,
                        std::integral_constant<bool, (Level < LayerNum - 1)>());
    return 0;
  }

  template <int Level, class NLocTupleT>
  inline void updateLowers(const NLocTupleT &nlocs, std::false_type) {}
  template <int Level, class NLocTupleT>
  void updateLowers(const NLocTupleT &nlocs, std::true_type) {
    auto &eles = internalElements<Level>();
    for (size_t i = 0; i < eles.size(); i++) {
      UpdateOldHandleContainer(std::get<Level - 1>(nlocs), eles[i].topo.lowers);
      RemoveInValidHandleFromContainer(eles[i].topo.lowers);
    }
  }

  template <int Level, class NLocTupleT>
  inline void updateUppers(const NLocTupleT &nlocs, std::false_type) {}
  template <int Level, class NLocTupleT>
  void updateUppers(const NLocTupleT &nlocs, std::true_type) {
    auto &eles = internalElements<Level>();
    for (size_t i = 0; i < eles.size(); i++) {
      UpdateOldHandleContainer(std::get<Level + 1>(nlocs), eles[i].topo.uppers);
      RemoveInValidHandleFromContainer(eles[i].topo.uppers);
    }
  }

  template <int... S> inline void clearUsingSequence(Sequence<S...>) {
    int dummy[] = {clearAtLevel<S>()...};
  }

  template <int Level> inline int clearAtLevel() {
    internalElements<Level>().clear();
    return 0;
  }

  template <int... S>
  inline bool hasGarbageUsingSequence(Sequence<S...>) const {
    for (bool r : {hasGarbageAtLevel<S>()...}) {
      if (r)
        return true;
    }
    return false;
  }

  template <int Level> inline bool hasGarbageAtLevel() const {
    for (const auto &t : internalElements<Level>()) {
      if (!t.exists)
        return true;
    }
    return false;
  }

public:
  template <class Archive> inline void serialize(Archive &ar) { ar(_contents); }

private:
  ContentsType _contents;
};

template <class VertDataT, class EdgeDataT>
using HomogeneousGraph02 =
    HomogeneousGraph<VertDataT, LayerConfig<EdgeDataT, 2>>;

template <class VertDataT, class EdgeDataT>
using HomogeneousGraph0x =
    HomogeneousGraph<VertDataT, LayerConfig<EdgeDataT, Dynamic>>;
}
}
