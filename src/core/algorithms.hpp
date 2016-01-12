#pragma once

#include "basic_types.hpp"
#include "utility.hpp"

namespace pano {
namespace core {

// generic algorithms
template <class IteratorT, class OutIteratorT, class IsCompatibleFunT>
void ForeachCompatibleWithLastElement(IteratorT begin, IteratorT end,
                                      OutIteratorT out,
                                      IsCompatibleFunT &&isCompWithLast) {
  if (begin == end)
    return;
  *out = *begin;
  ++out;
  IteratorT lastIter = begin;
  ++begin;
  while (begin != end) {
    if (isCompWithLast(*lastIter, *begin)) {
      *out = *begin;
      ++out;
      lastIter = begin;
    }
    ++begin;
  }
}

// merge, rearrange the input array
// DistanceFunctorT(a, b) -> DistanceT : compute the distance from a to b
// returns the begin iterators of merged groups
template <class IteratorT, class IterOutIteratorT, class DistanceT,
          class DistanceFunctorT = DefaultDistanceFunctor>
IterOutIteratorT
MergeNearNaive(IteratorT begin, IteratorT end, IterOutIteratorT itersOut,
               std::true_type, DistanceT thres,
               DistanceFunctorT &&distFun = DistanceFunctorT()) {
  if (begin == end)
    return itersOut;

  std::vector<IteratorT> gBegins(1, begin);
  for (auto i = std::next(begin); i != end; ++i) {
    DistanceT minDist = std::numeric_limits<DistanceT>::max();
    auto nearestGBeginIter = gBegins.end();
    for (auto giter = gBegins.begin(); giter != gBegins.end(); ++giter) {
      auto gBegin = *giter;
      DistanceT dist = distFun(*gBegin, *i);
      if (dist <= thres && dist < minDist) {
        minDist = dist;
        nearestGBeginIter = giter;
      }
    }
    if (nearestGBeginIter != gBegins.end()) { // found group
      if (std::next(nearestGBeginIter) != gBegins.end()) {
        auto nextGBegin = *std::next(nearestGBeginIter);
        std::rotate(nextGBegin, i, std::next(i));
        for (auto j = std::next(nearestGBeginIter); j != gBegins.end(); ++j)
          ++(*j);
      }
    } else { // add new group
      gBegins.push_back(i);
    }
  }
  return std::copy(gBegins.begin(), gBegins.end(), itersOut);
}

// merge, without rearrangement
// DistanceFunctorT(a, b) -> DistanceT : compute the distance from a to b
// returns the iterators pointing to group leaders
template <class IteratorT, class IterOutIteratorT, class DistanceT,
          class DistanceFunctorT = DefaultDistanceFunctor>
IterOutIteratorT
MergeNearNaive(IteratorT begin, IteratorT end, IterOutIteratorT itersOut,
               std::false_type, DistanceT thres,
               DistanceFunctorT &&distFun = DistanceFunctorT()) {
  if (begin == end)
    return itersOut;

  *(itersOut++) = begin;
  std::vector<IteratorT> gBegins(1, begin);
  for (auto i = std::next(begin); i != end; ++i) {
    auto giter = gBegins.begin();
    for (; giter != gBegins.end(); ++giter) {
      auto gBegin = *giter;
      auto dist = distFun(*gBegin, *i);
      if (dist <= thres) {
        break;
      }
    }
    if (giter == gBegins.end()) { // add new group
      gBegins.push_back(i);
      *(itersOut++) = i;
    }
  }

  return itersOut;
}

// merge using RTree, without rearrangement
// DistanceFunctorT(a, b) -> ? : compute the distance from a to b
// BoundingBoxFunctorT(a) -> Box<?,?> : compute the bounding box of a
// returns the iterators pointing to group leaders
template <class IteratorT, class IterOutIteratorT, class DistanceT,
          class DistanceFunctorT = DefaultDistanceFunctor,
          class BoundingBoxFunctorT = DefaultBoundingBoxFunctor>
IterOutIteratorT
MergeNearRTree(IteratorT begin, IteratorT end, IterOutIteratorT itersOut,
               std::false_type, DistanceT thres,
               DistanceFunctorT &&distFun = DistanceFunctorT(),
               BoundingBoxFunctorT &&getBoundingBox = BoundingBoxFunctorT()) {

  if (begin == end)
    return itersOut;

  using BoxType = decltype(getBoundingBox(*begin));
  using T = typename BoxType::Type;
  static const int N = BoxType::Dimension;

  third_party::RTree<IteratorT, T, N> rtree;
  for (auto i = begin; i != end; ++i) {
    Box<T, N> box = getBoundingBox(*i);
    for (int k = 0; k < N; k++) { // extend the box
      box.minCorner[k] -= thres * 2;
      box.maxCorner[k] += thres * 2;
    }
    // search in RTree
    int foundCount = 0;
    rtree.Search(box.minCorner.val, box.maxCorner.val,
                 [distFun, i, thres, &foundCount](IteratorT it) {
                   if (distFun(*i, *it) <= thres) {
                     foundCount++;
                     return false;
                   }
                   return true;
                 });
    if (foundCount == 0) {
      rtree.Insert(box.minCorner.val, box.maxCorner.val, i);
      *(itersOut)++ = i;
    }
  }

  return itersOut;
}

// Minimum Spanning Tree
// EdgeVertsGetterT(Edge e)->std::pair<Vert,Vert>
// EdgeCompareOnWeightT(Edge e1, Edge e2)->bool
//     determins whether weight of e1 is lower than weight of e2
// VertCompareT(Vert v1, Vert v2)->bool
//     used in std::map to register set id of vertices
template <class VertIteratorT, class EdgeIteratorT, class EdgeVertsGetterT,
          class EdgeOutputIteratorT, class EdgeCompareOnWeightT,
          class VertCompareT = std::less<
              typename std::iterator_traits<VertIteratorT>::value_type>>
void MinimumSpanningTree(VertIteratorT vertsBegin, VertIteratorT vertsEnd,
                         EdgeIteratorT edgesBegin, EdgeIteratorT edgesEnd,
                         EdgeOutputIteratorT MSTedges,
                         EdgeVertsGetterT &&vertsGetter,
                         EdgeCompareOnWeightT &&edgeCompareOnWeight,
                         VertCompareT &&vertCompare = VertCompareT()) {

  using Edge =
      typename std::iterator_traits<typename EdgeIteratorT>::value_type;
  using Vert =
      typename std::iterator_traits<typename VertIteratorT>::value_type;
  static_assert(
      std::is_same<
          std::decay_t<decltype(std::get<0>(vertsGetter(*edgesBegin)))>,
          Vert>::value &&
          std::is_same<
              std::decay_t<decltype(std::get<1>(vertsGetter(*edgesBegin)))>,
              Vert>::value,
      "result of EdgeVertsGetterT must be convertiable to std::tuple<Vert, "
      "Vert>!");

  std::vector<Edge> edges(edgesBegin, edgesEnd);
  std::sort(edges.begin(), edges.end(), edgeCompareOnWeight);

  std::map<Vert, int, VertCompareT> vertSetIds(vertCompare);
  int idx = 0;
  for (auto i = vertsBegin; i != vertsEnd; ++i)
    vertSetIds.insert(std::make_pair((*i), idx++));

  auto remainedEdgesBegin = edges.begin();
  while (remainedEdgesBegin != edges.end()) {
    Edge e = *remainedEdgesBegin;
    auto verts = vertsGetter(e);
    int fromid = vertSetIds[std::get<0>(verts)];
    int toid = vertSetIds[std::get<1>(verts)];
    if (fromid != toid) {
      *MSTedges++ = e;
      for (auto &vtoid : vertSetIds) {
        if (vtoid.second == toid) {
          vtoid.second = fromid;
        }
      }
    }
    ++remainedEdgesBegin;
  }
}

// DepthFirstSearch
template <class VertIteratorT, class NeighborVertsContainerGetterT,
          class VertCallbackT,
          class VertCompareT = std::less<
              typename std::iterator_traits<VertIteratorT>::value_type>>
void DepthFirstSearch(
    VertIteratorT vertsBegin, VertIteratorT vertsEnd,
    NeighborVertsContainerGetterT &&neighborVertsContainerGetter,
    VertCallbackT &&vertCallback, VertCompareT &&vertCompare = VertCompareT()) {

  using Vert =
      typename std::iterator_traits<typename VertIteratorT>::value_type;
  static_assert(
      std::is_same<
          Vert, std::decay_t<decltype(*std::begin(neighborVertsContainerGetter(
                    std::declval<Vert>())))>>::value,
      "NeighborVertsContainerGetterT should returns a container of Vert");
  static_assert(
      std::is_same<Vert,
                   std::decay_t<decltype(*std::end(neighborVertsContainerGetter(
                       std::declval<Vert>())))>>::value,
      "NeighborVertsContainerGetterT should returns a container of Vert");

  struct {
    bool operator()(Vert root, std::map<Vert, bool, VertCompareT> &vVisited,
                    NeighborVertsContainerGetterT vNeighborsGetter,
                    VertCallbackT vCallback) const {
      if (vVisited[root])
        return true;
      if (!vCallback(root))
        return false;

      vVisited[root] = true;
      auto vNeighborsContainer = vNeighborsGetter(root);
      for (const auto &v : vNeighborsContainer) {
        if (!(*this)(v, vVisited, vNeighborsGetter, vCallback))
          return false;
      }
      return true;
    }
  } depthFirstSearchOneTree;

  std::map<Vert, bool, VertCompareT> visited(vertCompare);
  for (auto i = vertsBegin; i != vertsEnd; ++i)
    visited[*i] = false;
  while (true) {
    auto rootIter = vertsBegin;
    while (rootIter != vertsEnd && visited[*rootIter]) {
      ++rootIter;
    }
    if (rootIter == vertsEnd)
      break;
    if (!depthFirstSearchOneTree(*rootIter, visited,
                                 neighborVertsContainerGetter, vertCallback))
      break;
  }
}

// BreadthFirstSearch
template <class VertIteratorT, class NeighborVertsContainerGetterT,
          class VertCallbackT,
          class VertCompareT = std::less<
              typename std::iterator_traits<VertIteratorT>::value_type>>
void BreadthFirstSearch(
    VertIteratorT vertsBegin, VertIteratorT vertsEnd,
    NeighborVertsContainerGetterT neighborVertsContainerGetter,
    VertCallbackT vertCallback, VertCompareT vertCompare = VertCompareT()) {

  using Vert =
      typename std::iterator_traits<typename VertIteratorT>::value_type;
  static_assert(
      std::is_same<
          Vert, std::decay_t<decltype(*std::begin(neighborVertsContainerGetter(
                    std::declval<Vert>())))>>::value,
      "NeighborVertsContainerGetterT should returns a container of Vert");
  static_assert(
      std::is_same<Vert,
                   std::decay_t<decltype(*std::end(neighborVertsContainerGetter(
                       std::declval<Vert>())))>>::value,
      "NeighborVertsContainerGetterT should returns a container of Vert");

  struct {
    bool operator()(const Vert &root,
                    std::map<Vert, bool, VertCompareT> &vVisited,
                    const NeighborVertsContainerGetterT &vNeighborsGetter,
                    VertCallbackT vCallback) {
      std::queue<Vert> Q;
      Q.push(root);
      vVisited[root] = true;
      while (!Q.empty()) {
        Vert v = Q.front();
        Q.pop();
        if (!vCallback(v))
          return false;
        auto vNeighborsContainer = vNeighborsGetter(v);
        for (const auto &vv : vNeighborsContainer) {
          if (vVisited.at(vv))
            continue;
          Q.push(vv);
          vVisited[vv] = true;
        }
      }
      return true;
    }
  } breadthFirstSearchOneTree;

  std::map<Vert, bool, VertCompareT> visited(vertCompare);
  for (auto i = vertsBegin; i != vertsEnd; ++i)
    visited[*i] = false;
  while (true) {
    auto rootIter = vertsBegin;
    while (rootIter != vertsEnd && visited[*rootIter]) {
      ++rootIter;
    }
    if (rootIter == vertsEnd)
      break;
    if (!breadthFirstSearchOneTree(*rootIter, visited,
                                   neighborVertsContainerGetter, vertCallback))
      break;
  }
}

// Topological Sort (using Depth First Search)
template <class VertIteratorT, class VertOutIteratorT,
          class PredecessorVertsContainerGetterT,
          class VertCompareT = std::less<
              typename std::iterator_traits<VertIteratorT>::value_type>>
void TopologicalSort(
    VertIteratorT vertsBegin, VertIteratorT vertsEnd,
    VertOutIteratorT sortedVertsBegin,
    PredecessorVertsContainerGetterT predecessorVertsContainerGetter,
    VertCompareT vertCompare = VertCompareT()) {

  using Vert =
      typename std::iterator_traits<typename VertIteratorT>::value_type;
  static_assert(
      std::is_same<Vert, std::decay_t<decltype(
                             *std::begin(predecessorVertsContainerGetter(
                                 std::declval<Vert>())))>>::value,
      "PredecessorVertsContainerGetterT should returns a container of Vert");
  static_assert(
      std::is_same<
          Vert, std::decay_t<decltype(*std::end(predecessorVertsContainerGetter(
                    std::declval<Vert>())))>>::value,
      "PredecessorVertsContainerGetterT should returns a container of Vert");

  struct {
    void operator()(
        Vert root, std::map<Vert, bool, VertCompareT> &vVisited,
        VertOutIteratorT sortedVertsOut,
        PredecessorVertsContainerGetterT predecessorVertsContainerGetter) {
      if (vVisited[root])
        return;

      vVisited[root] = true;
      for (const auto &v : predecessorVertsContainerGetter(root)) {
        (*this)(v, vVisited, sortedVertsOut, predecessorVertsContainerGetter);
      }

      *sortedVertsOut = root;
      ++sortedVertsOut;
    }
  } depthFirstSearchOneTree;

  std::map<Vert, bool, VertCompareT> visited(vertCompare);
  for (auto i = vertsBegin; i != vertsEnd; ++i)
    visited[*i] = false;
  while (true) {
    auto rootIter = vertsBegin;
    while (rootIter != vertsEnd && visited[*rootIter]) {
      ++rootIter;
    }
    if (rootIter == vertsEnd)
      break;
    depthFirstSearchOneTree(*rootIter, visited, sortedVertsBegin,
                            predecessorVertsContainerGetter);
  }
}

// Connected Components
template <class VertIteratorT, class NeighborVertsContainerGetterT,
          class VertexTypeRecorderT,
          class VertCompareT = std::less<
              typename std::iterator_traits<VertIteratorT>::value_type>>
int ConnectedComponents(
    VertIteratorT vertsBegin, VertIteratorT vertsEnd,
    NeighborVertsContainerGetterT neighborVertsContainerGetter,
    VertexTypeRecorderT vertTypeRecorder,
    VertCompareT vertCompare = VertCompareT()) {

  using Vert =
      typename std::iterator_traits<typename VertIteratorT>::value_type;
  static_assert(
      std::is_same<
          Vert, std::decay_t<decltype(*std::begin(neighborVertsContainerGetter(
                    std::declval<Vert>())))>>::value,
      "NeighborVertsContainerGetterT should returns a container of Vert");
  static_assert(
      std::is_same<Vert,
                   std::decay_t<decltype(*std::end(neighborVertsContainerGetter(
                       std::declval<Vert>())))>>::value,
      "NeighborVertsContainerGetterT should returns a container of Vert");

  struct {
    void operator()(const Vert &root,
                    std::map<Vert, bool, VertCompareT> &vVisited,
                    const NeighborVertsContainerGetterT &vNeighborsGetter,
                    const VertexTypeRecorderT &vTypeRecorder, int cid) {
      std::queue<Vert> Q;
      Q.push(root);
      vVisited[root] = true;
      while (!Q.empty()) {
        Vert v = Q.front();
        Q.pop();
        vTypeRecorder(v, cid);
        auto vNeighborsContainer = vNeighborsGetter(v);
        for (const auto &vv : vNeighborsContainer) {
          if (vVisited.at(vv))
            continue;
          Q.push(vv);
          vVisited[vv] = true;
        }
      }
    }
  } breadthFirstSearchOneTree;

  std::map<Vert, bool, VertCompareT> visited(vertCompare);
  for (auto i = vertsBegin; i != vertsEnd; ++i)
    visited[*i] = false;

  int cid = 0;
  while (true) {
    auto rootIter = vertsBegin;
    while (rootIter != vertsEnd && visited[*rootIter]) {
      ++rootIter;
    }
    if (rootIter == vertsEnd)
      break;
    breadthFirstSearchOneTree(*rootIter, visited, neighborVertsContainerGetter,
                              vertTypeRecorder, cid);
    cid++;
  }

  return cid;
}

// SimulatedAnnealing
// - EnergyFunT: (StateT)->Scalar
// - TemperatureFunT: (int iter)->Scalar
// - NeighborsFunT: (StateT, int iter, (StateT)->void )
template <class T> struct DefaultTemperatureFunctor {
  DefaultTemperatureFunctor(T t0) : t0(t0) {}
  inline T operator()(int iter) const { return t0 / iter; }
  T t0;
};
template <class StateT, class EnergyFunT, class TemperatureFunT,
          class NeighborsFunT, class RNG>
int SimulatedAnnealing(StateT &state, EnergyFunT energyFun,
                        TemperatureFunT temperatureFun,
                        NeighborsFunT neighborsFun,
                        RNG &&rng) {

  std::uniform_real_distribution<double> dist(0.0, 1.0);
  double energy = energyFun(state);
  int i = 0;

  while (true) {
    double temperature = temperatureFun(i);

    StateT newStateWithLowestEnergy;
    double lowestEnergy = std::numeric_limits<double>::infinity();
    bool hasNewState = false;
    neighborsFun(state, i, [&newStateWithLowestEnergy, &lowestEnergy,
                            &hasNewState, &energyFun](auto &&newState) {
      double newEnergy = energyFun(newState);
      if (newEnergy < lowestEnergy) {
        newStateWithLowestEnergy = newState;
        lowestEnergy = newEnergy;
        hasNewState = true;
      }
    });

    if (!hasNewState) {
      break;
    }

    double prob = 1.0;
    if (energy <= lowestEnergy) {
      prob = exp(-(lowestEnergy - energy) / temperature);
    }
    if (prob >= dist(rng)) {
      state = newStateWithLowestEnergy;
      energy = lowestEnergy;
    }

    ++i;
  }

  return i;
}
}
}
