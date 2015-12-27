#pragma once

#include "handle.hpp"

namespace pano {
namespace core {

// ForestTopo
struct ForestTopo {
  Handle<ForestTopo> hd;
  std::set<Handle<ForestTopo>> children;
  Handle<ForestTopo> parent;

  template <class Archive> inline void serialize(Archive &ar) {
    ar(hd, children, parent);
  }
};

// Forest
template <class T> class Forest {
public:
  Forest() {}
  Forest(const Forest &) = default;
  Forest &operator=(const Forest &) = default;
  Forest(Forest &&f) { _nodes = std::move(f._nodes); }
  Forest &operator=(Forest &&f) {
    _nodes = std::move(f._nodes);
    return *this;
  }

  using NodeHandle = Handle<ForestTopo>;
  using NodeExistsPred = TripletExistsPred<ForestTopo, T>;

  inline const T &data(NodeHandle h) const { return _nodes[h.id].data; }
  inline T &data(NodeHandle h) { return _nodes[h.id].data; }
  inline const ForestTopo &topo(NodeHandle h) const {
    return _nodes[h.id].topo;
  }
  inline NodeHandle parent(NodeHandle h) const {
    return _nodes[h.id].topo.parent;
  }

  inline ConstConditionalContainerWrapper<TripletArray<ForestTopo, T>,
                                          NodeExistsPred>
  nodes() const {
    return ConstConditionalContainerWrapper<TripletArray<ForestTopo, T>,
                                            NodeExistsPred>(&_nodes);
  }
  inline ConditionalContainerWrapper<TripletArray<ForestTopo, T>,
                                     NodeExistsPred>
  nodes() {
    return ConditionalContainerWrapper<TripletArray<ForestTopo, T>,
                                       NodeExistsPred>(&_nodes);
  }
  inline const TripletArray<ForestTopo, T> &internalNodes() const {
    return _nodes;
  }
  inline NodeHandle firstRoot() const {
    for (auto &n : _nodes) {
      if (n.topo.parent.invalid())
        return n.topo.hd;
    }
    return NodeHandle();
  }

  inline NodeHandle add(NodeHandle parent, const T &data) {
    ForestTopo topo;
    topo.hd = NodeHandle(_nodes.size());
    topo.parent = parent;
    _nodes.emplace_back(std::move(topo), data);
    if (parent.valid()) {
      _nodes[parent.id].topo.children.insert(topo.hd);
    }
    return topo.hd;
  }

  inline NodeHandle add(NodeHandle parent, T &&data) {
    ForestTopo topo;
    topo.hd = NodeHandle(_nodes.size());
    topo.parent = parent;
    _nodes.emplace_back(std::move(topo), std::move(data));
    if (parent.valid()) {
      _nodes[parent.id].topo.children.insert(topo.hd);
    }
    return topo.hd;
  }

  inline NodeHandle addRoot(const T &data) { return add(NodeHandle(), data); }
  inline NodeHandle addRoot(T &&data) {
    return add(NodeHandle(), std::move(data));
  }
  inline bool isRoot(NodeHandle nh) const {
    return _nodes[nh.id].topo.parent.invalid();
  }
  inline bool isLeaf(NodeHandle nh) const {
    auto &children = _nodes[nh.id].topo.children;
    for (auto &ch : children) {
      if (ch.valid())
        return false;
    }
    return true;
  }

  inline void remove(NodeHandle h) {
    _nodes[h.id].exists = false;
    for (auto &ch : _nodes[h.id].topo.children) {
      remove(ch);
    }
  }

  void clear() { _nodes.clear(); }

  template <class NodeHandlePtrContainerT = HandlePtrArray<ForestTopo>>
  void gc(const NodeHandlePtrContainerT &nhPtrs = NodeHandlePtrContainerT()) {
    std::vector<NodeHandle> nnlocs;
    RemoveAndMap(_nodes, nnlocs);
    for (auto &node : _nodes) {
      UpdateOldHandle(nnlocs, node.topo.hd);
      UpdateOldHandle(nnlocs, node.topo.parent);
      UpdateOldHandleContainer(nnlocs, node.topo.children);
      RemoveInValidHandleFromContainer(node.topo.children);
    }
    for (auto &nhPtr : nhPtrs) {
      UpdateOldHandle(nnlocs, *nhPtr);
    }
  }

  template <class NodeHandleCallbackFunT>
  bool depthFirstSearch(NodeHandle asRoot,
                        const NodeHandleCallbackFunT &callback) const {
    assert(_nodes[asRoot.id].exists);
    if (!callback(asRoot))
      return false;
    for (auto &ch : _nodes[asRoot.id].topo.children) {
      if (_nodes[ch.id].exists) {
        if (!depthFirstSearch(ch, callback))
          return false;
      }
    }
    return true;
  }

  template <class NodeHandleCallbackFunT>
  bool breadthFirstSearch(NodeHandle asRoot,
                          const NodeHandleCallbackFunT &callback) const {
    assert(_nodes[asRoot.id].exists);
    std::queue<NodeHandle> nhs;
    nhs.push(asRoot);
    while (!nhs.empty()) {
      NodeHandle nh = nhs.front();
      nhs.pop();
      if (!callback(nh))
        return false;
      for (auto &ch : _nodes[nh.id].topo.children) {
        if (_nodes[ch.id].exists)
          nhs.push(ch);
      }
    }
    return true;
  }

  template <class Archive> inline void serialize(Archive &ar) { ar(_nodes); }

private:
  TripletArray<ForestTopo, T> _nodes;
};
}
}
