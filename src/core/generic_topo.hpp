#ifndef PANORAMIX_CORE_GENERIC_TOPO_HPP
#define PANORAMIX_CORE_GENERIC_TOPO_HPP


#include "handle.hpp"

namespace panoramix {
    namespace core { 

        // the mesh class
        struct VertTopo;
        struct HalfTopo;
        struct FaceTopo;
        struct VertTopo {
            Handle<VertTopo> hd;
            HandleArray<HalfTopo> halfedges;
            template <class Archive> inline void serialize(Archive & ar) { ar(hd, halfedges); }
        };
        struct HalfTopo {
            Handle<HalfTopo> hd;
            Handle<VertTopo> endVertices[2];
            inline Handle<VertTopo> & from() { return endVertices[0]; }
            inline Handle<VertTopo> & to() { return endVertices[1]; }
            inline const Handle<VertTopo> & from() const { return endVertices[0]; }
            inline const Handle<VertTopo> & to() const { return endVertices[1]; }
            Handle<HalfTopo> opposite;
            Handle<FaceTopo> face;
            template <class Archive> inline void serialize(Archive & ar) { ar(hd, endVertices, opposite, face); }
        };
        struct FaceTopo {
            Handle<FaceTopo> hd;
            HandleArray<HalfTopo> halfedges;
            template <class Archive> inline void serialize(Archive & ar) { ar(hd, halfedges); }
        };

        template <class VertDataT, class HalfDataT = Dummy, class FaceDataT = Dummy>
        class Mesh {
        public:
            static const int LayerNum = 2;

            using VertData = VertDataT;
            using HalfData = HalfDataT;
            using FaceData = FaceDataT;

            using VertHandle = Handle<VertTopo>;
            using HalfHandle = Handle<HalfTopo>;
            using FaceHandle = Handle<FaceTopo>;

            using VertsTable = TripletArray<VertTopo, VertDataT>;
            using HalfsTable = TripletArray<HalfTopo, HalfDataT>;
            using FacesTable = TripletArray<FaceTopo, FaceDataT>;

            using VertExistsPred = TripletExistsPred<VertTopo, VertDataT>;
            using HalfExistsPred = TripletExistsPred<HalfTopo, HalfDataT>;
            using FaceExistsPred = TripletExistsPred<FaceTopo, FaceDataT>;

            using Vertex = typename VertsTable::value_type;
            using HalfEdge = typename HalfsTable::value_type;
            using Face = typename FacesTable::value_type;

            inline VertsTable & internalVertices() { return _verts; }
            inline HalfsTable & internalHalfEdges() { return _halfs; }
            inline FacesTable & internalFaces() { return _faces; }
            inline const VertsTable & internalVertices() const { return _verts; }
            inline const HalfsTable & internalHalfEdges() const { return _halfs; }
            inline const FacesTable & internalFaces() const { return _faces; }

            inline ConditionalContainerWrapper<VertsTable, VertExistsPred> vertices() {
                return ConditionalContainerWrapper<VertsTable, VertExistsPred>(&_verts);
            }
            inline ConditionalContainerWrapper<HalfsTable, HalfExistsPred> halfedges() {
                return ConditionalContainerWrapper<HalfsTable, HalfExistsPred>(&_halfs);
            }
            inline ConditionalContainerWrapper<FacesTable, FaceExistsPred> faces() {
                return ConditionalContainerWrapper<FacesTable, FaceExistsPred>(&_faces);
            }
            inline ConstConditionalContainerWrapper<VertsTable, VertExistsPred> vertices() const {
                return ConstConditionalContainerWrapper<VertsTable, VertExistsPred>(&_verts);
            }
            inline ConstConditionalContainerWrapper<HalfsTable, HalfExistsPred> halfedges() const {
                return ConstConditionalContainerWrapper<HalfsTable, HalfExistsPred>(&_halfs);
            }
            inline ConstConditionalContainerWrapper<FacesTable, FaceExistsPred> faces() const {
                return ConstConditionalContainerWrapper<FacesTable, FaceExistsPred>(&_faces);
            }

            inline VertTopo & topo(VertHandle v) { return _verts[v.id].topo; }
            inline HalfTopo & topo(HalfHandle h) { return _halfs[h.id].topo; }
            inline FaceTopo & topo(FaceHandle f) { return _faces[f.id].topo; }
            inline const VertTopo & topo(VertHandle v) const { return _verts[v.id].topo; }
            inline const HalfTopo & topo(HalfHandle h) const { return _halfs[h.id].topo; }
            inline const FaceTopo & topo(FaceHandle f) const { return _faces[f.id].topo; }

            inline VertDataT & data(VertHandle v) { return _verts[v.id].data; }
            inline HalfDataT & data(HalfHandle h) { return _halfs[h.id].data; }
            inline FaceDataT & data(FaceHandle f) { return _faces[f.id].data; }
            inline const VertDataT & data(VertHandle v) const { return _verts[v.id].data; }
            inline const HalfDataT & data(HalfHandle h) const { return _halfs[h.id].data; }
            inline const FaceDataT & data(FaceHandle f) const { return _faces[f.id].data; }

            VertHandle addVertex(const VertDataT & vd = VertDataT());
            HalfHandle addEdge(VertHandle from, VertHandle to,
                const HalfDataT & hd = HalfDataT(), const HalfDataT & hdrev = HalfDataT(), bool mergeDuplicateEdge = true);
            FaceHandle addFace(const HandleArray<HalfTopo> & halfedges,
                const FaceDataT & fd = FaceDataT());
            FaceHandle addFace(const HandleArray<VertTopo> & vertices, bool autoflip = true,
                const FaceDataT & fd = FaceDataT());
            template <class VertHandleIteratorT,
            class = std::enable_if_t<std::is_same<std::iterator_traits<VertHandleIteratorT>::value_type, VertHandle>::value >>
                FaceHandle addFace(VertHandleIteratorT vhBegin, VertHandleIteratorT vhEnd, bool autoflip = true, const FaceDataT & fd = FaceDataT());
            FaceHandle addFace(VertHandle v1, VertHandle v2, VertHandle v3, bool autoflip = true,
                const FaceDataT & fd = FaceDataT());
            FaceHandle addFace(VertHandle v1, VertHandle v2, VertHandle v3, VertHandle v4, bool autoflip = true,
                const FaceDataT & fd = FaceDataT());

            HalfHandle findEdge(VertHandle from, VertHandle to) const;

            inline bool removed(FaceHandle f) const { return !_faces[f.id].exists; }
            inline bool removed(HalfHandle e) const { return !_halfs[e.id].exists; }
            inline bool removed(VertHandle v) const { return !_verts[v.id].exists; }

            inline void remove(FaceHandle f);
            inline void remove(HalfHandle e);
            inline void remove(VertHandle v);

            Mesh & unite(const Mesh & m);

            // garbage collection
            template <class VertHandlePtrContainerT = HandlePtrArray<VertTopo>,
            class HalfHandlePtrContainerT = HandlePtrArray<HalfTopo>,
            class FaceHandlePtrContainerT = HandlePtrArray<FaceTopo>
            >
            void gc(const VertHandlePtrContainerT & vps = VertHandlePtrContainerT(),
            const HalfHandlePtrContainerT & hps = HalfHandlePtrContainerT(),
            const FaceHandlePtrContainerT & fps = FaceHandlePtrContainerT());

            void clear();

            template <class Archive> inline void serialize(Archive & ar) { ar(_verts, _halfs, _faces); }

        private:
            VertsTable _verts;
            HalfsTable _halfs;
            FacesTable _faces;
        };

        // implementation of class Mesh

        template <class VertDataT, class HalfDataT, class FaceDataT>
        typename Mesh<VertDataT, HalfDataT, FaceDataT>::VertHandle
            Mesh<VertDataT, HalfDataT, FaceDataT>::addVertex(const VertDataT & vd) {
                VertTopo topo;
                topo.hd.id = _verts.size();
                _verts.emplace_back(std::move(topo), vd, true);
                return _verts.back().topo.hd;
            }

        template <class VertDataT, class HalfDataT, class FaceDataT>
        typename Mesh<VertDataT, HalfDataT, FaceDataT>::HalfHandle
            Mesh<VertDataT, HalfDataT, FaceDataT>::addEdge(VertHandle from, VertHandle to, const HalfDataT & hd, const HalfDataT & hdrev, bool mergeDuplicateEdge) {
                if (from == to){
                    return HalfHandle();
                }
                // find existed halfedge
                if (mergeDuplicateEdge){
                    HalfHandle hh = findEdge(from, to);
                    if (hh.valid()){
                        _halfs[hh.id].data = hd;
                        _halfs[_halfs[hh.id].topo.opposite.id].data = hdrev;
                        return hh;
                    }
                }

                HalfHandle hh1(_halfs.size());
                Triplet<HalfTopo, HalfDataT> ht;
                ht.topo.hd.id = _halfs.size();
                ht.topo.from() = from;
                ht.topo.to() = to;
                ht.exists = true;
                ht.data = hd;
                //_halfs.push_back({ { { _halfs.size() }, { from, to }, { -1 }, { -1 } }, true, hd });
                _halfs.push_back(ht);
                HalfHandle hh2(_halfs.size());
                ht.topo.hd.id = _halfs.size();
                ht.topo.from() = to;
                ht.topo.to() = from;
                ht.exists = true;
                ht.data = hdrev;
                //_halfs.push_back({ { { _halfs.size() }, { to, from }, { -1 }, { -1 } }, true, hdrev });
                _halfs.push_back(ht);

                _halfs[hh1.id].topo.opposite = hh2;
                _halfs[hh2.id].topo.opposite = hh1;

                _verts[from.id].topo.halfedges.push_back(hh1);
                _verts[to.id].topo.halfedges.push_back(hh2);
                return hh1;
            }

        template <class VertDataT, class HalfDataT, class FaceDataT>
        typename Mesh<VertDataT, HalfDataT, FaceDataT>::FaceHandle
            Mesh<VertDataT, HalfDataT, FaceDataT>::addFace(const HandleArray<HalfTopo> & halfedges, const FaceDataT & fd) {
                Triplet<FaceTopo, FaceDataT> ft;
                ft.topo.hd.id = _faces.size();
                ft.topo.halfedges = halfedges;
                ft.exists = true;
                ft.data = fd;
                //_faces.push_back({ { { _faces.size() }, halfedges }, true, fd });
                _faces.push_back(ft);
                for (HalfHandle hh : halfedges){
                    _halfs[hh.id].topo.face = _faces.back().topo.hd;
                }
                return _faces.back().topo.hd;
            }

        template <class VertDataT, class HalfDataT, class FaceDataT>
        typename Mesh<VertDataT, HalfDataT, FaceDataT>::FaceHandle
            Mesh<VertDataT, HalfDataT, FaceDataT>::addFace(const HandleArray<VertTopo> & vertices, bool autoflip, const FaceDataT & fd) {
                HandleArray<HalfTopo> halfs;
                assert(vertices.size() >= 3);
                HalfHandle hh = findEdge(vertices.back(), vertices.front());
                auto verts = vertices;
                if (hh.valid() && _halfs[hh.id].topo.face.valid() && autoflip){
                    std::reverse(verts.begin(), verts.end());
                }

                for (size_t i = 0; i < verts.size(); i++){
                    size_t inext = (i + 1) % verts.size();
                    halfs.push_back(addEdge(verts[i], verts[inext]));
                }
                return addFace(halfs, fd);
            }

        template <class VertDataT, class HalfDataT, class FaceDataT>
        template <class VertHandleIteratorT, class>
        typename Mesh<VertDataT, HalfDataT, FaceDataT>::FaceHandle
            Mesh<VertDataT, HalfDataT, FaceDataT>::addFace(
            VertHandleIteratorT vhBegin, VertHandleIteratorT vhEnd, bool autoflip, const FaceDataT & fd) {
                HandleArray<HalfTopo> halfs;
                HalfHandle hh = findEdge(vertices.back(), vertices.front());
                HandleArray<VertTopo> verts(vhBegin, vhEnd);
                assert(verts.size() >= 3);
                if (hh.valid() && _halfs[hh.id].topo.face.valid() && autoflip){
                    std::reverse(verts.begin(), verts.end());
                }

                for (size_t i = 0; i < verts.size(); i++){
                    size_t inext = (i + 1) % verts.size();
                    halfs.push_back(addEdge(verts[i], verts[inext]));
                }
                return addFace(halfs, fd);
            }

        template <class VertDataT, class HalfDataT, class FaceDataT>
        typename Mesh<VertDataT, HalfDataT, FaceDataT>::FaceHandle
            Mesh<VertDataT, HalfDataT, FaceDataT>::addFace(VertHandle v1, VertHandle v2, VertHandle v3, bool autoflip, const FaceDataT & fd) {
                HalfHandle hh = findEdge(v3, v1);
                if (hh.valid() && _halfs[hh.id].topo.face.valid() && autoflip){
                    std::swap(v1, v3);
                }
                return addFace({
                    addEdge(v1, v2),
                    addEdge(v2, v3),
                    addEdge(v3, v1)
                }, fd);
            }

        template <class VertDataT, class HalfDataT, class FaceDataT>
        typename Mesh<VertDataT, HalfDataT, FaceDataT>::FaceHandle
            Mesh<VertDataT, HalfDataT, FaceDataT>::addFace(VertHandle v1, VertHandle v2, VertHandle v3, VertHandle v4, bool autoflip, const FaceDataT & fd) {
                HalfHandle hh = findEdge(v4, v1);
                if (hh.valid() && _halfs[hh.id].topo.face.valid() && autoflip){
                    std::swap(v1, v4);
                }
                return addFace({
                    addEdge(v1, v2),
                    addEdge(v2, v3),
                    addEdge(v3, v4),
                    addEdge(v4, v1)
                }, fd);
            }


        template <class VertDataT, class HalfDataT, class FaceDataT>
        typename Mesh<VertDataT, HalfDataT, FaceDataT>::HalfHandle
            Mesh<VertDataT, HalfDataT, FaceDataT>::findEdge(VertHandle from, VertHandle to) const {
                for (HalfHandle hh : _verts[from.id].topo.halfedges){
                    assert(_halfs[hh.id].topo.endVertices[0] == from);
                    if (_halfs[hh.id].topo.endVertices[1] == to){
                        return hh;
                    }
                }
                return HalfHandle();
            }

        template <class VertDataT, class HalfDataT, class FaceDataT>
        void Mesh<VertDataT, HalfDataT, FaceDataT>::remove(FaceHandle f) {
            if (f.invalid() || removed(f))
                return;
            _faces[f.id].exists = false;
            for (auto & hh : _faces[f.id].topo.halfedges){
                hh.reset();
            }
        }

        template <class VertDataT, class HalfDataT, class FaceDataT>
        void Mesh<VertDataT, HalfDataT, FaceDataT>::remove(HalfHandle h) {
            if (h.invalid() || removed(h))
                return;
            HalfHandle hop = _halfs[h.id].topo.opposite;
            _halfs[h.id].exists = false;
            _halfs[hop.id].exists = false;

            remove(_halfs[h.id].topo.face);
            remove(_halfs[hop.id].topo.face);

            _halfs[h.id].topo.from().reset();
            _halfs[hop.id].topo.to().reset();
            _halfs[h.id].topo.face.reset();
            _halfs[hop.id].topo.face.reset();
        }

        template <class VertDataT, class HalfDataT, class FaceDataT>
        void Mesh<VertDataT, HalfDataT, FaceDataT>::remove(VertHandle v) {
            if (v.invalid() || removed(v))
                return;
            _verts[v.id].exists = false;
            for (HalfHandle hh : _verts[v.id].topo.halfedges)
                remove(hh);
            _verts[v.id].topo.halfedges.clear();
        }

        template <class VertDataT, class HalfDataT, class FaceDataT>
        Mesh<VertDataT, HalfDataT, FaceDataT>& Mesh<VertDataT, HalfDataT, FaceDataT>::unite(const Mesh & m) {
            std::vector<VertHandle> vtable(m.Vertices().size());
            std::vector<HalfHandle> htable(m.HalfEdges().size());
            std::vector<FaceHandle> ftable(m.Faces().size());

            for (auto v : m.vertices()){
                vtable[v.topo.hd.id] = addVertex(v.data);
            }
            for (auto h : m.halfedges()){
                VertHandle oldfrom = h.topo.from();
                VertHandle oldto = h.topo.to();
                VertHandle newfrom = vtable[oldfrom.id];
                VertHandle newto = vtable[oldto.id];
                htable[h.topo.hd.id] = addEdge(newfrom, newto, h.data, m.data(h.topo.opposite));
            }
            for (auto f : m.faces()){
                HandleArray<HalfTopo> hs;
                hs.reserve(f.topo.halfedges.size());
                for (auto hh : f.topo.halfedges){
                    hs.push_back(htable[hh.id]);
                }
                ftable[f.topo.hd.id] = addFace(hs, f.data);
            }

            return *this;
        }



        template <class VertDataT, class HalfDataT, class FaceDataT>
        template <class VertHandlePtrContainerT, class HalfHandlePtrContainerT, class FaceHandlePtrContainerT>
        void Mesh<VertDataT, HalfDataT, FaceDataT>::gc(const VertHandlePtrContainerT & vps,
            const HalfHandlePtrContainerT & hps,
            const FaceHandlePtrContainerT & fps){
            std::vector<VertHandle> vnlocs;
            std::vector<HalfHandle> hnlocs;
            std::vector<FaceHandle> fnlocs;
            RemoveAndMap(_verts, vnlocs);
            RemoveAndMap(_halfs, hnlocs);
            RemoveAndMap(_faces, fnlocs);

            for (size_t i = 0; i < _verts.size(); i++){
                UpdateOldHandle(vnlocs, _verts[i].topo.hd);
                UpdateOldHandleContainer(hnlocs, _verts[i].topo.halfedges);
                RemoveInValidHandleFromContainer(_verts[i].topo.halfedges);
            }
            for (size_t i = 0; i < _halfs.size(); i++){
                UpdateOldHandle(hnlocs, _halfs[i].topo.hd);
                UpdateOldHandleContainer(vnlocs, _halfs[i].topo.endVertices);
                UpdateOldHandle(hnlocs, _halfs[i].topo.opposite);
                UpdateOldHandle(fnlocs, _halfs[i].topo.face);
            }
            for (size_t i = 0; i < _faces.size(); i++){
                UpdateOldHandle(fnlocs, _faces[i].topo.hd);
                UpdateOldHandleContainer(hnlocs, _faces[i].topo.halfedges);
                RemoveInValidHandleFromContainer(_faces[i].topo.halfedges);
            }
            for (auto vp : vps){
                UpdateOldHandle(vnlocs, *vp);
            }
            for (auto hp : hps){
                UpdateOldHandle(hnlocs, *hp);
            }
            for (auto fp : fps){
                UpdateOldHandle(fnlocs, *fp);
            }
        }

        template <class VertDataT, class HalfDataT, class FaceDataT>
        void Mesh<VertDataT, HalfDataT, FaceDataT>::clear() {
            _verts.clear();
            _halfs.clear();
            _faces.clear();
        }


        // Forest 
        struct ForestTopo {
            Handle<ForestTopo> hd;
            std::set<Handle<ForestTopo>> children;
            Handle<ForestTopo> parent;

            template <class Archive> inline void serialize(Archive & ar) { ar(hd, children, parent); }
        };

        template <class T>
        class Forest {
        public:
            using NodeHandle = Handle<ForestTopo>;
            using NodeExistsPred = TripletExistsPred<ForestTopo, T>;

            inline const T & data(NodeHandle h) const { return _nodes[h.id].data; }
            inline T & data(NodeHandle h) { return _nodes[h.id].data; }
            inline const ForestTopo & topo(NodeHandle h) const { return _nodes[h.id].topo; }
            inline NodeHandle parent(NodeHandle h) const { return _nodes[h.id].topo.parent; }

            inline ConstConditionalContainerWrapper<TripletArray<ForestTopo, T>, NodeExistsPred> nodes() const {
                return ConstConditionalContainerWrapper<TripletArray<ForestTopo, T>, NodeExistsPred>(&_nodes);
            }
            inline ConditionalContainerWrapper<TripletArray<ForestTopo, T>, NodeExistsPred> nodes() {
                return ConditionalContainerWrapper<TripletArray<ForestTopo, T>, NodeExistsPred>(&_nodes);
            }
            inline const TripletArray<ForestTopo, T> & internalNodes() const { return _nodes; }
            inline NodeHandle firstRoot() const {
                for (auto & n : _nodes){
                    if (n.topo.parent.invalid())
                        return n.topo.hd;
                }
                return NodeHandle();
            }

            inline NodeHandle add(NodeHandle parent, const T & data) {
                ForestTopo topo;
                topo.hd = NodeHandle(_nodes.size());
                topo.parent = parent;
                _nodes.emplace_back(std::move(topo), data);
                if (parent.valid()){
                    _nodes[parent.id].topo.children.insert(topo.hd);
                }
                return topo.hd;
            }

            inline NodeHandle add(NodeHandle parent, T && data) {
                ForestTopo topo;
                topo.hd = NodeHandle(_nodes.size());
                topo.parent = parent;
                _nodes.emplace_back(std::move(topo), std::move(data));
                if (parent.valid()){
                    _nodes[parent.id].topo.children.insert(topo.hd);
                }
                return topo.hd;
            }

            inline NodeHandle addRoot(const T & data){ return add(NodeHandle(), data); }
            inline NodeHandle addRoot(T && data) { return add(NodeHandle(), std::move(data)); }
            inline bool isRoot(NodeHandle nh) const { return _nodes[nh.id].topo.parent.invalid(); }
            inline bool isLeaf(NodeHandle nh) const {
                auto & children = _nodes[nh.id].topo.children;
                for (auto & ch : children){
                    if (ch.valid())
                        return false;
                }
                return true;
            }

            inline void remove(NodeHandle h) {
                _nodes[h.id].exists = false;
                for (auto & ch : _nodes[h.id].topo.children){
                    remove(ch);
                }
            }

            template <class NodeHandlePtrContainerT = HandlePtrArray<ForestTopo>>
            void gc(const NodeHandlePtrContainerT & nhPtrs = NodeHandlePtrContainerT()) {
                std::vector<NodeHandle> nnlocs;
                RemoveAndMap(_nodes, nnlocs);
                for (auto & node : _nodes){
                    UpdateOldHandle(nnlocs, node.topo.hd);
                    UpdateOldHandle(nnlocs, node.topo.parent);
                    UpdateOldHandleContainer(nnlocs, node.topo.children);
                    RemoveInValidHandleFromContainer(node.topo.children);
                }
                for (auto & nhPtr : nhPtrs){
                    UpdateOldHandle(nnlocs, *nhPtr);
                }
            }

            template <class NodeHandleCallbackFunT>
            bool depthFirstSearch(NodeHandle asRoot, const NodeHandleCallbackFunT & callback) const {
                assert(_nodes[asRoot.id].exists);
                if (!callback(asRoot))
                    return false;
                for (auto & ch : _nodes[asRoot.id].topo.children){
                    if (_nodes[ch.id].exists){
                        if (!depthFirstSearch(ch, callback))
                            return false;
                    }
                }
                return true;
            }


            template <class NodeHandleCallbackFunT>
            bool breadthFirstSearch(NodeHandle asRoot, const NodeHandleCallbackFunT & callback) const {
                assert(_nodes[asRoot.id].exists);
                std::queue<NodeHandle> nhs;
                nhs.push(asRoot);
                while (!nhs.empty()){
                    NodeHandle nh = nhs.front();
                    nhs.pop();
                    if (!callback(nh))
                        return false;
                    for (auto & ch : _nodes[nh.id].topo.children){
                        if (_nodes[ch.id].exists)
                            nhs.push(ch);
                    }
                }
                return true;
            }


            template <class Archive> inline void serialize(Archive & ar) { ar(_nodes); }

        private:
            TripletArray<ForestTopo, T> _nodes;
        };






        // helper classes for HomogeneousGraph
        template <class Tag, int L, int ChildN>
        struct Topo {
            static const int Level = L;
            std::array<HandleOfTypeAtLevel<Tag, Level - 1>, ChildN> lowers; // use std::array
            std::set<HandleOfTypeAtLevel<Tag, Level + 1>> uppers;
            HandleOfTypeAtLevel<Tag, Level> hd;
            explicit inline Topo(int id = -1) : hd(id){}
            explicit inline Topo(int id, std::initializer_list<HandleOfTypeAtLevel<Tag, Level - 1>> ls) : hd(id) {
                assert(ls.size() == lowers.size());
                std::copy(ls.begin(), ls.end(), lowers.begin());
            }
            template <class Archive> inline void serialize(Archive & ar) { ar(lowers, uppers, hd); }
        };

        // dynamic sized
        template <class Tag, int L>
        struct Topo<Tag, L, Dynamic> {
            static const int Level = L;
            std::vector<HandleOfTypeAtLevel<Tag, Level - 1>> lowers; // use std::array
            std::set<HandleOfTypeAtLevel<Tag, Level + 1>> uppers;
            HandleOfTypeAtLevel<Tag, Level> hd;
            explicit inline Topo(int id = -1) : hd(id){}
            explicit inline Topo(int id, std::initializer_list<HandleOfTypeAtLevel<Tag, Level - 1>> ls)
                : hd(id), lowers(ls) {}
            template <class Archive> inline void serialize(Archive & ar) { ar(lowers, uppers, hd); }
        };


        // zero sized
        template <class Tag, int L>
        struct Topo<Tag, L, 0> {
            static const int Level = L;
            std::set<HandleOfTypeAtLevel<Tag, Level + 1>> uppers;
            HandleOfTypeAtLevel<Tag, Level> hd;
            explicit inline Topo(int id = -1) : hd(id){}
            explicit inline Topo(int id, std::initializer_list<HandleOfTypeAtLevel<Tag, Level - 1>> ls)
                : hd(id) {
                assert(ls.size() == 0);
            }
            template <class Archive> inline void serialize(Archive & ar) { ar(uppers, hd); }
        };

        // configuration for each layer
        template <class DataT, int ChildN>
        struct LayerConfig {
            using DataType = DataT;
            enum { ChildNum = ChildN };
        };

        template <class T>
        struct IsLayerConfig : public std::false_type {};

        template <class DataT, int ChildN>
        struct IsLayerConfig<LayerConfig<DataT, ChildN>> : public std::true_type{};

        namespace {


            // content of each layer
            template <class Tag, int Level, class T>
            struct LayerContent {};

            template <class Tag, int Level, class DataT, int ChildN>
            struct LayerContent<Tag, Level, LayerConfig<DataT, ChildN>> {
                using TopoType = Topo<Tag, Level, ChildN>;
                using DataType = DataT;
                using TripletType = Triplet<TopoType, DataT>;
                using TableType = TripletArray<TopoType, DataT>;
                TableType contentTable;
                template <class Archive> inline void serialize(Archive & ar) { ar(contentTable); }
            };


            // make a layer content tuple with levels
            template <class Tag, class CConfTupleT, int Level>
            struct LayerContentFromConfigTuple {
                using type = LayerContent<Tag, Level, typename std::tuple_element<size_t(Level), CConfTupleT>::type>;
            };

            template <class Tag, class CConfTupleT, class SequenceT>
            struct LayerContentTuple {};

            template <class Tag, class CConfTupleT, int ...Levels>
            struct LayerContentTuple<Tag, CConfTupleT, Sequence<Levels...>> {
                using type = std::tuple<typename LayerContentFromConfigTuple<Tag, CConfTupleT, Levels>::type...>;
            };

        }


        // Homogeneous Graph
        template <class BaseDataT, class ...CConfs>
        class HomogeneousGraph {
            using Tag = HomogeneousGraph;

            // number of layers
            static const unsigned LayerNum = sizeof...(CConfs)+1;
            // tuple of all layers
            using ContentsType = typename LayerContentTuple<Tag, std::tuple<LayerConfig<BaseDataT, 0>, CConfs...>,
                typename SequenceGenerator<LayerNum>::type>::type;
            // layer content type at level Level
            template <int Level>
            struct LayerContentTypeStruct {
                using type = typename std::tuple_element<Level, ContentsType>::type;
            };
            // predicate of element existance at level Level
            template <int Level>
            using ElementExistsPred = TripletExistsPred<typename LayerContentTypeStruct<Level>::type::TopoType,
                typename LayerContentTypeStruct<Level>::type::DataType>;


        public:

            // element triplet type at level Level
            template <int Level>
            using TripletType = typename LayerContentTypeStruct<Level>::type::TableType::value_type;

            // internal elements
            template <int Level>
            inline const typename LayerContentTypeStruct<Level>::type::TableType & internalElements() const {
                return std::get<Level>(_contents).contentTable;
            }

            template <int Level>
            inline typename LayerContentTypeStruct<Level>::type::TableType & internalElements() {
                return std::get<Level>(_contents).contentTable;
            }


            // topo
            template <int Level>
            inline const typename LayerContentTypeStruct<Level>::type::TopoType & topo(HandleOfTypeAtLevel<Tag, Level> h) const {
                return internalElements<Level>()[h.id].topo;
            }

            template <int Level>
            inline typename LayerContentTypeStruct<Level>::type::TopoType & topo(HandleOfTypeAtLevel<Tag, Level> h) {
                return internalElements<Level>()[h.id].topo;
            }

            // data
            template <int Level>
            inline const typename LayerContentTypeStruct<Level>::type::DataType & data(HandleOfTypeAtLevel<Tag, Level> h) const {
                return internalElements<Level>()[h.id].data;
            }

            template <int Level>
            inline typename LayerContentTypeStruct<Level>::type::DataType & data(HandleOfTypeAtLevel<Tag, Level> h) {
                return internalElements<Level>()[h.id].data;
            }

            // traverse
            template <int Level>
            inline ConditionalContainerWrapper<typename LayerContentTypeStruct<Level>::type::TableType, ElementExistsPred<Level>> elements() {
                return ConditionalContainerWrapper<typename LayerContentTypeStruct<Level>::type::TableType, ElementExistsPred<Level>>(
                    &(internalElements<Level>()));
            }

            template <int Level>
            inline ConstConditionalContainerWrapper<typename LayerContentTypeStruct<Level>::type::TableType, ElementExistsPred<Level>> elements() const {
                return ConstConditionalContainerWrapper<typename LayerContentTypeStruct<Level>::type::TableType, ElementExistsPred<Level>>(
                    &(internalElements<Level>()));
            }


            // add element
            template <int Level>
            HandleOfTypeAtLevel<Tag, Level> add(std::initializer_list<HandleOfTypeAtLevel<Tag, Level - 1>> depends,
                const typename LayerContentTypeStruct<Level>::type::DataType & d = typename LayerContentTypeStruct<Level>::type::DataType()) {
                int id = static_cast<int>(internalElements<Level>().size());
                internalElements<Level>().emplace_back(typename LayerContentTypeStruct<Level>::type::TopoType(id, depends), d, true);
                for (const HandleOfTypeAtLevel<Tag, Level - 1> & lowh : depends) {
                    topo(lowh).uppers.insert(HandleOfTypeAtLevel<Tag, Level>(id));
                }
                return HandleOfTypeAtLevel<Tag, Level>(id);
            }

            template <int Level>
            HandleOfTypeAtLevel<Tag, Level> add(std::initializer_list<HandleOfTypeAtLevel<Tag, Level - 1>> depends,
                typename LayerContentTypeStruct<Level>::type::DataType && d) {
                int id = static_cast<int>(internalElements<Level>().size());
                internalElements<Level>().emplace_back(typename LayerContentTypeStruct<Level>::type::TopoType(id, depends),
                    std::forward<typename LayerContentTypeStruct<Level>::type::DataType>(d), true);
                for (const HandleOfTypeAtLevel<Tag, Level - 1> & lowh : depends) {
                    topo(lowh).uppers.insert(HandleOfTypeAtLevel<Tag, Level>(id));
                }
                return HandleOfTypeAtLevel<Tag, Level>(id);
            }

            // add element of the lowest level
            HandleOfTypeAtLevel<Tag, 0> add(const typename LayerContentTypeStruct<0>::type::DataType & d = LayerContentTypeStruct<0>::type::DataType()) {
                int id = static_cast<int>(internalElements<0>().size());
                internalElements<0>().emplace_back(typename LayerContentTypeStruct<0>::type::TopoType(id), d, true);
                return HandleOfTypeAtLevel<Tag, 0>(id);
            }

            HandleOfTypeAtLevel<Tag, 0> add(typename LayerContentTypeStruct<0>::type::DataType && d){
                int id = static_cast<int>(internalElements<0>().size());
                internalElements<0>().emplace_back(typename LayerContentTypeStruct<0>::type::TopoType(id),
                    std::forward<typename LayerContentTypeStruct<0>::type::DataType>(d), true);
                return HandleOfTypeAtLevel<Tag, 0>(id);
            }

            // removed
            template <int Level>
            inline bool removed(HandleOfTypeAtLevel<Tag, Level> h) const {
                return !internalElements<Level>()[h.id].exists;
            }

            // remove
            template <int Level>
            void remove(HandleOfTypeAtLevel<Tag, Level> h) {
                if (h.invalid() || removed(h))
                    return;
                cleanLowers<Level>(h, std::integral_constant<bool, (Level > 0)>());
                internalElements<Level>()[h.id].exists = false; // mark as deleted
                cleanUppers<Level>(h, std::integral_constant<bool, (Level < LayerNum - 1)>());
            }



        private:
            template <int Level>
            inline void cleanLowers(HandleOfTypeAtLevel<Tag, Level> h, std::false_type) {}
            template <int Level>
            void cleanLowers(HandleOfTypeAtLevel<Tag, Level> h, std::true_type) {
                auto & c = internalElements<Level>()[h.id];
                auto & clowerTable = internalElements<Level - 1>();

                for (auto & lowh : c.topo.lowers) { // remove this h from all lowers' uppers set
                    if (lowh.invalid() || removed(lowh))
                        continue;
                    auto & low = clowerTable[lowh.id];
                    low.topo.uppers.erase(h);
                }
            }

            template <int Level>
            inline void cleanUppers(HandleOfTypeAtLevel<Tag, Level> h, std::false_type) {}
            template <int Level>
            void cleanUppers(HandleOfTypeAtLevel<Tag, Level> h, std::true_type) {
                auto & c = internalElements<Level>()[h.id];
                for (auto & uph : c.topo.uppers) {
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
                return hasGarbageUsingSequence(typename SequenceGenerator<LayerNum>::type());
            }

            inline bool isDense() const { return !hasGarbage(); }

        private:
            template <int ...S>
            void gcUsingSequence(Sequence<S...>) {
                std::tuple<std::vector<HandleOfTypeAtLevel<Tag, S>>...> nlocs;
                int dummy[] = { RemoveAndMap(internalElements<S>(), std::get<S>(nlocs))... };
                int dummy2[] = { updateEachLayerHandles<S>(nlocs)... };
            }

            template <int Level, class NLocTupleT>
            int updateEachLayerHandles(const NLocTupleT & nlocs) {
                auto & eles = internalElements<Level>();
                for (size_t i = 0; i < eles.size(); i++){
                    UpdateOldHandle(std::get<Level>(nlocs), eles[i].topo.hd);
                }
                updateLowers<Level>(nlocs, std::integral_constant<bool, (Level > 0)>());
                updateUppers<Level>(nlocs, std::integral_constant<bool, (Level < LayerNum - 1)>());
                return 0;
            }

            template <int Level, class NLocTupleT>
            inline void updateLowers(const NLocTupleT & nlocs, std::false_type) {}
            template <int Level, class NLocTupleT>
            void updateLowers(const NLocTupleT & nlocs, std::true_type) {
                auto & eles = internalElements<Level>();
                for (size_t i = 0; i < eles.size(); i++){
                    UpdateOldHandleContainer(std::get<Level - 1>(nlocs), eles[i].topo.lowers);
                    RemoveInValidHandleFromContainer(eles[i].topo.lowers);
                }
            }

            template <int Level, class NLocTupleT>
            inline void updateUppers(const NLocTupleT & nlocs, std::false_type) {}
            template <int Level, class NLocTupleT>
            void updateUppers(const NLocTupleT & nlocs, std::true_type) {
                auto & eles = internalElements<Level>();
                for (size_t i = 0; i < eles.size(); i++){
                    UpdateOldHandleContainer(std::get<Level + 1>(nlocs), eles[i].topo.uppers);
                    RemoveInValidHandleFromContainer(eles[i].topo.uppers);
                }
            }

            template <int ...S>
            inline void clearUsingSequence(Sequence<S...>) {
                int dummy[] = { clearAtLevel<S>() ... };
            }

            template <int Level>
            inline int clearAtLevel() {
                internalElements<Level>().clear();
                return 0;
            }

            template <int ...S>
            inline bool hasGarbageUsingSequence(Sequence<S...>) const {
                for (bool r : { hasGarbageAtLevel<S>() ... }){
                    if (r)
                        return true;
                }
                return false;
            }

            template <int Level>
            inline bool hasGarbageAtLevel() const {
                for (const auto & t : internalElements<Level>()){
                    if (!t.exists)
                        return true;
                }
                return false;
            }

        public:
            template <class Archive> inline void serialize(Archive & ar) { ar(_contents); }

        private:
            ContentsType _contents;
        };

        template <class VertDataT, class EdgeDataT>
        using HomogeneousGraph02 = HomogeneousGraph<VertDataT, LayerConfig<EdgeDataT, 2>>;








    }

}
 
#endif