#ifndef PANORAMIX_CORE_MESH_HPP
#define PANORAMIX_CORE_MESH_HPP

#include "basic_types.hpp"
#include "misc.hpp"

namespace panoramix {
    namespace core {
        
        struct Dummy {};
        
        /**
         * @brief Handle struct
         */
        template <class Tag>
        struct Handle {
            int64_t id;
            inline Handle(int64_t id_ = -1) : id(id_){}
            inline bool operator == (Handle h) const { return id == h.id; }
            inline bool operator != (Handle h) const { return id != h.id; }
            inline void reset() { id = -1; }
            inline bool isValid() const { return id >= 0; }
            inline bool isInValid() const { return id < 0; }
        };
        template <class Tag>
        using HandleArray = std::vector<Handle<Tag>>;
        template <class Tag>
        using HandlePtrArray = std::vector<Handle<Tag>*>;
        
        /**
         * @brief Topology structs
         */
        struct VertTopo;
        struct HalfTopo;
        struct FaceTopo;
        struct VertTopo {
            Handle<VertTopo> hd;
            HandleArray<HalfTopo> halfedges;
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
        };
        struct FaceTopo {
            Handle<FaceTopo> hd;
            HandleArray<HalfTopo> halfedges;
        };
        
        /**
         * @brief Triplet struct
         */
        template <class TopoT, class DataT>
        struct Triplet {
            TopoT topo;
            uint16_t exists;
            DataT data;
        };
        template <class TopoT, class DataT>
        using TripletArray = std::vector<Triplet<TopoT, DataT>>;
        
        
        /**
         * @brief Helper functions for the Mesh class
         */
        template <class ComponentTableT, class UpdateHandleTableT>
        void RemoveAndMap(ComponentTableT & v, UpdateHandleTableT & newlocations) {
            // ComponentTableT : std::vector<Triplet<TopoT, DataT>>
            // UpdateHandleTableT: std::vector<Handle<TopoT>>
            newlocations.resize(v.size());
            int64_t index = 0;
            for (size_t i = 0; i < v.size(); i++){
                newlocations[i] = { v[i].exists == false ? -1 : (index++) };
            }
            for (int i = int(v.size() - 1); i >= 0; --i){
                if (!v[i].exists){
                    v.erase(v.begin() + i);
                }
            }
        }
        template <class UpdateHandleTableT, class TopoT>
        inline void UpdateOldHandle(const UpdateHandleTableT & newlocationTable, Handle<TopoT> & h) {
            // UpdateHandleTableT: std::vector<Handle<TopoT>>
            h = newlocationTable[h.id];
        }
        template <class UpdateHandleTableT, class ContainerT>
        inline void UpdateOldHandleContainer(const UpdateHandleTableT& newlocationTable, ContainerT & hs) {
            // UpdateHandleTableT: std::vector<Handle<TopoT>>
            for (auto & h : hs){
                h = newlocationTable[h.id];
            }
        }
        template <class ContainerT>
        inline void RemoveInValidHandleFromContainer(ContainerT & hs) {
            auto invalid = typename std::iterator_traits<decltype(std::begin(hs))>::value_type();
            invalid.reset();
            hs.erase(std::remove(std::begin(hs), std::end(hs), invalid), std::end(hs));
        }
        
        /**
         * @brief The Mesh class, halfedge structure
         */
        template <class VertDataT, class HalfDataT = Dummy, class FaceDataT = Dummy>
        class Mesh {
            
        public:
            inline Mesh(){}
            
            using VertHandle = Handle<VertTopo>;
            using HalfHandle = Handle<HalfTopo>;
            using FaceHandle = Handle<FaceTopo>;
            using VertsTable = TripletArray<VertTopo, VertDataT>;
            using HalfsTable = TripletArray<HalfTopo, HalfDataT>;
            using FacesTable = TripletArray<FaceTopo, FaceDataT>;
            using Vertex = typename VertsTable::value_type;
            using HalfEdge = typename HalfsTable::value_type;
            using Face = typename FacesTable::value_type;
            
            inline VertsTable & internalVertices() { return _verts; }
            inline HalfsTable & internalHalfEdges() { return _halfs; }
            inline FacesTable & internalFaces() { return _faces; }
            inline const VertsTable & internalVertices() const { return _verts; }
            inline const HalfsTable & internalHalfEdges() const { return _halfs; }
            inline const FacesTable & internalFaces() const { return _faces; }
            
            inline JumpContainerWrapper<VertsTable> vertices() { return JumpContainerWrapper<VertsTable>(&_verts); }
            inline JumpContainerWrapper<HalfsTable> halfedges() { return JumpContainerWrapper<HalfsTable>(&_halfs); }
            inline JumpContainerWrapper<FacesTable> faces() { return JumpContainerWrapper<FacesTable>(&_faces); }
            inline ConstJumpContainerWrapper<VertsTable> vertices() const { return ConstJumpContainerWrapper<VertsTable>(&_verts); }
            inline ConstJumpContainerWrapper<HalfsTable> halfedges() const { return ConstJumpContainerWrapper<HalfsTable>(&_halfs); }
            inline ConstJumpContainerWrapper<FacesTable> faces() const { return ConstJumpContainerWrapper<FacesTable>(&_faces); }
            
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
                               const HalfDataT & hd = HalfDataT(), const HalfDataT & hdrev = HalfDataT());
            FaceHandle addFace(const HandleArray<HalfTopo> & halfedges,
                               const FaceDataT & fd = FaceDataT());
            FaceHandle addFace(const HandleArray<VertTopo> & vertices, bool autoflip = true,
                               const FaceDataT & fd = FaceDataT());
            
            HalfHandle findEdge(VertHandle from, VertHandle to) const;
            
            inline bool removed(FaceHandle f) const { return !_faces[f.id].exists; }
            inline bool removed(HalfHandle e) const { return !_halfs[e.id].exists; }
            inline bool removed(VertHandle v) const { return !_verts[v.id].exists; }
            
            inline void remove(FaceHandle f);
            inline void remove(HalfHandle e);
            inline void remove(VertHandle v);
            
            Mesh & unite(const Mesh & m);
            
            /**
             * @brief garbage collection
             */
            template <class VertHandlePtrContainerT = HandlePtrArray<VertTopo>,
            class HalfHandlePtrContainerT = HandlePtrArray<HalfTopo>,
            class FaceHandlePtrContainerT = HandlePtrArray<FaceTopo>
            >
            void gc(const VertHandlePtrContainerT & vps = VertHandlePtrContainerT(),
                    const HalfHandlePtrContainerT & hps = HalfHandlePtrContainerT(),
                    const FaceHandlePtrContainerT & fps = FaceHandlePtrContainerT());
            
            void clear();
            
        private:
            VertsTable _verts;
            HalfsTable _halfs;
            FacesTable _faces;
        };
        
        // implementation of class Mesh
        
        template <class VertDataT, class HalfDataT, class FaceDataT>
        typename Mesh<VertDataT, HalfDataT, FaceDataT>::VertHandle
        Mesh<VertDataT, HalfDataT, FaceDataT>::addVertex(const VertDataT & vd) {
            Triplet<VertTopo, VertDataT> t;
            t.data = vd;
            t.topo.hd.id = _verts.size();
            t.exists = true;
            _verts.push_back(t);
            return _verts.back().topo.hd;
        }
        
        template <class VertDataT, class HalfDataT, class FaceDataT>
        typename Mesh<VertDataT, HalfDataT, FaceDataT>::HalfHandle
        Mesh<VertDataT, HalfDataT, FaceDataT>::addEdge(VertHandle from, VertHandle to, const HalfDataT & hd, const HalfDataT & hdrev) {
            if (from == to){
                return HalfHandle();
            }
            // find existed halfedge
            HalfHandle hh = findEdge(from, to);
            if (hh.isValid()){
                _halfs[hh.id].data = hd;
                _halfs[_halfs[hh.id].topo.opposite.id].data = hdrev;
                return hh;
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
            if (hh.isValid() && _halfs[hh.id].topo.face.isValid() && autoflip){
                std::reverse(verts.begin(), verts.end());
            }
            
            for (size_t i = 0; i < verts.size(); i++){
                size_t inext = (i + 1) % verts.size();
                halfs.push_back(addEdge(verts[i], verts[inext]));
            }
            return addFace(halfs, fd);
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
            if (f.isInValid() || removed(f))
                return;
            _faces[f.id].exists = false;
            for(auto & hh : _faces[f.id].topo.halfedges){
                hh.reset();
            }
        }
        
        template <class VertDataT, class HalfDataT, class FaceDataT>
        void Mesh<VertDataT, HalfDataT, FaceDataT>::remove(HalfHandle h) {
            if (h.isInValid() || removed(h))
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
            if (v.isInValid() || removed(v))
                return;
            _verts[v.id].exists = false;
            for (HalfHandle hh : _verts[v.id].topo.halfedges)
                remove(hh);
            _verts[v.id].topo.halfedges.clear();
        }
        
        template <class VertDataT, class HalfDataT, class FaceDataT>
        Mesh<VertDataT, HalfDataT, FaceDataT>& Mesh<VertDataT, HalfDataT, FaceDataT>::unite(const Mesh & m) {
            std::vector<VertHandle> vtable(m.internalVertices().size());
            std::vector<HalfHandle> htable(m.internalHalfEdges().size());
            std::vector<FaceHandle> ftable(m.internalFaces().size());
            
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
        
    }
}

#endif