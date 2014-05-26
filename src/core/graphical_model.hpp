#ifndef PANORAMIX_CORE_GRAPHICAL_MODEL_HPP
#define PANORAMIX_CORE_GRAPHICAL_MODEL_HPP

#include <cstdint>
#include <array>
#include <set>
#include <vector>

#include <Eigen/Core>

#include "template_utilities.hpp"
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
        template <class Tag>
        struct CompareHandle {
            inline bool operator()(Handle<Tag> a, Handle<Tag> b) const {
                return a.id < b.id;
            }
        };
        template <class Tag>
        struct HashHandle {
            inline uint64_t operator()(Handle<Tag> a) const {
                return static_cast<uint64_t>(a.id);
            }
        };

        template <int L>
        struct AtLevel {
            static const int Level = L;
        };
        template <int L>
        using HandleAtLevel = Handle<AtLevel<L>>;


        
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
            uint8_t exists;
            DataT data;
            inline Triplet(){}
            inline Triplet(const TopoT & t, const DataT & d, uint8_t e = true) 
                : topo(t), exists(e), data(d){}
        };
        template <class TopoT, class DataT>
        struct TripletExistsPred {
            inline uint8_t operator()(const Triplet<TopoT, DataT> & t) const { 
                return t.exists; 
            }
        };
        template <class TopoT, class DataT>
        using TripletArray = std::vector<Triplet<TopoT, DataT>>;
        
        
        /**
         * @brief Helper functions for the Mesh class
         */
        template <class ComponentTableT, class UpdateHandleTableT>
        int RemoveAndMap(ComponentTableT & v, UpdateHandleTableT & newlocations) {
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
            return 0;
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
            static const int LayerNum = 2;
            
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
                if (hh.isValid()){
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





        /**
         * @brief The ConstraintGraph class
         */
        struct ConstraintTopo;
        struct ComponentTopo {
            Handle<ComponentTopo> hd;
            HandleArray<ConstraintTopo> constraints;
        };
        struct ConstraintTopo {
            Handle<ConstraintTopo> hd;
            HandleArray<ComponentTopo> components;
        };

        template <class ComponentDataT, class ConstraintDataT>
        class ConstraintGraph {
        public:
            inline ConstraintGraph(){}

            using ComponentHandle = Handle<ComponentTopo>;
            using ConstraintHandle = Handle<ConstraintTopo>;

            using ComponentsTable = TripletArray<ComponentTopo, ComponentDataT>;
            using ConstraintsTable = TripletArray<ConstraintTopo, ConstraintDataT>;

            using ComponentExistsPred = TripletExistsPred<ComponentTopo, ComponentDataT>;
            using ConstraintExistsPred = TripletExistsPred<ConstraintTopo, ConstraintDataT>;

            using Component = typename ComponentsTable::value_type;
            using Constraint = typename ConstraintsTable::value_type;

            inline ComponentsTable & internalComponents() { return _components; }
            inline ConstraintsTable & internalConstraints() { return _constraints; }
            inline const ComponentsTable & internalComponents() const { return _components; }
            inline const ConstraintsTable & internalConstraints() const { return _constraints; }

            inline ConditionalContainerWrapper<ComponentsTable, ComponentExistsPred> components() {
                return ConditionalContainerWrapper<ComponentsTable, ComponentExistsPred>(&_components);
            }
            inline ConditionalContainerWrapper<ConstraintsTable, ConstraintExistsPred> constraints() {
                return ConditionalContainerWrapper<ConstraintsTable, ComponentExistsPred>(&_constraints);
            }
            inline ConstConditionalContainerWrapper<ComponentsTable, ComponentExistsPred> components() const {
                return ConstConditionalContainerWrapper<ComponentsTable, ComponentExistsPred>(&_components);
            }
            inline ConstConditionalContainerWrapper<ConstraintsTable, ConstraintExistsPred> constraints() const {
                return ConstConditionalContainerWrapper<ConstraintsTable, ComponentExistsPred>(&_constraints);
            }


            inline ComponentTopo & topo(ComponentHandle h) { return _components[h.id].topo; }
            inline ConstraintTopo & topo(ConstraintHandle h) { return _constraints[h.id].topo; }
            inline const ComponentTopo & topo(ComponentHandle h) const { return _components[h.id].topo; }
            inline const ConstraintTopo & topo(ConstraintHandle h) const { return _constraints[h.id].topo; }

            inline ComponentDataT & data(ComponentHandle h) { return _components[h.id].data; }
            inline ConstraintDataT & data(ConstraintHandle h) { return _constraints[h.id].data; }
            inline const ComponentDataT & data(ComponentHandle h) const { return _components[h.id].data; }
            inline const ConstraintDataT & data(ConstraintHandle h) const { return _constraints[h.id].data; }

            ComponentHandle addComponent(const ComponentDataT & compData = ComponentDataT());
            ConstraintHandle addConstraint(const HandleArray<ComponentTopo> & components, const ConstraintDataT & consData = ConstraintDataT());

            inline bool removed(ComponentHandle f) const { return !_components[f.id].exists; }
            inline bool removed(ConstraintHandle e) const { return !_constraints[e.id].exists; }

            inline void remove(ComponentHandle f);
            inline void remove(ConstraintHandle e);

            /**
            * @brief garbage collection
            */
            template <class ComponentHandlePtrContainerT = HandlePtrArray<ComponentTopo>,
                class ConstraintHandlePtrContainerT = HandlePtrArray<ConstraintTopo>>
            void gc(const ComponentHandlePtrContainerT & compps = ComponentHandlePtrContainerT(),
            const ConstraintHandlePtrContainerT & consps = ConstraintHandlePtrContainerT());

            void clear();

        private:
            ComponentsTable _components;
            ConstraintsTable _constraints;
        };


        // implementation of ConstraintGraph
        template <class ComponentDataT, class ConstraintDataT>
        typename ConstraintGraph<ComponentDataT, ConstraintDataT>::ComponentHandle 
            ConstraintGraph<ComponentDataT, ConstraintDataT>::addComponent(const ComponentDataT & compData) {
            ComponentTopo topo;
            topo.hd.id = _components.size();
            _components.emplace_back(std::move(topo), compData, true);
            return _components.back().topo.hd;
        }

        template <class ComponentDataT, class ConstraintDataT>
        typename ConstraintGraph<ComponentDataT, ConstraintDataT>::ConstraintHandle
            ConstraintGraph<ComponentDataT, ConstraintDataT>::addConstraint(const HandleArray<ComponentTopo> & components, const ConstraintDataT & consData) {
            ConstraintTopo topo;
            topo.hd.id = _constraints.size();
            topo.components = components;
            for (auto & component : components) {
                _components[component.id].topo.constraints.push_back(topo.hd);
            }
            _constraints.emplace_back(std::move(topo), consData, true);
            return _constraints.back().topo.hd;
        }
        
        template <class ComponentDataT, class ConstraintDataT>
        void ConstraintGraph<ComponentDataT, ConstraintDataT>::remove(ConstraintHandle h) {
            if (h.isInValid() || removed(h))
                return;
            _constraints[h.id].exists = false;
            for (auto & comp : _constraints[h.id].topo.components){
                comp.reset();
            }
        }

        template <class ComponentDataT, class ConstraintDataT>
        void ConstraintGraph<ComponentDataT, ConstraintDataT>::remove(ComponentHandle h) {
            if (h.isInValid() || removed(h))
                return;
            _components[h.id].exists = false;
            for (ConstraintHandle hh : _components[h.id].topo.constraints)
                remove(hh);
            _components[h.id].topo.constraints.clear();
        }


        template <class ComponentDataT, class ConstraintDataT>
        template <class ComponentHandlePtrContainerT, class ConstraintHandlePtrContainerT>
        void ConstraintGraph<ComponentDataT, ConstraintDataT>::gc(const ComponentHandlePtrContainerT & compps,
            const ConstraintHandlePtrContainerT & consps){
            std::vector<ComponentHandle> vnlocs;
            std::vector<ConstraintHandle> hnlocs;
            RemoveAndMap(_components, vnlocs);
            RemoveAndMap(_constraints, hnlocs);

            for (size_t i = 0; i < _components.size(); i++){
                UpdateOldHandle(vnlocs, _components[i].topo.hd);
                UpdateOldHandleContainer(hnlocs, _components[i].topo.constraints);
                RemoveInValidHandleFromContainer(_components[i].topo.constraints);
            }
            for (size_t i = 0; i < _constraints.size(); i++){
                UpdateOldHandle(hnlocs, _constraints[i].topo.hd);
                UpdateOldHandleContainer(vnlocs, _constraints[i].topo.components);
                RemoveInValidHandleFromContainer(_constraints[i].topo.components);
            }

            for (auto compp : compps){
                UpdateOldHandle(vnlocs, *compp);
            }
            for (auto consp : consps){
                UpdateOldHandle(hnlocs, *consp);
            }
        }

        template <class ComponentDataT, class ConstraintDataT>
        void ConstraintGraph<ComponentDataT, ConstraintDataT>::clear() {
            _constraints.clear();
            _components.clear();
        }






        using Eigen::Dynamic;

        // fix sized
        template <int L, int ChildN>
        struct Topo {
            static const int Level = L;
            std::array<Handle<AtLevel<Level - 1>>, ChildN> lowers; // use std::array
            std::set<Handle<AtLevel<Level + 1>>, CompareHandle<AtLevel<Level + 1>>> uppers;
            HandleAtLevel<Level>hd;
            explicit inline Topo(int id = -1) : hd(id){}
            explicit inline Topo(int id, std::initializer_list<Handle<AtLevel<Level - 1>>> ls) : hd(id) {
                assert(ls.size() == lowers.size());
                std::copy(ls.begin(), ls.end(), lowers.begin());
            }
        };

        // dynamic sized
        template <int L>
        struct Topo<L, Dynamic> {
            static const int Level = L;
            std::vector<Handle<AtLevel<Level - 1>>> lowers; // use std::array
            std::set<Handle<AtLevel<Level + 1>>, CompareHandle<AtLevel<Level + 1>>> uppers;
            HandleAtLevel<Level>hd;
            explicit inline Topo(int id = -1) : hd(id){}
            explicit inline Topo(int id, std::initializer_list<Handle<AtLevel<Level - 1>>> ls)
                : hd(id), lowers(ls) {}
        };


        // zero sized
        template <int L>
        struct Topo<L, 0> {
            static const int Level = L;
            std::set<Handle<AtLevel<Level + 1>>, CompareHandle<AtLevel<Level + 1>>> uppers;
            HandleAtLevel<Level>hd;
            explicit inline Topo(int id = -1) : hd(id){}
            explicit inline Topo(int id, std::initializer_list<Handle<AtLevel<Level - 1>>> ls)
                : hd(id) {
                assert(ls.size() == 0);
            }
        };

        // configuration for each layer
        template <class DataT, int ChildN>
        struct LayerConfig {
            using DataType = DataT;
            static const int ChildNum = ChildN;
        };

        template <class T>
        struct IsLayerConfig : public std::false_type {};

        template <class DataT, int ChildN>
        struct IsLayerConfig<LayerConfig<DataT, ChildN>> : public std::true_type{};

        // content of each layer
        template <int Level, class T>
        struct LayerContent {};

        template <int Level, class DataT, int ChildN>
        struct LayerContent<Level, LayerConfig<DataT, ChildN>> {
            using TopoType = Topo<Level, ChildN>;
            using DataType = DataT;
            using TripletType = Triplet<TopoType, DataT>;
            using TableType = TripletArray<TopoType, DataT>;
            TableType contentTable;
        };

        template <class CConfTupleT, int Level>
        struct LayerContentFromConfigTuple {
            using type = LayerContent<Level, typename std::tuple_element<size_t(Level), CConfTupleT>::type>;
        };

        // make a layer content tuple with levels
        template <class CConfTupleT, class SequenceT>
        struct LayerContentTuple {};

        template <class CConfTupleT, int ...Levels>
        struct LayerContentTuple<CConfTupleT, Sequence<Levels...>> {
            using type = std::tuple<typename LayerContentFromConfigTuple<CConfTupleT, Levels>::type...>;
        };



        // Graphical Model
        template <class BaseDataT, class ...CConfs>
        class GraphicalModel {
            static const unsigned LayerNum = sizeof...(CConfs)+1;

            using ContentsType =
                typename LayerContentTuple<std::tuple<LayerConfig<BaseDataT, 0>, CConfs...>,
                typename SequenceGenerator<LayerNum>::type>::type;

            template <int Level>
            struct LayerContentTypeStruct {
                using type = typename std::tuple_element<Level, ContentsType>::type;
            };

            template <int Level>
            using ElementExistsPred = TripletExistsPred<typename LayerContentTypeStruct<Level>::type::TopoType,
                typename LayerContentTypeStruct<Level>::type::DataType>;

        public:

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
            inline const typename LayerContentTypeStruct<Level>::type::TopoType & topo(HandleAtLevel<Level> h) const {
                return internalElements<Level>()[h.id].topo;
            }

            template <int Level>
            inline typename LayerContentTypeStruct<Level>::type::TopoType & topo(HandleAtLevel<Level> h) {
                return internalElements<Level>()[h.id].topo;
            }

            // data
            template <int Level>
            inline const typename LayerContentTypeStruct<Level>::type::DataType & data(HandleAtLevel<Level> h) const {
                return internalElements<Level>()[h.id].data;
            }

            template <int Level>
            inline typename LayerContentTypeStruct<Level>::type::DataType & data(HandleAtLevel<Level> h) {
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
            HandleAtLevel<Level> add(std::initializer_list<HandleAtLevel<Level - 1>> depends,
                const typename LayerContentTypeStruct<Level>::type::DataType & d = typename LayerContentTypeStruct<Level>::type::DataType()) {
                int id = internalElements<Level>().size();
                internalElements<Level>().emplace_back(typename LayerContentTypeStruct<Level>::type::TopoType(id, depends), d, true);
                for (const HandleAtLevel<Level - 1> & lowh : depends) {
                    topo(lowh).uppers.insert(HandleAtLevel<Level>(id));
                }
                return HandleAtLevel<Level>(id);
            }

            // add element of the lowest level
            HandleAtLevel<0> add(const typename LayerContentTypeStruct<0>::type::DataType & d = LayerContentTypeStruct<0>::type::DataType()) {
                int id = internalElements<0>().size();
                internalElements<0>().emplace_back(typename LayerContentTypeStruct<0>::type::TopoType(id), d, true);
                return HandleAtLevel<0>(id);
            }

            // removed
            template <int Level>
            inline bool removed(HandleAtLevel<Level> h) const {
                return !internalElements<Level>()[h.id].exists;
            }

            // remove
            template <int Level>
            void remove(HandleAtLevel<Level> h) {
                if (h.isInValid() || removed(h))
                    return;
                cleanLowers<Level>(h, std::integral_constant<bool, (Level>0) > ());
                internalElements<Level>()[h.id].exists = false; // mark as deleted
                cleanUppers<Level>(h, std::integral_constant < bool, (Level<LayerNum - 1)>());
            }

        private:
            template <int Level>
            inline void cleanLowers(HandleAtLevel<Level> h, std::false_type) {}
            template <int Level>
            void cleanLowers(HandleAtLevel<Level> h, std::true_type) {
                auto & c = internalElements<Level>()[h.id];
                auto & clowerTable = internalElements<Level - 1>();

                for (auto & lowh : c.topo.lowers) { // remove this h from all lowers' uppers set
                    if (lowh.isInValid() || removed(lowh))
                        continue;
                    auto & low = clowerTable[lowh.id];
                    low.topo.uppers.erase(h);
                }
            }

            template <int Level>
            inline void cleanUppers(HandleAtLevel<Level> h, std::false_type) {}
            template <int Level>
            void cleanUppers(HandleAtLevel<Level> h, std::true_type) {
                auto & c = internalElements<Level>()[h.id];
                for (auto & uph : c.topo.uppers) {
                    remove(uph);
                }
            }


        public:
            // garbage collection
            inline void gc() {
                gcUsingSequence(typename SequenceGenerator<LayerNum>::type());
            }

        private:
            template <int ...S>
            void gcUsingSequence(Sequence<S...>) {
                std::tuple<std::vector<HandleAtLevel<S>>...> nlocs;
                int dummy[] = { RemoveAndMap(internalElements<S>(), std::get<S>(nlocs))... };

            }

            template <int Level, class NLocTupleT>
            inline void updateEachLayerHandles(const NLocTupleT & nlocs) {
                auto & eles = internalElements<Level>();
                for (size_t i = 0; i < eles.size(); i++){
                    UpdateOldHandle(std::get<Level>(nlocs), eles[i].topo.hd);
                }
                updateLowers<Level>(nlocs, std::integral_constant<bool, (Level > 0)>());
                updateUppers<Level>(nlocs, std::integral_constant<bool, (Level<LayerNum - 1)>());
            }

            template <int Level, class NLocTupleT>
            inline void updateLowers(const NLocTupleT & nlocs, std::false_type) {}
            template <int Level, class NLocTupleT>
            inline void updateLowers(const NLocTupleT & nlocs, std::true_type) {
                auto & eles = internalElements<Level>();
                for (size_t i = 0; i < eles.size(); i++){
                    UpdateOldHandleContainer(std::get<Level - 1>(nlocs), eles[i].topo.lowers);
                    RemoveInValidHandleFromContainer(eles[i].topo.lowers);
                }
            }

            template <int Level, class NLocTupleT>
            inline void updateUppers(const NLocTupleT & nlocs, std::false_type) {}
            template <int Level, class NLocTupleT>
            inline void updateUppers(const NLocTupleT & nlocs, std::true_type) {
                auto & eles = internalElements<Level>();
                for (size_t i = 0; i < eles.size(); i++){
                    UpdateOldHandleContainer(std::get<Level + 1>(nlocs), eles[i].topo.uppers);
                    RemoveInValidHandleFromContainer(eles[i].topo.uppers);
                }
            }

        private:
            ContentsType _contents;
        };






    }
}

#endif