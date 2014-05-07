#ifndef PANORAMIX_CORE_SHAPE_HPP
#define PANORAMIX_CORE_SHAPE_HPP
 
namespace panoramix {
    namespace core {

    	static const int Dyn = 0;

        // shape
        template <int ... Ds>
        class Shape;

        // empty shape
        template <>
        class Shape<> {
        public:
            static const int StaticVolume = 1;
            static const int Rank = 0;
            static const int DynamicNum = 0;
            inline Shape(){}
            inline int volume() const { return 1; }
        };

        // shape with fixed dimension
        template <int D, int ... Ds>
        class Shape<D, Ds...> : private Shape<Ds...>{
        public:
            using BaseType = Shape<Ds...>;
            static const bool LengthIsDynamic = false;
            static const int StaticLength = D;
            static const int StaticVolume = BaseType::StaticVolume * D;
            static const int Rank = sizeof...(Ds)+1;
            static const int DynamicNum = BaseType::DynamicNum;

            // ignore fixed dimension in constructor
            template <class... Ints>
            inline explicit Shape(Ints ... ints)
                : BaseType(ints...){}

            // promote
            inline explicit Shape(const BaseType & s)
                : BaseType(s){}

            // base
            inline const BaseType& base() const {
                return (const BaseType&)(*this);
            }

            // first order size
            inline int length() const { return D; }

            // volume
            inline int volume() const { return ((const BaseType&)(*this)).volume() * length(); }

            // size
            template <int Ord> inline int size() const {
                return ((const BaseType&)(*this)).size<Ord - 1>();
            }
            template <> inline int size<0>() const { return length(); }

            // resize
            template <int Ord> inline void resize(int s) {
                return ((BaseType&)(*this)).resize<Ord - 1>(s);
            }
            template <> inline void resize<0>(int s) {
                static_assert(false, "unresizable on fixed dimension!");
            }

            // dynamic size
            template <int DOrd> inline int dsize() const {
                return ((const BaseType&)(*this)).dsize<DOrd>();
            }

            // dynamic resize
            template <int DOrd> inline void dresize(int ds) {
                ((BaseType&)(*this)).dresize<DOrd>(ds);
            }


            // static size
            template <int SOrd> inline int ssize() const {
                return ((const BaseType&)(*this)).ssize<SOrd - 1>();
            }
            template <> inline int ssize<0>() const {
                return D;
            }

        };

        // shape with dynamic dimension
        template <int ... Ds>
        class Shape<Dyn, Ds...> : private Shape<Ds...>{
        public:
            using BaseType = Shape<Ds...>;
            static const bool LengthIsDynamic = true;
            static const int StaticLength = Dyn;
            static const int StaticVolume = Dyn;
            static const int Rank = sizeof...(Ds)+1;
            static const int DynamicNum = BaseType::DynamicNum + 1;

            // set size for Dyn dimension
            template <class Int, class... Ints>
            inline explicit Shape(Int i, Ints ... ints)
                : BaseType(ints...), _length(i){}

            // promote
            inline Shape(int len, const Shape<Ds...> & s)
                : BaseType(s), _length(len){}



            // base
            inline const BaseType& base() const {
                return (const BaseType&)(*this);
            }

            // first order size
            inline int length() const { return _length; }

            // volume
            inline int volume() const { return ((const BaseType&)(*this)).volume() * length(); }

            // size
            template <int Ord> inline int size() const {
                return ((const BaseType&)(*this)).size<Ord - 1>();
            }
            template <> inline int size<0>() const { return _length; }

            // resize
            template <int Ord> inline void resize(int s) {
                return ((BaseType&)(*this)).resize<Ord - 1>(s);
            }
            template <> inline void resize<0>(int s) {
                _length = s;
            }

            // dynamic size
            template <int DOrd> inline int dsize() const {
                return ((const BaseType&)(*this)).dsize<DOrd - 1>();
            }
            template <> inline int dsize<0>() const {
                return _length;
            }

            // dynamic resize
            template <int DOrd> inline void dresize(int ds) {
                ((BaseType&)(*this)).dresize<DOrd - 1>(ds);
            }
            template <> inline void dresize<0>(int ds) {
                _length = ds;
            }

            // static size
            template <int SOrd> inline int ssize() const {
                return ((const BaseType&)(*this)).ssize<SOrd>();
            }


        private:
            int _length;
        };


    }
}
 
#endif