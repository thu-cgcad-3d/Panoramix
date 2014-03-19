#ifndef PANORAMIX_CORE_UTIL_HPP
#define PANORAMIX_CORE_UTIL_HPP
 
namespace panoramix {
	namespace core {



		template <class T>
		inline T Square(const T & v) {
			return v * v;
		}

		template <class Vec3T, class ValueT>
		inline void LongitudeLatitudeFromDirection(const Vec3T & d, ValueT & longi, ValueT & lati) {
			longi = atan2(d(1), d(0));
			lati = atan(d(2) / sqrt(Square(d(1)) + Square(d(0))));
		}

		template <class Vec3T, class ValueT>
		inline void DirectionFromLongitudeLatitude(const ValueT & longi, const ValueT & lati, Vec3T & d) {
			d = Vec3T(cos(longi)*cos(lati), sin(longi)*cos(lati), sin(lati));
		}

		template <class T, class K>
		inline T WrapBetween(const T& input, const K& low, const K& high) {
			if (low >= high)
				return input;
			if (low <= input && input < high)
				return input;
			const K sz = high - low;
			return input - int((input - low) / sz) * sz + (input < low ? sz : 0);
		}
 
		template <class Mat4T, class Vec3T>
		Mat4T Matrix4MakeLookAt(const Vec3T & eye, const Vec3T & center, const Vec3T & up, const Mat4T & base) {
			Vec3T zaxis = (center - eye).normalized();
			Vec3T xaxis = up.cross(zaxis).normalized();
			Vec3T yaxis = zaxis.cross(xaxis);
			Mat4T m;
			m <<
				xaxis(0), yaxis(0), zaxis(0), 0,
				xaxis(1), yaxis(1), zaxis(1), 0,
				xaxis(2), yaxis(2), zaxis(2), 0,
				-xaxis.dot(eye), -yaxis.dot(eye), -zaxis.dot(eye), 1;
			return m.transpose() * base;
		}

		template <class Mat4T, class ValueT>
		Mat4T Matrix4MakePerspective(const ValueT & fovyRadians, const ValueT & aspect, const ValueT & nearZ, const ValueT & farZ, const Mat4T & base) {
			//assert(fovyRadians > 0 && aspect != 0);
			ValueT cotan = ValueT(1.0) / std::tan(fovyRadians / 2.0);
			Mat4T m;
			m <<
				cotan / aspect, 0, 0, 0,
				0, cotan, 0, 0,
				0, 0, (farZ + nearZ) / (nearZ - farZ), -1,
				0, 0, (2 * farZ * nearZ) / (nearZ - farZ), 0;
			return m.transpose() * base;
		}

	}
}
 
#endif