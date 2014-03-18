#ifndef PANORAMIX_CORE_UTIL_HPP
#define PANORAMIX_CORE_UTIL_HPP
 
namespace panoramix {
	namespace core {
 
		template <class Mat4T, class Vec3T>
		Mat4T Matrix4MakeLookAt(const Vec3T & eye, const Vec3T & center, const Vec3T & up, const Mat4T & base) {
			Vec3T n = (eye - center).normalized();
			Vec3T u = up.cross(n).normalized();
			Vec3T v = n.cross(u);
			Mat4T m;
			m <<
				u(0), v(0), n(0), 0,
				u(1), v(1), n(1), 0,
				u(2), v(2), n(2), 0,
				(-u).dot(eye), (-v).dot(eye), (-n).dot(eye), 1;
			return m * base;
		}

		template <class Mat4T, class ValueT>
		Mat4T Matrix4MakePerspective(const ValueT & fovyRadians, const ValueT & aspect, const ValueT & nearZ, const ValueT & farZ, const Mat4T & base) {
			ValueT cotan = ValueT(1.0) / tan(fovyRadians / 2.0);
			Mat4T m;
			m <<
				cotan / aspect, 0, 0, 0,
				0, cotan, 0, 0,
				0, 0, (farZ + nearZ) / (nearZ - farZ), -1,
				0, 0, (2 * farZ * nearZ) / (nearZ - farZ), 0;
			return m * base;
		}
 
	}
}
 
#endif