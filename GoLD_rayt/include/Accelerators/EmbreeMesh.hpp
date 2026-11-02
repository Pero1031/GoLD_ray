/**
 * @file Accelerators/EmbreeMesh.hpp
 * @brief Triangle mesh intersection with Intel Embree 4 (header-only).
 *
 * Embree works in single precision and is used here only to FIND the hit: the closest triangle
 * and its barycentric coordinates. The hit record (position, normals, material) is then built by
 * the Triangle itself in Real precision (Triangle::fillHit), so shading code is unchanged and
 * this works in a double build as well as in a float build (Real = float in Core/Precision.hpp).
 *
 * Notes:
 *  - Rays are converted to float; with scene scales around 1 the 1e-5 ray offset used by
 *    SurfaceInteraction::SpawnRay is far above float round-off (about 1e-7 * distance).
 *  - The scene is built with the robust (watertight) flag and high build quality.
 *  - Intended as a drop-in for BuildMeshBVH(); the rest of the scene (spheres, quads) stays in an
 *    Aggregate / BVH of its own.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include <embree4/rtcore.h>

#include "Geometry/TriangleMesh.hpp"

namespace rayt {

    namespace detail {
        /// One process-wide Embree device (created on first use, released at exit).
        inline RTCDevice embreeDevice() {
            struct Holder {
                RTCDevice device = nullptr;
                Holder() { device = rtcNewDevice(nullptr); }
                ~Holder() { if (device) rtcReleaseDevice(device); }
            };
            static Holder holder;
            return holder.device;
        }
    } // namespace detail

    class EmbreeMesh : public Hittable {
    public:
        EmbreeMesh(const std::shared_ptr<const TriangleMesh>& mesh, const std::shared_ptr<Material>& material) {
            RTCDevice device = detail::embreeDevice();
            if (!device) throw std::runtime_error("Embree: rtcNewDevice failed (is embree4.dll next to the exe?)");

            for (size_t i = 0; i < mesh->numTriangles(); ++i)
                m_tris.push_back(std::make_shared<Triangle>(mesh, i, material));

            m_scene = rtcNewScene(device);
            rtcSetSceneFlags(m_scene, RTC_SCENE_FLAG_ROBUST);
            rtcSetSceneBuildQuality(m_scene, RTC_BUILD_QUALITY_HIGH);

            RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
            float* vb = static_cast<float*>(rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0,
                RTC_FORMAT_FLOAT3, 3 * sizeof(float), mesh->p.size()));
            for (size_t i = 0; i < mesh->p.size(); ++i) {
                vb[3 * i] = static_cast<float>(mesh->p[i].x);
                vb[3 * i + 1] = static_cast<float>(mesh->p[i].y);
                vb[3 * i + 2] = static_cast<float>(mesh->p[i].z);
            }
            unsigned* ib = static_cast<unsigned*>(rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0,
                RTC_FORMAT_UINT3, 3 * sizeof(unsigned), mesh->numTriangles()));
            for (size_t i = 0; i < mesh->indices.size(); ++i) ib[i] = static_cast<unsigned>(mesh->indices[i]);

            rtcCommitGeometry(geom);
            rtcAttachGeometry(m_scene, geom);
            rtcReleaseGeometry(geom);
            rtcCommitScene(m_scene);

            RTCBounds b;
            rtcGetSceneBounds(m_scene, &b);
            m_bounds = AABB(Vector3(b.lower_x, b.lower_y, b.lower_z), Vector3(b.upper_x, b.upper_y, b.upper_z));
        }

        ~EmbreeMesh() override { if (m_scene) rtcReleaseScene(m_scene); }

        EmbreeMesh(const EmbreeMesh&) = delete;
        EmbreeMesh& operator=(const EmbreeMesh&) = delete;

        bool hit(const Ray& r, SurfaceInteraction& rec) const override {
            RTCRayHit rh;
            setRay(rh.ray, r);
            rh.hit.geomID = RTC_INVALID_GEOMETRY_ID;
            rh.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
            rtcIntersect1(m_scene, &rh);
            if (rh.hit.geomID == RTC_INVALID_GEOMETRY_ID) return false;

            m_tris[rh.hit.primID]->fillHit(r, Real(rh.ray.tfar), Real(rh.hit.u), Real(rh.hit.v), rec);
            return true;
        }

        bool occluded(const Ray& r) const override {
            RTCRay ray;
            setRay(ray, r);
            rtcOccluded1(m_scene, &ray);
            return ray.tfar < 0.0f;   // Embree sets tfar = -inf on a hit
        }

        bool intersectT(const Ray& r, Real& t) const override {
            RTCRayHit rh;
            setRay(rh.ray, r);
            rh.hit.geomID = RTC_INVALID_GEOMETRY_ID;
            rh.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
            rtcIntersect1(m_scene, &rh);
            if (rh.hit.geomID == RTC_INVALID_GEOMETRY_ID) return false;
            t = Real(rh.ray.tfar);
            return true;
        }

        AABB bounds() const override { return m_bounds; }

        size_t numTriangles() const { return m_tris.size(); }
        const std::vector<std::shared_ptr<Triangle>>& triangles() const { return m_tris; }

    private:
        static void setRay(RTCRay& ray, const Ray& r) {
            ray.org_x = static_cast<float>(r.o.x);
            ray.org_y = static_cast<float>(r.o.y);
            ray.org_z = static_cast<float>(r.o.z);
            ray.dir_x = static_cast<float>(r.d.x);
            ray.dir_y = static_cast<float>(r.d.y);
            ray.dir_z = static_cast<float>(r.d.z);
            ray.tnear = static_cast<float>(std::max(r.tMin, Real(0)));
            ray.tfar = (r.tMax >= static_cast<Real>(std::numeric_limits<float>::max()))
                ? std::numeric_limits<float>::infinity() : static_cast<float>(r.tMax);
            ray.time = 0.0f;
            ray.mask = static_cast<unsigned>(-1);
            ray.id = 0;
            ray.flags = 0;
        }

        RTCScene m_scene = nullptr;
        std::vector<std::shared_ptr<Triangle>> m_tris;
        AABB m_bounds;
    };

} // namespace rayt
