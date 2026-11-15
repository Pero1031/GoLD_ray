/**
 * @file Scene/MovingObject.hpp
 * @brief An object that moves in a straight line during the exposure (motion blur).
 */

#pragma once

#include <memory>

#include "Scene/Hittable.hpp"

namespace rayt {

    /**
     * @brief Wraps a Hittable that moves with a constant velocity (metres per second).
     *
     * The object sits where it was built at ray time 0 (the opening of the shutter) and at `velocity * time` from there at
     * ray time `time` (seconds). A ray is moved into the frame of the object instead of the object being moved: the origin
     * is shifted back by the displacement, the inner object is tested, and the hit point is shifted forward again. The
     * direction and the distance t do not change, and neither do the normals and tangents (a translation does not turn
     * anything).
     *
     * The bounds cover the whole path up to `shutter` seconds, so the BVH above it stays valid for every ray time.
     * Lights do not move with it (the light samplers know the position at time 0 only): a moving emitter must not be a light.
     */
    class MovingObject final : public Hittable {
    public:
        MovingObject(std::shared_ptr<Hittable> inner, const Vector3& velocity, Real shutterSeconds)
            : m_inner(std::move(inner)), m_velocity(velocity), m_shutter(shutterSeconds) {}

        bool hit(const Ray& r, SurfaceInteraction& rec) const override {
            const Vector3 shift = m_velocity * r.time;
            Ray local = r;
            local.o = r.o - shift;
            if (!m_inner->hit(local, rec)) return false;
            rec.p += shift;
            return true;
        }

        bool occluded(const Ray& r) const override {
            Ray local = r;
            local.o = r.o - m_velocity * r.time;
            return m_inner->occluded(local);
        }

        bool intersectT(const Ray& r, Real& t) const override {
            Ray local = r;
            local.o = r.o - m_velocity * r.time;
            return m_inner->intersectT(local, t);
        }

        AABB bounds() const override {
            const AABB b = m_inner->bounds();
            const Vector3 end = m_velocity * m_shutter;
            return AABB::unite(b, AABB(b.min + end, b.max + end));
        }

    private:
        std::shared_ptr<Hittable> m_inner;
        Vector3 m_velocity;
        Real m_shutter;
    };

} // namespace rayt
