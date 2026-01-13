/**
 * @file Geometry/HittableList.hpp
 * @brief A container for multiple Hittable objects.
 * * Provides a simple aggregate structure to manage a collection of scene objects.
 * Performs linear intersection testing and computes the collective bounding box.
 */


#pragma once

#include <vector>
#include <memory>

#include "Geometry/Hittable.hpp"
#include "Core/AABB.hpp"

namespace rayt {

    /**
     * @brief A collection of Hittable objects.
     * * This class acts as a simple scene graph or world container. During intersection
     * tests, it iterates through all contained objects and ensures that 'rec'
     * always contains the information of the closest hit.
     */
    class Aggregate : public Hittable {
    public:
        Aggregate() {}

        /**
         * @brief Constructs a list with an initial object.
         */
        Aggregate(std::shared_ptr<Hittable> object) { add(object); }

        /**
         * @brief Adds a hittable object to the collection.
         */
        void add(std::shared_ptr<Hittable> object) { objects.push_back(object); }

        /**
         * @brief Clears all objects from the list.
         */
        void clear() { objects.clear(); }

        /**
        * @brief Intersects a ray with all objects in the list and finds the closest hit.
        * * The ray is treated as read-only. If a hit is found, the hit distance is written
        * to rec.t. This method keeps track of the closest hit using a local variable
        * (closestSoFar) and prunes farther tests by passing a reduced tMax via a copied Ray.
        * @param r   The incident ray.
        * @param rec Output structure to store the closest intersection details.
        * @return True if at least one object in the list was hit within the ray's interval.
        */
        virtual bool hit(const Ray& r, SurfaceInteraction& rec) const override {
            SurfaceInteraction tempRec;
            bool hitAnything = false;
            Real closestSoFar = r.tMax;

            for (const auto& object : objects) {
                Ray testRay = r;
                testRay.tMax = closestSoFar; 

                if (object->hit(testRay, tempRec)) {
                    hitAnything = true;
                    closestSoFar = tempRec.t;
                    rec = tempRec;
                }
            }

            return hitAnything;
        }

        /**
         * @brief Computes the Axis-Aligned Bounding Box (AABB) that encloses all objects.
         * * Used to construct higher-level acceleration structures (like BVH) using this
         * list as a node.
         * @return A consolidated AABB covering every object in the collection.
         */
        AABB bounds() const override {
            if (objects.empty()) {
                // Returns an empty/invalid AABB if the list is empty.
                return AABB(Point3(0), Point3(0));
            }

            AABB b = objects[0]->bounds();
            for (size_t i = 1; i < objects.size(); ++i) {
                b = AABB::unite(b, objects[i]->bounds());
            }
            return b;
        }

    public:
        // Container of polymorphic hittable objects.
        std::vector<std::shared_ptr<Hittable>> objects;
    };

} // namespace rayt