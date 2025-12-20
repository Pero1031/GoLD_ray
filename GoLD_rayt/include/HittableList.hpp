#pragma once

// #include "pch.h"
#include "Hittable.hpp"
#include <vector>
#include <memory>

namespace rayt {

    class HittableList : public Hittable {
    public:
        HittableList() {}
        HittableList(std::shared_ptr<Hittable> object) { add(object); }

        void add(std::shared_ptr<Hittable> object) { objects.push_back(object); }
        void clear() { objects.clear(); }

        // リスト内の全オブジェクトと当たり判定を行い、
        // "一番手前(closest)" の交差情報を rec に記録する
        // NOTE:
        // Implementations of hit() may shrink r.tMax (mutable) when a closer
        // intersection is found. r.tMax must only decrease.
        virtual bool hit(const Ray& r, SurfaceInteraction& rec) const override {
            SurfaceInteraction tempRec;
            bool hitAnything = false;

            for (const auto& object : objects) {
                if (object->hit(r, tempRec)) {
                    hitAnything = true;
                    r.tMax = tempRec.t;
                    rec = tempRec;
                }
            }

            return hitAnything;
        }

    public:
        std::vector<std::shared_ptr<Hittable>> objects;
    };

}