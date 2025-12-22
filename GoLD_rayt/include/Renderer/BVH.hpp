#pragma once

#include <memory>
#include <vector>
#include <algorithm>
#include <random>
#include <limits>

#include "Core/Types.hpp"
#include "Core/Ray.hpp"
#include "Core/Interaction.hpp"
#include "Core/AABB.hpp"
#include "Geometry/Hittable.hpp"   

namespace rayt {

    class BVHNode : public Hittable {
    public:
        std::shared_ptr<Hittable> left;
        std::shared_ptr<Hittable> right;
        AABB box;
        int splitAxis = 0;

        BVHNode() = default;

        BVHNode(std::vector<std::shared_ptr<Hittable>>& objects,
            size_t start, size_t end)
        {
            // ここは「とりあえず動く」median split（SAHは後でOK）
            const size_t span = end - start;

            // split axis: 最大extentの軸を使う（そこそこ良い）
            AABB centroidBox;
            {
                Vector3 cmin(std::numeric_limits<Real>::infinity());
                Vector3 cmax(-std::numeric_limits<Real>::infinity());
                for (size_t i = start; i < end; ++i) {
                    const AABB b = objects[i]->bounds();
                    const Vector3 c = (b.min + b.max) * Real(0.5);
                    cmin = glm::min(cmin, c);
                    cmax = glm::max(cmax, c);
                }
                centroidBox = AABB(cmin, cmax);
            }
            const Vector3 e = centroidBox.extent();
            int axis = 0;
            if (e.y > e.x) axis = 1;
            if (e.z > (axis == 0 ? e.x : e.y)) axis = 2;

            this->splitAxis = axis;

            auto cmp = [axis](const std::shared_ptr<Hittable>& a,
                const std::shared_ptr<Hittable>& b) {
                    const AABB ba = a->bounds();
                    const AABB bb = b->bounds();
                    const Vector3 ca = (ba.min + ba.max) * Real(0.5);
                    const Vector3 cb = (bb.min + bb.max) * Real(0.5);
                    return ca[axis] < cb[axis];
                };

            if (span == 1) {
                left = objects[start];
                right = nullptr;
            }
            else if (span == 2) {
                if (cmp(objects[start], objects[start + 1])) {
                    left = objects[start];
                    right = objects[start + 1];
                }
                else {
                    left = objects[start + 1];
                    right = objects[start];
                }
            }
            else {
                std::sort(objects.begin() + start, objects.begin() + end, cmp);
                const size_t mid = start + span / 2;
                left = std::make_shared<BVHNode>(objects, start, mid);
                right = std::make_shared<BVHNode>(objects, mid, end);
            }

            // AABBの結合 (葉ノードの場合は片方だけ)
            if (!right) {
                box = left->bounds();
            }
            else {
                box = AABB::unite(left->bounds(), right->bounds());
            }
        }

        bool hit(const Ray& r, SurfaceInteraction& rec) const override {
            // AABBヒット判定 (r.tMaxが短くなるのを利用)
            if (!box.intersect(r, r.tMin, r.tMax))
                return false;

            // 葉ノード（右がない）なら、左だけ判定して終わり
            if (!right) {
                return left->hit(r, rec);
            }

            // ここからは「内部ノード（左右あり）」の処理
            // レイの向きに合わせて手前(first)と奥(second)を決める
            const bool dirNeg = r.d[splitAxis] < 0;

            // rightはnullptrではないことが保証されているので安全
            auto first = dirNeg ? right : left;
            auto second = dirNeg ? left : right;

            bool hitAny = false;

            // 手前をチェック (ヒットすれば r.tMax が縮む)
            if (first->hit(r, rec)) hitAny = true;

            // 奥をチェック (縮んだ r.tMax で判定されるので高速)
            // ※ firstと同じノードでないことは保証されているのでそのまま呼ぶ
            if (second->hit(r, rec)) hitAny = true;

            return hitAny;
        }

        AABB bounds() const override { return box; }
    };

} // namespace rayt
