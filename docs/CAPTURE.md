# 撮った物を取り込む（スキャン、3D Gaussian Splatting）（2026-11-13）

**状態:** 段階 1（メッシュを経由して取り込み、質感を物理材質で差し替える）を実装済み。3DGS の学習は、この PC（GPU なし）では回さず、別の PC で回したメッシュを持ってくる前提。段階 2（ガウシアンそのものを周囲として描く）、段階 3（再照明可能な GS の推定値との比較）は、まだ。

## 1. ねらい

3DGS は、多視点写真から「その照明での見え方」を再現する表現で、反射率や散乱係数のような物理量を持たない。鏡面や半透明は崩れ、照明を変えられない。このレンダラーは逆で、物理量（スペクトル、表面下散乱、層）はあるが、形は手で作る。**3DGS から形を借り、質感はこちらで差し替える**のが、両方の長所をそのまま使う筋。

| 段階 | 内容 | 状態 |
|---|---|---|
| 1 | 撮影 → 3DGS → メッシュ → 整える → 読む → 物理材質で再照明 | 済み（この文書） |
| 2 | 3DGS の `.ply`（ガウシアン）を直接読み、発光する粒子として周囲に置く（NVIDIA の 3DGRT と同じ考え） | 未 |
| 3 | 再照明可能な GS（GS-IR など）が推定した albedo・粗さ・法線を材質に流し込み、同じカメラで 3DGS の新規視点と物理的な再描画を比べる | 未 |

## 2. 手順（段階 1）

**(1) 撮影 → 3DGS → メッシュ**（GPU のある PC で）。メッシュを出せる実装を使う: 2DGS（`mesh extraction` で PLY）、Gaussian Opacity Fields、SuGaR（テクスチャ付き OBJ）。GPU が弱いときは、画像を長辺 1000〜1600 px に縮小し、反復を 15k〜30k に抑える。2DGS は本家の 3DGS より軽く、メッシュも素直に出る。対象は静物（石、金属、ガラス、半透明の物）。人の顔や肌は動くので向かない。

**(2) 整える**（`tools/data/mesh_clean.py`。OBJ と PLY を読み書きする）:

```bat
python tools\data\mesh_clean.py mesh.ply out.ply --keep-uv --up=-y --drop-small 0.01 --height 0.3 --ground
```

- `--up=-y`: COLMAP / 3DGS の座標（y が下）を、レンダラーの y 上に回す。`--up=z` は z 上（Blender、Open3D）。
- `--drop-small 0.01`: 全体の 1% 未満の三角形しか持たない、離れた破片（floaters）を捨てる。
- `--height 0.3 --ground`: 高さ 0.3 m に合わせ、底を y = 0、足元の中心を x = z = 0 に置く（単位はメートル）。
- `--keep-uv`: 撮影テクスチャを貼るなら uv を残す（なければ捨てて小さくする）。`--cap`: 散乱媒質や `randomWalk` で中を満たすときに、穴を閉じる。
- 法線は、常に捨てて、幾何から作り直す（スキャンの法線は向きが信用できない）。符号付き体積が負なら `--flip`。

**(3) シーンに置く**（`docs/SCENE_FILE.md`）。`mesh` の `file` に `.ply` を書けばそのまま読める。撮影の見え方は `diffuse` の `texture`、差し替えは好きな材質。

```jsonc
"materials": { "captured": { "type": "diffuse", "albedo": 1.0, "texture": "assets/x/texture.jpg" },
               "gold":     { "type": "conductor", "metal": "Au", "roughness": 0.2 } },
"objects":   [ { "shape": "mesh", "file": "assets/x/out.ply", "material": "gold" } ]
```

**例:** `scenes/capture_erato.json`。Erato のスキャン（`assets/Meshes/erato/`）を `--keep-uv` で PLY にし、左は撮影テクスチャのまま、右は金にして、同じ 3 灯で描く。絵は `out\studies\capture\capture_erato.png`。

## 3. 確かめたこと

- `--scene mesh` の c2: PLY の書き出しと読み込みの往復（ascii の頂点 uv + 色の読み飛ばし + `TextureFile` のコメント、バイナリの面ごとの `texcoord` 列）。三角形数・法線・体積が一致し、uv の v が反転して入る（テクスチャの約束: `Core/Texture.hpp` は v が下向き。OBJ の読み込みも同じ反転にそろえた）。RGB ビルドで PASS。
- `mesh_clean.py` の往復: トーラス（`--up=-y --height 0.1 --ground` → PLY → OBJ）、Erato（OBJ 1,237,495 の角 → uv を保って溶接 227,503 頂点、12 秒、PLY 12.6 MB。OBJ の 107 MB から 1/8）。
- 描画: テクスチャは正しい向き（顔、竪琴）で貼れ、同じ光の下で金に差し替わる。900 × 600、128 spp で約 3 分（開発用ビルド）。

## 4. 限界

- **撮影テクスチャの色は、撮ったときの照明込み。** 反射率ではないので、`texture` のまま別の光で描いても正しい再照明ではない（「撮影の見え方」の見本として使う）。反射率を取り出すのは段階 3。
- PLY の頂点の色は読まない（テクスチャか材質で色を決める）。1 物体 1 材質（部位で分けるなら、MeshLab などでメッシュを分ける）。
- `--cap` は単純な穴だけ閉じる。3DGS のメッシュの質（穴、二重の面、法線の向き）は、対象と手法で変わり、試すまで分からない。
- 3DGS の学習と描画は CUDA が前提で、この PC では回らない。ガウシアン自体は、まだ読めない（段階 2）。
