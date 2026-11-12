# コードの構成（2026-11-13）

どこに何があるかの地図。機能ごとの理論と結果は、各 `docs/*.md`。全体で約 3.4 万行（エンジン約 2.2 万、CLI 約 7 千、自己検査約 5 千）。

## 1. 二層構造と、1 枚描くときの流れ

**エンジン**（`GoLD_rayt/include/` のヘッダが中心。本体が大きい物は `src/` の同じ名前の .cpp: `Materials/`（CoatedDiffuse、Dielectric、LayeredCoat、RoughConductor、ThinFilm、Hair、WrinkleField）、`Render/`（Integrator、Film、CausticPhotonMap）、`Geometry/`（ProceduralMesh、GaussianSplats）、`IO/`（PmxLoader、PlyLoader、ImageIO）、`Accelerators/`（BVH））と、その上に乗る**アプリ**（`include/Cli` と `src/Cli`。コマンドライン `main_cli.cpp`）の二層。GUI（`Application`、`UI`、`main_ui.cpp`）は別の薄い層で、同じエンジンを呼ぶ。**GUI は凍結中**（CLI が本線。GPU 実装ができたらリアルタイム性が生きる。将来は「描くものを決める操作画面」に転用の可能性があるので、`Cli/SceneDescription`（エンジンに依存しないデータ）はその土台として保つ）。

**新しい .cpp は `GoLD_rayt.vcxproj`（と `.filters`）に登録しないとリンクされない。** 先頭は `#include "pch.h"`。

1 枚描くときの流れ:

1. `main_cli` → `Cli/Options`（引数の解釈）
2. `Cli/Scenes/SceneRegistry`（`--scene` 名から組み立て関数を選ぶ。JSON のシーンは `Scenes/FileScene`）
3. `Cli/RenderPipeline`（共通の手順: 標準材質、環境、シーン全体の BVH、光源の登録、カメラ、積分器の選択、描画、保存）
4. `Render/Integrator`（経路追跡の本体 `PathIntegrator::Li`）
5. `Render/Film`（画素の蓄積、色域、PNG / PFM）→ `Render/Denoiser`

## 2. ディレクトリごとの役割

| ディレクトリ | 中身 |
|---|---|
| `Core` | 型、乱数と標本抽出、AABB、`Ray`、`SurfaceInteraction`（当たった点の情報。肌の機能が足した欄も多い）、`Spectrum`（RGB か波長標本かをマクロ `RAYT_SPECTRAL_MODE` で切り替える。2 種類のビルドがあるのはこのため）、テクスチャ |
| `Color` | 分光と RGB の変換、等色関数、色域（sRGB / P3 / Rec.2020）、CAT16 のホワイトバランス、照明のスペクトル、金属の複素屈折率の表、Sellmeier の分散 |
| `Materials` | `Material` が BSDF の共通インターフェイス（eval / sample / pdf に加え、発光、不透明度、媒質の境界、表面下散乱、法線の摂動などのフック）。実装: 拡散、粗い導体、誘電体、薄膜、コート層、髪、肌（`CoatedDiffuse` が肌の本体で、`SkinRelief`・`WrinkleField`・`SkinPigment`・`SkinFilm`・`SkinLotion`・`SkinScatter`・`CosmeticLayer`・`KubelkaMunkLayer` を束ねる）、表面下散乱（`Subsurface` の Burley 型と `SubsurfaceWalk` の体積歩行）、媒質（`Medium`、`DensityField`、`TabulatedPhase`）、褪色（`FadingReflectance`、`FadingDye`） |
| `Microfacet` | GGX と、その方向アルベド（多重散乱の補償） |
| `Lights` | 面光源、点光源、環境光、光源の選び方（`LightSampler`: 一様、放射パワー、空間的な重要度） |
| `Geometry` | 球、四角形、三角形メッシュ、髪の曲線と生成、回転体と海のメッシュ、ガウシアンの粒子（3DGS） |
| `Scene` | `Hittable`（交差の共通インターフェイス）、`Primitive`（形と材質と光源を束ねる）、`Aggregate`（リスト）、`Scene`（全体: BVH、光源、環境、粒子） |
| `Accelerators` | BVH（自前の SAH）と Embree |
| `Render` | 経路追跡（`Integrator.cpp`。媒質、粒子、コースティクスのフックもここ）、影の光線の透過率、表面下散乱の出口探索、コースティクスの光子マップ、多様体 NEE、LVC-BPT の試作、カメラ、フィルム、デノイザー |
| `IO` | OBJ、PLY、PMX（モデル、モーフ、ポーズ）、環境マップ、CSV のスペクトル、褪色の測定系列 |

## 3. CLI 側

| 場所 | 中身 |
|---|---|
| `Cli/Options` | 引数の既定値、構造体、解釈 |
| `Cli/Commands` | 描かないコマンド（検査、測定、モデルの一覧） |
| `Cli/Scenes/` | 組み込みシーンの組み立て: 金属、ガラスと媒質、スタジオ、水、PMX、粒子。`SceneRegistry` が表 |
| `Cli/SceneFile.cpp`、`SceneFileMaterials.cpp`、`SceneDescription.hpp`、`Scenes/FileScene.cpp` | JSON シーン。読む → 素の構造体 → エンジンの物を作る、の三段 |
| `Cli/RenderPipeline` | 共通の描画手順 |
| `Cli/Checks/` | 話題ごとの自己検査（`--scene furnace \| lights \| aov \| sss \| volume \| caustic \| manifoldtest \| metals \| mesh \| splats \| cosmetic \| relief \| ...`） |
| `Cli/Tools/` | 測定と診断（`--scene skinlayers`、`sssrw`、`bench`、`--trace`） |

## 4. 機能と場所の対応

| 機能 | 場所 | 文書 |
|---|---|---|
| 分光レンダリング、色域、ホワイトバランス | `Core/Spectrum`、`Color/*` | `METAL_ACCURACY.md`、`SPECTRUM_FROM_RGB.md` |
| 金属の精度、測定反射率からの n, k | `Materials/RoughConductor`、`Color/ComplexIorTable`、`ConductorFromReflectance` | `METAL_ACCURACY.md`、`MEASURED_GOLD.md` |
| 薄膜、分散、コート | `Materials/ThinFilm`、`Color/Sellmeier`、`Materials/LayeredCoat` | `THIN_FILM.md`、`GLOSS.md` |
| 表面下散乱 | `Materials/Subsurface*`、`Render/SubsurfaceProbe` | `SSS.md` |
| 肌（凹凸、色素、保湿、化粧） | `Materials/CoatedDiffuse` + `Skin*`、`WrinkleField`、`CosmeticLayer` | `SKIN.md`、`SKIN_MOISTURE.md`、`COSMETICS.md`、`WETNESS.md` |
| 髪 | `Materials/Hair`、`Geometry/HairCurves`、`HairGenerator` | `HAIR.md` |
| 水、海、雲と霧 | `Geometry/ProceduralMesh`、`Materials/Medium`、`DensityField`、`TabulatedPhase` | `WATER.md`、`CLOUD.md` |
| コースティクス、多様体 NEE、LVC-BPT | `Render/CausticPhotonMap`、`SpecularNee`、`VolumeBPT` | `LVC_BPT.md` |
| 撮った物の取り込み、3DGS | `IO/PlyLoader`、`Geometry/GaussianSplats`、`Cli/Scenes/FileScene` | `CAPTURE.md` |
| MMD | `IO/Pmx*`、`Cli/PmxSupport`、`Scenes/PmxScene` | — |
| シーンファイル | `Cli/SceneFile*`、`SceneDescription` | `SCENE_FILE.md` |

## 5. 整理の記録と候補

**済み（2026-11-13、分割）:** 本体が長く、よく直し、インターフェイスが安定している物を .hpp と .cpp に分けた: `WrinkleField`（構築は .cpp、レイごとの参照だけヘッダ）、`Hair`（ローブの数式 `hairmath` はヘッダ、材質は .cpp）、`CausticPhotonMap`、`PmxLoader`、`PlyLoader`、`GaussianSplats`。分けなかった物: `VolumeBPT`・`SpecularNee`（試作で、触らない）、`CosmeticLayer`・`SkinRelief` など（レイごとに呼ぶ小さな関数が多く、開発用ビルド（LTCG なし）で遅くなり得る）、テンプレート（`Spectrum`）。回帰は `tools/test/regress.py`（小さな描画の平均色 + 自己検査 15 本）で確かめた。**ビルド時間（開発用ビルド）:** 分割した物の中身（.cpp）を直したときの再ビルドが、全体（約 80 秒）から約 4 秒になった。ヘッダを触ったときは、分割の前後とも全体の再ビルド（`Cli/EngineHeaders.hpp` が全ヘッダを全 CLI ファイルに配るため。次の候補の 1 つ目）。

**候補:**
- **`Cli/EngineHeaders.hpp` が 63 個のヘッダを束ね、CLI の 41 ファイル全部が読む。** エンジンのヘッダを 1 つ触ると CLI 全体が再コンパイルされる。各ファイルが使う物だけ include する形に。
- **`Materials` に性質の違う物が同居**している: BSDF（`RoughConductor` など）、肌の「場」（`SkinRelief`、`SkinPigment`、`WrinkleField`、`DensityField`）、時間のモデル（`SkinLotion`、`Fading*`）。
- **`SurfaceInteraction` が肌の機能で太っている**（`ns`、`filmCover`、`surfScatter`、`nl`、`albedoTint`、`glossScale`、`roughScale` など）。すべての交差が運ぶ構造体なので、結合が強い。
- **`Integrator.cpp` に機能のフックが集まっている**（媒質、粒子、光子マップ、多様体 NEE、不透明度、無効境界）。
- **`SceneFileMaterials.cpp` が 734 行の分岐一本**で、材質を足すたびに伸びる。
- **試作のままの物**（`VolumeBPT`、`SpecularNee` は prototype と書いてある）と本線の区別が、ディレクトリからは見えない。
- **GUI（`Application`）が本線と同じ速度で追従しているか**は、未確認。
