# シーンファイル（JSON）

CLI は、シーンを JSON ファイルで指定できます。**ファイルは、先頭の引数**に置きます。

```bat
GoLD_rayt.exe scenes\gold_bowl_tungsten.json
GoLD_rayt.exe scenes\gold_bowl_tungsten.json --spp 16 --width 400   :: 後ろの指定が、ファイルの値より優先される
GoLD_rayt.exe --list-scenes                                          :: 組み込みシーンとコマンドの一覧
```

Visual Studio では、「プロジェクトのプロパティ → デバッグ → コマンド引数」に `scenes\gold_bowl_tungsten.json` と書き、シーンを変えるときは、このファイル名だけを変える。作業ディレクトリは `$(ProjectDir)`。

## 2 つの形式

| 形式 | 中身 | 使いどころ |
|---|---|---|
| 1. 組み込みシーン＋オプション | `"scene"` と `"options"` | コマンドラインの内容を、ファイルに残す（`scenes/presets/`） |
| 2. データで書いたシーン | `"objects"` ほか | 球・板・メッシュ・回転体・光源・材質を、再コンパイルなしで書く（`scenes/*.json`） |

形式 2 のファイルは、シーン名 `file` として描画される。`"scene"` と `"objects"` は、同時に書けない。

```jsonc
// 形式 1
{ "scene": "goldbowl", "options": { "light": "tungsten", "width": 800, "spp": 256 } }
```

```jsonc
// 形式 2（抜粋）
{
  "options": { "width": 800, "spp": 256 },
  "camera": { "from": [0.02, 0.27, 0.33], "at": [0.035, 0.035, 0], "fov": 30 },
  "materials": {
    "gold": { "type": "conductor", "metal": "Au", "roughness": 0.12 },
    "bulb": { "type": "emission", "spectrum": "tungsten", "radiance": 90 }
  },
  "objects": [
    { "shape": "sphere", "center": [0, 0.03, 0], "radius": 0.028, "material": "gold" },
    { "shape": "sphere", "center": [-0.18, 0.2, 0.1], "radius": 0.03, "material": "bulb" }
  ]
}
```

## 書き方の約束

- `//` のコメントが使える。キーの綴りの誤りは、**許されるキーの一覧つきのエラー**になる。
- ファイル内の相対パスは、プログラムを実行するフォルダ（`GoLD_rayt/`）からの相対パス（コマンドラインと同じ）。
- `"options"`: コマンドラインの名前から `--` を除いたもの。数値・文字列は値、`true` はフラグ。例: `{ "width": 800, "denoise": true, "out": "x.png" }`。出力の色域と PNG のビット数も同じ（`"gamut": "p3"`、`"bits": 16`）。
- 単位はメートル。ベクトルは `[x, y, z]`、色は `[r, g, b]`（数値 1 つなら灰色）。

## キーの一覧（形式 2）

| キー | 中身 |
|---|---|
| `camera` | `from`、`at`、`up`、`fov`（度）。省略すると、`[0, 0.5, 2.5]` から `[0, 0, -1]` を 35° で見る |
| `environment` | `"none"`（既定）、`"hdr"`、`"uniform"`（白）、または `{ "type", "file", "scale" }`。`hdr` の既定のファイルは `--env` |
| `materials` | 名前 → 材質（下） |
| `lights` | 点光源の配列 `{ "type": "point", "position", "color", "intensity" }` |
| `objects` | 物体の配列（下） |

**材質（`type`）**

| type | キー |
|---|---|
| `diffuse` | `albedo`（RGB）、または `reflectance`（測定した反射率 R(λ) の CSV。**どちらか 1 つ**） |
| `conductor` | `metal`（`Au`・`Ag`・`Cu`: Johnson & Christy の表）、`table`（RefractiveIndex.info の n, k の CSV）、`reflectance`（測定した垂直入射の反射率 R(λ) の CSV。Gulbrandsen の方法で n, k にする）の**どれか 1 つ**。ほかに `roughness`、`anisotropy`（−1〜1。> 0: 接線の向きに粗い、< 0: 接線に沿った溝のヘアライン仕上げ）、`rotation`（度。接線を法線の周りに回す）、`multipleScattering` |
| `dielectric` | `ior`、`roughness`、`abbe`（分散。0 = なし）、`sellmeier`（分散を Sellmeier の式で与える。ガラスの名前 `"diamond"`・`"sapphire"`・`"bk7"`・`"silica"`、または `{ "b": [b1, b2, b3], "c": [c1, c2, c3] }`（c は µm²）。`abbe` とは両立しない。RGB ビルドは d 線の屈折率）、`distanceRoughness`（`{ "roughness", "start", "end" }`: y 軸からの距離が start から end の間に、粗さがこの値まで増える。メッシュで解像できない遠くの小さな波の傾きを、粗さとして返す。海に使う）、`medium`（内側を散乱媒質で満たす。下） |
| `thinFilm` | 滑らかな面の上の、薄い透明な膜の干渉（シャボン玉、油膜、コーティングされた金属）。`ior`（膜）、`thickness`（nm の数値、または勾配 `from`/`to`/`axis`/`range`/`exponent`、`warp`、`noise`（`stretch` あり））、`coverage`（膜のある場所のマスク）、`substrate`（`"free"`、`{ "dielectric": n }`、`{ "metal": "Au" }`）。理論と例は `docs/THIN_FILM.md` |
| `layered` | 透明なコートが、下地の上にある材質（漆を塗った金、車の塗装）。`coat`: `{ "ior", "roughness", "tint", "samples" }`（コートの屈折率、粗さ、法線方向の透過率 RGB（1 = 透明）、1 回の評価での追跡の数）。`base`: `{ "albedo" }`（拡散）または `{ "metal" | "table" | "reflectance", "roughness", "anisotropy", "rotation", "multipleScattering" }`（金属）のどちらか 1 つ。理論と例は `docs/GLOSS.md` |
| `skin` | 肌（`CoatedDiffuse`）。`albedo`（既定 0.80, 0.52, 0.42）、`roughness`（既定 0.5）、`specular`（皮脂のコートの重み、既定 0.4）。`scatterDistanceMm`（RGB、mm。例 `[3.67, 1.37, 0.68]`）を書くと**表面下散乱**（Burley のプロファイル、`sssScale` で倍率）。`relief`（省略可。毛穴と皮溝の凹凸: `cellMm`、`furrowDepthMm`、`furrowWidth`、`stretch`、`coarseScale`、`coarseWeight`、`poreDepthMm`、`poreRadiusMm`、`poreFraction`、`poreSizeVariation`、`poreDensityVariation`、`poreDensityFreq`、`poreDarkness`（毛穴を暗くする、0〜1）、`seed`。法線だけを曲げる）。`relief` の中の `wrinkles`（**皮溝と毛穴をグラフ合成で作る**: `tileMm`（既定 24）、`resolution`（2048）、`steps`（300）、`depthScale`、`seed`、`levels`（省略可。層ごとの `spacingMm`、`angleDeg`、`base`、`distance`、`continuation`、`similarity`、`crossing`、`deposit`、`width`、`poreWidth`、`poreMin`、`poreBlend`、`poreSkew`、`wobble`、`wobbleFreq`、`depthUm`、`bulge`、`darkPores`）。これを書くと、Worley の皮溝と毛穴は使われない）。`gloss`（省略可。皮脂の多い所と少ない所: `variation`（0〜0.95）、`frequency`（1 m あたりの斑の数）、`seed`）。`pigment`（省略可。色の揺らぎ: `freckleCellMm`、`freckleFraction`、`freckleRadiusMm`、`freckleStrength`、`clusterFreq`、`blotchFreq`、`blotchStrength`、`rednessFreq`、`rednessStrength`、`seed`）。`docs/SKIN.md` |
| `hair` | 髪の毛 1 本の散乱（メラニンの色づき、R・TT・TRT のハイライト。`docs/HAIR.md`）。吸収は、`melanin`（`[真性, 亜性]`。黒 `[8, 0]`、茶 `[1.3, 0.1]`、赤 `[0.4, 1.4]`、金髪 `[0.25, 0.15]`、白髪 `[0, 0]`）、`color`（見た目の色 RGB、各 (0, 1]）、`sigmaA`（毛の半径あたりの吸収 RGB）の**どれか 1 つ**（既定は `melanin` [1.3, 0]）。`eta`（既定 1.55）、`betaM`・`betaN`（縦・方位の粗さ、0〜1、既定 0.3）、`tilt`（キューティクルの傾き、度、既定 2）、`lobes`（`["R", "TT", "TRT", "rest"]` の一部だけを描く。説明用） |
| `emission` | `spectrum`、`radiance`、`twoSided` |

`dielectric` の `medium`: 閉じた形（球、向きの正しい OBJ）の内側を、散乱媒質で満たす。ガラスの境界から入った光が、内部で散乱してから出てくる（蝋、ミルク、翡翠）。

```jsonc
"jade": { "type": "dielectric", "roughness": 0.2,
          "medium": { "sigmaT": [250, 150, 250], "albedo": [0.85, 0.998, 0.95], "g": 0.3 } }
```

- `sigmaT`: 消散係数（1/m）。大きいほど濃い。`albedo`: 散乱で光が残る割合（1 = 吸収なし）。どちらも、数値 1 つ（灰色）か `[r, g, b]`（色つき）。`g`: 位相関数の非等方性（−1〜1。正で前方散乱）。`g2` と `weight` を足すと、2 ローブ（前方＋後方）の位相関数 `weight·HG(g) + (1−weight)·HG(g2)` になる。平均コサインは `weight·g + (1−weight)·g2`。同じ σs′ = σs(1−平均コサイン) なら、厚い媒質では見た目がほぼ同じで、薄い部分・逆光・表面近くの光で差が出る。
- `phaseTable`: 位相関数を表で与える（`g`・`g2`・`weight` とは同時に使えない）。CSV は「角度（度）, 値」の行で、0° から 180° まで等間隔、`#` はコメント。値は 1 ステラジアンあたりの密度（大きさは自動で正規化）。雲の水滴の Mie 散乱は `assets/Cloud/mie_reff10.csv`（`tools/studies/mie_phase.py` で作る。`docs/CLOUD.md`）。
- `--media-nee`（`options` では `"media-nee": true`）: 媒質の衝突ごとに光源を標本抽出し、屈折率 1 の境界を素通りする。太陽のような小さな光源で雲や霧を照らすときに要る（`scenes/cloud_mie.json`）。
- 色つきの媒質は、1 本の経路で代表の波長（RGB 版は 1 チャンネル）だけで距離を選び、ほかの波長は確率密度の比で重みを付ける（`--scene volume` が白炉と灰色との一致を確認する）。色つきは、灰色より、同じ spp でノイズが多い。
- `noise`（省略可）を足すと、密度が場所で変わる媒質（雲、霧）になる。`sigmaT` は、密度が最大の場所の値。

```jsonc
"cloud": { "type": "dielectric", "ior": 1.0,   // ior 1: 境界が見えない
           "medium": { "sigmaT": 420, "albedo": 0.9995, "g": 0.6,
                       "noise": { "frequency": 22, "octaves": 5, "lo": 0.42, "hi": 0.62, "seed": 7,
                                  "maskCenter": [0, 0, 0], "maskRadius": 0.1, "maskStart": 0.45 } } }
```

  `frequency`: 最初の層の 1 m あたりの模様の数（雲の大きさは約 1/frequency）。`octaves`: 層の数。`lo`/`hi`: ノイズが lo 以下で密度 0、hi 以上で密度 1（幅が狭いと輪郭がはっきり、広いと霞）。`seed`: 別の雲。`maskCenter`/`maskRadius`/`maskStart`: 球の表面に向けて密度を 0 にする（球の中の雲の輪郭を、球から外す）。密度は世界座標で決まる。
  仕組みは、ヌル衝突（密度の最大値で距離を選び、確率 1−q で何も起こさない）。
- 経路の長さの上限は、`--volume-cap`。

分光放射計の測定データは、`tools/data/sr5100.py` で、2 列の CSV（`波長_nm, 値`）に変換できる（下）。

`spectrum`: プリセット名 `"d65"`・`"tungsten"`・`"candle"`・`"led"`、または `{ "blackbody": 2856 }`、`{ "rgb": [1, 1, 1] }`、`{ "file": "x.csv" }`（波長と値の表）。

**物体（`shape`）**

| shape | キー |
|---|---|
| `sphere` | `center`、`radius` |
| `quad` | `corner`、`edge1`、`edge2`（`cross(edge1, edge2)` の向きに面する） |
| `mesh` | `file`（OBJ）、`scale`、`rotate`（`{ "axis", "degrees" }`）、`translate`。順は、拡大 → 回転 → 移動 |
| `lathe` | `profile`（`[半径, 高さ]` の配列、または `{ "shell": { "radius", "wall", "endDegrees", "steps" } }`: 椀の壁）、`segments`、`translate` |
| `ocean` | 海面（波のスペクトルから作る。`docs/WATER.md`）。原点の周りの円形のメッシュ。`windSpeed`（m/s、1〜40、既定 8）、`windDirection`（度、波が進む向き `atan2(z, x)`）、`spread`（向きの広がり s）、`peakEnhance`（JONSWAP の γ、1 = Pierson–Moskowitz）、`chop`（0〜1、波頭のとがり）、`time`（秒）、`minWavelength`（m）、`seed`、`rInner`・`rOuter`（m）、`rings`・`sectors`（メッシュの細かさ） |
| `hair` | 手続き的に作った毛（`docs/HAIR.md`）。`style`: `"lock"`（上から垂れる毛束）または `"head"`（球の頭に流れる髪）。共通: `strands`、`segments`、`length`（m）、`rootRadius`・`tipRadius`（毛の半径、m。実物は約 0.00005）、`clumps`（房の数）、`seed`。`lock`: `root`（根元の位置）、`width`、`tipWidth`（先端の幅 / 根元の幅）、`waveAmp`・`waveLength`、`twist`、`frizz`、`flyaway`（束から外れる毛の割合）。`head`: `center`、`radius`（頭の球）、`droop`、`outward`、`sway`・`swayFreq`、`headFrizz`、`clumpStrength`、`maxPolar`（度）、`hairline`。材質は `hair` |

すべての物体に `material`（必須）。`emission` の材質を持つ球・板は、既定で**光源**（次点光源推定の対象）になる。`"light": false` なら、光源にしない（全周の淡い光、のような使い方）。

## 確認したこと

`scenes/` の 3 つ（`gold_spheres`、`gold_bowl_tungsten`、`studio`）は、コード版の `gold`、`goldbowl --light tungsten`、`studio` を書き起こしたもの。RGB 版・スペクトル版の両方で、4 回ずつ描いた平均色が、コード版と標準誤差の 3 倍以内で一致し、画像も区別がつかない。

## 材質や形状を足すとき

3 か所を直す。

1. `include/Cli/SceneDescription.hpp`: 型に項目を足す。
2. `src/Cli/SceneFile.cpp`: JSON を読んで検証する（キーの一覧に足す）。
3. `src/Cli/Scenes/FileScene.cpp`: その型から、エンジンの物を作る。

## 分光放射計（SR-5100）のデータを使う

`tools/data/sr5100.py` が、SR-5100 の CSV の書き出し（スポットと画像の 2 形式）を読んで、このレンダラが読む 2 列の CSV に変換する。

```bat
python tools\data\sr5100.py info 測定.csv                                       :: 中身の一覧（点の数、ピーク、測定距離、積分時間）
python tools\data\sr5100.py spots 測定.csv --out 出力フォルダ                     :: スポットを 1 点ずつ「波長_nm, radiance」の CSV に
python tools\data\sr5100.py roi 画像.csv --peak 2 --out 光源.csv                  :: 画像の最も明るい画素の周り（5×5）の平均スペクトル
python tools\data\sr5100.py reflectance 試料.csv 白.csv --out 反射率.csv          :: 試料 ÷ 白色（白色が弱い波長は出力しない）
```

- LED のスペクトルは、`"spectrum": { "file": "光源.csv" }` で指定する（値の単位は何でもよい。明るさは 1 にそろえて使う）。
- 測定した反射率は、`"type": "conductor", "reflectance": "反射率.csv"` で金属として、`"type": "diffuse", "reflectance": "反射率.csv"` で拡散面（布など）として使える。拡散面では、スペクトル版は曲線をそのまま使い、RGB 版は D65 の下の色に直して使う（測定範囲の外は、端の値で延ばす。1 を超える値は警告）。
- 測定データは実験データなので、リポジトリには入れない（`assets/Metal_measured/` は Git の対象外）。
