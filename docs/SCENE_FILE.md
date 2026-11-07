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
- `"options"`: コマンドラインの名前から `--` を除いたもの。数値・文字列は値、`true` はフラグ。例: `{ "width": 800, "denoise": true, "out": "x.png" }`。
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
| `diffuse` | `albedo` |
| `conductor` | `metal`（`Au`・`Ag`・`Cu`: Johnson & Christy の表）、`table`（RefractiveIndex.info の n, k の CSV）、`reflectance`（測定した垂直入射の反射率 R(λ) の CSV。Gulbrandsen の方法で n, k にする）の**どれか 1 つ**。ほかに `roughness`、`multipleScattering` |
| `dielectric` | `ior`、`roughness`、`abbe`（分散。0 = なし） |
| `emission` | `spectrum`、`radiance`、`twoSided` |

分光放射計の測定データは、`tools/sr5100.py` で、2 列の CSV（`波長_nm, 値`）に変換できる（下）。

`spectrum`: プリセット名 `"d65"`・`"tungsten"`・`"candle"`・`"led"`、または `{ "blackbody": 2856 }`、`{ "rgb": [1, 1, 1] }`、`{ "file": "x.csv" }`（波長と値の表）。

**物体（`shape`）**

| shape | キー |
|---|---|
| `sphere` | `center`、`radius` |
| `quad` | `corner`、`edge1`、`edge2`（`cross(edge1, edge2)` の向きに面する） |
| `mesh` | `file`（OBJ）、`scale`、`rotate`（`{ "axis", "degrees" }`）、`translate`。順は、拡大 → 回転 → 移動 |
| `lathe` | `profile`（`[半径, 高さ]` の配列、または `{ "shell": { "radius", "wall", "endDegrees", "steps" } }`: 椀の壁）、`segments` |

すべての物体に `material`（必須）。`emission` の材質を持つ球・板は、既定で**光源**（次点光源推定の対象）になる。`"light": false` なら、光源にしない（全周の淡い光、のような使い方）。

## 確認したこと

`scenes/` の 3 つ（`gold_spheres`、`gold_bowl_tungsten`、`studio`）は、コード版の `gold`、`goldbowl --light tungsten`、`studio` を書き起こしたもの。RGB 版・スペクトル版の両方で、4 回ずつ描いた平均色が、コード版と標準誤差の 3 倍以内で一致し、画像も区別がつかない。

## 材質や形状を足すとき

3 か所を直す。

1. `include/Cli/SceneDescription.hpp`: 型に項目を足す。
2. `src/Cli/SceneFile.cpp`: JSON を読んで検証する（キーの一覧に足す）。
3. `src/Cli/Scenes/FileScene.cpp`: その型から、エンジンの物を作る。

## 分光放射計（SR-5100）のデータを使う

`tools/sr5100.py` が、SR-5100 の CSV の書き出し（スポットと画像の 2 形式）を読んで、このレンダラが読む 2 列の CSV に変換する。

```bat
python tools\sr5100.py info 測定.csv                                       :: 中身の一覧（点の数、ピーク、測定距離、積分時間）
python tools\sr5100.py spots 測定.csv --out 出力フォルダ                     :: スポットを 1 点ずつ「波長_nm, radiance」の CSV に
python tools\sr5100.py roi 画像.csv --peak 2 --out 光源.csv                  :: 画像の最も明るい画素の周り（5×5）の平均スペクトル
python tools\sr5100.py reflectance 試料.csv 白.csv --out 反射率.csv          :: 試料 ÷ 白色（白色が弱い波長は出力しない）
```

- LED のスペクトルは、`"spectrum": { "file": "光源.csv" }` で指定する（値の単位は何でもよい。明るさは 1 にそろえて使う）。
- 測定した反射率は、`"type": "conductor", "reflectance": "反射率.csv"` で金属として使える。
- 測定データは実験データなので、リポジトリには入れない（`assets/Metal_measured/` は Git の対象外）。
