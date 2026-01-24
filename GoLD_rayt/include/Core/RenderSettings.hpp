/*
* @file include/Core/RenderSettings.hpp
*/

#pragma once

// ビューモード用（必要に応じて増やす）
enum class ViewMode {
    Beauty = 0,
    Normal,
    Depth,
    Albedo
};

struct RenderSettings {
    // レンダリング制御
    bool paused = false;
    bool requestReset = false;
    bool requestSave = false;

    // 画質・設定
    int sppPerFrame = 1;
    int targetSpp = 1024;
    int maxDepth = 8;

    // カメラ・表示
    float exposure = 0.0f;
    ViewMode viewMode = ViewMode::Beauty;

    // 解像度（初期値）
    int width = 1280;
    int height = 800;
};