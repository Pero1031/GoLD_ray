// src/ImageIO.cpp

#include "pch.h"

// 実装を有効にするためのマクロ定義 (このファイルだけで行う)
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

// 警告抑制 (stbは少し古い書き方をしている箇所があるため)
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4996) // 廃止された関数への警告などを無視
#endif

#include "stb_image.h"
#include "stb_image_write.h"

#ifdef _MSC_VER
#pragma warning(pop)
#endif