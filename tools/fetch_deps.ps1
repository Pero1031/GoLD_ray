<#
  Downloads the third-party libraries that this repository does not contain, into the places the Visual Studio project expects.

    powershell -ExecutionPolicy Bypass -File tools\fetch_deps.ps1            download what is missing
    powershell -ExecutionPolicy Bypass -File tools\fetch_deps.ps1 -Force     download everything again

  What it fetches (each keeps its own license; see README section 9):

    Embree 4.4.0            librarys\embree-4.4.0.x64.windows\     ray tracing kernels (Apache-2.0)
    Open Image Denoise 2.5.1   librarys\oidn-2.5.1.x64.windows\    denoiser (Apache-2.0)
    GLM                     GoLD_rayt\external\glm\                 math (MIT)
    nlohmann/json 3.11.3    GoLD_rayt\external\nlohmann\json.hpp    scene files (MIT)
    stb_image, stb_image_write   GoLD_rayt\include\stb\             image reading / writing (public domain / MIT)
    GLFW 3.4 (Windows binaries)  GoLD_rayt\external\glfw-3.4.bin.WIN64\   window for the GUI (zlib)
    Dear ImGui (docking branch)  external\imgui\                    GUI (MIT)

  glad (the OpenGL loader, 3 small generated files) is part of the repository: GoLD_rayt\external\glad\.
  The project compiles the GUI sources in every build, so ImGui and GLFW are needed even for the command-line renderer.
#>
param([switch]$Force)
$ErrorActionPreference = "Stop"
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
$root = Split-Path $PSScriptRoot -Parent
$proj = Join-Path $root "GoLD_rayt"
$tmp = Join-Path $env:TEMP ("raytdeps_" + [Guid]::NewGuid().ToString("N").Substring(0, 8))
New-Item -ItemType Directory -Force $tmp | Out-Null

$free = (Get-PSDrive -Name $root.Substring(0, 1)).Free / 1GB
if ($free -lt 1.5) { throw ("Only {0:N1} GB are free on drive {1}: the downloads need about 1 GB. Free some space and run it again." -f $free, $root.Substring(0, 1)) }

function Get-File($url, $dest) {
    Write-Host "  downloading $url"
    $ProgressPreference = "SilentlyContinue"   # the progress bar makes Invoke-WebRequest very slow
    Invoke-WebRequest -Uri $url -OutFile $dest -UseBasicParsing
}

function Get-Zip($url, $name) {
    $zip = Join-Path $tmp "$name.zip"
    Get-File $url $zip
    $out = Join-Path $tmp $name
    Expand-Archive -Path $zip -DestinationPath $out -Force
    Remove-Item $zip
    return $out
}

function Step($label, $target, [scriptblock]$action) {
    if ((Test-Path $target) -and -not $Force) { Write-Host ("[ok]      {0}: already there ({1})" -f $label, $target); return }
    Write-Host ("[fetch]   {0}" -f $label)
    & $action
    if (-not (Test-Path $target)) { throw "${label}: expected ${target} after the download, but it is not there" }
}

# Embree and Open Image Denoise: the zip contains a folder with the name below
Step "Embree 4.4.0" (Join-Path $root "librarys\embree-4.4.0.x64.windows\include") {
    $d = Get-Zip "https://github.com/RenderKit/embree/releases/download/v4.4.0/embree-4.4.0.x64.windows.zip" "embree"
    New-Item -ItemType Directory -Force (Join-Path $root "librarys") | Out-Null
    $src = Get-ChildItem $d -Directory | Select-Object -First 1
    $dst = Join-Path $root "librarys\embree-4.4.0.x64.windows"
    if (Test-Path $dst) { Remove-Item $dst -Recurse -Force }
    Move-Item $src.FullName $dst
}
Step "Open Image Denoise 2.5.1" (Join-Path $root "librarys\oidn-2.5.1.x64.windows\include") {
    $d = Get-Zip "https://github.com/RenderKit/oidn/releases/download/v2.5.1/oidn-2.5.1.x64.windows.zip" "oidn"
    New-Item -ItemType Directory -Force (Join-Path $root "librarys") | Out-Null
    $src = Get-ChildItem $d -Directory | Select-Object -First 1
    $dst = Join-Path $root "librarys\oidn-2.5.1.x64.windows"
    if (Test-Path $dst) { Remove-Item $dst -Recurse -Force }
    Move-Item $src.FullName $dst
}

# GLM: the code includes <glm/glm.hpp>, so GoLD_rayt\external\glm must be the folder that holds glm.hpp
Step "GLM 1.0.1" (Join-Path $proj "external\glm\glm.hpp") {
    $d = Get-Zip "https://github.com/g-truc/glm/releases/download/1.0.1/glm-1.0.1-light.zip" "glm"
    $hpp = Get-ChildItem $d -Recurse -Filter glm.hpp | Where-Object { $_.Directory.Name -eq "glm" } | Select-Object -First 1
    if (-not $hpp) { throw "glm.hpp was not found in the GLM download" }
    $dst = Join-Path $proj "external\glm"
    if (Test-Path $dst) { Remove-Item $dst -Recurse -Force }
    New-Item -ItemType Directory -Force (Split-Path $dst) | Out-Null
    Move-Item $hpp.Directory.FullName $dst
}

# single headers
Step "nlohmann/json 3.11.3" (Join-Path $proj "external\nlohmann\json.hpp") {
    New-Item -ItemType Directory -Force (Join-Path $proj "external\nlohmann") | Out-Null
    Get-File "https://github.com/nlohmann/json/releases/download/v3.11.3/json.hpp" (Join-Path $proj "external\nlohmann\json.hpp")
}
Step "stb_image.h and stb_image_write.h" (Join-Path $proj "include\stb\stb_image_write.h") {
    New-Item -ItemType Directory -Force (Join-Path $proj "include\stb") | Out-Null
    Get-File "https://raw.githubusercontent.com/nothings/stb/master/stb_image.h" (Join-Path $proj "include\stb\stb_image.h")
    Get-File "https://raw.githubusercontent.com/nothings/stb/master/stb_image_write.h" (Join-Path $proj "include\stb\stb_image_write.h")
}

# GLFW (the zip contains glfw-3.4.bin.WIN64)
Step "GLFW 3.4" (Join-Path $proj "external\glfw-3.4.bin.WIN64\include") {
    $d = Get-Zip "https://github.com/glfw/glfw/releases/download/3.4/glfw-3.4.bin.WIN64.zip" "glfw"
    $src = Get-ChildItem $d -Directory | Select-Object -First 1
    $dst = Join-Path $proj "external\glfw-3.4.bin.WIN64"
    if (Test-Path $dst) { Remove-Item $dst -Recurse -Force }
    New-Item -ItemType Directory -Force (Split-Path $dst) | Out-Null
    Move-Item $src.FullName $dst
}

# Dear ImGui: the layout uses docking, so the docking branch (a snapshot; the version the author used was 1.92.6 WIP)
Step "Dear ImGui (docking branch)" (Join-Path $root "external\imgui\imgui.cpp") {
    $d = Get-Zip "https://github.com/ocornut/imgui/archive/refs/heads/docking.zip" "imgui"
    $src = Get-ChildItem $d -Directory | Select-Object -First 1
    $dst = Join-Path $root "external\imgui"
    if (Test-Path $dst) { Remove-Item $dst -Recurse -Force }
    New-Item -ItemType Directory -Force (Split-Path $dst) | Out-Null
    Move-Item $src.FullName $dst
}

Remove-Item $tmp -Recurse -Force -ErrorAction SilentlyContinue
Write-Host ""
Write-Host "Done. Build with:  powershell -ExecutionPolicy Bypass -File tools\build_cli.ps1"
