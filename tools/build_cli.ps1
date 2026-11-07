<#
  Builds the command-line renderer (Release build: about 5x faster than Debug) into out\cli and copies the DLLs it needs.

    powershell -ExecutionPolicy Bypass -File tools\build_cli.ps1              RGB build        -> out\cli
    powershell -ExecutionPolicy Bypass -File tools\build_cli.ps1 -Spectral    spectral build   -> out\cli-spectral
    powershell -ExecutionPolicy Bypass -File tools\build_cli.ps1 -Dev         development build -> out\cli-dev (-Spectral: out\cli-dev-spectral)

  -Dev switches off link-time code generation: changing one .cpp then rebuilds in about 5 s instead of 20 s (changing a
  header that many files include takes about 45 s instead of 25 s). The renderer is a few percent slower: use it for trying
  changes, and the normal build for timing and for the final images.

  Works from any PowerShell (MSBuild is found with vswhere when it is not on the PATH). Do not stop it while it runs.
#>
param([switch]$Spectral, [switch]$Dev)
$ErrorActionPreference = "Stop"
$root = Split-Path $PSScriptRoot -Parent
$name = if ($Dev) { if ($Spectral) { "cli-dev-spectral" } else { "cli-dev" } } else { if ($Spectral) { "cli-spectral" } else { "cli" } }
$out = Join-Path $root "out\$name"
$obj = Join-Path $root "out\$name-obj"

# MSBuild: on the PATH (Developer PowerShell), or found with vswhere
$msbuild = (Get-Command msbuild -ErrorAction SilentlyContinue).Source
if (-not $msbuild) {
    $vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio\Installer\vswhere.exe"
    if (Test-Path $vswhere) {
        $msbuild = & $vswhere -latest -requires Microsoft.Component.MSBuild -find "MSBuild\**\Bin\MSBuild.exe" | Select-Object -First 1
    }
}
if (-not $msbuild) { throw "MSBuild was not found. Install Visual Studio 2022 with the 'Desktop development with C++' workload, or open a Developer PowerShell." }

# A full build writes about 400 MB; the compiler fails with C1088 when the disk is full
$free = (Get-PSDrive -Name $root.Substring(0, 1)).Free / 1GB
$needed = if (Test-Path (Join-Path $out "GoLD_rayt.exe")) { 0.3 } else { 1.0 }   # a rebuild writes little; the first build writes about 400 MB
if ($free -lt $needed) { throw ("Only {0:N1} GB are free on drive {1}: the build needs about {2:N1} GB. Free some space and run it again." -f $free, $root.Substring(0, 1), $needed) }

if ($Spectral) { $env:CL = "/DRAYT_SPECTRAL_MODE" } else { Remove-Item Env:CL -ErrorAction SilentlyContinue }
New-Item -ItemType Directory -Force $out, $obj | Out-Null
Write-Host "MSBuild: $msbuild"
Write-Host ("Building the {0}{1} CLI into {2} ..." -f $(if ($Spectral) { "spectral" } else { "RGB" }), $(if ($Dev) { " (development, no link-time code generation)" } else { "" }), $out)
& $msbuild (Join-Path $root "GoLD_rayt\GoLD_rayt.vcxproj") /nologo /m /v:minimal /p:Configuration=Release /p:Platform=x64 /p:RaytCli=true "/p:WholeProgramOptimization=$(if ($Dev) { 'false' } else { 'true' })" "/p:OutDir=$out\" "/p:IntDir=$obj\"
$code = $LASTEXITCODE
Remove-Item Env:CL -ErrorAction SilentlyContinue
if ($code -ne 0) { throw "The build failed (exit code $code). The first 'error' line above says why." }

# the DLLs of Embree (with TBB) and Intel Open Image Denoise (CPU device only: the GPU devices are large and not needed)
foreach ($dll in "embree4.dll", "tbb12.dll", "tbbmalloc.dll") {
    Copy-Item (Join-Path $root "librarys\embree-4.4.0.x64.windows\bin\$dll") $out -Force
}
foreach ($dll in "OpenImageDenoise.dll", "OpenImageDenoise_core.dll", "OpenImageDenoise_device_cpu.dll") {
    Copy-Item (Join-Path $root "librarys\oidn-2.5.1.x64.windows\bin\$dll") $out -Force
}

Write-Host ""
Write-Host "Done: $out\GoLD_rayt.exe"
Write-Host "Run it from the project folder:"
Write-Host "    cd GoLD_rayt"
Write-Host "    ..\out\$name\GoLD_rayt.exe scenes\gold_bowl_tungsten.json --spp 16 --width 400"
