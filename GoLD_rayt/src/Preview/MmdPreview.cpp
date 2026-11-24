// A real-time OpenGL preview of an MMD model with the toon recipe; see Preview/MmdPreview.hpp.
#include "pch.h"
#include "Preview/MmdPreview.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Cli/PmxSupport.hpp"
#include "IO/PmxLoader.hpp"
#include "IO/PmxMorph.hpp"
#include "IO/PmxPose.hpp"
#include "IO/VpdLoader.hpp"

namespace rayt::preview {

    namespace {

        // ---- shaders (GLSL 3.30) ----

        const char* kMainVs = R"(#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aUv;
uniform mat4 uViewProj;
uniform mat4 uLightViewProj;
out vec3 vPos;
out vec3 vNormal;
out vec2 vUv;
out vec4 vLightPos;
void main() {
    vPos = aPos;
    vNormal = aNormal;
    vUv = aUv;
    vLightPos = uLightViewProj * vec4(aPos, 1.0);
    gl_Position = uViewProj * vec4(aPos, 1.0);
}
)";

        // MMD's self-shadow formula in display values, as ToonIntegrator's shadeLook
        const char* kMainFs = R"(#version 330 core
in vec3 vPos;
in vec3 vNormal;
in vec2 vUv;
in vec4 vLightPos;
uniform vec3 uDiffuse;
uniform vec3 uAmbient;
uniform vec3 uSpecular;
uniform float uShininess;
uniform float uAlpha;
uniform float uLightColor;
uniform sampler2D uTex;
uniform sampler2D uSphere;
uniform sampler2D uToon;
uniform sampler2D uShadowMap;
uniform bool uHasTex;
uniform bool uHasToon;
uniform int uSphereMode;
uniform bool uShadows;
uniform bool uFaceNormals;
uniform vec3 uFaceCenter;
uniform vec3 uFaceRadii;
uniform vec3 uToonShadow;
uniform vec3 uLightDir;
uniform vec3 uEye;
uniform vec3 uCamRight;
uniform vec3 uCamUp;
out vec4 fragColor;

float visibility(vec3 n) {
    vec3 p = vLightPos.xyz / vLightPos.w * 0.5 + 0.5;
    if (p.z > 1.0 || p.x < 0.0 || p.x > 1.0 || p.y < 0.0 || p.y > 1.0) return 1.0;
    float bias = 0.0008 + 0.002 * (1.0 - clamp(dot(n, uLightDir), 0.0, 1.0));
    vec2 texel = 1.0 / vec2(textureSize(uShadowMap, 0));
    float lit = 0.0;
    for (int x = -1; x <= 1; ++x)
        for (int y = -1; y <= 1; ++y)
            lit += (p.z - bias <= texture(uShadowMap, p.xy + vec2(x, y) * texel).r) ? 1.0 : 0.0;
    return lit / 9.0;
}

void main() {
    vec3 n = normalize(vNormal);
    vec3 wo = normalize(uEye - vPos);
    if (uFaceNormals) n = normalize((vPos - uFaceCenter) / (uFaceRadii * uFaceRadii));
    else if (dot(n, wo) < 0.0) n = -n;

    vec3 c = clamp(uDiffuse * uLightColor + uAmbient, 0.0, 1.0);
    float a = uAlpha;
    if (uHasTex) {
        vec4 t = texture(uTex, vUv);
        c *= t.rgb;
        a *= t.a;
    }
    if (a < 0.004) discard;
    if (uSphereMode == 1 || uSphereMode == 2) {
        vec2 su = vec2(dot(n, uCamRight) * 0.5 + 0.5, 0.5 - dot(n, uCamUp) * 0.5);
        vec3 sp = texture(uSphere, su).rgb;
        c = (uSphereMode == 1) ? c * sp : c + sp;
    }
    float ndl = dot(n, uLightDir);
    float vis = uShadows ? visibility(n) : 1.0;
    float comp = min(clamp(3.0 * ndl, 0.0, 1.0), vis);
    vec3 shade = uHasToon ? texture(uToon, vec2(0.5, 1.0 - 0.5 / float(textureSize(uToon, 0).y))).rgb : uToonShadow;
    c *= mix(shade, vec3(1.0), comp);
    if (uShininess > 0.0 && dot(uSpecular, vec3(1.0)) > 0.0) {
        vec3 h = normalize(uLightDir + wo);
        c += uSpecular * uLightColor * pow(max(0.0, dot(n, h)), uShininess) * comp;
    }
    fragColor = vec4(clamp(c, 0.0, 1.0), a);
}
)";

        // The inverted hull: pushed out along the normal by a width in pixels (in clip space, so the width does not depend on the distance)
        const char* kEdgeVs = R"(#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aUv;
uniform mat4 uViewProj;
uniform float uEdgePx;
uniform vec2 uViewport;
out vec2 vUv;
void main() {
    vec4 c = uViewProj * vec4(aPos, 1.0);
    vec4 cn = uViewProj * vec4(aPos + aNormal * 0.01, 1.0);
    vec2 d = cn.xy / cn.w - c.xy / c.w;
    float l = length(d);
    vec2 dir = (l > 1e-7) ? d / l : vec2(0.0);
    c.xy += dir * uEdgePx * 2.0 / uViewport * c.w;
    vUv = aUv;
    gl_Position = c;
}
)";
        const char* kEdgeFs = R"(#version 330 core
in vec2 vUv;
uniform vec3 uEdgeColor;
uniform sampler2D uTex;
uniform bool uHasTex;
out vec4 fragColor;
void main() {
    if (uHasTex && texture(uTex, vUv).a < 0.5) discard;
    fragColor = vec4(uEdgeColor, 1.0);
}
)";

        const char* kShadowVs = R"(#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec2 aUv;
uniform mat4 uLightViewProj;
out vec2 vUv;
void main() {
    vUv = aUv;
    gl_Position = uLightViewProj * vec4(aPos, 1.0);
}
)";
        const char* kShadowFs = R"(#version 330 core
in vec2 vUv;
uniform sampler2D uTex;
uniform bool uHasTex;
void main() {
    if (uHasTex && texture(uTex, vUv).a < 0.5) discard;
}
)";

        GLuint Compile(GLenum type, const char* src, std::string* error) {
            const GLuint s = glCreateShader(type);
            glShaderSource(s, 1, &src, nullptr);
            glCompileShader(s);
            GLint ok = 0;
            glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
            if (!ok) {
                char log[2048];
                glGetShaderInfoLog(s, sizeof(log), nullptr, log);
                if (error) *error = std::string("shader compile: ") + log;
                glDeleteShader(s);
                return 0;
            }
            return s;
        }

        GLuint Link(const char* vs, const char* fs, std::string* error) {
            const GLuint v = Compile(GL_VERTEX_SHADER, vs, error);
            if (!v) return 0;
            const GLuint f = Compile(GL_FRAGMENT_SHADER, fs, error);
            if (!f) { glDeleteShader(v); return 0; }
            const GLuint p = glCreateProgram();
            glAttachShader(p, v);
            glAttachShader(p, f);
            glLinkProgram(p);
            glDeleteShader(v);
            glDeleteShader(f);
            GLint ok = 0;
            glGetProgramiv(p, GL_LINK_STATUS, &ok);
            if (!ok) {
                char log[2048];
                glGetProgramInfoLog(p, sizeof(log), nullptr, log);
                if (error) *error = std::string("shader link: ") + log;
                glDeleteProgram(p);
                return 0;
            }
            return p;
        }

        inline glm::vec3 F(const Vector3& v) { return glm::vec3(float(v.x), float(v.y), float(v.z)); }
        inline GLint U(GLuint program, const char* name) { return glGetUniformLocation(program, name); }

        Vector3 Direction(Real yawDeg, Real pitchDeg) {
            const Real y = glm::radians(yawDeg), p = glm::radians(pitchDeg);
            return Vector3(std::sin(y) * std::cos(p), std::sin(p), std::cos(y) * std::cos(p));
        }

    } // namespace

    Vector3 OrbitCamera::eye() const { return target + Direction(yaw, pitch) * distance; }
    Vector3 OrbitCamera::forward() const { return -Direction(yaw, pitch); }

    MmdPreview::~MmdPreview() {
        release();
        if (m_mainProgram) glDeleteProgram(m_mainProgram);
        if (m_edgeProgram) glDeleteProgram(m_edgeProgram);
        if (m_shadowProgram) glDeleteProgram(m_shadowProgram);
        if (m_msaaFbo) glDeleteFramebuffers(1, &m_msaaFbo);
        if (m_msaaColor) glDeleteRenderbuffers(1, &m_msaaColor);
        if (m_msaaDepth) glDeleteRenderbuffers(1, &m_msaaDepth);
        if (m_resolveFbo) glDeleteFramebuffers(1, &m_resolveFbo);
        if (m_colorTexture) glDeleteTextures(1, &m_colorTexture);
        if (m_shadowFbo) glDeleteFramebuffers(1, &m_shadowFbo);
        if (m_shadowTexture) glDeleteTextures(1, &m_shadowTexture);
        if (m_floorVao) glDeleteVertexArrays(1, &m_floorVao);
        if (m_floorVbo) glDeleteBuffers(1, &m_floorVbo);
    }

    void MmdPreview::release() {
        for (const auto& [tex, id] : m_glTextures) glDeleteTextures(1, &id);
        m_glTextures.clear();
        m_keepTextures.clear();
        if (m_vao) glDeleteVertexArrays(1, &m_vao);
        if (m_vbo) glDeleteBuffers(1, &m_vbo);
        if (m_ebo) glDeleteBuffers(1, &m_ebo);
        m_vao = m_vbo = m_ebo = 0;
        m_indexCount = 0;
        m_ranges.clear();
    }

    bool MmdPreview::ensurePrograms(std::string* error) {
        if (m_mainProgram) return true;
        m_mainProgram = Link(kMainVs, kMainFs, error);
        m_edgeProgram = m_mainProgram ? Link(kEdgeVs, kEdgeFs, error) : 0;
        m_shadowProgram = m_edgeProgram ? Link(kShadowVs, kShadowFs, error) : 0;
        if (!m_shadowProgram) return false;
        for (GLuint p : { m_mainProgram, m_edgeProgram, m_shadowProgram }) {
            glUseProgram(p);
            glUniform1i(U(p, "uTex"), 0);
            if (p == m_mainProgram) {
                glUniform1i(U(p, "uSphere"), 1);
                glUniform1i(U(p, "uToon"), 2);
                glUniform1i(U(p, "uShadowMap"), 3);
            }
        }
        glUseProgram(0);

        // the floor: a 12 m square on y = 0 (position, normal, uv)
        const float s = 6.0f;
        const float floorVerts[] = {
            -s, 0, -s, 0, 1, 0, 0, 0,   s, 0, -s, 0, 1, 0, 1, 0,   s, 0, s, 0, 1, 0, 1, 1,
            -s, 0, -s, 0, 1, 0, 0, 0,   s, 0, s, 0, 1, 0, 1, 1,   -s, 0, s, 0, 1, 0, 0, 1 };
        glGenVertexArrays(1, &m_floorVao);
        glGenBuffers(1, &m_floorVbo);
        glBindVertexArray(m_floorVao);
        glBindBuffer(GL_ARRAY_BUFFER, m_floorVbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(floorVerts), floorVerts, GL_STATIC_DRAW);
        for (int a = 0; a < 3; ++a) glEnableVertexAttribArray(a);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
        glBindVertexArray(0);

        // the shadow map
        glGenTextures(1, &m_shadowTexture);
        glBindTexture(GL_TEXTURE_2D, m_shadowTexture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, kShadowSize, kShadowSize, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glGenFramebuffers(1, &m_shadowFbo);
        glBindFramebuffer(GL_FRAMEBUFFER, m_shadowFbo);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, m_shadowTexture, 0);
        glDrawBuffer(GL_NONE);
        glReadBuffer(GL_NONE);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        return true;
    }

    void MmdPreview::ensureTargets(int width, int height) {
        if (width == m_w && height == m_h && m_msaaFbo) return;
        m_w = width;
        m_h = height;
        if (!m_msaaFbo) {
            glGenFramebuffers(1, &m_msaaFbo);
            glGenRenderbuffers(1, &m_msaaColor);
            glGenRenderbuffers(1, &m_msaaDepth);
            glGenFramebuffers(1, &m_resolveFbo);
            glGenTextures(1, &m_colorTexture);
        }
        GLint maxSamples = 4;
        glGetIntegerv(GL_MAX_SAMPLES, &maxSamples);
        const int samples = std::min(8, static_cast<int>(maxSamples));
        glBindRenderbuffer(GL_RENDERBUFFER, m_msaaColor);
        glRenderbufferStorageMultisample(GL_RENDERBUFFER, samples, GL_RGBA8, width, height);
        glBindRenderbuffer(GL_RENDERBUFFER, m_msaaDepth);
        glRenderbufferStorageMultisample(GL_RENDERBUFFER, samples, GL_DEPTH_COMPONENT24, width, height);
        glBindFramebuffer(GL_FRAMEBUFFER, m_msaaFbo);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, m_msaaColor);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, m_msaaDepth);

        glBindTexture(GL_TEXTURE_2D, m_colorTexture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glBindFramebuffer(GL_FRAMEBUFFER, m_resolveFbo);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_colorTexture, 0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glBindRenderbuffer(GL_RENDERBUFFER, 0);
    }

    GLuint MmdPreview::uploadTexture(const std::shared_ptr<const Texture2D>& tex, bool clampToEdge) {
        if (!tex) return 0;
        if (const auto it = m_glTextures.find(tex.get()); it != m_glTextures.end()) return it->second;
        GLuint id = 0;
        glGenTextures(1, &id);
        glBindTexture(GL_TEXTURE_2D, id);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, tex->width(), tex->height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, tex->rgba8().data());
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        const GLint wrap = clampToEdge ? GL_CLAMP_TO_EDGE : GL_REPEAT;
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrap);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrap);
        m_glTextures[tex.get()] = id;
        m_keepTextures.push_back(tex);
        return id;
    }

    bool MmdPreview::load(const std::filesystem::path& pmx, const std::filesystem::path& pose, std::string* error) {
        const auto t0 = std::chrono::steady_clock::now();
        if (!ensurePrograms(error)) return false;

        io::PmxModel model;
        if (!io::loadPmx(pmx, model, error)) return false;
        if (!pose.empty()) {
            io::VpdPose vpd;
            std::string poseErr;
            if (!io::loadVpd(pose, vpd, &poseErr)) { if (error) *error = "pose: " + poseErr; return false; }
            // the expression entries of the pose file, then the bones (as the ray tracer's PMX scene)
            std::map<int, Real> weights;
            for (const auto& [name, w] : vpd.morphs)
                for (size_t i = 0; i < model.morphs.size(); ++i)
                    if (model.morphs[i].name == name) weights[static_cast<int>(i)] = w;
            std::vector<io::VpdBone> morphBones;
            if (!weights.empty()) model = io::applyMorphs(model, weights, &morphBones);
            io::mergePoseBones(vpd, morphBones);
            model = io::applyPose(model, vpd);
        }

        Vector3 size(1.0);
        const auto mesh = cli::PmxToMesh(model, Real(0), &m_height, &size);
        m_name = model.name;

        std::map<int, std::shared_ptr<Texture2D>> cache;
        auto loadTexture = [&](int index) -> std::shared_ptr<Texture2D> {
            if (index < 0 || index >= static_cast<int>(model.textures.size())) return nullptr;
            if (const auto it = cache.find(index); it != cache.end()) return it->second;
            std::string rel = model.textures[static_cast<size_t>(index)];
            std::replace(rel.begin(), rel.end(), '\\', '/');
            auto tex = Texture2D::load(pmx.parent_path() / cli::PathFromUtf8(rel));
            if (!tex) std::cout << "[Preview] could not load texture " << rel << "\n";
            cache[index] = tex;
            return tex;
        };
        const std::vector<ToonLook> looks = cli::MakePmxToonLooks(model, *mesh, loadTexture, false);

        // vertices: position, normal, uv (float)
        std::vector<float> verts;
        verts.reserve(mesh->p.size() * 8);
        m_lo = Vector3(1e30);
        m_hi = Vector3(-1e30);
        for (size_t i = 0; i < mesh->p.size(); ++i) {
            const Vector3 p = mesh->p[i];
            const Vector3 n = i < mesh->n.size() ? mesh->n[i] : Vector3(0, 1, 0);
            const Point2 uv = i < mesh->uv.size() ? mesh->uv[i] : Point2(0, 0);
            verts.insert(verts.end(), { float(p.x), float(p.y), float(p.z), float(n.x), float(n.y), float(n.z), float(uv.x), float(uv.y) });
            m_lo = glm::min(m_lo, p);
            m_hi = glm::max(m_hi, p);
        }
        std::vector<unsigned int> indices(mesh->indices.begin(), mesh->indices.end());

        // Which winding do the front faces have after the conversion (mirrored z)? The one whose face normal agrees with the vertex normals.
        double agree = 0;
        for (size_t t = 0; t + 2 < indices.size(); t += 3) {
            const Vector3& a = mesh->p[indices[t]];
            const Vector3 fn = glm::cross(mesh->p[indices[t + 1]] - a, mesh->p[indices[t + 2]] - a);
            agree += glm::dot(fn, mesh->n[indices[t]] + mesh->n[indices[t + 1]] + mesh->n[indices[t + 2]]);
        }
        m_ccwFront = agree >= 0;

        release();
        glGenVertexArrays(1, &m_vao);
        glGenBuffers(1, &m_vbo);
        glGenBuffers(1, &m_ebo);
        glBindVertexArray(m_vao);
        glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(verts.size() * sizeof(float)), verts.data(), GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLsizeiptr>(indices.size() * sizeof(unsigned int)), indices.data(), GL_STATIC_DRAW);
        for (int a = 0; a < 3; ++a) glEnableVertexAttribArray(a);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
        glBindVertexArray(0);
        m_indexCount = static_cast<int>(indices.size());

        for (size_t k = 0; k < model.materials.size(); ++k) {
            const auto& pm = model.materials[k];
            DrawRange r;
            r.first = pm.firstIndex;
            r.count = pm.indexCount;
            r.look = looks[k];
            r.doubleSided = pm.doubleSided();
            r.hidden = pm.diffuse.w <= Real(0.001);
            r.texture = uploadTexture(r.look.texture, false);
            r.sphere = uploadTexture(r.look.sphere, true);
            r.toon = uploadTexture(r.look.toon, true);
            m_ranges.push_back(r);
        }
        m_loadMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        std::cout << "[Preview] " << m_name << ": " << mesh->p.size() << " vertices, " << m_indexCount / 3 << " triangles, "
            << m_ranges.size() << " materials, " << m_glTextures.size() << " textures, height " << m_height << " m; "
            << (m_ccwFront ? "counter-clockwise" : "clockwise") << " front faces; " << m_loadMs << " ms\n";
        return true;
    }

    OrbitCamera MmdPreview::framingCamera() const {
        OrbitCamera c;
        const Vector3 size = m_hi - m_lo;
        const Real extent = std::max(m_height, size.x / Real(0.75));
        c.target = Vector3(0, m_height * Real(0.5), 0);
        c.vfov = 35;
        c.distance = extent * Real(1.22) / (Real(2) * std::tan(glm::radians(Real(17.5))));
        c.pitch = glm::degrees(std::atan2(m_height * Real(0.05), c.distance));
        c.yaw = 0;
        return c;
    }

    void MmdPreview::render(int width, int height, const OrbitCamera& camera, const PreviewOptions& o) {
        if (width < 1 || height < 1 || !m_mainProgram) return;
        ensureTargets(width, height);

        GLint prevFbo = 0, prevViewport[4];
        glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
        glGetIntegerv(GL_VIEWPORT, prevViewport);

        const glm::vec3 eye = F(camera.eye()), target = F(camera.target);
        const glm::mat4 view = glm::lookAt(eye, target, glm::vec3(0, 1, 0));
        const float nearZ = std::max(0.01f, float(camera.distance) * 0.02f);
        const glm::mat4 proj = glm::perspective(glm::radians(float(camera.vfov)), float(width) / float(height), nearZ, float(camera.distance) * 20.0f);
        const glm::mat4 viewProj = proj * view;
        const glm::vec3 camRight = glm::normalize(glm::vec3(view[0][0], view[1][0], view[2][0]));
        const glm::vec3 camUp = glm::normalize(glm::vec3(view[0][1], view[1][1], view[2][1]));

        // the key light (toward the light) and its shadow map's orthographic view round the model
        const glm::vec3 lightDir = F(Direction(o.lightAzimuth, o.lightElevation));
        const glm::vec3 center = F((m_lo + m_hi) * Real(0.5));
        const float radius = float(std::max(glm::length(m_hi - m_lo) * Real(0.55), m_height * Real(1.1)));
        const glm::vec3 lightUp = std::abs(lightDir.y) > 0.99f ? glm::vec3(0, 0, 1) : glm::vec3(0, 1, 0);
        const glm::mat4 lightView = glm::lookAt(center + lightDir * radius * 2.0f, center, lightUp);
        const glm::mat4 lightProj = glm::ortho(-radius, radius, -radius, radius, 0.01f, radius * 4.0f);
        const glm::mat4 lightViewProj = lightProj * lightView;

        auto bindTex = [](int unit, GLuint id) { glActiveTexture(GL_TEXTURE0 + unit); glBindTexture(GL_TEXTURE_2D, id); };

        // ---- shadow map ----
        if (o.shadows && loaded()) {
            glBindFramebuffer(GL_FRAMEBUFFER, m_shadowFbo);
            glViewport(0, 0, kShadowSize, kShadowSize);
            glClear(GL_DEPTH_BUFFER_BIT);
            glEnable(GL_DEPTH_TEST);
            glDepthFunc(GL_LESS);
            glDisable(GL_CULL_FACE);
            glDisable(GL_BLEND);
            glUseProgram(m_shadowProgram);
            glUniformMatrix4fv(U(m_shadowProgram, "uLightViewProj"), 1, GL_FALSE, glm::value_ptr(lightViewProj));
            glBindVertexArray(m_vao);
            for (const auto& r : m_ranges) {
                if (r.hidden || r.look.overlay || r.count == 0) continue;
                glUniform1i(U(m_shadowProgram, "uHasTex"), r.texture != 0);
                bindTex(0, r.texture);
                glDrawElements(GL_TRIANGLES, r.count, GL_UNSIGNED_INT, (void*)(static_cast<size_t>(r.first) * sizeof(unsigned int)));
            }
        }

        // ---- the picture ----
        glBindFramebuffer(GL_FRAMEBUFFER, m_msaaFbo);
        glViewport(0, 0, width, height);
        glClearColor(float(o.background.x), float(o.background.y), float(o.background.z), 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
        glEnable(GL_BLEND);
        glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
        glFrontFace(m_ccwFront ? GL_CCW : GL_CW);

        glUseProgram(m_mainProgram);
        const GLuint P = m_mainProgram;
        glUniformMatrix4fv(U(P, "uViewProj"), 1, GL_FALSE, glm::value_ptr(viewProj));
        glUniformMatrix4fv(U(P, "uLightViewProj"), 1, GL_FALSE, glm::value_ptr(lightViewProj));
        glUniform3fv(U(P, "uLightDir"), 1, glm::value_ptr(lightDir));
        glUniform3fv(U(P, "uEye"), 1, glm::value_ptr(eye));
        glUniform3fv(U(P, "uCamRight"), 1, glm::value_ptr(camRight));
        glUniform3fv(U(P, "uCamUp"), 1, glm::value_ptr(camUp));
        glUniform1f(U(P, "uLightColor"), 154.0f / 255.0f);
        glUniform1i(U(P, "uShadows"), o.shadows ? 1 : 0);
        bindTex(3, m_shadowTexture);

        auto setLook = [&](const ToonLook& L, GLuint tex, GLuint sphere, GLuint toon) {
            glUniform3fv(U(P, "uDiffuse"), 1, glm::value_ptr(F(L.diffuse)));
            glUniform3fv(U(P, "uAmbient"), 1, glm::value_ptr(F(L.ambient)));
            glUniform3fv(U(P, "uSpecular"), 1, glm::value_ptr(F(L.specular)));
            glUniform1f(U(P, "uShininess"), float(L.shininess));
            glUniform1f(U(P, "uAlpha"), float(L.alpha));
            glUniform1i(U(P, "uHasTex"), tex != 0);
            glUniform1i(U(P, "uHasToon"), toon != 0);
            glUniform1i(U(P, "uSphereMode"), sphere ? L.sphereMode : 0);
            glUniform3fv(U(P, "uToonShadow"), 1, glm::value_ptr(F(L.toonShadow)));
            glUniform1i(U(P, "uFaceNormals"), L.faceNormals ? 1 : 0);
            glUniform3fv(U(P, "uFaceCenter"), 1, glm::value_ptr(F(L.faceCenter)));
            glUniform3fv(U(P, "uFaceRadii"), 1, glm::value_ptr(F(L.faceRadii)));
            bindTex(0, tex);
            bindTex(1, sphere);
            bindTex(2, toon);
        };

        if (o.floor) {
            ToonLook floorLook;
            floorLook.diffuse = Vector3(0.75);
            floorLook.ambient = Vector3(0.35);
            floorLook.toonShadow = Vector3(0.68, 0.70, 0.80);
            setLook(floorLook, 0, 0, 0);
            glDisable(GL_CULL_FACE);
            glDepthMask(GL_TRUE);
            glBindVertexArray(m_floorVao);
            glDrawArrays(GL_TRIANGLES, 0, 6);
        }

        if (loaded()) {
            glBindVertexArray(m_vao);
            // in the model's order, as MMD draws (blending on, depth writes on); the overlay layers afterwards, slightly in front, without depth writes
            for (int pass = 0; pass < 2; ++pass) {
                for (const auto& r : m_ranges) {
                    if (r.hidden || r.count == 0 || r.look.overlay != (pass == 1)) continue;
                    if (r.doubleSided) glDisable(GL_CULL_FACE);
                    else { glEnable(GL_CULL_FACE); glCullFace(GL_BACK); }
                    if (pass == 1) { glEnable(GL_POLYGON_OFFSET_FILL); glPolygonOffset(-1.0f, -2.0f); glDepthMask(GL_FALSE); }
                    setLook(r.look, r.texture, r.sphere, r.toon);
                    glDrawElements(GL_TRIANGLES, r.count, GL_UNSIGNED_INT, (void*)(static_cast<size_t>(r.first) * sizeof(unsigned int)));
                }
            }
            glDisable(GL_POLYGON_OFFSET_FILL);
            glDepthMask(GL_TRUE);

            // the inverted hull: the edge-flagged materials again, their front faces culled, pushed out in screen space
            if (o.edges) {
                glUseProgram(m_edgeProgram);
                const GLuint E = m_edgeProgram;
                glUniformMatrix4fv(U(E, "uViewProj"), 1, GL_FALSE, glm::value_ptr(viewProj));
                glUniform2f(U(E, "uViewport"), float(width), float(height));
                glEnable(GL_CULL_FACE);
                glCullFace(GL_FRONT);
                glDisable(GL_BLEND);
                const float pxPerUnit = float(height) / 400.0f * float(o.edgeScale);   // as the toon render: 1 pixel for edge 0.5 at 800 pixels
                for (const auto& r : m_ranges) {
                    if (r.hidden || r.count == 0 || r.look.overlay || !r.look.edge || r.look.edgeSize <= Real(0)) continue;
                    glUniform1f(U(E, "uEdgePx"), float(r.look.edgeSize) * pxPerUnit);
                    glUniform3fv(U(E, "uEdgeColor"), 1, glm::value_ptr(F(r.look.edgeColor)));
                    glUniform1i(U(E, "uHasTex"), r.texture != 0);
                    bindTex(0, r.texture);
                    glDrawElements(GL_TRIANGLES, r.count, GL_UNSIGNED_INT, (void*)(static_cast<size_t>(r.first) * sizeof(unsigned int)));
                }
                glCullFace(GL_BACK);
            }
        }

        glBindVertexArray(0);
        glUseProgram(0);
        glDisable(GL_CULL_FACE);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_BLEND);
        glActiveTexture(GL_TEXTURE0);

        // resolve the multisampled picture into the texture ImGui shows
        glBindFramebuffer(GL_READ_FRAMEBUFFER, m_msaaFbo);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, m_resolveFbo);
        glBlitFramebuffer(0, 0, width, height, 0, 0, width, height, GL_COLOR_BUFFER_BIT, GL_NEAREST);

        glBindFramebuffer(GL_FRAMEBUFFER, static_cast<GLuint>(prevFbo));
        glViewport(prevViewport[0], prevViewport[1], prevViewport[2], prevViewport[3]);
    }

} // namespace rayt::preview
