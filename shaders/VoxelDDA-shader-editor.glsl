#version 300 es
#ifdef GL_ES
precision highp float;
precision highp int;
#endif

// ---------- SETTINGS ----------
// QUALITY: 0 = cheap, 1 = slightly low (mirip medium), 2 = medium, 3 = high
// 0:) Bentuk blok acak random, pohon, caves, clouds otomatis off
// 1:) Terrain, kalau aktif daun nya berbentuk: 3x3 2 lapis daun, lalu pada bagian atas berbentuk plus
// 2:) Daun pada pohon lebih detail: 4x4 dua lapis, bentuk plus seperti kualitas 1, kawah aktif
// 3:) Pohon lebih lebat dan lebih acak beserta daunnya, lebih berat
#define QUALITY 0

// Terrain stuff
#define ENABLE_TREES 1
#define ENABLE_CAVES 1
#define ENABLE_CLOUDS 1
#define ENABLE_ATMOSPHERE  1
#define ENABLE_SHADOW 1

// Bagian bawah akan terus dirender sejauh distance tanpa ditutup blok
#define ENDLESS_LOWER_HEIGHT 1
#define USE_TOUCH 1

// Karena belum ditambahkan bouncing light, anggap saja ini pajangan
#define BLOCK_LIGHT 0

// Fog stuff
#define FOG_DENSITY 0.0065
#define FOG_START 24.0
#define MAX_STEPS 40

#define T_LIMIT 180.0

// Ketinggian awan
#define CLOUD_HEIGHT 31

// Pengaturan kamera
#define CAMERA_HEIGHT 19
#define CAMERA_FOV 55 
#define CAMERA_SENSIVITY 2

// Pengaturan sunPos
#define SUN_PATH_ROTATION 10
#define SUN_ROT_SPEED 0.02

uniform vec2 resolution;
uniform vec2 touch;
uniform float time;

/* Nama sampler (customisasi sendiri, 
saran untuk sampler grass top dan leaves berwarna abu abu, 
karena nanti ada pencampuran warna hijau nanti) */
uniform sampler2D grass_top;
uniform sampler2D oak_leaves;
uniform sampler2D stone;
uniform sampler2D grass_side;
uniform sampler2D dirt;
uniform sampler2D oak_log_side;
uniform sampler2D oak_log_top;
/* Custom-able */

#define samplerGrass    grass_top
#define samplerLeaves   oak_leaves
#define samplerDirt     dirt
#define samplerStone    stone
#define samplerWood     oak_log_top
#define samplerWoodSide  oak_log_side
#define samplerGrassSide grass_side

out vec4 outp;
const float PI        = 3.141592653589793;
const float tileScale = 1.0;
const float atlasLod  = 0.0;
vec3 ro = vec3(0, 6372.0e3, 0);
float psCat = 35.0;
float rPl = 6372e3;
float rAt = 6471e3;
vec3 rlh = vec3(5.5e-6, 13e-6, 22.4e-6);
float kMi = 37e-6;
float rshl = 10e3;
float rlm  = 1.25e3;
float attn = 0.58;

struct TraceResult {
    bool hit;
    float t;
    ivec3 cell;
    vec3 normal;
};

const int BLOCK_EMISSIVE = 6;
const int BLOCK_CLOUD    = 7;

// --------- Random / Hash utilities ---------
// Purpose: cheap, deterministic pseudo-random values per coordinate.
// Note: fract(sin(dot(...))*C) is fast and common in shaders (see StackOverflow link above),
// but may have statistical shortcomings — consider integer/bitwise hashes (murmur) for
// better distribution in high-quality renderers.
float hash12(vec2 p) {
    return fract(sin(dot(p, vec2(127.1,311.7))) * 43758.5453123);
}

float hash13(vec3 p) {
    return fract(sin(dot(p, vec3(127.1,311.7,74.7))) * 43758.5453123);
}

// Small utility: floor vec3 -> ivec3 (used for lattice indexing)
ivec3 ivec3_floor(vec3 v) {
    return ivec3(int(floor(v.x)), int(floor(v.y)), int(floor(v.z)));
}

// Alternative rand one-liner variant used elsewhere in the codebase
float rnd(in vec2 co) {
    return fract(sin(dot(co.xy, vec2(12.9898,78.233))) * 43758.5453);
}

// Returns two t-values (entry, exit) or big values if no intersection.
vec2 raysi(vec3 r0, vec3 rd, float sr) {
    float a = dot(rd, rd);
    float b = 1.0 * dot(rd, r0);
    float c = dot(r0, r0) - (sr * sr);
    float d = (b*b) - 4.0*a*c;
    if (d < 0.0) return vec2(1e2, -1e5);
    return vec2(
        (-b - sqrt(d)) / (2.0 * a),
        (-b + sqrt(d)) / (2.0 * a)
    );
}

// Implementation: value-noise on unit lattice with Hermite smoothing f*f*(3-2f)
// and trilinear interpolation between the 8 corners. See Catlike Coding value-noise tutorial.
float noise3(vec3 p) {
    ivec3 i = ivec3(floor(p));
    vec3 f = fract(p);
    f = f * f * (3.0 - 2.0 * f); // Hermite smoothing (smoothstep-like)

    float n000 = hash13(vec3(i.x+0, i.y+0, i.z+0));
    float n100 = hash13(vec3(i.x+1, i.y+0, i.z+0));
    float n010 = hash13(vec3(i.x+0, i.y+1, i.z+0));
    float n110 = hash13(vec3(i.x+1, i.y+1, i.z+0));
    float n001 = hash13(vec3(i.x+0, i.y+0, i.z+1));
    float n101 = hash13(vec3(i.x+1, i.y+0, i.z+1));
    float n011 = hash13(vec3(i.x+0, i.y+1, i.z+1));
    float n111 = hash13(vec3(i.x+1, i.y+1, i.z+1));

    float nx00 = mix(n000, n100, f.x);
    float nx10 = mix(n010, n110, f.x);
    float nx01 = mix(n001, n101, f.x);
    float nx11 = mix(n011, n111, f.x);
    float nxy0 = mix(nx00, nx10, f.y);
    float nxy1 = mix(nx01, nx11, f.y);
    return mix(nxy0, nxy1, f.z);
}

// !Fbm ? Yes : No
float fbm3(vec3 p) {
    return 0.5 * noise3(p);
}


// faceUV: choose UV coords per face normal (top/side/front) for block atlasing
vec2 faceUV(vec3 localPos, vec3 faceNormal) {
    vec2 uv;
    if (abs(faceNormal.x) > 0.5)     uv = vec2(localPos.z, localPos.y) * tileScale;
    else if (abs(faceNormal.y) > 0.5) uv = vec2(localPos.x, localPos.z) * tileScale;
    else                              uv = vec2(localPos.x, localPos.y) * tileScale;
    return uv;
}

// Fast 'fetch' for a hardcoded grass green sample from an atlas (micro-optimisation)
vec3 fetchGrassSideGreenPixelFast() {
    return texture(samplerGrassSide, vec2(110.0 / 256.0, 222.0 / 256.0)).rgb;
}

// Material samplers: add per-cell tint/brightness variation using hash to avoid repetition
// Biar bervariasi aja sih brightnessnya
vec3 sampleLeaves(ivec3 cell, vec2 uv, vec3 faceNormal) {
    vec3 base = textureLod(samplerLeaves, uv, atlasLod).rgb;
    float v = hash12(vec2(float(cell.x) * 0.17, float(cell.y) * 0.29));
    base *= (0.82 + 0.28 * v);
    if (abs(faceNormal.y) > 0.5 && faceNormal.y < 0.0) base *= 0.70;
    vec3 grassGreen = fetchGrassSideGreenPixelFast();
    base = mix(base, grassGreen * length(base), 0.72);
    return base;
}

vec3 sampleWood(ivec3 cell, vec2 uv, vec3 faceNormal) {
	if (faceNormal.y > 0.0 || -faceNormal.y > 0.0){
    return textureLod(samplerWood, uv, atlasLod).rgb;
    } else {
    	return textureLod(samplerWoodSide, uv, atlasLod).rgb;
}    
}

// Top/side mix for grass: top uses grass texture tinted slightly by grass green sample
vec3 samplerGrassTopSide(vec3 faceNormal, vec2 uv, float atlasLodIn) {
    if (faceNormal.y > 0.0) {
        vec3 top = textureLod(samplerGrass, uv, atlasLodIn).rgb;
        vec3 tint = fetchGrassSideGreenPixelFast();
        float tintStrength = 0.92;
        vec3 result = mix(top, top * tint, tintStrength);
        return clamp(result, 0.0, 1.0);
    } else if (-faceNormal.y > 0.0) {
        return textureLod(samplerDirt, uv, atlasLodIn).rgb;
    } else {
        return textureLod(samplerGrassSide, uv, atlasLodIn).rgb;
    }
}

// Emissive color generation per cell using hash to decide tint and color choice
vec3 getEmissiveColor(ivec3 cell) {
    float r = hash13(vec3(float(cell.x) * 0.127, float(cell.y) * 0.223, float(cell.z) * 0.311));
    float tint = 0.6 + 0.8 * hash12(vec2(float(cell.x) * 0.13, float(cell.z) * 0.17));
    if (r < 0.6) return vec3(1.0, 0.86, 0.5) * tint;
    return vec3(0.5, 0.75, 1.0) * (tint * 0.9);
}

float getEmissiveIntensity(ivec3 cell) {
    float r = hash13(vec3(float(cell.x) * 0.37, float(cell.y) * 0.19, float(cell.z) * 0.29));
    return mix(1.2, 3.0, r);
}

// proceduralTerrainHeight: single-octave FBM scaling per quality level
float proceduralTerrainHeight(float fx, float fz, int q) {
    float s = (q == 1) ? 0.045 : (q == 2) ? 0.055 : 0.06;
    float h = fbm3(vec3(fx * s, 0.0, fz * s));
    float verticalScale = (q == 1) ? 10.0 : (q == 2) ? 11.0 : 12.0;
    float base = (q == 1) ? 3.5 : 4.0;
    return h * verticalScale + base;
}

// hasTreeAtColumn: probabilistic placement using hash thresholds (quality-dependent)
bool hasTreeAtColumn(int cx, int cz, int quality) {
#if ENABLE_TREES
    float r = hash12(vec2(float(cx) * 0.374, float(cz) * 0.912));
    if (quality == 1) return (r > 0.995);
    if (quality == 2) return (r > 0.985);
    return (r > 0.97);
#else
    return false;
#endif
}

// getBlockType: central world rule function — decides block content per cell (terrain, caves, trees, emissives)
int getBlockType(ivec3 cell) {
    const int q = QUALITY;
    #if ENDLESS_LOWER_HEIGHT
    #else
    if (cell.y <= -50) return 3;
    #endif
    if (q == 0) {
        
        float r = hash13(vec3(float(cell.x) * 0.12, float(cell.y) * 0.08, float(cell.z) * 0.15));
        float density = 0.30 + 0.360 * r + 0.15;

        if (density > ((cell.y >= 15) ? 0.8069 : 0.6)) {
            if (r > 0.999) return BLOCK_EMISSIVE;
            float pick = hash12(vec2(float(cell.x) * 0.19, float(cell.z) * 0.23));
            if (pick < 0.33) return 1;
            else if (pick < 0.4) return 4;
            else if (pick < 0.66) return 2;
            else return 3;
        }
        if (cell.y >= 15) return 0;

        return 0;
    }
    else {
        const int cloudY = CLOUD_HEIGHT;

#if ENABLE_CLOUDS
        if (cell.y == cloudY) {
            float ch = hash12(floor(vec2(float(cell.x) * 0.21, float(cell.z) * 0.13) * 2.0));
            if (ch > 0.74) return BLOCK_CLOUD;
            return 0;
        }
#else
        if (cell.y == cloudY) return 0;
#endif

        if (cell.y >= 19) return 0;

        float fx = float(cell.x);
        float fz = float(cell.z);

        float terrainH_f = proceduralTerrainHeight(fx, fz, q);
        int terrainH = int(floor(terrainH_f));

#if ENABLE_CAVES
        bool allowCaves = true;
#else
        bool allowCaves = false;
#endif

        bool isCave = false;
        if (allowCaves) {
            float caveNoise = 2.0 * fbm3(vec3(float(cell.x) * 0.10, float(cell.y) * 0.10, float(cell.z) * 0.10));
            isCave = (caveNoise > (q == 2 ? 0.68 : 0.66));
        }

        float rr = hash13(vec3(float(cell.x) * 0.23, float(cell.y) * 0.17, float(cell.z) * 0.19));
#if BLOCK_LIGHT
        if (isCave) {
            float caveThresh = (q == 1) ? 0.9975 : (q == 2) ? 0.9965 : 0.995;
            if (rr > caveThresh) return BLOCK_EMISSIVE;
        }

        if (cell.y == terrainH) {
            float surfThresh = (q == 1) ? 0.9995 : (q == 2) ? 0.999 : 0.998;
            if (rr > surfThresh) return BLOCK_EMISSIVE;
        }
        #endif

        if (cell.y <= terrainH && !isCave) {
            if (cell.y == terrainH) return 1;
            if (cell.y >= terrainH - 3) return 2;
            return 3;
        }

#if ENABLE_TREES
        if (hasTreeAtColumn(cell.x, cell.z, q)) {
            int trunkH = 3 + int(floor(hash12(vec2(float(cell.x) * 0.11, float(cell.z) * 0.08)) * 3.0));
            int baseY = terrainH + 1;
            int topY = baseY + trunkH;
            if (cell.y >= baseY && cell.y < topY) return 4;
        }
#endif

        int SEARCH_MAX = (q == 1) ? 1 : 3;
        int maxOffset = q == 1 ? 1 : (q == 2) ? 2 : 3;

// Leaves placement heuristics: limited region under tree top, quality-dependent search radius
        for (int ox = -SEARCH_MAX; ox <= SEARCH_MAX; ++ox) {
            for (int oz = -SEARCH_MAX; oz <= SEARCH_MAX; ++oz) {
                if (abs(ox) > maxOffset || abs(oz) > maxOffset) continue;
#if ENABLE_TREES
                int colX = cell.x + ox;
                int colZ = cell.z + oz;
                if (!hasTreeAtColumn(colX, colZ, q)) continue;
                float th_f = proceduralTerrainHeight(float(colX), float(colZ), q);
                int th = int(floor(th_f));
                int trunkH2 = 3 + int(floor(hash12(vec2(float(colX) * 0.11, float(colZ) * 0.08)) * 3.0));
                int baseY2 = th + 1;
                int topY2 = baseY2 + trunkH2;
                int dx = cell.x - colX;
                int dz = cell.z - colZ;
                int dy = cell.y - topY2;
// --
                if (dy == 0 || dy == 1) {
                    if (dx >= -1 && dx <= 2 && dz >= -1 && dz <= 2) return 5;
                    if (dx == -2 && dz >= -1 && dz <= 2) return 5;
                }
                if (dy == 2) {
                    if (dx == 0 && dz == 0) return 5;
                    if ((abs(dx) <= 1 && dz == 0) || (abs(dz) <= 1 && dx == 0)) return 5;
                }
#endif
            }
        }
    }

    return 0;
}

// Map block type -> sampled albedo color. Manual mapping
vec3 sampleBlockAlbedoTex_manual(ivec3 cell, vec3 localPos, vec3 faceNormal) {
    int b = getBlockType(cell);
    vec2 uv = faceUV(localPos, faceNormal);
    uv = fract(uv);
    if (b == 1)  return samplerGrassTopSide(faceNormal, uv, atlasLod);
    if (b == 2)  return textureLod(samplerDirt, uv, atlasLod).rgb;
    if (b == 3)  return textureLod(samplerStone, uv, atlasLod).rgb;
    if (b == 4)  return sampleWood(cell, uv, faceNormal);
    if (b == 5)  return sampleLeaves(cell, uv, faceNormal);
    if (b == BLOCK_EMISSIVE) return getEmissiveColor(cell);
    if (b == BLOCK_CLOUD)    return vec3(0.96, 0.96, 0.98);

    return vec3(0.0);
}

// init_voxel_dda: initialize stepping variables tMax and tDelta and starting voxel coordinates.
// See Amanatides & Woo: https://www.cs.yorku.ca/~amana/research/grid.pdf
void init_voxel_dda(vec3 rayOrigin, vec3 rayDir,
                    out ivec3 cellCoord, out ivec3 cellSign,
                    out vec3 tMaxOut, out vec3 tDeltaOut)
{
    cellCoord = ivec3_floor(rayOrigin);
    cellSign  = ivec3(rayDir.x > 0.0 ? 1 : -1,
                      rayDir.y > 0.0 ? 1 : -1,
                      rayDir.z > 0.0 ? 1 : -1);

    const float EPS = 1e-6;

    if (abs(rayDir.x) < EPS) { tDeltaOut.x = 1e20; tMaxOut.x = 1e20; }
    else {
        tDeltaOut.x = abs(1.0 / rayDir.x);
        float nextBoundaryX = float(cellCoord.x) + (rayDir.x > 0.0 ? 1.0 : 0.0);
        tMaxOut.x = (nextBoundaryX - rayOrigin.x) / rayDir.x;
        if (tMaxOut.x < 0.0) tMaxOut.x = 0.0;
    }

    if (abs(rayDir.y) < EPS) { tDeltaOut.y = 1e20; tMaxOut.y = 1e20; }
    else {
        tDeltaOut.y = abs(1.0 / rayDir.y);
        float nextBoundaryY = float(cellCoord.y) + (rayDir.y > 0.0 ? 1.0 : 0.0);
        tMaxOut.y = (nextBoundaryY - rayOrigin.y) / rayDir.y;
        if (tMaxOut.y < 0.0) tMaxOut.y = 0.0;
    }

    if (abs(rayDir.z) < EPS) { tDeltaOut.z = 1e20; tMaxOut.z = 1e20; }
    else {
        tDeltaOut.z = abs(1.0 / rayDir.z);
        float nextBoundaryZ = float(cellCoord.z) + (rayDir.z > 0.0 ? 1.0 : 0.0);
        tMaxOut.z = (nextBoundaryZ - rayOrigin.z) / rayDir.z;
        if (tMaxOut.z < 0.0) tMaxOut.z = 0.0;
    }
}

// - Amanatides & Woo, "A Fast Voxel Traversal Algorithm for Ray Tracing" (grid.pdf)
//   https://www.cs.yorku.ca/~amana/research/grid.pdf
TraceResult traceVoxel(vec3 rayOrigin, vec3 rayDir, float tLimit, int maxSteps) {
    TraceResult res;
    res.hit = false;
    res.t = 0.0;
    res.cell = ivec3(0);
    res.normal = vec3(0.0);

    ivec3 cell; ivec3 cellSign; vec3 tMax; vec3 tDelta;
    init_voxel_dda(rayOrigin, rayDir, cell, cellSign, tMax, tDelta);

    int startBlock = getBlockType(cell);
    if (startBlock != 0) {
        res.hit = true;
        res.t = 0.0;
        res.cell = cell;
        res.normal = vec3(0.0, 1.0, 0.0);
        return res;
    }

    float t = 0.0;

    for (int i = 0; i < maxSteps; ++i) {
        if (tMax.x <= tMax.y && tMax.x <= tMax.z) {
            t = tMax.x;
            cell.x += cellSign.x;
            tMax.x += tDelta.x;
            res.normal = vec3(-float(cellSign.x), 0.0, 0.0);
        }
        else if (tMax.y <= tMax.x && tMax.y <= tMax.z) {
            t = tMax.y;
            cell.y += cellSign.y;
            tMax.y += tDelta.y;
            res.normal = vec3(0.0, -float(cellSign.y), 0.0);
        }
        else {
            t = tMax.z;
            cell.z += cellSign.z;
            tMax.z += tDelta.z;
            res.normal = vec3(0.0, 0.0, -float(cellSign.z));
        }

        if (t > tLimit) break;

        int btype = getBlockType(cell);
        if (btype != 0) {
            res.hit = true;
            res.t = t;
            res.cell = cell;
            return res;
        }
    }

    return res;
}

// Modified code with removing loop
// https://github.com/wwwtyro/glsl-atmosphere/blob/master/index.glsl
vec3 atmosphere(vec3 r, vec3 r0, vec3 pSun, float iSun,
                float rPlanet, float rAtmos,
                vec3 kRlh, float kMie, float shRlh,
                float shMie, float g)
{
    vec2 psi = raysi(r0, r, rAtmos);
    if (psi.x > psi.y) return vec3(0.0);
    psi.y = min(psi.y, raysi(r0, r, rPlanet).x);

    float iStepSize = (psi.y - psi.x) / 18.0;
    float iTime = 0.0;

    vec3 totalRlh = vec3(0.0);
    vec3 totalMie = vec3(0.0);

    float iOdRlh = 0.0;
    float iOdMie = 0.0;

    float mu = dot(r, pSun);
    float mumu = mu * mu;
    float gg = g * g;

    float pRlh = clamp(3.0 / (12.0 * PI) * (1.0 + mumu), 0.0, 1.0);
    float pMie = clamp(3.0 / (22.0 * PI) * ((1.0 - gg) * (mumu + 1.0)) /
                       (pow(1.0 + gg - 2.0 * mu * g, 1.5) * (2.0 + gg)), 0.0, 1.0);

    vec3 iPos = r0 + r * (iTime + iStepSize * 0.5);
    float iHeight = length(iPos) - rPlanet;

    float odStepRlh = exp2(-iHeight / shRlh) * iStepSize;
    float odStepMie = pRlh * exp2(-iHeight / shMie) * iStepSize;

    iOdRlh += odStepRlh;
    iOdMie += odStepMie;

    float jStepSize = raysi(iPos, pSun, rAtmos).y / 5.0;
    float jTime = 0.0;

    vec3 jPos = iPos + pSun * (jTime + jStepSize * 0.5);
    float jHeight = length(jPos) - rPlanet;

    float jOdRlh = 0.0;
    float jOdMie = 0.0;

    jOdRlh += exp(-jHeight / shRlh) * jStepSize;
    jOdMie += exp(-jHeight / shMie) * jStepSize;

    vec3 attn2 = exp(-(kMie * (iOdMie + jOdMie) + kRlh * (iOdRlh + jOdRlh)));

    totalRlh += odStepRlh * attn2;
    totalMie += odStepMie * attn2;

    return iSun * (pRlh * kRlh * totalRlh + pMie * totalRlh * 5e-6 + (pMie * kMie * totalMie));
}

// I took the total rayleigh code, I was really lazy to make a new code, so I just cloned the code then deleted it and left the total rayleigh, the position used is vec3(0, 1, 0)
vec3 absorp(vec3 r, vec3 r0, vec3 pSun, float iSun,
            float rPlanet, float rAtmos, vec3 kRlh, float shRlh)
{
    vec2 psi = raysi(r0, r, rAtmos);
    if (psi.x > psi.y) return vec3(0.0);
    psi.y = min(psi.y, raysi(r0, r, rPlanet).x);

    float iStepSize = (psi.y - psi.x) / 18.0;
    float iTime = 1.0;

    vec3 iPos = r0 + r * (iTime + iStepSize * 0.5);
    float iHeight = length(iPos) - rPlanet;

    float odStepRlh = exp2(-iHeight / shRlh) * iStepSize;

    float jStepSize = raysi(iPos, pSun, rAtmos).y / 5.0;

    vec3 jPos = iPos + pSun * (0.0 + jStepSize * 0.5);
    float jHeight = length(jPos) - rPlanet;

    float jOdRlh = 0.0;
    jOdRlh += exp(-jHeight / shRlh) * jStepSize;

    vec3 attn2 = exp(-(kRlh * (odStepRlh + jOdRlh)));
    vec3 totalRlh = odStepRlh * attn2;

    return iSun * (totalRlh * 5e-6);
}

/* Sun direction: t drives elevation (sin/cos). SUN_PATH_ROTATION offsets azimuth. */
// - Scratchapixel — camera & ray generation tutorial
//   https://scratchapixel.com/lessons/3d-basic-rendering/perspective-and-camera
vec3 sunVec(float t) {
    return normalize(vec3(radians(float(SUN_PATH_ROTATION)), sin(t), cos(t)));
}

/* Camera position over time:
   - forward along Z (rail) + simple sway for motion feel.
   - CAMERA_HEIGHT is constant Y. */
vec3 camera_position_from_time(float t) {
    float speed = 3.0;
    float startZ = 48.0;
    float swayAmp = 1.2;
    float sway = sin(t * 0.8) * 0.5 + 0.2 * sin(t * 0.7);
    return vec3(sway * swayAmp, float(CAMERA_HEIGHT), startZ + t * speed);
}

/* Camera forward vector from touch + small procedural motion:
   - map touch -> [-1,1], apply deadzone, non-linear response, scale sensitivity,
   - compute yaw/pitch, return forward vector. */
vec3 camera_forward_vector(float t) {
    // Normalize touch to [0..1] independent of resolution, then center to [-1..1]
    vec2 touchUV = touch / max(resolution, vec2(1.0)); // avoid /0
    vec2 tn = (touchUV - vec2(0.5)) * 2.0;

    float aspect = resolution.x / max(1.0, resolution.y);
    tn.x *= aspect;

    tn.y = -tn.y;
    tn.y += 0.5;


    #if USE_TOUCH
        //
    #else
        // 
        tn = vec2(0.0, 1.0);
    #endif

    float baseSensitivity = float(CAMERA_SENSIVITY);
    float sensScale = baseSensitivity; // removed screen-dependent scaling

    // Input shaping
    float deadzone = 0.02;
    float nonLinearPow = 1.15;

    float len = length(tn);
    if (len < deadzone) {
        tn = vec2(0.0);
    } else {
        float m = (len - deadzone) / (1.0 - deadzone);
        tn *= (m / len);
    }

    tn = sign(tn) * pow(abs(tn), vec2(nonLinearPow));
    tn *= sensScale;

    // mapping to yaw/pitch limits
    float maxYaw = 1.2;
    float maxPitch = 0.8;
    float yawFromTouch = tn.x * maxYaw;
    float pitchFromTouch = tn.y * maxPitch;

    // subtle built-in motion
    float yaw = yawFromTouch + sin(t * 0.9) * 0.08 + sin(t * 0.28) * 0.12;
    float pitch = -0.45 + pitchFromTouch + sin(t * 1.1) * 0.05;

    float cx = cos(pitch);
    vec3 fwd;
    fwd.x = sin(yaw) * cx;
    fwd.y = sin(pitch);
    fwd.z = -cos(yaw) * cx;

    return normalize(-fwd);
}

/* Main: build ray, trace voxel, shade hit or sky, apply fog. */
void main() {
    // NDC and aspect correction
    vec2 ndc = (gl_FragCoord.xy / resolution.xy) * 2.0 - 1.0;
    ndc.x *= resolution.x / resolution.y;

    // camera basis
    vec3 camPos = camera_position_from_time(time);
    vec3 camFwd = camera_forward_vector(time);
    vec3 camRight = normalize(cross(camFwd, vec3(0.0, 1.0, 0.0)));
    vec3 camUp = normalize(cross(camRight, camFwd));

    // pinhole ray (FOV)
    float fov = float(CAMERA_FOV) * PI / 180.0;
    float px = ndc.x * tan(fov * 0.5);
    float py = ndc.y * tan(fov * 0.5);
    vec3 rayDir = normalize(camFwd + camRight * px + camUp * py);

    // trace voxels (Amanatides & Woo DDA)
    int maxSteps = MAX_STEPS * (QUALITY == 0 ? 1 : QUALITY == 1 ? 2 : QUALITY == 2 ? 4 : QUALITY == 3 ? 4 : 0);
    float tLimit = T_LIMIT;
    TraceResult tr = traceVoxel(camPos, rayDir, tLimit, maxSteps);

    bool hit = tr.hit;
    vec3 hitPosition = vec3(0.0);
    ivec3 hitCell = ivec3(0);
    vec3 hitFaceNormal = vec3(0.0);

    if (hit) {
        if (tr.t <= 0.0) hitPosition = camPos;
        else            hitPosition = camPos + rayDir * (tr.t + 1e-5);
        hitCell = tr.cell;
        hitFaceNormal = tr.normal;
    }

    // shading prep
    vec3 finalColor = vec3(0.0);
    vec3 N = normalize(hitFaceNormal);
    vec3 L = normalize(sunVec(time * SUN_ROT_SPEED));
    vec3 V = normalize(rayDir);

vec3 absorpColor = 
        	absorp(vec3(0,1,0), ro, L, psCat, rPl, rAt, rlh, rshl);

#if ENABLE_ATMOSPHERE
    vec3 color = atmosphere(rayDir, ro, L, psCat, rPl, rAt, rlh, kMi, rshl, rlm, attn); 
    vec3 sky = 1.0 - exp(-color); // Reducing so it doesn't overbright
#else
    vec3 sky = absorpColor; // Just absorb color
#endif

    if (hit) {
        // local coords & albedo mipmapping
        vec3 local = hitPosition - vec3(float(hitCell.x), float(hitCell.y), float(hitCell.z));
        vec3 localPos = fract(local);
        int b = getBlockType(hitCell);
        vec3 albedo = sampleBlockAlbedoTex_manual(hitCell, localPos, hitFaceNormal);

        // basic lighting
        float nDotL = max(0.0, dot(N, L));
        float lambert = nDotL;
        float rim = pow(1.0 - max(0.0, dot(N, V)), 2.2) * 0.12;

        vec3 lightColor = vec3(1.0, 0.96, 0.85);
        vec3 ambient = vec3(0.12, 0.14, 0.18);

        // single hard shadow, copying voxel tracing
        float shadow = 1.0;
        #if ENABLE_SHADOW
            TraceResult sh = traceVoxel(hitPosition + N * 0.01, L, tLimit, maxSteps);
            if (sh.hit) shadow = 0.0;
        #endif
        
        vec3 directLight   = lightColor * lambert;
        vec3 indirectLight = 0.25 * albedo;

        if (b == BLOCK_EMISSIVE) {
            vec3 emission = getEmissiveColor(hitCell) * getEmissiveIntensity(hitCell);
            indirectLight += 0.35 * emission;
            directLight += emission * 0.5;
        }

        vec3 diffAlbedo = albedo * absorpColor;
        vec3 lit = (shadow * directLight + ambient) * diffAlbedo + (indirectLight * absorpColor) + rim * 0.6 * absorpColor;

        // compact microfacet-like specular (GGX-like D + Schlick F)
        // https://github.com/glslify/glsl-ggx/blob/master/index.glsl
// - Walter et al., "Microfacet Models for Refraction through Rough Surfaces" (GGX/BRDF theory)
//   https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
        float roughness = 0.54;
        float F0 = 0.14;
        float alpha = roughness * roughness;

        vec3 H = normalize(L - V);
        float dotLH = max(0.0, dot(L, H));
        float dotNH = max(0.0, dot(N, H));
        float dotNL = max(0.0, dot(N, L));

        float alphaSqr = alpha * alpha;
        float denom = dotNH * dotNH * (alphaSqr - 1.0) + 1.0;
        float D = alphaSqr / (PI * denom * denom + 1e-6);
        float F = F0 + (1.0 - F0) * pow(1.0 - dotLH, 5.0);
        float k = 0.5 * alpha;
        float k2 = k * k;

        lit += shadow * (dotNL * D * F / (max(1e-6, dotLH * dotLH * (1.0 - k2)) + k2)) * absorpColor * lit;

        lit *= max(0.7, shadow);

        // simple post process
        finalColor = 1.0 - exp(-2.0 * lit);
    }
    else {
        finalColor = sky;
    }

    // exponential fog (Beer–Lambert style)
// - Beer–Lambert law (fog/absorption reference)
//   https://en.wikipedia.org/wiki/Beer–Lambert_law
    float viewDist = (tr.hit) ? max(0.0, tr.t) : tLimit;
    float fogRaw = 1.0 - exp(-FOG_DENSITY * max(0.0, viewDist - FOG_START));
    float fogFact = clamp(fogRaw, 0.0, 0.97);
    vec3 fogColor = sky;

    finalColor = mix(finalColor, fogColor, fogFact);

    outp = vec4(finalColor, 1.0);
}
