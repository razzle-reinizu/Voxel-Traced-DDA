# GLSL Voxel World Shader — Concise Technical README

## Summary

A fragment shader that renders a procedural voxel world using per-fragment **Amanatides–Woo DDA** traversal, hash/value-noise terrain, deterministic tree placement, simple material sampling, and an approximate atmosphere model.

## Features (at-a-glance)

* Voxel grid ray traversal (DDA)
* Heightmap-driven terrain + single-octave FBM
* Deterministic tree trunks & leaf volumes
* Block material mapping to texture samplers
* Single-sample hard shadows via re-tracing
* Atmospheric scattering approximation and fog
* Configurable quality modes (QUALITY 0..3) and compile-time feature flags

## Main rendering flow

1. Build camera (position, forward, right, up) from `time` and `touch`.
2. Compose pinhole ray (FOV and aspect-corrected NDC).
3. Trace voxels with `traceVoxel(rayOrigin, rayDir, tLimit, maxSteps)`.
4. If hit: compute local hit coords, sample albedo, evaluate lighting (direct + ambient + emissive + GGX-like specular), optionally test shadow by re-tracing toward sun.
5. If miss: compute sky color via `atmosphere` or `absorp`.
6. Apply exponential tone mapping and Beer–Lambert fog; output final color.

## Voxel traversal (DDA)

* `init_voxel_dda` computes `cellCoord`, `cellSign`, `tMax`, `tDelta`.
* `traceVoxel` advances along smallest `tMax` axis, updates `cell` and face normal, queries `getBlockType(cell)`, and returns first non-empty cell or runs out of steps.
* Immediate starting-voxel hit returns `t=0` with normal set to `(0,1,0)`.

## World rules (`getBlockType`)

* QUALITY==0: stochastic sparse blocks using `hash13` and density threshold; limited vertical range.
* QUALITY>0: compute `terrainH = floor(proceduralTerrainHeight(x,z))`; place grass/dirt/stone for `y <= terrainH` unless cave noise indicates void.
* Trees: `hasTreeAtColumn(cx,cz,q)` uses `hash12` thresholds; trunk and leaves generated deterministically with quality-dependent radii.
* Special IDs: `BLOCK_EMISSIVE`, `BLOCK_CLOUD`.

## Noise & terrain

* `hash12/hash13`: fast pseudo-random hash (fract(sin(dot(...))).
* `noise3`: value-noise on unit lattice with Hermite smoothing and trilinear interpolation.
* `fbm3`: `0.5 * noise3` (single-octave wrapper).
* `proceduralTerrainHeight` adjusts frequency/scale by `QUALITY`.

## Materials & UVs

* `faceUV` selects planar UV per dominant normal axis.
* `sampleBlockAlbedoTex_manual` maps block ID → sampler lookup:

  * 1: grass top/side mixing, 2: dirt, 3: stone, 4: wood, 5: leaves, emissives and cloud color special-cased.
* `sampleLeaves` applies per-cell tint and mixes a fixed grass-side texel for consistent tinting.

## Atmosphere & absorption

* `raysi` computes ray–sphere intersections used by `atmosphere` and `absorp`.
* `atmosphere` performs a condensed Rayleigh/Mie accumulation (loop removed; single-step integration) and returns scattering contribution.
* `absorp` computes a scaled Rayleigh absorption fallback.

## Lighting & shading

* Direct: Lambert (`max(dot(N,L),0)`), `lightColor`, ambient floor.
* Shadow: single hard occlusion via `traceVoxel` from hit toward `L`.
* Emissive blocks add emission to direct/indirect terms.
* Specular: GGX-like terms with `H = normalize(L - V)` (nonstandard half-vector) and Schlick F approximation.
* Tone mapping: `final = 1 - exp(-2 * lit)`.

## Camera & input

* `camera_position_from_time(t)`: forward motion on Z + horizontal sway; Y fixed (`CAMERA_HEIGHT`).
* `camera_forward_vector(t)`: touch → yaw/pitch with deadzone, non-linear shaping, and time-based wobble; returns `-fwd` normalized.

## Quality & flags

* `QUALITY` influences noise frequency, vertical scale, tree thresholds, `SEARCH_MAX`, and `maxSteps`.
* Compile-time flags: `ENABLE_TREES`, `ENABLE_CAVES`, `ENABLE_CLOUDS`, `ENABLE_ATMOSPHERE`, `ENABLE_SHADOW`, `BLOCK_LIGHT`, `ENDLESS_LOWER_HEIGHT`, `USE_TOUCH`.

## Limitations (code-observable)

* `atlasLod` fixed to `0.0` (no mip selection).
* `fbm3` is single-octave; caves use single-layer FBM.
* Atmosphere is an approximation (loop removed) not full multi-sample integration.
* `H = normalize(L - V)` differs from standard half-vector `(L+V)`.
* Start-voxel normal fallback `(0,1,0)` may be inaccurate for interior origins.

## Uniforms (required)

* `resolution: vec2`, `touch: vec2`, `time: float`, samplers: `grass_top, oak_leaves, stone, grass_side, dirt, oak_log_side, oak_log_top`.

## References (in-code)

* https://www.cs.yorku.ca/~amana/research/grid.pdf — Amanatides & Woo — fast voxel traversal 
* https://scratchapixel.com/lessons/3d-basic-rendering/perspective-and-camera — camera/ray generation
* https://github.com/wwwtyro/glsl-atmosphere/blob/master/index.glsl — referenced for atmospheric implementation
* https://github.com/glslify/glsl-ggx/blob/master/index.glsl GGX implemented for world shading


![1000374452](https://github.com/user-attachments/assets/cdb7a7e6-4672-422b-809d-9927d176123d)
![1000374453](https://github.com/user-attachments/assets/21cd91d4-0965-483b-b96f-9c1dcd046913)
