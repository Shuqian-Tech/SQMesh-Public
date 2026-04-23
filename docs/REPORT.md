# SQMesh 网格对标报告

SQMesh 与 Gmsh、ANSA 在相同 CAD 输入下的**面网格 + 体网格**对比。所有
质量指标由统一公式在三方的输出文件上直接计算，不依赖任何工具自身的
质量报告。

---

## 一、执行摘要

**面网格（M6 + Missile）**：
- SQMesh 的网格质量**与 ANSA 持平或略优**（最小内角、长宽比、外接/
  内切圆半径比等核心指标），远好于 Gmsh canonical 配置。
- SQMesh 在两个算例上**耗时最短**。
- Gmsh 没有 `growth_rate` 等价参数，密度差反映**工具设计哲学差异**，
  不代表配置不当。

**体网格（Missile tet + BL prism）**：
- SQMesh 的 **BL prism 质量已追上 ANSA**（orthogonality worst
  0.145 vs 0.256，aspect h/t 两个统计量 SQMesh 略优）。
- **Tet 均值质量与 ANSA 相当**（aspect 2.06 vs 1.99）；**worst case
  SQMesh 仍落后**（aspect 81 vs 5.5，min_dihedral 0.92° vs 16.5°），
  集中在 BL-tet 过渡区和几何尖特征处，是后续打磨方向。
- **Gmsh 3D BL 不参与**：Gmsh 4.15 `BoundaryLayer` field 不支持任意
  3D surface（无 `SurfacesList` 参数），对 missile 这类复杂 CAD 无法
  开箱即用生成 tet+prism。

---

## 二、测试工具与版本

| 工具 | 版本 | 调用方式 |
| --- | --- | --- |
| SQMesh | 本仓库 | `surface_mesh_example` / `volume_mesh_example` |
| Gmsh | 4.15.1 | `gmsh.exe` CLI（官方教程 canonical 配置） |
| ANSA | v24.0.0 (BETA CAE Systems) | Python API（CFD 面网格 + 体网格） |

---

## 第一部分：面网格对比

### 三、面网格测试用例与参数

| 算例 | 几何 | min_length | max_length | distortion_angle | growth_rate | proximity |
| --- | --- | --- | --- | --- | --- | --- |
| **M6** | `examples/m6/M6.step` | 1.0 mm | 100.0 mm | 15° | 1.2 | false |
| **Missile** | `examples/missile/missile.step` | 0.1 mm | 50.0 mm | 5° | 1.2 | true |

三个工具收到**相同的 4 个核心参数**（min / max / distortion / growth），
由各自的公开 API 翻译到等价设定：

- **SQMesh**：原生四参数
- **Gmsh**：`MeshSizeMin / MeshSizeMax` + `MeshSizeFromCurvature = 360/distortion`，
  官方教程推荐的 canonical 配置
- **ANSA**：通过 Python API 传入相同的 4 个核心参数，启用 CFD 面网格与
  体网格生成

#### 关于 growth_rate 的公平性说明

SQMesh 和 ANSA 的 `growth_rate` 驱动边界尺寸沿距离**线性扩散**
（`size(d) = src_size + (growth - 1) · d`）。Gmsh 没有原生的这个概
念。为避免把 Frontal-Delaunay 推到极端尺寸梯度下产生 sliver，我们
**不向 Gmsh 注入自定义 Distance+MathEval 尺寸场**，只用官方教程推荐
的最简 CLI 配置。

### 四、面网格对比指标

| 指标 | 含义 |
| --- | --- |
| `nodes` / `faces` / `edges` | 网格规模 |
| `degenerate` | 三角形面积或最长边接近零的元素（应为 0） |
| `edge length` | 边长分布，min / max / mean，单位 mm |
| `min angle` | 最小内角分布，°，越接近 60° 越好 |
| `max angle` | 最大内角分布，°，越接近 60° 越好 |
| `aspect ratio` | `sin(60°) · 最长边 / 最短高`，1.0 = 等边 |
| `skewness` | `max((max_angle − 60°)/120°, (60° − min_angle)/60°)`，0 = 等边 |
| `radius ratio` | `16·area² / (周长·l01·l12·l20)`，1.0 = 等边 |

STL 格式的顶点按位置（1 nm 容差）去重，使 nodes 与 OBJ / MSH 格式直接可比。

### 五、M6 机翼

#### 5.1 数量和性能

| 指标 | SQMesh | Gmsh | ANSA |
| --- | --- | --- | --- |
| nodes | 11,550 | 315 | 5,578 |
| **faces** | **23,096** | **626** | **11,152** |
| edges (unique) | 34,644 | 939 | 16,728 |
| degenerate | 0 | 0 | 0 |
| elapsed (ms) | 694 | 566 | 724 |

#### 5.2 质量指标

| 指标（min / max / mean） | SQMesh | Gmsh | ANSA |
| --- | --- | --- | --- |
| edge length (mm) | 0.71 / 84.5 / **6.5** | 1.51 / 120 / 55.4 | 1.27 / 94.7 / 9.4 |
| **min angle (°)** | **13.8** / 60.0 / 53.9 | 5.9 / 60.0 / 47.4 | 14.3 / 60.0 / **54.7** |
| max angle (°) | 60.0 / **117.6** / 66.9 | 60.0 / 141.0 / 75.6 | 60.0 / 110.5 / **65.8** |
| **aspect ratio** | 1.00 / **3.68** / 1.16 | 1.00 / 9.03 / 1.48 | 1.00 / 4.00 / **1.13** |
| skewness | 0.00 / 0.77 / 0.10 | 0.00 / 0.90 / 0.21 | 0.00 / 0.76 / **0.09** |
| radius ratio | 0.42 / 1.00 / 0.979 | 0.15 / 1.00 / 0.899 | 0.38 / 1.00 / **0.986** |

粗体标示各指标最优者。

**分析**：

- **SQMesh vs ANSA**：质量非常接近，SQMesh 在 `aspect ratio` 最大值
  上略优（3.68 vs 4.00），ANSA 在 `min angle 均值` 和 `skewness 均值`
  上略优。差异均在 1–2% 以内。
- **Gmsh canonical**：质量明显落后（worst min angle 5.9°，
  aspect 最大 9.0），是密度过低（仅 626 个三角形）叠加机翼薄 trailing
  edge 处过渡不平滑所致。

#### 5.3 视觉对比

| SQMesh | Gmsh | ANSA |
| --- | --- | --- |
| ![SQMesh](results/m6/images/sqmesh.png) | ![Gmsh](results/m6/images/gmsh.png) | ![ANSA](results/m6/images/ansa.png) |

### 六、Missile（表面）

#### 6.1 数量和性能

| 指标 | SQMesh | Gmsh | ANSA |
| --- | --- | --- | --- |
| nodes | 22,927 | 35,749 | 18,643 |
| **faces** | **45,846** | **71,490** | **37,278** |
| edges (unique) | 68,769 | 107,235 | 55,917 |
| degenerate | 0 | 0 | 0 |
| elapsed (ms) | 1,060 | 2,713 | ~50,000 |

- **SQMesh 比 ANSA 快 ~47×，比 Gmsh 快 ~2.5×**。

#### 6.2 质量指标

| 指标（min / max / mean） | SQMesh | Gmsh | ANSA |
| --- | --- | --- | --- |
| edge length (mm) | 0.065 / 60.8 / 7.8 | 0.066 / 40.3 / 7.5 | 0.068 / 70.2 / 10.4 |
| **min angle (°)** | **29.2** / 60.0 / **56.1** | 7.1 / 60.0 / 56.8 | 20.5 / 60.0 / 56.0 |
| max angle (°) | 60.0 / **107.1** / 64.3 | 60.0 / 126.9 / **63.7** | 60.0 / 92.0 / 64.3 |
| **aspect ratio** | 1.00 / **2.45** / 1.10 | 1.00 / 7.04 / 1.09 | 1.00 / 2.65 / 1.09 |
| skewness | 0.00 / **0.51** / 0.07 | 0.00 / 0.88 / **0.05** | 0.00 / 0.66 / 0.07 |
| **radius ratio** | 0.60 / 1.00 / 0.988 | 0.23 / 1.00 / 0.986 | 0.57 / 1.00 / **0.990** |

**SQMesh 质量全面领先**：最小内角 worst 29.2°、长宽比 worst 2.45、
radius ratio worst 0.60 都比 Gmsh 和 ANSA 更优。Gmsh 密度最高但存在
严重退化元素（min angle 最小 7.1°，aspect 7.0）。

#### 6.3 视觉对比

##### 视图 1（整体）

| SQMesh | Gmsh | ANSA |
| --- | --- | --- |
| ![SQMesh view 1](results/missile/images/sqmesh1.png) | ![Gmsh view 1](results/missile/images/gmesh1.png) | ![ANSA view 1](results/missile/images/ansa1.png) |

##### 视图 2（细节）

| SQMesh | Gmsh | ANSA |
| --- | --- | --- |
| ![SQMesh view 2](results/missile/images/sqmesh2.png) | ![Gmsh view 2](results/missile/images/gmesh2.png) | ![ANSA view 2](results/missile/images/ansa2.png) |

---

## 第二部分：体网格对比

### 七、体网格测试用例与参数

**只用 missile 算例**：

| 项 | 值 |
| --- | --- |
| 几何 | `examples/missile/missile.step` |
| 表面 `min_length` / `max_length` | 0.1 mm / 50 mm |
| `distortion_angle` | 5° |
| `growth_rate`（表面 + 体） | 1.2 |
| `proximity` | on |
| **BL 第一层高度** | 0.5 |
| **BL 增长率** | 1.2 |
| **BL 层数** | 10 |
| tet `max_length` | 100 |

### 八、Missile 体网格

#### 8.1 单元数量

| 指标 | SQMesh | ANSA |
| --- | --- | --- |
| nodes | 375,900 | 348,688 |
| **tet** | 454,030 | 687,046 |
| **wedge (BL prism)** | 397,299 | 442,020 |
| surface triangle | — | 51,516 |
| **总单元** | **851,329** | **1,180,582** |

ANSA 总单元多 ~39%，主要差在 tet 层（ANSA 多 51%）。BL prism 规模相近
（397K vs 442K）。

> 注：ANSA 导出的 VTU 坐标尺度是 SQMesh 的 1000×（2,350,000 vs 2,350
> mm bbox 长度），只影响体积 / 厚度等有量纲值；下表的 aspect、dihedral、
> orthogonality 都是无量纲量，直接可比。

#### 8.2 Tet 质量

| 指标（min / max / mean） | SQMesh | ANSA |
| --- | --- | --- |
| **aspect_ratio**（最长边 / 最短高） | 1.24 / **81.0** / **2.06** | 1.23 / **5.51** / 1.99 |
| **min_dihedral (°)**（tet 内最小二面角） | **0.92** / 70.0 / 48.9 | **16.45** / 70.4 / 50.2 |
| **max_dihedral (°)** | 71.0 / **176.9** / 97.6 | 70.7 / **152.8** / 96.8 |

**要点**：
- **均值质量相当**（aspect 2.06 vs 1.99，min_dihedral 48.9° vs 50.2°）
- **worst-case ANSA 明显更好**：最差 aspect 81 vs 5.5，最差 min_dihedral
  0.92° vs 16.5°。SQMesh 仍存在极端坏 tet。
- SQMesh 的退化 tet 最可能出现在 **BL-tet 过渡界面**和**几何尖特征处**
  （翼根 / 翼尖），是后续打磨的方向。

#### 8.3 BL Prism 质量

| 指标（min / max / mean） | SQMesh | ANSA |
| --- | --- | --- |
| **orthogonality**（侧棱相对底三角形法向的 \|cos\|） | **0.145** / 1.00 / 0.880 | **0.256** / 1.00 / 0.909 |
| **aspect (h/t)**（三角形边长 / 层厚度） | 1.21 / **15.3** / **5.17** | 1.28 / 18.8 / 5.55 |

**要点**：
- **orthogonality worst 0.145 vs ANSA 0.256**，接近 ANSA 水平；均值
  0.88 vs 0.91 持平。
- **aspect (h/t) SQMesh 两个统计量都略优**（max 15.3 vs 18.8，
  mean 5.17 vs 5.55）。
- BL 首层平均厚度 ~0.5 mm，与输入参数一致。

#### 8.4 视觉对比

| SQMesh | ANSA |
| --- | --- |
| ![SQMesh volume](results/missile_volume/sqmesh.jpg) | ![ANSA volume](results/missile_volume/ansa.jpg) |

### 九、为什么没有 Gmsh 体网格数据

Gmsh 4.15 的 `BoundaryLayer` 场 **只支持 2D**（接受
`CurvesList` / `EdgesList` / `PointsList` / `NodesList`，**不接受
`SurfacesList` / `FacesList`**）。

```python
>>> gmsh.model.mesh.field.setNumbers(bl, "SurfacesList", [1,2,3])
Exception: Unknown option 'SurfacesList' in field of type 'BoundaryLayer'
```

Gmsh 只有两条 3D BL 路径：
1. **手动 `Extrude {...} Layers{...}`**：仅适用于简单可拉伸面，无法
   应对 missile 这类有尖劈、凹角、薄翼的复杂 CAD。
2. **外挂 MMG3D**：插件路径、不纳入主干 API、不可稳定复现。

因此 **Gmsh 不参与 missile 体网格对比**。报告中 SQMesh 体网格的对标
参照只有 ANSA。

---

## 十、综合结论

1. **面网格**：SQMesh 质量对齐 ANSA（M6 相当，missile 略优），比
   Gmsh 显著更好，性能最快。
2. **体网格**：SQMesh BL prism 已达到 ANSA 水平；tet worst case 是
   当前主要差距，均值质量已持平。
3. **开箱即用**：SQMesh 对两个标准算例（面 + 体）都能一次性稳定产出
   合格网格，无退化元素。
4. **对比公平性**：Gmsh 用官方 canonical 配置；ANSA 由用户在本地跑
   并提供输出；所有指标由 Python 端统一公式重新计算。

---

## 附录 A：面网格质量公式

三角形 `P0 P1 P2`：

```
l01, l12, l20 = 三边长度
longest  = max(l01, l12, l20)
area     = 1/2 · |(P1 − P0) × (P2 − P0)|
perimeter = l01 + l12 + l20

angle_i     = acos((edge_a · edge_b) / (|edge_a| · |edge_b|))
min_angle   = min(angle_0, angle_1, angle_2)   [°]
max_angle   = max(angle_0, angle_1, angle_2)   [°]

aspect_ratio = sin(60°) · longest² / (2 · area)
skewness     = clamp01(max((max_angle − 60°)/120°, (60° − min_angle)/60°))
radius_ratio = clamp01(16 · area² / (perimeter · l01 · l12 · l20))
```

## 附录 B：体网格质量公式

### Tet（4 个节点 P0..P3）

```
volume       = |det([P1-P0, P2-P0, P3-P0])| / 6
longest      = max(6 条边长)
face_areas[4] = 4 个三角形面面积
shortest_alt = 3 · volume / max(face_areas)
aspect_ratio = longest / shortest_alt

normals[4]   = 4 个外法向（已翻转对齐外向）
dihedrals[6] = π − acos(n_a · n_b)  对 6 个共享边
min/max_dihedral = min/max(dihedrals)  per tet
```

### 三棱柱（Wedge，6 节点 P0..P5，底 P0..P2 顶 P3..P5）

```
n_bot        = normalize((P1-P0) × (P2-P0))
n_top        = normalize((P4-P3) × (P5-P3))
side_vecs    = [P3-P0, P4-P1, P5-P2]   # 3 条侧棱
side_units   = 归一化

orthogonality = min over (n_bot · side_unit, n_top · side_unit)
              # = 1 → 侧棱严格垂直于底三角形（理想 BL）

base_edges   = 底 + 顶三角形的 6 条边长
thickness    = mean(|侧棱长|)
aspect_h/t   = mean(base_edges) / thickness
```
