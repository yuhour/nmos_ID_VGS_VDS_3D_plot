<p align="center">
  <img src="https://github.com/user-attachments/assets/e448a317-334e-4ca7-91ef-4c0a0a00bc4a" width="720" alt="NMOS 3D Surface Preview">
</p>

<h1 align="center">NMOS I<sub>D</sub>–V<sub>GS</sub>–V<sub>DS</sub> 3D Characteristic Surface</h1>

<p align="center">
  <strong>基于 Level 1 MOSFET 模型 · Gismondi 平滑 · C¹ 连续过渡 · 本征增益分析</strong>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/python-3.8+-blue" alt="Python">
  <img src="https://img.shields.io/badge/numpy-≥1.20-4dabf7" alt="NumPy">
  <img src="https://img.shields.io/badge/matplotlib-≥3.4-ff9800" alt="Matplotlib">
  <img src="https://img.shields.io/badge/license-MIT-green" alt="License">
</p>

<p align="center">
  <a href="https://www.bilibili.com/video/BV1BTjg6KE3z/">🎥 B站演示视频</a> &nbsp;·&nbsp;
  <a href="#-功能特性">✨ 功能</a> &nbsp;·&nbsp;
  <a href="#-快速开始">🚀 快速开始</a> &nbsp;·&nbsp;
  <a href="#-物理模型">📐 模型</a>
</p>

---

## 📖 概述

一个面向**模拟/数模混合集成电路设计初学者和工程师**的 NMOS 特性可视化工具。

输入器件物理参数（栅宽、栅长、氧化层厚度、衬底掺杂浓度、费米势等），程序自动计算阈值电压等导出量，生成**交互式 3D 曲面图 + 2D 俯视热力图**。

基于 Level 1 MOSFET 模型（Shichman–Hodges），使用 **Gismondi 平滑 min 函数**实现线性区→饱和区的 C¹ 连续过渡，消除传统分段模型在边界的导数不连续。

> **v3 新增**：饱和区**本征增益 A<sub>v</sub> = g<sub>m</sub> / g<sub>ds</sub>** 的完整计算与显示——悬停浮窗、公式弹窗、控制台摘要三处同步展示。

---

## ✨ 功能特性

| # | 功能 | 说明 |
|---|------|------|
| 🎯 | **三维曲面图** | I<sub>D</sub> = f(V<sub>GS</sub>, V<sub>DS</sub>)，可旋转/缩放视角 |
| 🔥 | **二维热力图** | 俯视 I<sub>D</sub> 热力图，色阶直观 |
| 🖱️ | **鼠标悬停联动** | 2D 图上移动→3D 图红色高亮→浮窗显示 V<sub>GS</sub>, V<sub>DS</sub>, I<sub>D</sub>, g<sub>m</sub>, g<sub>ds</sub>, **A<sub>v</sub>**, R<sub>on</sub>, 工作区 |
| 🖊️ | **点击公式弹窗** | 左键点击→弹出深色公式窗口，逐步展示代入数值和完整计算推导 |
| 📊 | **本征增益分析** | 饱和区自动计算 **A<sub>v</sub> = g<sub>m</sub> / g<sub>ds</sub>**（V/V + dB），含 g<sub>ds</sub> 完整推导 |
| 🏷️ | **工作区自动标注** | 截止区、亚阈值区、深线性区、线性区、饱和区 + 边界线（按 `h` 键切换） |
| 🔄 | **C¹ 连续过渡** | Gismondi 平滑 min 函数，导数在 V<sub>DS</sub> = V<sub>ov</sub> 处连续 |
| ♻️ | **交互式重算** | 关闭窗口后修改参数并重绘，无需重启 |
| 🔍 | **两点查询** | 给定 V<sub>GS</sub>+I<sub>D</sub> 或 V<sub>DS</sub>+I<sub>D</sub>，网格搜索第三量 |
| 💾 | **图像保存** | 支持 PNG / JPG / PDF |

---

## 🚀 快速开始

### 依赖

```bash
pip install numpy matplotlib
```

| Package | 最低版本 |
|---------|----------|
| `numpy` | ≥ 1.20 |
| `matplotlib` | ≥ 3.4 |

### 三种使用方式

#### ① 命令行传参（一键出图）

```bash
python nmos_3d.py --W_um 10 --L_um 1 --t_ox_nm 20 --N_sub 1e17 --phi_F 0.35
```

全部核心参数通过命令行指定，跳过交互提示。

#### ② 交互式输入（逐步引导）

```bash
python nmos_3d.py
```

按提示逐一输入参数。输入 `d` 可快速加载演示默认值（W=10μm, L=1μm, t<sub>ox</sub>=20nm, N<sub>sub</sub>=1×10¹⁷, φ<sub>F</sub>=0.35V）。

#### ③ 命令行 + 交互式混合

```bash
python nmos_3d.py --W_um 20 --L_um 1.5
```

先在命令行指定部分参数，缺失的核心参数自动进入交互式提示补全。

### 完整参数列表

```bash
python nmos_3d.py \
  --W_um 10       # 栅宽 [μm]
  --L_um 1        # 沟道长度 [μm]
  --t_ox_nm 20    # 氧化层厚度 [nm]
  --N_sub 1e17    # 衬底掺杂浓度 [cm⁻³]
  --phi_F 0.35    # 费米势 [V]
  --C_ox 1.7e-7   # 栅氧电容 [F/cm²]（可选，留空自动从 t_ox 推算）
  --mu_n 300      # 电子迁移率 [cm²/(V·s)]（默认 300）
  --VFB -0.9      # 平带电压 [V]（默认从功函数差推算）
  --lam 0.02      # 沟道长度调制系数 λ [V⁻¹]（默认 0.02）
  --VGS_max 3.3   # VGS 扫描上限 [V]（默认 3.3）
  --VDS_max 6.0   # VDS 扫描上限 [V]（默认 6.0）
  --n_pts 100     # 每轴采样点数（默认 100）
  --save out.png  # 保存图像到指定路径
  --no-plot       # 仅打印控制台摘要，不显示图形
  --no-regions    # 不显示工作区标注和分割线
```

---

## 🖱️ 交互操作

### 图形窗口内

| 操作 | 效果 |
|------|------|
| 鼠标在 2D 热力图上移动 | 3D 图红色圆球联动高亮，浮窗实时显示该点全部参数 |
| 左键点击 2D 热力图 | 弹出**深色公式窗口**，逐步展示 I<sub>D</sub> 计算推导（含 Gismondi 中间量） |
| 按 `h` 键 | 切换工作区标注线/文字的显示与隐藏 |

### 关闭图形窗口后

进入修改/查询循环：

```
修改: key=value, key=value    →  W=20, L=1.5, VDS=8
查询: ? key=value, key=value  →  ? VGS=1.5, ID=0.8
留空:  重算当前参数
q:     退出
```

---

## 📊 本征增益（v3 新增）

在饱和区（V<sub>DS</sub> ≥ V<sub>ov</sub>），程序自动计算并展示**本征增益**（Intrinsic Gain）：

$$A_v = \frac{g_m}{g_{ds}}$$

其中：

$$g_m = \frac{\partial I_D}{\partial V_{GS}} \approx \beta \cdot V_{ov} \cdot (1 + \lambda V_{DS})$$

$$g_{ds} = \frac{\partial I_D}{\partial V_{DS}} \approx \frac{1}{2} \lambda \beta V_{ov}^2$$

$$A_v \approx \frac{2(1 + \lambda V_{DS})}{\lambda V_{ov}} \quad \text{[V/V]}$$

**展示位置：**

| 位置 | 内容 |
|------|------|
| 控制台摘要 | 参考工作点的 g<sub>m</sub>、g<sub>ds</sub>、A<sub>v</sub>（V/V + dB） |
| 悬停浮窗 | 鼠标指向饱和区任意点时，显示 g<sub>m</sub>、g<sub>ds</sub>、A<sub>v</sub> |
| 公式弹窗 | 饱和区点击后，完整展示 g<sub>ds</sub> 推导和 A<sub>v</sub> 计算 |

---

## 📐 物理模型

### 核心方程

**Level 1 MOSFET（Shichman–Hodges）：**

$$I_D = \beta \left[ V_{ov} \cdot V_{DS,\text{eff}} - \frac{1}{2} V_{DS,\text{eff}}^2 \right] (1 + \lambda V_{DS})$$

其中 β = μ<sub>n</sub>·C<sub>ox</sub>·(W/L)，V<sub>ov</sub> = V<sub>GS</sub> − V<sub>TH</sub>。

**Gismondi 平滑 min 函数：**

$$V_{DS,\text{eff}} = \frac{1}{2}\left(V_{DS} + V_{ov} - \sqrt{(V_{DS} - V_{ov})^2 + 4\delta^2}\right)$$

- δ → 0：精确还原 Level 1 分段模型（`V_DS,eff = min(V_DS, V_ov)`）
- δ > 0：在 V<sub>DS</sub> ≈ V<sub>ov</sub> 附近 C¹ 连续过渡

### 导出量

| 物理量 | 公式 |
|--------|------|
| 阈值电压 | V<sub>TH</sub> = V<sub>FB</sub> + 2φ<sub>F</sub> + γ√(2φ<sub>F</sub>) |
| 体效应系数 | γ = √(2q·ε<sub>si</sub>·N<sub>sub</sub>) / C<sub>ox</sub> |
| 平带电压 | V<sub>FB</sub> = −E<sub>g</sub>/2 − φ<sub>F</sub>（n⁺ poly / p-Si） |
| 工艺跨导 | k'<sub>n</sub> = μ<sub>n</sub>·C<sub>ox</sub> |

### 物理常数

| 符号 | 值 | 说明 |
|------|-----|------|
| q | 1.602×10⁻¹⁹ C | 电子电荷 |
| k | 1.381×10⁻²³ J/K | 玻尔兹曼常数 |
| T | 300 K | 温度 |
| ε₀ | 8.854×10⁻¹⁴ F/cm | 真空介电常数 |
| ε<sub>SiO₂</sub> | 3.9 | SiO₂ 相对介电常数 |
| ε<sub>Si</sub> | 11.7 | 硅相对介电常数 |
| E<sub>g</sub> | 1.12 eV | 硅禁带宽度 |

---

## 🧪 输出示例

### 控制台摘要

```
============================================================
  NMOS 器件物理参数摘要
============================================================
  栅宽 W          = 10.00 μm
  沟道长度 L      = 1.00 μm
  宽长比 W/L       = 10.00
  氧化层厚度 t_ox  = 20.00 nm
  C_ox             = 1.7265e-07 F/cm²
                   = 1.727 fF/μm²
  衬底掺杂 N_sub   = 1.00e+17 cm⁻³
  费米势 φ_F       = 0.350 V
  平带电压 V_FB    = -0.910 V
  体效应系数 γ     = 1.0552 V^(1/2)
  阈值电压 V_TH    = 0.673 V
  工艺跨导 k_n'    = 5.1796e-05 A/V²
  迁移率 μ_n       = 300.0 cm²/(V·s)
  沟长调制 λ       = 0.02 V⁻¹
------------------------------------------------------------
  VGS 扫描范围     = 0 – 3.3 V
  VDS 扫描范围     = 0 – 6.0 V
  采样点数         = 100 × 100
============================================================

  参考工作点 (VGS=2.17V, VDS=1.50V, 饱和):
    I_D ≈ 0.600 mA  (600.2 μA)
    g_m ≈ 0.800 mA/V
    g_ds ≈ 1.1654e-05 S  (11.65 μS)
    本征增益 A_v = g_m/g_ds ≈ 68.7 V/V  (36.7 dB)
```

---

## 📂 项目结构

```
NMOS_ID_VGS_VDS_3D_plot/
├── nmos_3d.py          # 主程序（当前版本 v3）
├── README.md           # 项目文档
└── 旧版本/             # 历史版本备份
    ├── v1.py           # 原始 v1 版本
    └── nmos_3d_v2_backup.py  # v2 版本（本征增益功能添加前）
```

---

## 📋 版本历史

| 版本 | 日期 | 更新 |
|------|------|------|
| **v3** | 2026.06 | 饱和区**本征增益** A<sub>v</sub> = g<sub>m</sub>/g<sub>ds</sub> 计算与显示（悬停/弹窗/摘要三处） |
| **v2** | 2026.06 | 点击公式弹窗（逐步代入计算）、点击蓝色标记、深线性区 R<sub>on</sub> 公式 |
| **v1** | 2026.06 | 基础 3D 曲面 + 2D 热力图、鼠标悬停联动、工作区标注、交互式重算、查询功能 |

---

## 📄 License

MIT

---

<p align="center">
  <sub>适合模拟集成电路设计初学者和工程师 · For analog/mixed-signal IC beginners and engineers</sub>
</p>
