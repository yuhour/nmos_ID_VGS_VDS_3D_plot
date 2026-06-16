# NMOS 器件 ID–VGS–VDS 三维特性曲面图

**NMOS Device ID–VGS–VDS 3D Characteristic Surface Plotter**
<img width="3088" height="1846" alt="5e9aedc21b234b3d0cead8ec19a5392b" src="https://github.com/user-attachments/assets/e448a317-334e-4ca7-91ef-4c0a0a00bc4a" />

---

基于 Level 1 MOSFET 模型（Shichman–Hodges），使用 Gismondi 平滑函数实现线性区到饱和区的 C¹ 连续过渡。生成交互式 3D 曲面图 + 2D 俯视热力图，支持鼠标悬停联动高亮和实时参数修改。

*Based on the Level 1 MOSFET model (Shichman–Hodges), with Gismondi smoothing for C¹-continuous transition between linear and saturation regions. Generates an interactive 3D surface + 2D top-down heatmap with hover-linked highlighting and real-time parameter modification.*

---

## 目录 / Table of Contents

- [功能特性 / Features](#功能特性--features)
- [安装 / Installation](#安装--installation)
- [使用方法 / Usage](#使用方法--usage)
  - [1. 命令行传参 / CLI Arguments](#1-命令行传参--cli-arguments)
  - [2. 交互式输入 / Interactive Input](#2-交互式输入--interactive-input)
  - [3. 对话中提供参数 / In-Session Parameters](#3-对话中提供参数--in-session-parameters)
- [参数说明 / Parameters](#参数说明--parameters)
- [交互功能 / Interactive Features](#交互功能--interactive-features)
- [物理模型 / Physical Model](#物理模型--physical-model)
- [输出示例 / Output Example](#输出示例--output-example)

---

## 功能特性 / Features

- **三维曲面图**：ID 作为 VGS 和 VDS 的函数，可旋转视角
- **二维热力图**：俯视视角，鼠标悬停显示 ID、gm、Ron 及工作区名称
- **工作区标注**：截止区、亚阈值区、深线性区、线性区、饱和区（按 `h` 键切换显隐）
- **C¹ 连续过渡**：Gismondi 平滑 min 函数消除分段模型的导数不连续
- **交互式修改**：关闭图形窗口后可修改参数并重绘，无需重启程序
- **查询功能**：给定两个已知量（如 VGS 和 ID），查找第三个量（VDS）
- **图像保存**：支持输出 PNG / JPG / PDF

---

## 安装 / Installation

**依赖 / Dependencies：**

```bash
pip install numpy matplotlib
```

| Package    | 最低版本 / Min Version |
|------------|------------------------|
| `numpy`    | ≥ 1.20                 |
| `matplotlib` | ≥ 3.4               |

---

## 使用方法 / Usage

### 1. 命令行传参 / CLI Arguments

一次性提供所有核心参数：

```bash
python nmos_3d.py --W_um 10 --L_um 1 --t_ox_nm 20 --N_sub 1e17 --phi_F 0.35
```

部分参数留空——缺失的核心参数将通过交互式提示补全：

```bash
python nmos_3d.py --W_um 20 --L_um 1.5
```

**完整参数列表：**

```bash
python nmos_3d.py \
  --W_um 10 \          # 栅宽 [μm]
  --L_um 1 \           # 沟道长度 [μm]
  --t_ox_nm 20 \       # 氧化层厚度 [nm]
  --N_sub 1e17 \       # 衬底掺杂浓度 [cm⁻³]
  --phi_F 0.35 \       # 费米势 [V]
  --mu_n 300 \         # 电子迁移率 [cm²/(V·s)]，可选
  --lam 0.02 \         # 沟道长度调制系数 [V⁻¹]，可选
  --VGS_max 3.3 \      # VGS 扫描上限 [V]，可选
  --VDS_max 6.0 \      # VDS 扫描上限 [V]，可选
  --n_pts 100 \        # 每轴采样点数，可选
  --save output.png \  # 保存图像，可选
  --no-plot            # 仅打印摘要不绘图
```

### 2. 交互式输入 / Interactive Input

直接运行，按提示逐一输入参数：

```bash
python nmos_3d.py
```

交互流程 / Interactive flow：

```
Enter 'd' for demo defaults [Enter=manual input]  →  输入 d 使用演示默认值
  W [μm] *required                                     →  逐一输入核心参数
  L [μm] *required
  ...
  Save to file [Enter to skip]                          →  可选保存路径
```

**演示默认值 / Demo defaults：**

| 参数 / Parameter | 值 / Value |
|------------------|------------|
| W                | 10 μm      |
| L                | 1 μm       |
| t_ox             | 20 nm      |
| N_sub            | 1×10¹⁷ cm⁻³ |
| φ_F              | 0.35 V     |

---

## 参数说明 / Parameters

### 核心器件参数 / Core Device Parameters（必须提供 / Required）

| 参数 / Parameter | 单位 / Unit | 说明 / Description |
|------------------|-------------|---------------------|
| `W_um`           | μm          | 栅宽 / Gate width |
| `L_um`           | μm          | 有效沟道长度 / Effective channel length |
| `t_ox_nm`        | nm          | 栅极氧化层厚度 / Gate oxide thickness |
| `C_ox`           | F/cm²       | 单位面积栅氧化层电容（留空则从 t_ox 推算） / Oxide capacitance per unit area (auto-calculated if blank) |
| `N_sub`          | cm⁻³        | 衬底掺杂浓度 / Substrate doping concentration |
| `phi_F`          | V           | 费米势 / Fermi potential |

### 辅助参数 / Auxiliary Parameters（可选，有默认值 / Optional, with defaults）

| 参数 / Parameter | 默认值 / Default | 说明 / Description |
|------------------|------------------|---------------------|
| `mu_n`           | 300              | 电子迁移率 [cm²/(V·s)] / Electron mobility |
| `VFB`            | auto             | 平带电压 [V]（留空则由功函数差推算） / Flat-band voltage (auto from work function difference) |
| `lam`            | 0.02             | 沟道长度调制系数 [V⁻¹] / Channel length modulation coefficient |
| `VGS_max`        | 3.3              | VGS 扫描上限 [V] / Max gate-source sweep voltage |
| `VDS_max`        | 6.0              | VDS 扫描上限 [V] / Max drain-source sweep voltage |
| `n_pts`          | 100              | 每轴采样点数 / Samples per axis |

---

## 交互功能 / Interactive Features

### 图形窗口 / Plot Window

| 操作 / Action | 说明 / Description |
|---------------|---------------------|
| 鼠标悬停 / Hover | 在 2D 热力图上移动，3D 图联动高亮，浮窗显示 VGS、VDS、ID、gm、Ron 及工作区名称 |
| 按 `h` 键 / Press `h` | 切换工作区标注和分割线的显示/隐藏 / Toggle region annotations |

### 关闭图形窗口后 / After Closing the Plot

进入修改/查询循环 / Enter the modify/query loop：

```
Modify: key=value (comma-sep)         |  修改：key=value（逗号分隔）
Query : ? key=value, key=value        |  查询：? key=value, key=value
Params: W, L, t_ox, N_sub, phi_F, C_ox, VGS, VDS, mu_n, lam, VFB, n_pts
e.g.: W=20, VDS=8        e.g.: ? VGS=1.5, ID=0.8
Enter=replot, q=quit                  |  留空重算，q 退出
```

**修改示例 / Modify example：**

```
> W=20, L=1.5
  W_um: 10.0 → 20.0
  L_um: 1.0 → 1.5
```

**查询示例 / Query example：**

```
> ? VGS=1.5, ID=0.8
  查询结果：VGS=1.500V  VDS=0.546V  ID=0.8012mA  [Saturation]
```

可用参数短名 / Short names：`W`, `L`, `t_ox` / `tox`, `VGS`, `VDS`

---

## 物理模型 / Physical Model

### 模型方程 / Model Equations

**Level 1 MOSFET（Shichman–Hodges）：**

$$I_D = \mu_n C_{ox} \frac{W}{L} \left[ (V_{GS} - V_{TH}) V_{DS,\text{eff}} - \frac{1}{2} V_{DS,\text{eff}}^2 \right] (1 + \lambda V_{DS})$$

**Gismondi 平滑 min 函数 / Gismondi Smooth Min：**

$$V_{DS,\text{eff}} = \frac{1}{2}\left(V_{DS} + V_{ov} - \sqrt{(V_{DS} - V_{ov})^2 + 4\delta^2}\right)$$

其中 $V_{ov} = V_{GS} - V_{TH}$，$\delta$ 为平滑参数。

- 当 $\delta \to 0$ 时，精确还原 Level 1 分段模型
- $\delta > 0$ 时，线性区 ↔ 饱和区 C¹ 连续

### 导出量 / Derived Quantities

| 量 / Quantity | 公式 / Formula |
|---------------|----------------|
| 阈值电压 / Threshold voltage | $V_{TH} = V_{FB} + 2\phi_F + \gamma\sqrt{2\phi_F}$ |
| 体效应系数 / Body effect coefficient | $\gamma = \frac{\sqrt{2q\varepsilon_{si} N_{sub}}}{C_{ox}}$ |
| 平带电压 / Flat-band voltage | $V_{FB} = -\frac{E_g}{2} - \phi_F$（n⁺ poly gate / p-Si） |
| 工艺跨导 / Process transconductance | $k'_n = \mu_n C_{ox}$ |

### 物理常数 / Physical Constants

| 符号 / Symbol | 值 / Value | 说明 / Description |
|---------------|------------|---------------------|
| $q$           | 1.602×10⁻¹⁹ C | 电子电荷 / Elementary charge |
| $k$           | 1.381×10⁻²³ J/K | 玻尔兹曼常数 / Boltzmann constant |
| $T$           | 300 K       | 温度 / Temperature |
| $\varepsilon_0$ | 8.854×10⁻¹⁴ F/cm | 真空介电常数 / Vacuum permittivity |
| $\varepsilon_{SiO_2}$ | 3.9 | SiO₂ 相对介电常数 |
| $\varepsilon_{Si}$ | 11.7 | 硅相对介电常数 |
| $E_g$         | 1.12 eV     | 硅禁带宽度 / Silicon bandgap |

---

## 输出示例 / Output Example

控制台摘要（demo 默认值）：

```
============================================================
  NMOS 器件物理参数摘要
============================================================
  栅宽 W          = 10.00 μm
  沟道长度 L      = 1.00 μm
  宽长比 W/L       = 100.00
  氧化层厚度 t_ox  = 20.00 nm
  C_ox             = 1.7265e-07 F/cm²
                   = 1.727 fF/μm²
  衬底掺杂 N_sub   = 1.00e+17 cm⁻³
  费米势 φ_F       = 0.350 V
  平带电压 V_FB    = -0.910 V
  体效应系数 γ     = 1.0545 V^(1/2)
  阈值电压 V_TH    = 0.673 V
  工艺跨导 k_n'    = 5.1796e-05 A/V²
  迁移率 μ_n       = 300 cm²/(V·s)
  沟长调制 λ       = 0.02 V⁻¹
------------------------------------------------------------
  VGS 扫描范围     = 0 – 3.3 V
  VDS 扫描范围     = 0 – 6.0 V
  采样点数         = 100 × 100
============================================================
## 图/ plot
<img width="892" height="576" alt="1d8e2817f57661ef10bf84a8bf7705e6" src="https://github.com/user-attachments/assets/117e7271-b912-495f-b25b-2b3189902688" />
<img width="1454" height="598" alt="53bb2b290d74a7f77abd53d543fb86aa" src="https://github.com/user-attachments/assets/01dc7dcf-bc0d-4c01-87e5-45bca8099f21" />

