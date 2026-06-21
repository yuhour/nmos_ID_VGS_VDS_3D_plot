#!/usr/bin/env python3
"""
NMOS 器件 ID–VGS–VDS 三维特性曲面图 (v2)
===========================================

基于 Level 1 MOSFET 模型 (Shichman–Hodges)。所有器件物理参数由用户提供，
代码内不预存默认值。

v2 新增：
  - 3D 视图内显示鼠标悬停点对应的计算公式（含代入数值）
  - 深线性区额外显示 R_on 计算公式
  - 各工作区显示对应的经典公式及 Gismondi 平滑中间量

必须提供的 6 个器件参数：
  - W      : 栅宽 [μm]
  - L      : 有效沟道长度 [μm]
  - t_ox   : 栅极氧化层厚度 [nm]
  - C_ox   : 单位面积栅氧化层电容 [F/cm²]（可留空，由 t_ox 推算）
  - N_sub  : 衬底掺杂浓度 [cm⁻³]
  - phi_F  : 费米势 [V]

三种使用方式（按优先级）：
  1. 命令行传参（适合自动化）：
     python v2.py --W_um 10 --L_um 1 --t_ox_nm 20 --N_sub 1e16 --phi_F 0.35
  2. 交互式输入（运行后按提示逐一输入）：
     python v2.py
  3. 对话中提供参数（在 CodeWhale 会话中直接告诉我数值）

依赖：numpy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker
import argparse
import sys
import os

# Windows GBK → UTF-8 兼容
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

# ============================================================================
# 物理常数 (CGS 单位制)
# ============================================================================
Q        = 1.602e-19      # 电子电荷 [C]
K        = 1.381e-23      # 玻尔兹曼常数 [J/K]
T        = 300.0          # 温度 [K]
EPS0     = 8.854e-14      # 真空介电常数 [F/cm]
EPS_SIO2 = 3.9            # SiO₂ 相对介电常数
EPS_SI   = 11.7           # 硅相对介电常数
EG       = 1.12           # 硅禁带宽度 [eV]

# ============================================================================
# 辅助参数的默认值（非器件物理参数，可覆盖）
# ============================================================================
AUX_DEFAULTS = {
    "mu_n":     300.0,    # 电子迁移率 [cm²/(V·s)]
    "VFB":      None,     # 平带电压 [V]；None → 由功函数差推算
    "lam":      0.02,     # 沟道长度调制系数 [V⁻¹]
    "VGS_max":  3.3,      # VGS 扫描上限 [V]
    "VDS_max":  6.0,      # VDS 扫描上限 [V]
    "n_pts":    100,      # 每轴采样点数
}

# ---- 演示默认值（快捷测试用） ----
DEMO_DEFAULTS = {
    "W_um":     10.0,
    "L_um":      1.0,
    "t_ox_nm":  20.0,
    "C_ox":     None,       # None → 从 t_ox 推算
    "N_sub":    1e17,
    "phi_F":    0.35,
    "mu_n":     300.0,
    "VFB":      None,       # None → 由功函数差推算
    "lam":      0.02,
    "VGS_max":   3.3,
    "VDS_max":   6.0,
    "n_pts":    100,
}

# ---- 核心器件参数的提示文本 ----
PARAM_PROMPTS = {
    "W_um":    ("栅宽 W",                            "μm"),
    "L_um":    ("有效沟道长度 L",                     "μm"),
    "t_ox_nm": ("栅极氧化层厚度 t_ox",                "nm"),
    "N_sub":   ("衬底掺杂浓度 N_sub",                  "cm⁻³"),
    "phi_F":   ("费米势 φ_F",                         "V"),
    "C_ox":    ("单位面积栅氧化层电容 C_ox（留空则从 t_ox 推算）", "F/cm²"),
}


# ============================================================================
# 交互式输入
# ============================================================================

def prompt_param(label, unit):
    """提示用户输入一个浮点数参数，留空返回 None。"""
    prompt_text = f"  {label} [{unit}]: "
    while True:
        try:
            raw = input(prompt_text).strip()
            if raw == "":
                return None
            return float(raw)
        except ValueError:
            print(f"    [!] 无法解析 '{raw}'，请输入数值（如 10、1.5、1e16）")
        except (EOFError, KeyboardInterrupt):
            print("\n[中断] 用户取消输入。")
            sys.exit(1)


def prompt_required_param(label, unit):
    """提示用户输入一个必须提供的浮点数参数，留空时重新提示。"""
    prompt_text = f"  {label} [{unit}] *required | *必须: "
    while True:
        try:
            raw = input(prompt_text).strip()
            if raw == "":
                print(f"    [!] 此参数必须提供，请输入数值")
                continue
            return float(raw)
        except ValueError:
            print(f"    [!] 无法解析 '{raw}'，请输入数值")
        except (EOFError, KeyboardInterrupt):
            print("\n[中断] 用户取消输入。")
            sys.exit(1)


def prompt_optional_param(label, unit, default_val):
    """
    提示用户输入一个浮点数参数，直接回车使用默认值。
    返回 float 或 None。
    """
    prompt = f"  {label} [{unit}] (default {default_val}): "
    while True:
        try:
            raw = input(prompt).strip()
            if raw == "":
                return default_val
            return float(raw)
        except ValueError:
            print(f"    [!] 无法解析 '{raw}'，请输入数值或直接回车用默认值")
        except (EOFError, KeyboardInterrupt):
            print("\n[中断] 用户取消输入。")
            sys.exit(1)


def collect_params_interactive(args):
    """
    交互式收集所有必要参数。
    命令行已有的参数跳过询问；缺少的核心参数逐一提示。
    返回 (params, save_path) 元组。
    """
    print("\n" + "=" * 60)
    print("  NMOS Device Parameter Input  |  NMOS 器件参数输入")
    print("  * Required — cannot be empty  |  * 标记为必须参数，不可留空")
    print("  Press Enter to skip optional  |  直接回车跳过可选参数（使用默认值）")
    print("=" * 60)

    # ---- 快捷选项：使用演示默认值 ----
    quick = input("  Enter 'd' for demo defaults [Enter=manual input]  |  输入 'd' 使用演示默认值 [留空=逐一输入]: ").strip().lower()
    if quick in ('d', 'default', 'y', 'yes'):
        params = dict(DEMO_DEFAULTS)  # 浅拷贝
        # 命令行已有的参数仍以命令行优先
        for key in DEMO_DEFAULTS:
            cli_val = getattr(args, key, None)
            if cli_val is not None:
                params[key] = cli_val
        print("  ✓ Demo defaults loaded (W=10μm, L=1μm, t_ox=20nm, N_sub=1e17, φ_F=0.35V)")
        print()

        # 仍询问保存路径
        save_input = input("  Save to file [Enter to skip, .png/.jpg/.pdf]  |  保存到文件 [留空跳过]: ").strip()
        save_path = save_input if save_input != "" else None
        print("-" * 55)
        return params, save_path

    params = {}

    # ---- 必需的核心器件参数（不可留空） ----
    required_keys = ["W_um", "L_um", "t_ox_nm", "N_sub", "phi_F"]
    for key in required_keys:
        cli_val = getattr(args, key, None)
        if cli_val is not None:
            params[key] = cli_val
            label, unit = PARAM_PROMPTS[key]
            print(f"  {label} [{unit}]: {cli_val}  (from CLI | 来自命令行)")
        else:
            label, unit = PARAM_PROMPTS[key]
            params[key] = prompt_required_param(label, unit)

    # ---- 可选的核心参数 C_ox（留空则从 t_ox 推算） ----
    cli_val = getattr(args, "C_ox", None)
    if cli_val is not None:
        params["C_ox"] = cli_val
        label, unit = PARAM_PROMPTS["C_ox"]
        print(f"  {label} [{unit}]: {cli_val}  (来自命令行)")
    else:
        label, unit = PARAM_PROMPTS["C_ox"]
        val = prompt_param(label, unit)
        params["C_ox"] = val  # None 表示从 t_ox 推算

    # ---- 辅助参数（有默认值） ----
    aux_keys = ["mu_n", "VFB", "lam", "VGS_max", "VDS_max", "n_pts"]
    for key in aux_keys:
        cli_val = getattr(args, key, None)
        default = AUX_DEFAULTS[key]
        if cli_val is not None:
            params[key] = cli_val
        else:
            if isinstance(default, float):
                params[key] = prompt_optional_param(
                    f"{key}（辅助参数）",
                    "—" if key in ("n_pts",) else
                    ("V" if key in ("VFB", "VGS_max", "VDS_max") else
                     "cm²/(V·s)" if key == "mu_n" else "V⁻¹"),
                    default,
                )
            else:
                # default 为 None 的情况（如 VFB）
                params[key] = default

    # ---- 询问保存路径 ----
    print()
    save_input = input("  Save to file [Enter to skip, .png/.jpg/.pdf]  |  保存到文件 [留空跳过]: ").strip()
    save_path = save_input if save_input != "" else None

    print("-" * 60)
    return params, save_path


# ============================================================================
# 物理量计算
# ============================================================================

def compute_derived(params):
    """
    由原始器件参数计算 C_ox、VFB、γ、VTH 等导出量。
    返回 derived 字典。
    """
    p = dict(params)

    # 单位转换
    W_cm    = p["W_um"]    * 1e-4
    L_cm    = p["L_um"]    * 1e-4
    t_ox_cm = p["t_ox_nm"] * 1e-7

    # C_ox：用户提供则用，否则从 t_ox 推算
    if p.get("C_ox") is not None:
        C_ox = p["C_ox"]
    else:
        C_ox = EPS0 * EPS_SIO2 / t_ox_cm

    # VFB：用户提供则用，否则按 n⁺-poly / p-Si 功函数差推算
    if p.get("VFB") is not None:
        VFB = p["VFB"]
    else:
        VFB = -(EG / 2.0) - p["phi_F"]       # ≈ -(0.56 + φ_F) V

    # 体效应系数 γ
    eps_si = EPS0 * EPS_SI
    gamma = np.sqrt(2.0 * Q * eps_si * p["N_sub"]) / C_ox

    # 阈值电压
    VTH = VFB + 2.0 * p["phi_F"] + gamma * np.sqrt(2.0 * p["phi_F"])

    # 工艺跨导 k_n' = μ_n · C_ox
    kp = p["mu_n"] * C_ox

    # 宽长比
    WL = W_cm / L_cm

    # Gismondi 平滑参数 δ（使线性区→饱和区 C¹ 连续）
    Vov_max_d = max(p["VGS_max"] - VTH, 0.1)
    delta = 0.02 * Vov_max_d

    return {
        "C_ox_F_cm2":   C_ox,
        "C_ox_fF_um2":  C_ox * 1e7,
        "VFB_V":        VFB,
        "gamma_V12":    gamma,
        "VTH_V":        VTH,
        "kp_A_V2":      kp,
        "WL":           WL,
        "W_cm":         W_cm,
        "L_cm":         L_cm,
        "t_ox_cm":      t_ox_cm,
        "delta_V":      delta,
    }


# ============================================================================
# ID 网格计算
# ============================================================================

def compute_id_grid(params, derived):
    """
    在 VGS × VDS 网格上计算漏极电流 ID [A]。
    返回 (VGS_grid, VDS_grid, ID_grid)。
    """
    vgs = np.linspace(0.0, params["VGS_max"], params["n_pts"])
    vds = np.linspace(0.0, params["VDS_max"], params["n_pts"])

    VGS, VDS = np.meshgrid(vgs, vds)          # shape: (n_VDS, n_VGS)
    VTH  = derived["VTH_V"]
    mu_n = params["mu_n"]
    C_ox = derived["C_ox_F_cm2"]
    WL   = derived["WL"]
    lam  = params["lam"]

    ID = np.zeros_like(VGS)
    Vov = VGS - VTH                           # 过驱动电压
    beta = mu_n * C_ox * WL                  # 增益因子 [A/V²]

    # Gismondi 平滑 min 函数：消除线性区 ↔ 饱和区导数不连续
    # VDS_eff ≈ VDS（线性区） → Vov（饱和区），C¹ 连续
    # 当 δ→0 时精确还原 Level 1 分段模型
    delta = derived["delta_V"]
    VDS_eff = 0.5 * (VDS + Vov - np.sqrt((VDS - Vov)**2 + 4.0 * delta**2))

    mask_active = (VGS > VTH) & (VDS >= 0)
    ID[mask_active] = (
        beta
        * (Vov[mask_active] * VDS_eff[mask_active]
           - 0.5 * VDS_eff[mask_active] ** 2)
        * (1.0 + lam * VDS[mask_active])
    )

    return VGS, VDS, ID


# ============================================================================
# 三维绘图 (v2: 新增 3D 公式显示)
# ============================================================================

def plot_id_3d(VGS, VDS, ID_A, derived, params, query_point=None):
    """
    绘制 ID–VGS–VDS 三维曲面图 + 俯视热力图，鼠标悬停联动高亮。
    v2 新增：3D 视图内公式显示框，展示计算公式及代入数值。
    query_point: (vgs, vds, id_mA) 或 None — 查询高亮点。
    """
    ID_mA = ID_A * 1e3
    VTH = derived["VTH_V"]
    VGS_max = params["VGS_max"]
    VDS_max = params["VDS_max"]
    mu_n = params["mu_n"]
    C_ox = derived["C_ox_F_cm2"]
    WL   = derived["WL"]
    lam  = params["lam"]
    n_pts = params["n_pts"]
    delta = derived["delta_V"]

    beta_A  = mu_n * C_ox * WL                # [A/V²]
    beta_mA = beta_A * 1e3                    # [mA/V²]

    # ---- 1D 网格（用于交互式查找） ----
    vgs_1d = np.linspace(0.0, VGS_max, n_pts)
    vds_1d = np.linspace(0.0, VDS_max, n_pts)

    # ========================================================================
    #  图面布局：左 3D + 右 2D 俯视
    # ========================================================================
    fig = plt.figure(figsize=(22, 10))
    ax3d = fig.add_subplot(1, 2, 1, projection="3d")
    ax2d = fig.add_subplot(1, 2, 2)

    # ========================================================================
    #  左侧 — 3D 曲面
    # ========================================================================
    surf = ax3d.plot_surface(
        VGS, VDS, ID_mA,
        cmap=cm.viridis,
        linewidth=0,
        antialiased=False,
        alpha=0.88,
        rstride=max(1, n_pts // 15),
        cstride=max(1, n_pts // 15),
    )

    ax3d.view_init(elev=28, azim=-58)

    ax3d.set_xlabel("Gate–Source Voltage $V_{GS}$ [V]", fontsize=11, labelpad=10)
    ax3d.set_ylabel("Drain–Source Voltage $V_{DS}$ [V]", fontsize=11, labelpad=10)
    ax3d.set_zlabel("Drain Current $I_D$ [mA]", fontsize=11, labelpad=10)
    ax3d.set_xlim(0.0, VGS_max)
    ax3d.set_ylim(0.0, VDS_max)
    ax3d.set_box_aspect([VGS_max, VDS_max, max(ID_mA.max(), 1e-6)])

    # 3D 标题
    title_lines = [
        "NMOS $I_D$–$V_{GS}$–$V_{DS}$ Characteristic Surface  (v2: formula display)",
        (
            f"$W={params['W_um']}\\,\\mu\\mathrm{{m}},\\;"
            f"L={params['L_um']}\\,\\mu\\mathrm{{m}},\\;"
            f"t_\\mathrm{{ox}}={params['t_ox_nm']}\\,\\mathrm{{nm}},\\;"
            f"V_\\mathrm{{TH}}={VTH:.3f}\\,\\mathrm{{V}},\\;"
            f"\\mu_n={mu_n}\\,\\mathrm{{cm}}^2\\!/\\mathrm{{V·s}}$"
        ),
    ]
    ax3d.set_title("\n".join(title_lines), fontsize=10, pad=14, linespacing=1.3)

    # ========================================================================
    #  3D 区域边界线 & 文字标签（可选 — 按 h 键切换）
    # ========================================================================
    region_artists = []

    if params.get("show_regions", True):
        delta_plot = derived["delta_V"]

        def id_at(vgs_val, vds_val):
            if vgs_val <= VTH or vds_val < 0:
                return 0.0
            vov = vgs_val - VTH
            vds_eff = 0.5 * (vds_val + vov
                              - np.sqrt((vds_val - vov)**2 + 4.0 * delta_plot**2))
            return (mu_n * C_ox * WL
                    * (vov * vds_eff - 0.5 * vds_eff**2)
                    * (1.0 + lam * vds_val)) * 1e3

        # 边界线 1: VGS = VTH
        if 0 < VTH < VGS_max:
            nb = 80
            vds_b1 = np.linspace(0, VDS_max, nb)
            line1 = ax3d.plot([VTH] * nb, vds_b1, np.zeros(nb),
                              color="crimson", linewidth=2.0, linestyle="--", zorder=10)
            region_artists.extend(line1)
            txt1 = ax3d.text(VTH - 0.03, VDS_max * 0.55, 0,
                             "$V_{GS}=V_{TH}$", color="crimson", fontsize=9,
                             fontweight="bold", ha="right", va="bottom")
            region_artists.append(txt1)

        # 边界线 2: VDS = VGS − VTH（饱和边界，浅色虚线）
        vgs_sat_max = min(VGS_max, VTH + VDS_max)
        if vgs_sat_max > VTH + 0.01:
            vgs_b2 = np.linspace(VTH + 0.002, vgs_sat_max, 200)
            vds_b2 = vgs_b2 - VTH
            id_b2  = np.array([id_at(vg, vd) for vg, vd in zip(vgs_b2, vds_b2)])
            line2 = ax3d.plot(vgs_b2, vds_b2, id_b2,
                              color="lightsalmon", linewidth=1.8, linestyle="--",
                              alpha=0.9, zorder=10)
            region_artists.extend(line2)
            mid = len(vgs_b2) // 2 + len(vgs_b2) // 6
            txt2 = ax3d.text(vgs_b2[mid], vds_b2[mid], id_b2[mid] * 1.05,
                             "$V_{DS}=V_{GS}{-}V_{TH}$", color="lightsalmon",
                             fontsize=9, fontweight="bold", ha="left", va="bottom")
            region_artists.append(txt2)

            line2p = ax3d.plot(vgs_b2, vds_b2, np.zeros_like(vgs_b2),
                               color="lightsalmon", linewidth=1.0, linestyle=":",
                               alpha=0.35, zorder=6)
            region_artists.extend(line2p)

        # 边界线 3: VDS = 0.15×(VGS−VTH)（深线性上限）
        if vgs_sat_max > VTH + 0.01:
            vgs_b3 = np.linspace(VTH + 0.005, vgs_sat_max, 200)
            vds_b3 = 0.15 * (vgs_b3 - VTH)
            id_b3  = np.array([id_at(vg, vd) for vg, vd in zip(vgs_b3, vds_b3)])
            line3 = ax3d.plot(vgs_b3, vds_b3, id_b3,
                              color="darkorange", linewidth=1.5, linestyle="--",
                              alpha=0.85, zorder=8)
            region_artists.extend(line3)

        # 文字标注
        has_active = (VTH < VGS_max)
        Vov_range = VGS_max - VTH if has_active else 0.0
        if has_active and Vov_range > 0.1:
            vgs_active = VTH + Vov_range * 0.60
            vov_active = vgs_active - VTH
            vds_deep = vov_active * 0.06
            vds_lin  = vov_active * 0.50
            vds_sat  = min(vov_active * 1.40, VDS_max)
            z_deep = id_at(vgs_active, vds_deep)
            z_lin  = id_at(vgs_active, vds_lin)
            z_sat  = id_at(vgs_active, vds_sat)

            if vds_deep < VDS_max and z_deep > 0:
                t_deep = ax3d.text(vgs_active, vds_deep,
                                   z_deep + max(ID_mA.max() * 0.04, 0.02),
                                   "Deep Triode", color="#c44e00", fontsize=10,
                                   fontweight="bold", ha="center", va="bottom",
                                   bbox=dict(boxstyle="round,pad=0.3",
                                             facecolor="lightyellow",
                                             edgecolor="#c44e00", alpha=0.82))
                region_artists.append(t_deep)
            if vds_lin < VDS_max and z_lin > 0:
                t_lin = ax3d.text(vgs_active, vds_lin,
                                  z_lin + max(ID_mA.max() * 0.04, 0.02),
                                  "Triode", color="#1a6fb5", fontsize=10,
                                  fontweight="bold", ha="center", va="bottom",
                                  bbox=dict(boxstyle="round,pad=0.3",
                                            facecolor="lightyellow",
                                            edgecolor="#1a6fb5", alpha=0.82))
                region_artists.append(t_lin)
            if vds_sat <= VDS_max and z_sat > 0:
                t_sat = ax3d.text(vgs_active, vds_sat,
                                  z_sat + max(ID_mA.max() * 0.04, 0.02),
                                  "Saturation", color="#1b6c1b", fontsize=10,
                                  fontweight="bold", ha="center", va="bottom",
                                  bbox=dict(boxstyle="round,pad=0.3",
                                            facecolor="lightyellow",
                                            edgecolor="#1b6c1b", alpha=0.82))
                region_artists.append(t_sat)

        # 截止区
        if VTH > 0.15:
            vgs_cutoff = VTH * 0.45
        else:
            vgs_cutoff = max(VTH - 0.20, 0.02)
        if vgs_cutoff > 0 and vgs_cutoff < VGS_max * 0.9:
            t_cut = ax3d.text(vgs_cutoff, VDS_max * 0.52,
                              max(ID_mA.max() * 0.015, 0.005),
                              "Cutoff", color="#777777", fontsize=10,
                              fontweight="bold", ha="center", va="bottom",
                              bbox=dict(boxstyle="round,pad=0.3",
                                        facecolor="whitesmoke",
                                        edgecolor="#999999", alpha=0.82))
            region_artists.append(t_cut)

        # 亚阈值区
        if VTH > 0.08:
            vgs_sub = VTH - 0.06
        else:
            vgs_sub = max(VTH - 0.04, 0.01)
        if vgs_sub > 0.005 and vgs_sub < VGS_max * 0.95:
            t_sub = ax3d.text(vgs_sub, VDS_max * 0.52,
                              max(ID_mA.max() * 0.015, 0.005),
                              "Subthreshold", color="#b07c2e", fontsize=10,
                              fontweight="bold", ha="center", va="bottom",
                              bbox=dict(boxstyle="round,pad=0.3",
                                        facecolor="lemonchiffon",
                                        edgecolor="#b07c2e", alpha=0.82))
            region_artists.append(t_sub)

        # 图例
        lines_proxy = [
            plt.Line2D([0], [0], color="lightsalmon", linewidth=1.8, linestyle="--",
                       label="VDS=VGS-VTH (sat.)"),
            plt.Line2D([0], [0], color="crimson", linewidth=2.0, linestyle="--",
                       label="VGS=VTH (cutoff)"),
            plt.Line2D([0], [0], color="darkorange", linewidth=1.5, linestyle="--",
                       label="VDS=0.15Vov (deep tri.)"),
        ]
        leg = ax3d.legend(handles=lines_proxy, loc="upper left", fontsize=8,
                          framealpha=0.85, ncol=1)
        region_artists.append(leg)

    # 按键事件
    if region_artists:
        def _toggle_regions(event):
            if event.key == 'h':
                for art in region_artists:
                    art.set_visible(not art.get_visible())
                fig.canvas.draw_idle()

        fig.canvas.mpl_connect('key_press_event', _toggle_regions)
        fig.text(0.22, 0.005, "Press 'h' to toggle region annotations",
                 ha='center', fontsize=9, color='gray', style='italic',
                 transform=fig.transFigure)

    # 3D 网格
    ax3d.xaxis.pane.fill = False
    ax3d.yaxis.pane.fill = False
    ax3d.zaxis.pane.fill = False
    ax3d.grid(True, alpha=0.4)

    # ========================================================================
    #  右侧 — 2D 俯视热力图
    # ========================================================================
    heatmap = ax2d.pcolormesh(VGS, VDS, ID_mA, cmap=cm.viridis, shading="auto")
    cbar = fig.colorbar(heatmap, ax=ax2d, shrink=0.85, aspect=20, pad=0.04)
    cbar.set_label("Drain Current $I_D$ [mA]", fontsize=11, labelpad=8)

    ax2d.set_xlabel("Gate–Source Voltage $V_{GS}$ [V]", fontsize=12)
    ax2d.set_ylabel("Drain–Source Voltage $V_{DS}$ [V]", fontsize=12)
    ax2d.set_title("Top-Down View  —  hover for details", fontsize=12)
    ax2d.set_xlim(0.0, VGS_max)
    ax2d.set_ylim(0.0, VDS_max)
    ax2d.set_aspect("auto")

    # ---- 交互元素 ----
    cursor_2d, = ax2d.plot([], [], "r+", markersize=14, markeredgewidth=2,
                           zorder=100)
    highlight_3d, = ax3d.plot([], [], [], "ro", markersize=10,
                              markeredgecolor="darkred", markeredgewidth=1.5,
                              zorder=200)

    # v2: 点击蓝色标记
    click_marker_2d, = ax2d.plot([], [], "bs", markersize=12, markeredgewidth=2,
                                 zorder=150)
    click_marker_3d, = ax3d.plot([], [], [], "bD", markersize=10,
                                 markeredgecolor="darkblue", markeredgewidth=1.5,
                                 zorder=210)
    click_annot = None   # 点击标注 (将在 _on_click 中创建)

    # 信息浮窗（固定在 2D 面板左上角，避免遮挡深线性区）
    info_box = ax2d.text(
        0.02, 0.98, "", transform=ax2d.transAxes,
        fontsize=10, verticalalignment="top", horizontalalignment="left",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white",
                  edgecolor="#555555", alpha=0.92),
    )

    # ========================================================================
    #  v2 新增：公式弹窗引用（点击触发）
    # ========================================================================
    popup_fig_ref = [None]   # mutable container for nonlocal access in _on_click

    # ---- 工作区判断 ----
    def _region_name(vgs, vds):
        if vgs < VTH - 0.06:
            return "Cutoff"
        if vgs < VTH:
            return "Subthreshold"
        vov = vgs - VTH
        if vov <= 0:
            return "Cutoff"
        if vds < 0.15 * vov:
            return "Deep Triode"
        if vds < vov:
            return "Triode"
        return "Saturation"

    # ========================================================================
    #  v2 新增：构建公式文本（各工作区不同）
    # ========================================================================

    def _build_formula_text(vgs_val, vds_val, id_val, region):
        """
        根据工作区生成公式显示文本，含代入数值。
        返回 ANSI 风格的多行字符串（matplotlib 用 $...$ 数学模式）。
        """
        vov_val = vgs_val - VTH

        # 计算 Gismondi VDS_eff
        if vgs_val > VTH and vds_val >= 0:
            vds_eff_val = 0.5 * (vds_val + vov_val
                                 - np.sqrt((vds_val - vov_val)**2 + 4.0 * delta**2))
        else:
            vds_eff_val = 0.0

        region_colors = {
            "Cutoff":       "#aaaaaa",
            "Subthreshold": "#d4a843",
            "Deep Triode":  "#ff6b35",
            "Triode":       "#4da6ff",
            "Saturation":   "#4ec94e",
        }
        color = region_colors.get(region, "#ffffff")

        header = f"═══ {region} ═══"

        lines = [header, ""]

        # ── 共通量 ──
        lines.append(
            f"$\\beta = \\mu_n C_{{ox}} (W/L)"
            f" = {mu_n} \\times {C_ox:.4e} \\times {WL:.2f}"
            f" = {beta_A:.4e}$ A/V$^2$"
        )

        if region in ("Cutoff", "Subthreshold"):
            # 截止 / 亚阈值区
            lines.append("")
            if region == "Cutoff":
                lines.append(
                    f"$V_{{GS}} = {vgs_val:.4f} < V_{{TH}} = {VTH:.4f}$ "
                    f"$\\Rightarrow$ channel not formed"
                )
            else:
                lines.append(
                    f"$V_{{GS}} = {vgs_val:.4f} \\lesssim V_{{TH}} = {VTH:.4f}$ "
                    f"(weak inversion)"
                )
            lines.append(
                f"$I_D \\approx 0$  (Level 1 model)"
            )
            if id_val > 0:
                lines.append(
                    f"$I_D$ (Gismondi) = {id_val:.6f} mA"
                )

        elif region == "Deep Triode":
            # 深线性区：给出 Ron 公式
            r_on = vds_val / (id_val * 1e-3) if id_val > 1e-12 else float('inf')

            lines.append(
                f"$V_{{ov}} = V_{{GS}} - V_{{TH}}"
                f" = {vgs_val:.4f} - {VTH:.4f} = {vov_val:.4f}$ V"
            )
            lines.append(
                f"$V_{{DS,eff}} \\approx V_{{DS}}$"
                f" (Gismondi: {vds_eff_val:.4f} V)"
            )
            lines.append("")
            lines.append(
                "$\\mathbf{Deep\\ Triode\\ (V_{DS} \\ll V_{ov}):}$"
            )
            lines.append(
                f"$I_D \\approx \\beta \\cdot V_{{ov}} \\cdot V_{{DS}}$"
            )
            lines.append(
                f"$\\;\\;\\; = {beta_A:.4e} \\times {vov_val:.4f}"
                f" \\times {vds_val:.4f}$"
            )
            id_approx_dt = beta_mA * vov_val * vds_val
            lines.append(
                f"$\\;\\;\\; = {id_approx_dt:.4f}$ mA"
            )
            lines.append("")
            lines.append(
                "$\\mathbf{R_{on}\\ (on-resistance):}$"
            )
            lines.append(
                f"$R_{{\\mathrm{{on}}}}"
                f" = \\frac{{V_{{DS}}}}{{I_D}}"
                f" \\approx \\frac{{1}}{{\\beta \\cdot V_{{ov}}}}$"
            )
            lines.append(
                f"$\\;\\;\\; = \\frac{{1}}"
                f"{{{beta_A:.4e} \\times {vov_val:.4f}}}$"
            )
            ron_approx = 1.0 / (beta_A * vov_val)
            if ron_approx < 1e3:
                lines.append(
                    f"$\\;\\;\\; \\approx {ron_approx:.2f}$ Ω"
                )
            else:
                lines.append(
                    f"$\\;\\;\\; \\approx {ron_approx/1e3:.3f}$ kΩ"
                )
            lines.append("")
            if r_on < 1e3:
                lines.append(
                    f"$R_{{\\mathrm{{on}}}}$ (from $I_D$ exact)"
                    f" = {r_on:.2f} Ω"
                )
            else:
                lines.append(
                    f"$R_{{\\mathrm{{on}}}}$ (from $I_D$ exact)"
                    f" = {r_on/1e3:.3f} kΩ"
                )
            lines.append("")
            gm = beta_mA * vds_val
            lines.append(f"$g_m \\approx \\beta \\cdot V_{{DS}}$ = {gm:.4f} mA/V")

        elif region == "Triode":
            # 线性区
            lines.append(
                f"$V_{{ov}} = V_{{GS}} - V_{{TH}}"
                f" = {vgs_val:.4f} - {VTH:.4f} = {vov_val:.4f}$ V"
            )
            lines.append(
                f"$V_{{DS,eff}}$ (Gismondi) = {vds_eff_val:.4f} V"
            )
            lines.append("")
            lines.append(
                "$\\mathbf{Triode\\ (Linear)\\;Region:}$"
            )
            lines.append(
                f"$I_D = \\beta \\left["
                f"V_{{ov}} V_{{DS,eff}} - \\frac{{1}}{{2}} V_{{DS,eff}}^2"
                f"\\right](1 + \\lambda V_{{DS}})$"
            )
            lines.append("")
            lines.append(" ── Step-by-step ──")
            term1 = vov_val * vds_eff_val
            term2 = 0.5 * vds_eff_val**2
            paren = term1 - term2
            mod   = 1.0 + lam * vds_val
            lines.append(
                f"$V_{{ov}} \\cdot V_{{DS,eff}}"
                f" = {vov_val:.4f} \\times {vds_eff_val:.4f}"
                f" = {term1:.4f}$"
            )
            lines.append(
                f"$\\frac{{1}}{{2}} V_{{DS,eff}}^2"
                f" = 0.5 \\times {vds_eff_val:.4f}^2"
                f" = {term2:.4f}$"
            )
            lines.append(
                f"$[\\;\\;\\;\\;\\;\\;]"
                f" = {term1:.4f} - {term2:.4f}"
                f" = {paren:.4f}$ V$^2$"
            )
            lines.append(
                f"$(1 + \\lambda V_{{DS}})"
                f" = 1 + {lam} \\times {vds_val:.4f}"
                f" = {mod:.4f}$"
            )
            lines.append(
                f"$I_D = {beta_A:.4e} \\times {paren:.4f}"
                f" \\times {mod:.4f}$"
            )
            id_exact = beta_mA * paren * mod
            lines.append(
                f"$I_D = {id_exact:.4f}$ mA"
            )
            lines.append("")
            gm = beta_mA * vds_val * (1.0 + lam * vds_val)
            lines.append(f"$g_m = \\beta \\cdot V_{{DS}} \\cdot (1+\\lambda V_{{DS}})$"
                         f" = {gm:.4f} mA/V")

        elif region == "Saturation":
            # 饱和区
            lines.append(
                f"$V_{{ov}} = V_{{GS}} - V_{{TH}}"
                f" = {vgs_val:.4f} - {VTH:.4f} = {vov_val:.4f}$ V"
            )
            lines.append(
                f"$V_{{DS,eff}} \\approx V_{{ov}}$"
                f" (Gismondi: {vds_eff_val:.4f} V)"
            )
            lines.append("")
            lines.append(
                "$\\mathbf{Saturation\\;Region}\\;(V_{DS} \\geq V_{ov}):$"
            )
            lines.append(
                f"$I_D \\approx \\frac{{1}}{{2}}"
                f"\\beta \\cdot V_{{ov}}^2 "
                f"\\cdot (1 + \\lambda V_{{DS}})$"
            )
            lines.append("")
            lines.append(" ── Step-by-step ──")
            term_sq = 0.5 * vov_val**2
            mod_sat = 1.0 + lam * vds_val
            lines.append(
                f"$\\frac{{1}}{{2}} V_{{ov}}^2"
                f" = 0.5 \\times {vov_val:.4f}^2"
                f" = {term_sq:.4f}$ V$^2$"
            )
            lines.append(
                f"$(1 + \\lambda V_{{DS}})"
                f" = 1 + {lam} \\times {vds_val:.4f}"
                f" = {mod_sat:.4f}$"
            )
            lines.append(
                f"$I_D \\approx {beta_A:.4e} \\times {term_sq:.4f}"
                f" \\times {mod_sat:.4f}$"
            )
            id_approx_sat = beta_mA * term_sq * mod_sat
            lines.append(
                f"$I_D \\approx {id_approx_sat:.4f}$ mA"
            )
            lines.append("")
            lines.append(" ── Gismondi exact ──")
            # Use the actual Gismondi formula for comparison
            vds_eff_sat = vds_eff_val
            term1_s = vov_val * vds_eff_sat
            term2_s = 0.5 * vds_eff_sat**2
            paren_s = term1_s - term2_s
            id_gis = beta_mA * paren_s * mod_sat
            lines.append(
                f"$I_D$ (Gismondi) = {id_gis:.4f} mA"
            )
            lines.append("")
            gm = beta_mA * vov_val * mod_sat
            lines.append(
                f"$g_m = \\beta \\cdot V_{{ov}} \\cdot (1+\\lambda V_{{DS}})$"
                f" = {gm:.4f} mA/V"
            )

            # ---- 本征增益 (Intrinsic Gain) ----
            # gds = ∂ID/∂VDS = (1/2)*β*Vov²*λ  (饱和区近似)
            gds_S = 0.5 * beta_A * vov_val**2 * lam
            Av = gm * 1e-3 / gds_S   # gm [mA/V] → [A/V]，Av 无量纲 V/V
            lines.append("")
            lines.append("$\\mathbf{Intrinsic\\ Gain\\ (A_v = g_m / g_{ds}):}$")
            lines.append(
                f"$g_{{\\mathrm{{ds}}}} \\approx "
                f"\\frac{{\\partial I_D}}{{\\partial V_{{\\mathrm{{DS}}}}}}"
                f" = \\frac{{1}}{{2}} \\lambda \\beta V_{{\\mathrm{{ov}}}}^2$"
            )
            lines.append(
                f"$\\;\\;\\; = 0.5 \\times {lam} \\times {beta_A:.4e}"
                f" \\times {vov_val:.4f}^2$"
            )
            if gds_S < 1e-3:
                lines.append(
                    f"$\\;\\;\\; = {gds_S:.4e}$ S"
                    f"  (${gds_S*1e6:.2f}$ µS)"
                )
            else:
                lines.append(
                    f"$\\;\\;\\; = {gds_S:.4e}$ S"
                )
            lines.append("")
            lines.append(
                f"$A_v = \\frac{{g_m}}{{g_{{\\mathrm{{ds}}}}}}"
                f" = \\frac{{{gm:.4f}\\ \\mathrm{{mA/V}}}}"
                f"{{{gds_S:.4e}\\ \\mathrm{{S}}}}$"
            )
            lines.append(
                f"$\\;\\;\\; = {Av:.1f}$ V/V"
                f"  (${20 * np.log10(max(Av, 1e-12)):.1f}$ dB)"
            )

        return "\n".join(lines)

    # ========================================================================
    #  鼠标悬停联动（含 v2 公式显示）
    # ========================================================================

    def _on_hover(event):
        if event.inaxes != ax2d:
            return
        if event.xdata is None or event.ydata is None:
            return

        # 最近网格点
        iv = np.argmin(np.abs(vgs_1d - event.xdata))
        id_ = np.argmin(np.abs(vds_1d - event.ydata))
        vgs_val = vgs_1d[iv]
        vds_val = vds_1d[id_]
        id_val  = ID_mA[id_, iv]
        region  = _region_name(vgs_val, vds_val)

        # 更新 2D 十字光标
        cursor_2d.set_data([vgs_val], [vds_val])

        # 更新 2D 信息浮窗（含 gm / R_on）
        vov_val = vgs_val - VTH

        lines = [
            f"$V_{{GS}}$ = {vgs_val:.3f} V",
            f"$V_{{DS}}$ = {vds_val:.3f} V",
            f"$I_D$ = {id_val:.4f} mA",
        ]

        if region == "Deep Triode":
            if id_val > 1e-12:
                r_on = vds_val / (id_val * 1e-3)
                lines.append(f"$R_{{\\mathrm{{on}}}}$ = {r_on:.2f} Ω" if r_on < 1e3
                             else f"$R_{{\\mathrm{{on}}}}$ = {r_on/1e3:.3f} kΩ")
            gm = beta_mA * vds_val
            lines.append(f"$g_m$ = {gm:.4f} mA/V")
        elif region == "Triode":
            gm = beta_mA * vds_val * (1.0 + lam * vds_val)
            lines.append(f"$g_m$ = {gm:.4f} mA/V")
        elif region == "Saturation":
            gm = beta_mA * vov_val * (1.0 + lam * vds_val)
            gds_S = 0.5 * beta_A * vov_val**2 * lam
            Av = gm * 1e-3 / gds_S if gds_S > 1e-30 else float('inf')
            lines.append(f"$g_m$ = {gm:.4f} mA/V")
            lines.append(f"$g_{{\\mathrm{{ds}}}}$ = {gds_S:.4e} S")
            lines.append(f"$A_v$ = {Av:.1f} V/V  ({20*np.log10(max(Av, 1e-12)):.1f} dB)")
        elif region == "Subthreshold":
            lines.append(f"$g_m$ ≈ 0 mA/V")

        lines.append(f"Region: {region}")
        info_box.set_text("\n".join(lines))

        # 更新 3D 红色高亮点
        highlight_3d.set_data([vgs_val], [vds_val])
        highlight_3d.set_3d_properties([id_val])

        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", _on_hover)

    # ========================================================================
    #  v2 新增：左键点击 → 公式弹窗
    # ========================================================================

    def _show_formula_popup(vgs_val, vds_val, id_val, region, rcolor):
        """弹出独立公式窗口，显示计算过程。重复点击关闭旧窗重开。"""
        # 关闭旧弹窗
        if popup_fig_ref[0] is not None:
            try:
                plt.close(popup_fig_ref[0])
            except Exception:
                pass

        formula_text = _build_formula_text(vgs_val, vds_val, id_val, region)
        raw_lines = formula_text.split("\n")

        # 估算窗口高度
        n_lines = len(raw_lines)
        fig_h = max(4.5, n_lines * 0.32 + 1.2)

        popup_fig = plt.figure(figsize=(8.0, fig_h))
        popup_fig_ref[0] = popup_fig

        # 深色背景
        popup_fig.patch.set_facecolor("#0b1220")

        # 逐行渲染：各行独立判定是否为 mathtext（以 $ 包裹为数学行）
        spacing = 0.94 / max(n_lines, 1)
        for i, line in enumerate(raw_lines):
            y_pos = 0.96 - i * spacing
            if y_pos < 0.02:
                break
            stripped = line.strip()
            if not stripped:
                continue
            # header 行（纯文本，带颜色和加粗）
            if stripped.startswith("═══") and stripped.endswith("═══"):
                popup_fig.text(
                    0.06, y_pos, stripped,
                    fontsize=12, fontweight="bold",
                    color=rcolor, family="monospace",
                    va="top",
                )
            elif stripped.startswith("──") or stripped.startswith("──"):
                # 分隔线
                popup_fig.text(
                    0.06, y_pos, stripped,
                    fontsize=9, color="#888888", family="monospace",
                    va="top",
                )
            elif stripped.startswith("$"):
                # mathtext 数学行
                popup_fig.text(
                    0.06, y_pos, stripped,
                    fontsize=10, color="#e0e0e0",
                    va="top",
                )
            else:
                # 普通文本行
                popup_fig.text(
                    0.06, y_pos, stripped,
                    fontsize=9.5, color="#c8c8c8", family="monospace",
                    va="top",
                )

        # 标题
        popup_fig.suptitle(
            f"Formula  —  {region}",
            fontsize=13, color=rcolor, fontweight="bold",
            x=0.06, ha="left", y=0.995,
        )

        popup_fig.subplots_adjust(top=0.93, bottom=0.03, left=0.02, right=0.98)
        popup_fig.show()

    def _on_click(event):
        """左键点击 2D 热力图 → 蓝色标记 + 弹出公式窗口"""
        nonlocal click_annot

        if event.inaxes != ax2d:
            return
        if event.xdata is None or event.ydata is None:
            return
        if event.button != 1:   # 仅左键
            return

        iv = np.argmin(np.abs(vgs_1d - event.xdata))
        id_ = np.argmin(np.abs(vds_1d - event.ydata))
        vgs_val = vgs_1d[iv]
        vds_val = vds_1d[id_]
        id_val = ID_mA[id_, iv]
        region = _region_name(vgs_val, vds_val)

        # 更新 2D 蓝色方块标记
        click_marker_2d.set_data([vgs_val], [vds_val])
        # 更新 3D 蓝色菱形标记
        click_marker_3d.set_data([vgs_val], [vds_val])
        click_marker_3d.set_3d_properties([id_val])

        # 更新 / 创建点击标注
        if click_annot is not None:
            try:
                click_annot.remove()
            except Exception:
                pass
        click_annot = ax2d.annotate(
            f"VGS={vgs_val:.3f}\nVDS={vds_val:.3f}\nID={id_val:.4f}",
            xy=(vgs_val, vds_val), xytext=(14, 14), textcoords="offset points",
            fontsize=8.5, color="#1a3a8a", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#e8f0ff",
                      edgecolor="#3366cc", alpha=0.90),
            zorder=202,
        )

        rcolor = {
            "Cutoff":       "#aaaaaa",
            "Subthreshold": "#d4a843",
            "Deep Triode":  "#ff6b35",
            "Triode":       "#4da6ff",
            "Saturation":   "#4ec94e",
        }.get(region, "#ffffff")

        _show_formula_popup(vgs_val, vds_val, id_val, region, rcolor)

        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("button_press_event", _on_click)

    # ---- 2D 面板上的区域边界投影 ----
    if params.get("show_regions", True) and 0 < VTH < VGS_max:
        # VGS = VTH
        ax2d.axvline(VTH, color="crimson", linewidth=1.2, linestyle="--", alpha=0.6)
        # VDS = VGS - VTH
        vgs_b = np.linspace(VTH, vgs_sat_max, 100)
        ax2d.plot(vgs_b, vgs_b - VTH, color="lightsalmon", linewidth=1.2,
                  linestyle="--", alpha=0.7)

    # ---- 查询高亮（蓝色标记） ----
    if query_point is not None:
        q_vgs, q_vds, q_id = query_point
        # 2D 蓝色星形
        ax2d.plot(q_vgs, q_vds, "b*", markersize=16, markeredgewidth=1.5,
                  zorder=200, label="Query")
        # 3D 蓝色菱形
        ax3d.plot([q_vgs], [q_vds], [q_id], "bD", markersize=11,
                  markeredgecolor="darkblue", markeredgewidth=1.5, zorder=200)
        # 简要标注
        ax2d.annotate(
            f"VGS={q_vgs:.3f}\nVDS={q_vds:.3f}\nID={q_id:.3f}",
            xy=(q_vgs, q_vds), xytext=(12, 12), textcoords="offset points",
            fontsize=9, color="darkblue", fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightcyan",
                      edgecolor="darkblue", alpha=0.88),
            zorder=201,
        )

    plt.tight_layout()
    return fig


# ============================================================================
# 查询辅助：给定两个量查第三个
# ============================================================================

def find_point(params, derived, target):
    """
    target: {"VGS": val, "ID": val} 或 {"VDS": val, "ID": val} 等
    返回 (vgs_found, vds_found, id_mA_found) 或 None
    """
    VGS_grid, VDS_grid, ID_A = compute_id_grid(params, derived)
    ID_mA = ID_A * 1e3
    vgs_1d = VGS_grid[0, :]
    vds_1d = VDS_grid[:, 0]

    if "VGS" in target and "VDS" in target:
        iv = np.argmin(np.abs(vgs_1d - target["VGS"]))
        id_ = np.argmin(np.abs(vds_1d - target["VDS"]))
        return vgs_1d[iv], vds_1d[id_], ID_mA[id_, iv]

    if "VGS" in target and "ID" in target:
        iv = np.argmin(np.abs(vgs_1d - target["VGS"]))
        col = ID_mA[:, iv]
        id_ = np.argmin(np.abs(col - target["ID"]))
        return vgs_1d[iv], vds_1d[id_], ID_mA[id_, iv]

    if "VDS" in target and "ID" in target:
        id_ = np.argmin(np.abs(vds_1d - target["VDS"]))
        row = ID_mA[id_, :]
        iv = np.argmin(np.abs(row - target["ID"]))
        return vgs_1d[iv], vds_1d[id_], ID_mA[id_, iv]

    return None


# ============================================================================
# 控制台摘要
# ============================================================================

def print_summary(params, derived):
    """打印所有器件物理参数到控制台。"""
    print()
    print("=" * 60)
    print("  NMOS 器件物理参数摘要")
    print("=" * 60)
    print(f"  栅宽 W          = {params['W_um']:.2f} μm")
    print(f"  沟道长度 L      = {params['L_um']:.2f} μm")
    print(f"  宽长比 W/L       = {derived['WL']:.2f}")
    print(f"  氧化层厚度 t_ox  = {params['t_ox_nm']:.2f} nm")
    print(f"  C_ox             = {derived['C_ox_F_cm2']:.4e} F/cm²")
    print(f"                   = {derived['C_ox_fF_um2']:.3f} fF/μm²")
    print(f"  衬底掺杂 N_sub   = {params['N_sub']:.2e} cm⁻³")
    print(f"  费米势 φ_F       = {params['phi_F']:.3f} V")
    print(f"  平带电压 V_FB    = {derived['VFB_V']:.3f} V")
    print(f"  体效应系数 γ     = {derived['gamma_V12']:.4f} V^(1/2)")
    print(f"  阈值电压 V_TH    = {derived['VTH_V']:.3f} V")
    print(f"  工艺跨导 k_n'    = {derived['kp_A_V2']:.4e} A/V²")
    print(f"  迁移率 μ_n       = {params['mu_n']} cm²/(V·s)")
    print(f"  沟长调制 λ       = {params['lam']} V⁻¹")
    print("-" * 60)
    print(f"  VGS 扫描范围     = 0 – {params['VGS_max']} V")
    print(f"  VDS 扫描范围     = 0 – {params['VDS_max']} V")
    print(f"  采样点数         = {params['n_pts']} × {params['n_pts']}")
    print("=" * 60)

    # 一个参考工作点（饱和区）
    VTH = derived["VTH_V"]
    if VTH < params["VGS_max"]:
        vgs_sample = min(VTH + 1.5, params["VGS_max"])
        vds_sample = min(vgs_sample - VTH, params["VDS_max"])
        Vov = vgs_sample - VTH
        id_sat = 0.5 * params["mu_n"] * derived["C_ox_F_cm2"] * derived["WL"] * Vov**2
        mod_sat = 1.0 + params["lam"] * vds_sample
        id_sat_mod = id_sat * mod_sat
        gm_sat = params["mu_n"] * derived["C_ox_F_cm2"] * derived["WL"] * Vov * mod_sat
        gds_sat = 0.5 * params["mu_n"] * derived["C_ox_F_cm2"] * derived["WL"] * Vov**2 * params["lam"]
        Av_sat = gm_sat / gds_sat if gds_sat > 1e-30 else float('inf')
        print(f"\n  参考工作点 (VGS={vgs_sample:.2f}V, VDS={vds_sample:.2f}V, 饱和):")
        print(f"    I_D ≈ {id_sat_mod*1e3:.3f} mA  ({id_sat_mod*1e6:.1f} μA)")
        print(f"    g_m ≈ {gm_sat*1e3:.3f} mA/V")
        print(f"    g_ds ≈ {gds_sat:.4e} S  ({gds_sat*1e6:.2f} μS)")
        print(f"    本征增益 A_v = g_m/g_ds ≈ {Av_sat:.1f} V/V  ({20*np.log10(max(Av_sat, 1e-12)):.1f} dB)")
    print()


# ============================================================================
# 命令行接口
# ============================================================================

def build_parser():
    p = argparse.ArgumentParser(
        description="NMOS ID–VGS–VDS 三维特性曲面图生成器 (v2: formula display)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
 示例:
   # 命令行一次性传参
   python v2.py --W_um 10 --L_um 1 --t_ox_nm 20 --N_sub 1e16 --phi_F 0.35

   # 交互式输入（无命令行参数时自动进入）
   python v2.py

   # 部分命令行 + 其余交互式补充
   python v2.py --W_um 20 --L_um 1.5
         """,
    )

    g_dev = p.add_argument_group("核心器件参数（必须提供，无默认值）")
    g_dev.add_argument("--W_um",    type=float, default=None, help="栅宽 [μm]")
    g_dev.add_argument("--L_um",    type=float, default=None, help="有效沟道长度 [μm]")
    g_dev.add_argument("--t_ox_nm", type=float, default=None, help="栅极氧化层厚度 [nm]")
    g_dev.add_argument("--C_ox",    type=float, default=None, help="单位面积栅氧化层电容 [F/cm²]（留空则从 t_ox 推算）")
    g_dev.add_argument("--N_sub",   type=float, default=None, help="衬底掺杂浓度 [cm⁻³]")
    g_dev.add_argument("--phi_F",   type=float, default=None, help="费米势 [V]")

    g_aux = p.add_argument_group("辅助参数（可选，有默认值）")
    g_aux.add_argument("--mu_n",    type=float, default=None, help=f"电子迁移率 [cm²/(V·s)]（默认 {AUX_DEFAULTS['mu_n']}）")
    g_aux.add_argument("--VFB",     type=float, default=None, help=f"平带电压 [V]（默认由功函数差推算）")
    g_aux.add_argument("--lam",     type=float, default=None, help=f"沟道长度调制系数 [V⁻¹]（默认 {AUX_DEFAULTS['lam']}）")

    g_sweep = p.add_argument_group("扫描 / 输出设置")
    g_sweep.add_argument("--VGS_max", type=float, default=None, help=f"VGS 扫描上限 [V]（默认 {AUX_DEFAULTS['VGS_max']}）")
    g_sweep.add_argument("--VDS_max", type=float, default=None, help=f"VDS 扫描上限 [V]（默认 {AUX_DEFAULTS['VDS_max']}）")
    g_sweep.add_argument("--n_pts",   type=int,   default=None, help=f"每轴采样点数（默认 {AUX_DEFAULTS['n_pts']}）")
    g_sweep.add_argument("--no-plot", action="store_true",        help="不显示图形（仅打印摘要）")
    g_sweep.add_argument("--save",       type=str,   default=None,   help="保存图像到指定路径（如 output.png）")
    g_sweep.add_argument("--no-regions", action="store_true",         help="不显示工作区域标注和分割线（图形中仍可按 h 键切换）")

    return p


# ============================================================================
# 主入口
# ============================================================================

def main():
    parser = build_parser()
    args = parser.parse_args()

    # ---- 确定是需要交互输入还是全部来自命令行 ----
    core_keys = ["W_um", "L_um", "t_ox_nm", "N_sub", "phi_F"]
    all_core_provided = all(getattr(args, k, None) is not None for k in core_keys)

    if all_core_provided:
        # 全部核心参数已通过命令行提供 —— 跳过交互输入
        params = {
            "W_um":    args.W_um,
            "L_um":    args.L_um,
            "t_ox_nm": args.t_ox_nm,
            "C_ox":    args.C_ox,
            "N_sub":   args.N_sub,
            "phi_F":   args.phi_F,
            "mu_n":    args.mu_n    if args.mu_n    is not None else AUX_DEFAULTS["mu_n"],
            "VFB":     args.VFB     if args.VFB     is not None else AUX_DEFAULTS["VFB"],
            "lam":     args.lam     if args.lam     is not None else AUX_DEFAULTS["lam"],
            "VGS_max": args.VGS_max if args.VGS_max is not None else AUX_DEFAULTS["VGS_max"],
            "VDS_max": args.VDS_max if args.VDS_max is not None else AUX_DEFAULTS["VDS_max"],
            "n_pts":        args.n_pts   if args.n_pts   is not None else AUX_DEFAULTS["n_pts"],
            "show_regions": not args.no_regions,
        }
    else:
        # 缺少核心参数 → 交互式补充
        params, interactive_save = collect_params_interactive(args)
        # 合并交互式提供的保存路径（命令行 --save 优先）
        if interactive_save and not args.save:
            args.save = interactive_save
        # 补充命令行未指定的辅助参数
        for key, default in AUX_DEFAULTS.items():
            if key not in params or params[key] is None:
                params[key] = default
        params["show_regions"] = not args.no_regions

    # ---- 参数短名映射（命令提示符修改用） ----
    KEY_MAP = {
        "W": "W_um", "L": "L_um", "tox": "t_ox_nm", "t_ox": "t_ox_nm",
        "VGS": "VGS_max", "VDS": "VDS_max",
    }
    INT_KEYS = {"n_pts"}

    # ====================================================================
    #  显示 → 关闭 → 修改 → 重算 循环
    # ====================================================================
    _first = True
    query_point = None   # 查询高亮点（跨轮保持）

    def _region_static(vgs, vds, vth):
        if vgs < vth - 0.06:
            return "Cutoff"
        if vgs < vth:
            return "Subthreshold"
        vov = vgs - vth
        if vov <= 0:
            return "Cutoff"
        if vds < 0.15 * vov:
            return "Deep Triode"
        if vds < vov:
            return "Triode"
        return "Saturation"

    while True:
        # ---- 计算 ----
        try:
            derived = compute_derived(params)
        except Exception as e:
            print(f"\n[错误] 物理量计算失败: {e}", file=sys.stderr)
            break

        VGS, VDS, ID = compute_id_grid(params, derived)
        print_summary(params, derived)

        if args.no_plot:
            break   # --no-plot 模式只跑一次

        # ---- 显示图形 ----
        fig = plot_id_3d(VGS, VDS, ID, derived, params, query_point)
        if _first and args.save:
            fig.savefig(args.save, dpi=150, bbox_inches="tight")
            print(f"[信息] 图像已保存至: {args.save}")
        _first = False

        print("\n  \033[33m\033[1m触控提示\033[0m  图形窗口中按 \033[1mh\033[0m 键可切换区域标注的显示/隐藏")
        print("  \033[33m\033[1mv2 功能\033[0m   鼠标悬停 → 3D 视图左上角显示计算公式与代入数值")
        print("  \033[33m\033[1m修改提示\033[0m  关闭图形窗口后可修改参数并重新生成\n")
        plt.show()

        # ---- 修改 / 查询循环 ----
        print("-" * 58)
        print("  Modify: key=value (comma-sep)         |  修改：key=value（逗号分隔）")
        print("  Query : ? key=value, key=value        |  查询：? key=value, key=value")
        print("  Params: W, L, t_ox, N_sub, phi_F, C_ox, VGS, VDS, mu_n, lam, VFB, n_pts")
        print("  e.g.: W=20, VDS=8        e.g.: ? VGS=1.5, ID=0.8")
        print("  Enter=replot, q=quit                  |  留空重算，q 退出")
        print("-" * 58)

        query_point = None   # 本轮查询结果

        cmd = input("> ").strip()
        if cmd.lower() in ("q", "quit", "exit"):
            break
        if cmd == "":
            continue

        # ---- 查询模式（以 ? 或 find 开头） ----
        if cmd.startswith("?") or cmd.lower().startswith("find"):
            raw = cmd[1:].strip() if cmd.startswith("?") else cmd[4:].strip()
            target = {}
            for part in raw.split(","):
                part = part.strip()
                if not part or "=" not in part:
                    continue
                k, v = part.split("=", 1)
                k = k.strip().upper()
                try:
                    target[k] = float(v.strip())
                except ValueError:
                    print(f"  [!] '{v.strip()}' 不是有效数值")

            if len(target) != 2:
                print("  [!] 查询需要恰好两个已知量，如 ? VGS=1.5, ID=0.8")
                continue

            result = find_point(params, derived, target)
            if result is None:
                print("  [!] 未找到匹配点")
                continue

            q_vgs, q_vds, q_id = result
            region = _region_static(q_vgs, q_vds, derived["VTH_V"])
            print(f"  查询结果：VGS={q_vgs:.3f}V  VDS={q_vds:.3f}V  ID={q_id:.4f}mA  [{region}]")
            query_point = (q_vgs, q_vds, q_id)
            continue

        # ---- 修改模式 ----
        updated = 0
        for part in cmd.split(","):
            part = part.strip()
            if not part:
                continue
            if "=" not in part:
                print(f"  [!] 忽略 '{part}'（需要 key=value 格式）")
                continue

            key, val_str = part.split("=", 1)
            key = key.strip()
            val_str = val_str.strip()

            key = KEY_MAP.get(key, key)           # 短名 → 全名

            if key in params:
                try:
                    old = params[key]
                    if key in INT_KEYS:
                        params[key] = int(float(val_str))
                    else:
                        new_val = float(val_str)
                        params[key] = new_val if new_val != 0.0 or key not in ("C_ox", "VFB") else None
                    print(f"  {key}: {old} → {params[key]}")
                    updated += 1
                except ValueError:
                    print(f"  [!] '{val_str}' 不是有效数值")
            else:
                print(f"  [!] 未知参数 '{key}'（可用：{', '.join(params.keys())}）")

        if updated == 0:
            print("  [!] 未检测到有效修改，将使用原参数重算\n")


if __name__ == "__main__":
    main()
