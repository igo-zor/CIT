"""
第十章實驗介面用：論文符號與中文語義對照（對齊《約束世界論》§10.1.1 與各節參數表）。

僅供 Streamlit 顯示，不影響數值計算。
"""

from __future__ import annotations

from typing import Any

# ---------------------------------------------------------------------------
# 表單標籤：「論文記號 — 中文（節次）」
# ---------------------------------------------------------------------------

L10_CAND_ENUM = "候選超邊母集大小 M（枚舉預估用；§10.2 候選母集）"

L10_N = "$n:=|V|$ — 節點數（§10.1.1）"
L10_K_MAX = "$k_{max}$ — 單條超邊最大階數（§10.1.1）"
L10_K_MIN = "$k_{min}$ — 單條超邊最小階數（§10.1.1）"
L10_M_MAX = "$m_{max}$ — 超邊數上限（程式）；論文常記 $m=|E|$（§10.1.1）"
L10_SAMPLE_LIMIT = "$N_{cand}$ — 候選抽樣數／近似 $|\\mathcal{C}_{cand}|$（§10.2.6）"
L10_N_CFG = (
    "$N_{cfg}$ — 自可採用域抽樣之觀測集大小 $|\\mathcal C_{obs}|$（§10.3.2–10.3.4；"
    "§10.5.10 表亦用此符號表示過濾後用於解析統計之合法樣本數；若可採用數不足則以全集為觀測集並於輸出標示）"
)
L10_N_REP = "$N_{rep}$ — 重複次數／觀測集重抽次數（§10.3.3；靜態重抽穩健性）"
# seed：RNG 基底整數；論文參數表之「種子批次數」為 $N_{seed}$（僅 §10.2 候選抽樣），兩者不同。
L10_SEED = (
    "seed — 偽隨機基底種子（整數；固定 RNG 以重現結果；"
    "非論文 $N_{seed}$，後者為 §10.2.6「種子批次數」）"
)
L10_N_SEED_BATCHES = (
    "$N_{seed}$ — 種子批次數（候選抽樣輪數；§10.2.6；"
    "程式用 seed+i·1009 派生各輪 RNG；與左欄基底 seed 不同）"
)

L10_D_MAX = "$d_{max}$ — 頂點最大允許度數（域型約束之一；§10.2.6）"
L10_CONNECTED = "連通 — 2-section 圖連通（域型 $\Lambda_{dom}$ 條件；§10.2）"
L10_FORBID_TRI = "禁二元△ — 禁止僅由二元超邊構成之三角形（域型／forbidden motif；§10.2）"

L10_DELTA = "$\delta$ — 解析簽名距離閾值，定義 $\\approx_\delta$（§10.1.1、§10.3）"
L10_SIG_WEAK = "簽名層次 — 對應論文 $\mathrm{Sig}_\Lambda$ 之 weak / medium / strong 操作化（§10.3）"

L10_CFG_ADMISSIBLE = "|Cfg_Λ| — 可採用配置數（合法域；§10.1.1）"
L10_S_LAMBDA = "|S_{Λ,δ}| — 解析單元數（穩定化後；§10.1.1）"
L10_OVERLAP = "重疊率 — $R_{overlap}$，相容鄰域交疊程度（§10.3.5）"
L10_TABLE_T3_DELTA = "表 10-3 簽名比較掃描之 $\\delta$（§10.3）"

# §10.4
L10_PI = "π — 粗細投影 π_{Λ′→Λ}（§10.1.1、§10.4）"
L10_SIG_COARSE = "粗層 $\\mathrm{Sig}_\\Lambda$（weak／medium／strong；§10.4）"
L10_SIG_FINE = "細層 $\\mathrm{Sig}_{\\Lambda'}$（§10.4）"
L10_DELTA_COARSE = "粗層閾值 $\\delta$（§10.4）"
L10_DELTA_FINE = "細層閾值 $\\delta'$（§10.4）"
L10_SIG_A = "A 鏈簽名（motif 細分側；§10.4 進階順序比較）"
L10_SIG_B = "B 鏈簽名（度數細分側；§10.4 進階順序比較）"
L10_KERNEL = "條件核 — 纖維上 $K_s$（uniform / proportional；§10.4）"
L10_PUSH_ERR = "推前最大誤差 — $\varepsilon_{push}$ 型度量（§10.4.2）"
L10_H_COARSE = "H(p_Λ) — 粗層熵 $H(p)$（§10.4.2）"
L10_H_FINE = "H(p_{Λ′}) — 細層熵（§10.4.2）"

# §10.5（《約束世界論 30》小節編號）
L10_NA_NB = "n_A / n_B — 二分 $\mathcal{V}=\mathcal{V}_A\\sqcup\mathcal{V}_B$（§10.5.2（一））"
L10_M_EDGES = "$m=|E|$ — 超邊總數（§10.1.1、§10.5.4（三））"
L10_ALPHA_CROSS = "$\\alpha_{cross}$ — 跨區塊邊生成傾向（§10.5.4（一））"
L10_DELTA_ENT = (
    "δ_ent — 整體解析映射 $q_\\Lambda^{ent}$ 之距離閾 $\delta$（程式操作化；§10.5.3、§10.5.6（三））"
)
L10_N_CFG_10_5 = (
    "$N_{cfg}$ — 域型過濾後用於解析統計之合法樣本數（§10.5.4（五）、§10.5.10）"
)
L10_RHO_IRRED = "ρ_irred — 解析不可分解單元比例（§10.5.5（一））"
L10_D_SEP = "D_sep — 與 $p_A\\otimes p_B$ 之偏離（§10.5.5（三））"
L10_I_AB = "I(A;B) — 互資訊（§10.5.5（四））"

# §10.6
L10_N_CTX = "$N_{ctx}$ — 上下文族樣本數（§10.6.4（五）、§10.6.6 建議表）"
L10_NODES_BIT = "$n$ — 全域位元變數個數 $|V|$（§10.6.6 建議表；parity 系統）"
L10_M_CTX = "$M$ — 上下文數／局部視窗個數（§10.6.4（一））"
L10_W_CTX = (
    "$w$ — 每個上下文視窗所含節點數（§10.6.4（二）；"
    "與 §10.9 時間聚合視窗之 $w$ 不同）"
)
L10_ETA_CTX = (
    "$\\eta$ — 相鄰視窗交疊節點數（§10.6.4（三）；"
    "非 §10.8 微擾強度 $\\eta$）"
)
L10_T_LOC = "$T_{\\mathrm{loc}}$ — 局部條件型別數（§10.6.4（四）；程式現以 parity 為主）"
L10_N_SEARCH_106 = "$N_{\\mathrm{search}}$ — 全域搜尋步數上限（§10.6.3、§10.6.6 建議表；parity 合併實為多項式時間）"
L10_CHI_GLUE = r"$\chi_{\mathrm{glue}}$ — 不可全域延拓指標（§10.6.5（一））"
L10_RHO_GLUE = r"$\rho_{\mathrm{glue}}$ — 拼合障礙比例（§10.6.5（二））"
L10_R_MIN = r"$r_{\min}$ — 最小衝突核大小（§10.6.5（三）；程式以「上下文」子族計）"
L10_RHO_VAL = r"$\rho_{\mathrm{val}}$ — 共同賦值可行率（§10.6.5（四））"
L10_N_VAL_106 = r"$N_{\mathrm{val}}$ — 共同賦值可行樣本數（§10.6.5（四））"
L10_MODE_CTX = "mode — obstruction／satisfiable（§10.6.6；parity 循環視窗）"

# §10.7
L10_T_STEPS = "$T$ — 演化步數（§10.1.1、§10.7.4）"
L10_N0_RUNS = "$N_0$ — 軌道／初始樣本數（§10.7.6）"
L10_JS_TERM = "JS — 終端分布 Jensen–Shannon（bit；§10.4.1、§10.7）"
L10_OBS_SIG = (
    "$\\Lambda_{\\mathrm{obs}}$（程式鍵 signature／sig_obs）— weak／medium／strong 解析觀測層（§10.7.4、§10.8、§10.9）"
)
L10_M_TRIAL = (
    "$M_{\\mathrm{trial}}$ — 每步候選更新數（§10.7.4、§10.7.6、§10.9.6；§10.7 全實驗接線；§10.9 軌道已接線）"
)
L10_W_H = (
    "$w_H$ — 熵計算滑動窗口寬度（§10.7.4；本版僅逐步 $H_\\Lambda^{(\\ell)}$，尚未接線）"
)
L10_W_A = (
    "$w_A$ — 吸引子／週期檢測窗口（§10.7.4；本版核心尚未接線）"
)
L10_P_MAX = "$P_{\\max}$ — 最大週期搜尋長度（§10.7.6；本版核心尚未接線）"
L10_N_SEED_107 = (
    "$N_{\\mathrm{seed}}$ — §10.7.6 建議之種子重跑次數（穩健性；"
    "本版單次實驗僅用基底 seed，尚未以迴圈接線）"
)

# §10.8
L10_ETA = "$\eta$ — 微擾強度（§10.1.1、§10.8.4）"
L10_T_SB = "$T_{sb}$ — 破缺／後續重組步數（§10.8.4）"
L10_R_DEPTH = "$r$ — 局部鄰域深度（incidence 圖；§10.1.1）"
L10_N_TYPE = "$N_{type}$ — 局部型別數（§10.8.5）"
L10_INIT_FAMILY = "init_family — 初態族 sym／pert／rand（§10.8）"
L10_N_SAMPLES = "初態抽取樣本數（§10.8）"

# §10.9
L10_W_WIN = "$w$ — 時間聚合視窗寬度（§10.1.1、§10.9.2）"
L10_H_MACRO = "$H_{macro}^{(w)}$ — 宏觀型別熵（§10.9.2）"
L10_W_LIST = "window_list — 多個 $w$ 以逗號分隔（§10.9.2）"
L10_DELTA_T = "$\\Delta t$ — 宏觀觀測時間步進（聚合起點間隔；§10.9.4）"
L10_N_HIST = "$N_{\\mathrm{hist}}$ — 微觀歷史條數（§10.9.6）"
L10_N_SEED_RUNS_109 = (
    "$N_{\\mathrm{seed}}$ — §10.9.6 每組參數重跑次數（穩健性；"
    "程式鍵 n_seed_runs；**非** §10.2 候選採樣之種子批次數）"
)
L10_EPS_PLAT_109 = (
    r"$\varepsilon_{\mathrm{plat}}^{(w)}$ — 宏觀平台判定閾（§10.9.4、§10.9.6；"
    r"本實作施於 $C_{\mathrm{eff}}$ 代理序列）"
)


MD_10_1_1_CORE = r"""
**§10.1.1 統一符號（與《約束世界論》正文對齊）**

| 記號 | 中文／語義 |
|------|------------|
| $c=(V,E)$ | 有限超圖配置；$V$ 節點集，$E$ 超邊集 |
| $n:=|V|$, $m:=|E|$ | 節點數、超邊數 |
| $k_{min},k_{max}$ | 允許之最小／最大超邊階數 |
| $\mathcal{C}_{cand}$ | 候選母集 |
| $\Lambda_{dom}$, $\mathrm{Cfg}_\Lambda$ | 域型約束、合法配置域 |
| $\mathrm{Sig}_\Lambda$, $d_{\mathrm{sig}}$ | 解析簽名、簽名空間差異量 |
| $\delta$ | 閾值；$c_i\approx_\delta c_j \Leftrightarrow d_{\mathrm{sig}}\le\delta$ |
| $\sim_\delta$, $S_{\Lambda,\delta}$, $q_{\Lambda,\delta}$ | 穩定化等價、解析單元空間、解析投影 |
| $\Lambda\sqsubseteq\Lambda'$ | 細化；$\pi_{\Lambda'\to\Lambda}$ 粗細投影 |
| $p_\Lambda$, $H(p)$ | 粗層分布、熵 |
| $T$ | 演化序列長度 $c^{(0)}\to\cdots\to c^{(T)}$ |
| $\eta$, $T_{sb}$ | 微擾強度、破缺步數（§10.8） |
| $r$, $w$ | 局部鄰域深度、時間聚合視窗（§10.9） |
"""


def render_section_glossary(streamlit_module: Any, section: str) -> None:
    """
    在頁面中插入「論文符號對照」摺疊區。

    Args:
        streamlit_module: 已 import 之 ``streamlit``。
        section: ``"10.1"`` | ``"10.2"`` | … | ``"10.9"``。
    """
    st = streamlit_module
    if section == "10.1":
        with st.expander("論文符號對照（§10.1.1　統一符號與操作型定義）", expanded=False):
            st.markdown(MD_10_1_1_CORE)
        return

    bodies = {
        "10.2": r"""
| 程式／輸出 | 論文記號 | 中文 |
|------------|----------|------|
| n, k_max, m_max | $n$, $k_{max}$, 上限 $m$ | 節點數、最大超邊階數、超邊數上限（候選生成） |
| sample_limit | $N_{cand}$ | 候選樣本數 |
| seed（整數） | （基底種子） | 偽隨機基底種子；固定後可重現同一抽樣／動力學 |
| num_seeds | $N_{seed}$ | 種子批次數：候選抽樣輪數（每輪派生 seed+i·1009）；**不是** seed 本體 |
| \|Cfg\| 各層 | $\|\mathrm{Cfg}_{\Lambda^{(i)}}\|$ | 各级域型下合法配置數 |
| 保留比例 | $r_{keep}^{(i)}$ | 佔候選之比例（§10.2.5） |
""",
        "10.3": r"""
| 程式 | 論文記號 | 中文 |
|------|----------|------|
| sample_limit | $N_{cand}$ | 候選生成上限；先篩得 $\mathrm{Cfg}_\Lambda$ 再抽觀測集 |
| n_cfg | $N_{cfg}$ | 自可採用域抽樣之觀測配置數 $|\mathcal C_{obs}|$（不足時用全集並提示） |
| n_rep | $N_{rep}$ | 重抽次數；固定其餘參數後重複抽取觀測集以檢驗穩健性 |
| δ | $\delta$ | 解析相容閾值，定義 $\approx_\delta$ |
| seed | （基底種子） | 偽隨機基底種子；**非** §10.2.6 之 $N_{seed}$（種子批次數） |
| 簽名 weak/medium/strong | $\sigma_{res}$ | 解析簽名版本（§10.3.4） |
| 重疊率 | $R_{overlap}$ | 相容鄰域交疊（§10.3.5） |
| 傳遞違反率 | $R_{trans\text{-}viol}$ | $\approx_\delta$ 非傳遞程度 |
| \|S_Λ\| | $\|S_{\Lambda,\delta}\|$ | 解析單元數 |
""",
        "10.4": r"""
| 程式 | 論文記號 | 中文 |
|------|----------|------|
| seed | （基底種子） | 偽隨機基底種子；**非** $N_{seed}$ 種子批次數（後者僅 §10.2） |
| 粗/細簽名+δ | $\Lambda$, $\Lambda'$ | 粗細解析層（$\Lambda\sqsubseteq\Lambda'$） |
| A/B 鏈簽名 | $A,B$ | 進階順序比較：A→B 與 B→A（A 可對應 motif 細分，B 可對應度數細分） |
| π | $\pi_{\Lambda'\to\Lambda}$ | 細類投影到粗類 |
| 條件核 | $K_s$ | 纖維上條件分布（§10.4） |
| 推前誤差 | $\varepsilon_{push}$ | 粗細推前一致性 |
| JS | $\mathrm{JS}(\cdot,\cdot)$ | 終端分布差異（bit；§10.4.1） |
""",
        "10.5": r"""
| 程式 | 論文記號 | 中文 |
|------|----------|------|
| seed | （基底種子） | 偽隨機基底種子；**非** $N_{seed}$ 種子批次數 |
| n_A, n_B | $n_A$, $n_B$ | 二分兩側節點數（§10.5.4（二）） |
| k_min, k_max | $k_{min}$, $k_{max}$ | 超邊階數範圍（§10.5.4（四）、§10.5.10） |
| m_edges | $m=|E|$ | 超邊總數（§10.5.4（三）） |
| sample_limit | $N_{cand}$ | 候選樣本數（§10.5.4（五）、§10.5.10） |
| n_cfg | $N_{cfg}$ | 過濾後解析統計用合法樣本數（§10.5.10） |
| α_cross | $\alpha_{cross}$ | 跨區塊耦合傾向（§10.5.4（一）） |
| δ_ent | （程式） | $q_\Lambda^{ent}$ 分類距離閾（§10.5.6（三）） |
| ρ_irred | $\rho_{irred}$ | 不可分解單元比例（§10.5.5（一）） |
| D_sep, I(A;B) | $D_{sep}$, $I(A;B)$ | 與局部分離之偏離、互資訊（§10.5.5） |
""",
        "10.6": r"""
| 程式 | 論文記號 | 中文 |
|------|----------|------|
| seed | （基底種子） | 偽隨機基底種子；**非** §10.2 之 $N_{seed}$ 種子批次數 |
| n_nodes | $n$ | 節點總數（§10.6.6 建議表） |
| n_ctx | $N_{ctx}$ | 隨機上下文族樣本數（§10.6.4（五）） |
| M | $M$ | 上下文數（§10.6.4（一）） |
| w_ctx | $w$ | 視窗大小（§10.6.4（二）） |
| eta_ctx | $\eta$ | 交疊大小（§10.6.4（三）；非 §10.8 微擾 $\eta$） |
| T_loc | $T_{\mathrm{loc}}$ | 局部型別數（§10.6.4（四）；現僅 parity） |
| n_search | $N_{\mathrm{search}}$ | 搜尋上限（§10.6.3；GF(2) 實際極快） |
| n_val | $N_{\mathrm{val}}$ | 可行樣本數（§10.6.5（四）） |
| χ_glue | $\chi_{\mathrm{glue}}$ | 不可全域延拓指標（§10.6.5（一）） |
| ρ_glue | $\rho_{\mathrm{glue}}$ | 拼合障礙比例（§10.6.5（二）） |
| r_min | $r_{\min}$ | 最小衝突核（§10.6.5（三）；程式：上下文子族） |
| ρ_val | $\rho_{\mathrm{val}}$ | 共同賦值可行率（§10.6.5（四）） |
""",
        "10.7": r"""
| 程式 | 論文記號 | 中文 |
|------|----------|------|
| seed | （基底種子） | 偽隨機基底種子；**非** §10.2 之種子批次數，亦**異於** §10.7.6 表之 $N_{\mathrm{seed}}$（重跑次數） |
| steps | $T$ | 演化步數（§10.7.4） |
| runs | $N_0$ | 初始軌道數 |
| signature | $\Lambda_{\mathrm{obs}}$ | 解析觀測層 weak／medium／strong |
| m_trial, w_h, w_a, p_max, n_seed_107 | $M_{\mathrm{trial}},w_H,w_A,P_{\max},N_{\mathrm{seed}}$ | 論文建議主線已於介面記錄；軌道核心尚未全部接線者見參數表備註 |
| 熵序列 | $H_\Lambda^{(\ell)}$ | 解析觀測層上之熵（§10.7.5） |
| n_reach_mean 等 | $N_{\mathrm{reach}}$ | 軌跡上相異解析單元數之彙總（§10.7.5） |
| JS（實驗 B） | $\mathrm{JS}$ | 雙路徑終端分布差異（擴充；非 §10.7.5 主表） |
""",
        "10.8": r"""
| 程式 | 論文記號 | 中文 |
|------|----------|------|
| seed | （基底種子） | 偽隨機基底種子；**非** $N_{seed}$ 種子批次數 |
| η | $\eta$ | 微擾強度 |
| T_sb | $T_{sb}$ | 破缺步數 |
| r | $r$ | 局部鄰域深度 |
| N_type, H_type | $N_{type}$, $H_{type}$ | 型別數、型別熵 |
""",
        "10.9": r"""
| 程式 | 論文記號 | 中文 |
|------|----------|------|
| n, k_max, m_max, d_max, sample_limit | $n$, $k_{max}$, 上限 $m$, $d_{max}$, $N_{cand}$ | 域型與候選採樣（§10.9.3；與 §10.7 域型約定類似） |
| connected, forbid_pair_triangles | （域型） | 2-section 連通、禁二元三角 |
| steps | $T$ | 每條微觀歷史演化步數 |
| m_trial | $M_{\mathrm{trial}}$ | 每步候選更新數 |
| sig_obs | $\sigma_{\mathrm{obs}}$ | weak／medium／strong 解析簽名（宏觀標籤） |
| delta_t | $\Delta t$ | 宏觀觀測沿時間前進之步進 |
| epsilon_plat | $\varepsilon_{\mathrm{plat}}^{(w)}$ | 平台判定閾（施於 $C_{\mathrm{eff}}$ 代理序列） |
| n_hist | $N_{\mathrm{hist}}$ | 微觀歷史條數（輸出對各歷史再平均） |
| n_seed_runs | $N_{\mathrm{seed}}$ | §10.9.6 穩健性重跑次數（**非** §10.2 種子批次數） |
| seed | （基底種子） | 偽隨機基底種子 |
| window_list | $w$ 集合 | 時間聚合視窗寬度（建議含 1 作參照） |
| r | $r$ | 空間鄰域深度（宏觀標籤） |
| R_edge_bar | $\bar R_{\mathrm{edge}}$ | 微觀邊周轉率平均 |
| H_macro | $H_{\mathrm{macro}}^{(w)}$ | 宏觀型別熵 |
| JS_vs_w1 | $\mathrm{JS}_w$ | 相對 $w=1$ 之分布差異（bit；§10.9.2） |
""",
    }
    body = bodies.get(section, "")
    if not body:
        return
    with st.expander(f"論文符號對照（§{section}）", expanded=False):
        st.markdown(body.strip())


def render_full_ch10_sidebar_note(streamlit_module: Any) -> None:
    """側欄最下方簡短提示（可選）。"""
    st = streamlit_module
    st.sidebar.caption(
        "符號以論文 §10.1.1 為準；各頁摺疊「論文符號對照」可對照節次。"
        " seed 為偽隨機基底種子；§10.2 另欄論文 N_seed 為種子批次數（兩者不同）。"
    )
