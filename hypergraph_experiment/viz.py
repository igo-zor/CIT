"""
圖表與超圖視覺化：供 Streamlit 與封存 PNG 使用（相依 matplotlib、networkx）。
"""

from __future__ import annotations

import io
import platform
import statistics
from itertools import combinations
from typing import Any, List, Sequence, Tuple

import matplotlib.pyplot as plt
import networkx as nx
from matplotlib import font_manager
from matplotlib.patches import Circle, Polygon

# 模組層級：確保全圖使用中文字型（Microsoft JhengHei / YaHei 等），避免 DejaVu 缺字警告
_CONFIGURED_CJK: bool = False


def _configure_matplotlib_chinese_font() -> None:
    """
    依作業系統挑選可覆蓋 CJK 的字型，並寫入 rcParams。

    Windows 通常內建「微軟正黑體」「微軟雅黑」；若皆不可用則仍可能出現警告。
    """
    global _CONFIGURED_CJK
    if _CONFIGURED_CJK:
        return
    _CONFIGURED_CJK = True
    plt.rcParams["axes.unicode_minus"] = False

    if platform.system() == "Windows":
        candidates = (
            "Microsoft JhengHei",
            "Microsoft YaHei",
            "Microsoft JhengHei UI",
            "Microsoft YaHei UI",
            "SimHei",
        )
    elif platform.system() == "Darwin":
        candidates = ("PingFang TC", "Heiti TC", "STHeiti", "Arial Unicode MS")
    else:
        candidates = (
            "Noto Sans CJK TC",
            "Noto Serif CJK TC",
            "WenQuanYi Zen Hei",
            "Droid Sans Fallback",
        )

    # 取第一個實際解析為非 DejaVu 的 family（findfont 找不到時會落到 DejaVu）
    primary: list[str] = []
    for family in candidates:
        path = font_manager.findfont(font_manager.FontProperties(family=family))
        if path and "dejavu" not in path.lower():
            primary = [family]
            break

    if primary:
        tail = [f for f in plt.rcParams.get("font.sans-serif", []) if f not in primary]
        plt.rcParams["font.sans-serif"] = primary + tail
    elif platform.system() == "Windows":
        # 後備：直接寫入常見中文族名（若above檢測失敗仍嘗試）
        plt.rcParams["font.sans-serif"] = ["Microsoft JhengHei", "Microsoft YaHei", "sans-serif"]


_configure_matplotlib_chinese_font()


def _tab20_color(index: int) -> Any:
    """自 tab20 取一色（相容新舊 matplotlib 之 colormap API）。"""
    import matplotlib as mpl

    try:
        cmap = mpl.colormaps["tab20"]
    except (AttributeError, KeyError):
        cmap = mpl.cm.get_cmap("tab20")
    colors = getattr(cmap, "colors", None)
    if colors is not None:
        return colors[index % len(colors)]
    return cmap((index % 20) / 19.0)


def fig_to_png_bytes() -> bytes:
    """將當前 matplotlib figure 轉成 PNG 位元組並關閉 figure。"""
    _configure_matplotlib_chinese_font()
    buf = io.BytesIO()
    plt.savefig(buf, format="png", bbox_inches="tight", dpi=120)
    plt.close()
    buf.seek(0)
    return buf.read()


def matplotlib_figure_to_png_bytes(fig: Any) -> bytes:
    """
    將指定 matplotlib Figure 光栅化為 PNG 位元組（**不**關閉 figure）。

    典型流程：繪圖後呼叫本函式取得位元組，再以 ``st.pyplot(fig)`` 顯示，最後由呼叫端
    ``plt.close(fig)`` 釋放資源；並可搭配 ``st.download_button`` 提供檔案下載（與 §10.3 頁面輔圖一致）。

    Args:
        fig: matplotlib 之 ``Figure`` 實例（``matplotlib.figure.Figure``）。

    Returns:
        PNG 圖檔之原始位元組。
    """
    _configure_matplotlib_chinese_font()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
    buf.seek(0)
    return buf.read()


def plot_hypergraph_incidence(config: dict[str, Any]) -> bytes:
    """
    以二部圖呈現超圖：一側為頂點、一側為超邊節點。

    Args:
        config: 含 ``vertices``（列表）、``hyperedges``（列表之列表）之字典。

    Returns:
        PNG 圖檔位元組。
    """
    verts = list(config.get("vertices") or [])
    edges_raw = config.get("hyperedges") or []

    G = nx.Graph()
    for v in verts:
        G.add_node(f"v{v}", bipartite=0, label=str(v))
    for i, e in enumerate(edges_raw):
        ename = f"e{i}"
        G.add_node(ename, bipartite=1, label=f"e{i}|{len(e)}")
        for v in e:
            G.add_edge(f"v{v}", ename)

    pos = {}
    n_v = len(verts)
    for j, v in enumerate(sorted(verts)):
        pos[f"v{v}"] = (0.0, j / max(n_v - 1, 1))
    n_e = len(edges_raw)
    for j in range(n_e):
        pos[f"e{j}"] = (1.0, j / max(n_e - 1, 1) if n_e > 1 else 0.5)

    plt.figure(figsize=(6, 4))
    v_nodes = [n for n, d in G.nodes(data=True) if d.get("bipartite") == 0]
    e_nodes = [n for n, d in G.nodes(data=True) if d.get("bipartite") == 1]
    nx.draw_networkx_nodes(G, pos, nodelist=v_nodes, node_color="#4ecdc4", node_size=500, label="vertices")
    nx.draw_networkx_nodes(G, pos, nodelist=e_nodes, node_color="#ff6b6b", node_size=400, label="hyperedges")
    nx.draw_networkx_edges(G, pos, alpha=0.7)
    labels = {n: n.replace("v", "") if n.startswith("v") else n for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)
    plt.axis("off")
    plt.title("超圖（頂點–超邊 二部圖）")
    return fig_to_png_bytes()


def _convex_hull_2d(points: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """
    平面點集之凸包（單調鏈，無 SciPy 相依）。
    共線時可能僅剩兩點，上層繪製應改畫線段。
    """
    pts = sorted(set(points))
    if len(pts) <= 2:
        return pts

    def cross(o: Tuple[float, float], a: Tuple[float, float], b: Tuple[float, float]) -> float:
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: List[Tuple[float, float]] = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper: List[Tuple[float, float]] = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def plot_hypergraph_spatial(
    config: dict[str, Any],
    *,
    layout_seed: int = 42,
    scale: float = 1.2,
) -> bytes:
    """
    以接近 Wolfram Physics Project 常見之**空間嵌入**呈現超圖：

    - 僅顯示頂點節點（無二部圖之「超邊方塊」）。
    - 佈局：在 **2-section 圖**（同超邊內頂點互鄰）上做 spring layout。
    - 每條超邊：二元畫線段；三元以上畫**半透明凸包**（著色區域）。

    與 Wolfram ``HypergraphPlot`` 仍可能有細節差異（其內部為專用排版與動畫），
    但語意上同為「多點共屬一高階關係」之平面展示。

    Args:
        config: ``vertices``、``hyperedges`` 與 ``plot_hypergraph_incidence`` 相同。
        layout_seed: ``spring_layout`` 隨機種子，便於重現。
        scale: 座標縮放，避免貼邊。

    Returns:
        PNG 位元組。
    """
    verts = sorted(config.get("vertices") or [])
    edges_raw = config.get("hyperedges") or []
    vset = set(verts)

    G = nx.Graph()
    for v in verts:
        G.add_node(v)
    for e in edges_raw:
        el = [x for x in e if x in vset]
        for a, b in combinations(sorted(el), 2):
            G.add_edge(a, b)

    if not verts:
        plt.figure(figsize=(5, 3))
        plt.text(0.5, 0.5, "無頂點", ha="center", va="center")
        plt.axis("off")
        return fig_to_png_bytes()

    pos = nx.spring_layout(G, seed=layout_seed, k=scale / max(len(verts), 1) ** 0.5, iterations=50)
    pos = {v: (scale * x, scale * y) for v, (x, y) in pos.items()}

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect("equal")
    ax.axis("off")

    for i, e in enumerate(edges_raw):
        el = sorted({v for v in e if v in pos})
        color = _tab20_color(i)
        pts = [pos[v] for v in el]
        if len(el) == 0:
            continue
        if len(el) == 1:
            x, y = pts[0]
            circ = Circle((x, y), 0.06, facecolor=color, edgecolor=color, linewidth=1.5, alpha=0.85, zorder=3)
            ax.add_patch(circ)
        elif len(el) == 2:
            ax.plot(
                [pts[0][0], pts[1][0]],
                [pts[0][1], pts[1][1]],
                color=color,
                solid_capstyle="round",
                linewidth=3.5,
                alpha=0.75,
                zorder=2,
            )
        else:
            hull = _convex_hull_2d(list(pts))
            if len(hull) == 2:
                ax.plot(
                    [hull[0][0], hull[1][0]],
                    [hull[0][1], hull[1][1]],
                    color=color,
                    linewidth=3.5,
                    alpha=0.75,
                    zorder=2,
                )
            else:
                poly = Polygon(
                    hull,
                    facecolor=color,
                    edgecolor=color,
                    linewidth=1.4,
                    alpha=0.33,
                    zorder=1,
                )
                ax.add_patch(poly)

    for v in verts:
        x, y = pos[v]
        ax.scatter([x], [y], s=220, c="#1a1a1a", edgecolors="white", linewidths=1.2, zorder=5)

    ax.set_title("超圖（2-section 彈簧佈局 + 超邊凸包，Wolfram 風）")
    plt.tight_layout()
    return fig_to_png_bytes()


def plot_hypergraph(
    config: dict[str, Any],
    *,
    style: str = "spatial",
    layout_seed: int = 42,
) -> bytes:
    """
    依樣式繪製超圖。

    Args:
        config: 配置字典。
        style: ``\"incidence\"`` 二部圖；``\"spatial\"`` Wolfram 風空間嵌入。
        layout_seed: 僅 ``spatial`` 使用。

    Returns:
        PNG 位元組。
    """
    if style == "incidence":
        return plot_hypergraph_incidence(config)
    return plot_hypergraph_spatial(config, layout_seed=layout_seed)


def plot_largest_class_bars(analysis: dict[str, Any]) -> bytes:
    """
    靜態實驗：依 ``largest_classes`` 繪製前若干等價類大小長條圖。
    """
    sizes: Sequence[int] = analysis.get("largest_classes") or []
    if not sizes:
        plt.figure(figsize=(5, 3))
        plt.text(0.5, 0.5, "無類別大小資料", ha="center", va="center")
        plt.axis("off")
        return fig_to_png_bytes()

    plt.figure(figsize=(6, 4))
    xs = list(range(len(sizes)))
    plt.bar(xs, sizes, color="#45b7d1")
    plt.xlabel("排序後之等價類索引（由大至小）")
    plt.ylabel("類別大小 |C|")
    plt.title("前 10 大等價類大小")
    return fig_to_png_bytes()


def plot_entropy_time_series(
    series: Sequence[float],
    title: str = "解析熵 H（比特）沿時間步",
    *,
    plateau_epsilon: float | None = None,
) -> bytes:
    """
    動力學實驗：熵序列折線圖。

    Args:
        series: 各步之熵（bit）。
        title: 圖標題。
        plateau_epsilon: 若設定，於序列平均高度繪製 ±ε 之淺色帶，對應平台判定量級（§10.7）。
    """
    plt.figure(figsize=(7, 4))
    ys = list(series)
    plt.plot(range(len(ys)), ys, marker="o", markersize=2, linewidth=1)
    if plateau_epsilon is not None and ys:
        m = statistics.mean(ys)
        e = float(plateau_epsilon)
        plt.axhspan(m - e, m + e, alpha=0.2, color="0.5", label=f"±ε_plat={e:g}")
        plt.legend(loc="upper right", fontsize=8)
    plt.xlabel("時間步 t")
    plt.ylabel("H (bits)")
    plt.title(title)
    plt.grid(True, alpha=0.3)
    return fig_to_png_bytes()
