# CIT — 約束世界論 · 第十章超圖實驗

本專案將《約束世界論》**第十章**之有限**超圖 toy model** 實作為可重現程式：域型約束篩選、解析簽名與等價類、靜態指標、簡易動力學，以及 **§10.4** 之解析**細化鏈**、投影 **π<sub>Λ′→Λ</sub>**、纖維**條件核 K<sub>s</sub>** 與推前一致性。目標是**實作約束、觀察收束、比較湧現結構**，而非擬合連續物理方程。

## 專案結構

| 路徑 | 說明 |
|------|------|
| [`experiment.py`](experiment.py) | 命令列入口，輸出 JSON |
| [`hypergraph_experiment/core.py`](hypergraph_experiment/core.py) | 超圖配置、約束、簽名、等價類、動力學、`run_full_experiment` |
| [`hypergraph_experiment/refinement.py`](hypergraph_experiment/refinement.py) | §10.4：細化、π、K_s、雙路徑比較 |
| [`hypergraph_experiment/storage.py`](hypergraph_experiment/storage.py) | 實驗封存（`manifest.json`、圖檔路徑） |
| [`hypergraph_experiment/viz.py`](hypergraph_experiment/viz.py) | 圖表與超圖二部圖視覺化 |
| [`streamlit_app.py`](streamlit_app.py) | 網頁控制台（表單、圖表、封存） |
| [`pages/2_Experiment_library.py`](pages/2_Experiment_library.py) | 第二頁：已封存實驗列表、比對表、熵曲線重疊、CSV／ZIP |
| [`requirements.txt`](requirements.txt) | Streamlit／繪圖／NetworkX／pandas |
| [`約束世界論 29.md`](約束世界論%2029.md) | 論文稿本參考（含第十章參數與定義；檔名視您本機版本而定） |
| [`README_CH10_FIELDS.md`](README_CH10_FIELDS.md) | §10.2–10.9 表格欄位與論文語義對照、CSV 表頭說明 |

多頁應用之分頁腳本位於 [`pages/`](pages/)（檔名以數字前綴排序）；實際側邊欄標題以 Streamlit 顯示為準。

## 環境與安裝

建議使用專案內虛擬環境（Windows PowerShell 範例）：

```powershell
cd D:\code\CIT
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

核心運算（`experiment.py`、僅 import `hypergraph_experiment.core`）以**標準函式庫**為主；GUI 與圖表額外需要 `requirements.txt` 中的套件。

## 命令列實驗

靜態分析（預設）：

```powershell
python experiment.py --n 5 --max-edge-size 3 --max-edges 4 --sample-limit 2000
```

動力學：

```powershell
python experiment.py --mode dynamics --n 5 --steps 30 --runs 50
```

**§10.4 細化分析**（附加於 `static` 結果中的 `refinement_10_4` 欄位）：

```powershell
python experiment.py --refinement --refine-coarse-sig weak --refine-coarse-delta 3 `
  --refine-fine-sig medium --refine-fine-delta 0 --refine-kernel uniform
```

可加上 `--no-refine-compare-chains` 跳過 §10.4.1 雙路徑 JS／熵差。完整旗標說明請見 `experiment.py` 內 `build_parser()`。

## Streamlit 網頁介面

建議啟動方式（語意最清楚）：

```powershell
streamlit run streamlit_app.py
```

若誤用 `python streamlit_app.py`，腳本會嘗試自動改為以 `streamlit run` 啟動。更細的說明、倉儲路徑與字體提示見 **[README_STREAMLIT.md](README_STREAMLIT.md)**。

主控台另內建（見 [`streamlit_app.py`](streamlit_app.py)）：

- **枚舉規模預估**：`sample_limit=0` 時之候選配置總數（避免誤觸過大枚舉）。
- **表 10-2**：域型約束梯子；欄位為中文無符號（合法配置數、佔候選比例、相對候選收縮率、相對上一層收縮率／保留率、累計排除數、本層新增排除數、子集與鏈式關係等）。詳見 [README_CH10_FIELDS.md](README_CH10_FIELDS.md)。
- **表 10-3**：同一組可採用配置上比較 weak／medium／strong；欄位含解析單元數、解析壓縮比、解析熵位元等（中文表頭）。
- **匯入 JSON**：上傳 `result.json` 或 `manifest.json` 以檢視摘要並與當前實驗參數差異比對。

批次邏輯實作於 [`hypergraph_experiment/paper_tables.py`](hypergraph_experiment/paper_tables.py)；`sample_candidates_and_filter` 見 [`hypergraph_experiment/core.py`](hypergraph_experiment/core.py)。

## 實驗資料儲存

- 預設目錄：`experiments_data/<UTC時間>_<短ID>/`
- 內容包含 `manifest.json`、`sample_hypergraphs.json`、可選 `figures/*.png`
- 可設定環境變數 **`CIT_EXPERIMENTS_ROOT`** 指向自訂根目錄

## 與論文之對應（摘要）

- **表 10-1**：節點數 n、超邊限制、d_max、δ、軌道長度、偽隨機基底種子（seed）等 — 對應 CLI／Streamlit 表單；§10.2.6 之 **種子批次數**（\(N_{seed}\)）僅用於表 10-2 候選抽樣多輪，與基底種子不同。
- **§10.2–10.3**：配置域、簽名距離、等價類、壓縮比與熵 — 對應 `analyze_static`／`analyze_dynamics`；表 10-2／10-3 之資料與欄位定義見 `paper_tables.py` 與 [README_CH10_FIELDS.md](README_CH10_FIELDS.md)；Streamlit 匯出之 CSV 表頭與畫面一致。
- **§10.4**：Λ→Λ′、π、K_s、推前一致性、雙路徑比較 — 對應 `refinement.py` 與 `--refinement`。

## 授權與引用

若於學術脈絡使用，請依您持有的《約束世界論》版本自行註明出處；本程式碼為對應第十章之實作範例，不取代論文證明。

---

**硬體備註**：實驗以 CPU 枚舉／抽樣為主，無需 GPU。
