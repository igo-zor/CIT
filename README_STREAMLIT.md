# 超圖實驗 Streamlit 控制台

本目錄實作《約束世界論》第十章之有限超圖 toy model：**命令列**維持 [`experiment.py`](experiment.py)；**網頁介面**以 Streamlit 提供參數輸入、圖表、超圖視覺化、實驗封存與多 run 比對。

## 環境（建議虛擬環境）

於專案根目錄 `d:\code\CIT`（或您本機之對應路徑）：

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

## 啟動（務必使用虛擬環境內的 streamlit）

啟動後瀏覽器會開啟主控台；側邊欄可切換至 **「實驗庫與比對」** 分頁（`pages/2_Experiment_library.py`）。

```powershell
.\.venv\Scripts\Activate.ps1
streamlit run streamlit_app.py
```

若不想啟用 shell 整合，可直接：

```powershell
.\.venv\Scripts\streamlit.exe run streamlit_app.py
```

命令列實驗亦請使用同一解譯器，以免裝在系統 Python 與 venv 不一致：

```powershell
.\.venv\Scripts\python.exe experiment.py --n 5 --sample-limit 500
```

## 實驗資料儲存位置

- 預設目錄：`experiments_data/<UTC時間>_<短UUID>/`
- 內容：`manifest.json`（參數 + 分析結果）、`sample_hypergraphs.json`、`figures/*.png`
- 可藉環境變數 **`CIT_EXPERIMENTS_ROOT`** 指定其他根路徑

## 與 CLI 對照

同一組參數下，`run_full_experiment` 之 JSON 結構應與下列指令一致：

```powershell
python experiment.py --n 5 --sample-limit 2000 --mode static
```

大型 `n` 或 `sample_limit` 可能阻塞 Streamlit 執行緒，建議改用 CLI 批次實驗後，將產出之 JSON 手動放入自訂目錄（進階）或縮小參數於網頁試算。

## Quadro ／ 顯示卡

本實驗以 CPU 枚舉與抽樣為主，**無需** GPU；Matplotlib 僅用於靜態圖輸出。

## 表格匯出與欄位命名（§10.2–10.9）

- 多數分頁使用 `render_table_with_copy_csv`（定義於 [`hypergraph_experiment/streamlit_common.py`](hypergraph_experiment/streamlit_common.py)）：可**一鍵複製 Markdown** 與**下載 CSV**，兩者表頭相同。
- 扁平原表欄位（如批次實驗之技術鍵名）會經 `build_ch10_column_name_map` 轉成**繁體中文、不含符號**之顯示欄名；內部 JSON／封存檔仍維持程式原有鍵名。
- 表 10-2、10-3 之資料欄位由 `paper_tables.py` 直接以中文鍵產出，並與論文 §10.2 觀察重點（配置域收縮、鏈式子集等）對齊說明。

完整欄位語義表請見專案根目錄 **[README_CH10_FIELDS.md](README_CH10_FIELDS.md)**。
