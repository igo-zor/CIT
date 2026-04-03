# 第十章實驗欄位與論文對照

本文件說明 Streamlit 分頁與 CSV 匯出之**顯示欄名**如何對照《約束世界論》第十章敘事，以及程式內資料來源。

## 偽隨機基底種子與論文 \(N_{seed}\)

- **偽隨機基底種子**（程式鍵 `seed`，介面常標為 `seed — …`）：單一整數，用於固定 RNG 狀態、重現抽樣與動力學。
- **論文 \(N_{seed}\)（種子批次數）**：僅在 **§10.2** 表 10-2 候選生成脈絡使用（程式鍵 `num_seeds`）；表示以 `seed + i·1009` 派生之多輪抽樣，與基底種子**不是**同一參數。
- 批次表欄名「偽隨機基底種子」若缺列，舊版 session 之「隨機種子」仍可由 `batch_cell` 自動相容。

## §10.3：候選上限 \(N_{cand}\) 與觀測集 \(N_{cfg}\)

- **`sample_limit`（介面：候選採樣上限）**：論文 **\(N_{cand}\)** 之操作化——候選母集生成／抽樣之上限；經域型篩選後得 **\(|\mathrm{Cfg}_\Lambda|\)**（程式鍵 `num_admissible_configs`）。
- **`n_cfg`（介面：輸入配置數（§10.3）；合法樣本數（§10.5））**：論文 **\(N_{cfg}\)**——自 **\(\mathrm{Cfg}_\Lambda\)** **不重複**抽樣得出觀測集 **\(\mathcal C_{\mathrm{obs}}\)**（程式鍵 `num_obs_configs`）；靜態／動力學分析僅在此集合上運算。若可採用配置數小於請求值，程式採用**全部**可採用集並於 `n_cfg_notice` 提示。
- **§10.5** 之子樣本使用專用 RNG 鏈節偏移（`OBS_SUBSAMPLE_RNG_OFFSET_10_5`），與 §10.3 之 `subsample_obs_configs` 預設偏移解耦，避免同一 `seed` 下兩章觀測集序列重疊。
- **`n_rep`（介面：重複次數）**：論文 **\(N_{rep}\)**——在固定其餘參數下重抽觀測集之次數（僅 static）；輸出 `analysis_repetitions`（逐次）與 `analysis_rep_summary`（平均/標準差）供穩健性檢視。
- **壓縮比**：`analyze_static` 之 **`compression_ratio_U`** 對齊正文 **\(U_\Lambda = |S_{\Lambda,\delta}| / N_{cfg}\)**（此處分母為觀測集大小，即 `len(configs)` 或 `num_obs_configs`）。舊版實作曾為 **\(|\mathcal C_{\mathrm{obs}}|/|S|\)**（平均壓縮倍數），已改為與論文同向。

## 表 10-2（配置域梯子）

資料由 [`hypergraph_experiment/paper_tables.py`](hypergraph_experiment/paper_tables.py) 的 `table_10_2_domain_ladder` 產生；欄名已採**繁體中文、不含符號**（與論文「合法配置數、保留比例、被排除構型類型數、包含關係、配置域收縮率」之語義對齊；程式輔以額外衍生欄）。

| 顯示欄名 | 語義說明 |
| ---------- | ---------- |
| 域型約束層級 | 該列對應之約束梯子層（標籤可由使用者自訂） |
| 合法配置數 | 該層域型約束下，自候選母集篩出之可採用配置個數 |
| 佔候選比例 | 合法配置數 ÷ 候選母集大小（對應論文保留比例之操作化） |
| 相對候選收縮率 | 1 減去佔候選比例（相對整體候選母集之收縮） |
| 累計排除數 | 候選母集大小減去該層合法配置數 |
| 本層新增排除數 | 相對上一層合法集合，於該層多排除之個數 |
| 是否為上一層子集 | 該層合法配置集合是否為上一層之子集 |
| 鏈式子集成立 | 自弱至強逐層之子集關係是否全程成立 |
| 相對上一層保留率 | 該層合法配置數 ÷ 上一層合法配置數（第一層無上一層時為空值） |
| 相對上一層收縮率 | 1 減去相對上一層保留率（論文配置域收縮率之層間版；第一層為空值） |
| 該層排除主因…筆數 | 對**整個候選母集**中本層遭排除之配置，依「首個未通過條件」分類計數（超邊數、超邊大小、度數、二階連通、禁二元三角） |
| 本層新增排除主因…筆數 | 僅上一層仍合法、本層始遭排除之配置之主因計數（對應 §10.2.7 第二點之「結構性收縮」操作化起點） |
| 違規主因集中度／本層新增違規主因集中度 | 各主因筆數中最大者所占比例（0–1；無排除時為空） |

**觀察重點**（對照論文 §10.2.7）：合法配置數與保留比例是否隨約束增強單調下降；鏈式子集是否穩定成立；參數小幅變動時趨勢是否仍大致保持。第二點（排除是否集中於特定違規模式而非均勻散落）可搭配**主因分布與集中度**解讀；本實作以域型梯子已定義之五類條件為「主因」標籤，**未**細分論文中更細的 motif 語彙（局部過密等需擴充標籤後才能對應）。

## 表 10-3

同檔 `table_10_3_signature_comparison`；欄名包含：**解析簽名**、**解析單元數**、**平均單元大小**、**解析壓縮比**、**重疊率**、**相容孤立率**（對應論文 **R_iso** 之操作化）、**傳遞違反率**、**解析熵位元**。

產生表 10-3 時應傳入與 §10.3 單次實驗相同之 **`n_cfg`** 與 **`seed`**：先對可採用配置全集做觀測集子抽樣，再於**同一** \(\mathcal C_{\mathrm{obs}}\) 上對 weak／medium／strong 呼叫 `analyze_static`。

### §10.3 靜態 JSON 鍵（`analyze_static`）

| 鍵名 | 說明 |
| ------ | ------ |
| `delta` | **整數**簽名距離閾值；與 `signature_distance` 之離散刻度一致。論文若以較細小數標度敘述，與本實作僅為尺度表述差異，以本專案數值為準。 |
| `s_min` | 重疊率計算時僅納入 **\|T(c)\| ≥ s_min** 之配置參與配對平均；`0` 表示不過濾（預設、向後相容）。 |
| `isol_rate_compat_graph` | **R_iso**：在「不同配置且簽名距離 ≤ δ」之相容關係圖上，無此类鄰邊之節點占比。 |
| `overlap_rate` | 鄰域對 Jaccard 風格重疊之平均（受 `s_min` 影響）。 |
| `compression_ratio_U` | **\(U_\Lambda = \|S\|/N_{cfg}\)**（觀測集大小為分母）；`avg_class_size` 為算術平均類別大小，恒滿足 \(\sum_i \|C_i\| = N\)，故平均為 \(N/\|S\|\)。 |
| `num_configs` | 本列分析所使用之配置筆數（即為觀測集 \(\|\mathcal C_{\mathrm{obs}}\|\)）。 |

### §10.7 動力學（`analyze_dynamics`）

| 鍵名 | 說明 |
| ------ | ------ |
| `r_adm_mean`／`r_adm_time_series_mean` 等 | 論文 §10.7.2（四）之合法更新率：每步提出 $M_{\mathrm{trial}}$ 候選，通過合法性之比例 $M_{\mathrm{adm}}/M_{\mathrm{trial}}$，並對步／軌跡彙總。 |
| `legal_update_step_fraction_mean` | 各軌跡相鄰步狀態是否相異之比例，再對軌跡平均；保留作為語義相近之代理量（向後相容）。 |
| `epsilon_plat` | 熵時間序列平台判定：相鄰 \|ΔH\| ≤ ε 視為同段。 |
| `entropy_summary.plateau_*` | 逐步聚合熵 $H_\Lambda^{(\ell)}$ 平台段摘要。 |
| `entropy_time_series_wH`／`entropy_summary_wH` | 窗口熵 $H_{\Lambda,w_H}^{(\ell)}$ 與其平台段摘要（§10.7.2（六））。 |
| `p_cycle_summary`／`ell_a_summary` | 週期長度 $P_{\mathrm{cycle}}$ 與吸引子進入時刻 $\ell_A$（以解析類別標籤序列與熵差分操作化；§10.7.2（八））。 |
| `n_reach_mean`／`n_reach_per_run` 等 | §10.7.5（五）$N_{\mathrm{reach}}$：各軌跡上相異解析單元（類別索引）個數，再對軌跡取平均／最小／最大。 |

Streamlit 頁 `pages/10_07_Dynamics_10_7.py` 已拆分**實驗 A**（單路徑，§10.7.5 主線）與**實驗 B**（雙路徑終端比較，擴充）。實驗 A 結果區除扁平 CSV 外，另有 **§10.7.5 論文輸出參數對照表**（[`hypergraph_experiment/ch10_section_10_7_tables.py`](hypergraph_experiment/ch10_section_10_7_tables.py)）；其中 $\bar r_{\mathrm{adm}}$ 已對應 `r_adm_mean`（論文定義），並保留 `legal_update_step_fraction_mean` 作為向後相容之代理量。介面欄位 $M_{\mathrm{trial}}$、$w_H$、$w_A$、$P_{\max}$、$N_{\mathrm{seed}}$（程式鍵 `m_trial`、`w_h` 等）已接線至動力學模擬與彙總，確保參數掃描能得到一致且可重現之結果。

## §10.3–10.9 單次與批次扁平原表

扁平化鍵名（如 `r3s_parameters_n`）經 [`hypergraph_experiment/streamlit_common.py`](hypergraph_experiment/streamlit_common.py) 的 `build_ch10_column_name_map` 對應為**中文顯示欄名**；[`render_table_with_copy_csv`](hypergraph_experiment/streamlit_common.py) 會將同一組映射同時套用到畫面、一鍵複製 Markdown 與下載 CSV，確保表頭一致。

自本版起，各實驗結果區固定追加 **輸入實驗參數表（固定參數＋變數）**，並與論文符號對照：

- 單次執行：由 `result.parameters`（或頁面等價參數來源）產生參數表，欄位含 `參數鍵 / 論文記號 / 論文語義 / 參數值`，並同樣支援一鍵複製 Markdown 與 CSV 下載（`render_parameters_table`）。
- 批次逐組：`run_batch_per_run_rows` 會保留每列原始參數（`param_row`），`render_batch_per_run_tables` 會在每組結果前自動渲染該組參數表，便於逐組實驗追溯與記錄。

內部 JSON/API 鍵名（例如 `run_full_experiment` 回傳結構）不因顯示層改名而改動；比對程式碼或封存檔時請以 JSON 鍵為準。

## 批次參數表欄名（§10.02–10.09）

各頁 `st.data_editor` 預設採**繁體中文、無符號**欄名（與論文建議參數語義對齊），列執行時經 `batch_cell` 對應回程式引數；若 session 中仍保留舊版英文欄名，會自動回退讀取。

## 論文建議 preset 對齊（§10.2–§10.9）

頁面單次預設與「載入論文建議批次模板」已對齊《約束世界論 29–30.md》第十章主線建議；共用基線常數集中於
[`hypergraph_experiment/ch10_paper_presets.py`](hypergraph_experiment/ch10_paper_presets.py)。

- §10.3：`N_cfg=300`、`s_min=2`（並保留 weak/medium/strong 與 `delta` 掃描）。
- §10.4：`N_cfg^Λ/N_cfg^Λ'=2000/2000`、`ε_push*=0.01`、`ε_JS*=0.01`；批次模板第一列為主線，第二列保留論文建議鏈比較。
- §10.5（v30 §10.5.10）：`N_{cand}=5000`、`N_{cfg}=2000`、`k_{min}=2`、`k_{max}=3`、`m=16`、`n_A=n_B=6`、`alpha_cross=0.30`；頁面可一鍵載入 **§10.5.4（一）** 之離散 \(\alpha\) 掃描列。輸出參數小節編號為 **§10.5.5**。
- §10.6：對齊 §10.6.6 建議表之 **`n=8,M=4,w=2,η=1,N_ctx=1000,T_loc=2,N_search=5000`**（程式鍵 `w_ctx`／`eta_ctx`）；模式 obstruction／satisfiable；輸出含 **`n_val`**（§10.6.5（四）分子）。
- §10.7：`CH10_7_BASELINE` 含 §10.7.6 主線（含 `m_trial`、`w_h`、`w_a`、`p_max`、`n_seed_107` 等記錄用欄位）與 `epsilon_plat=0.01`；批次掃描歸屬實驗 A；`CH10_7_T_CHOICES`／`CH10_7_N0_CHOICES` 對應 §10.7.4 建議離散集合。
- §10.8：主線 `eta=0.10`、`T_sb=20`，並保留 `eta` / `T_sb` 掃描列。
- §10.9：主線視窗清單 `1,2,4,8,16`，批次模板保留 `{1,2,4,8,16}` 與保守版 `{1,4,8}`。

名詞映射提醒：

- `seed`：偽隨機**基底**種子（固定 RNG）。
- `N_seed`：種子批次數（僅 §10.2 候選抽樣輪數）；兩者不同。

## §10.8／§10.9 補充

- **§10.8**：`metrics` 含 `N_iso_mean`（球面規範簽名去重之代理）、`A_reach_mean`（2-section 可達對之最短路平均）、`A_reach_pair_fraction_mean`。
- **§10.9**：軌跡層 `R_edge_bar`（邊周轉率）；各視窗列含 `tau_unit_*`（宏觀標籤游程）、`L_plat_*`（有效邊數序列之平台摘要）。

## §10.4 纖維核

`refinement.apply_refinement_step` 之 `detail.kernel_stability`：各纖維上條件核之平均非零比例與平均熵（bit），供觀察核之尖銳度。

§10.4 頁面可調實驗層控制（雙實驗獨立）：

- `refine_coarse_sample_size`、`refine_fine_sample_size`：粗/細層解析所用樣本數上限（細層取自粗層子集）。
- `refine_fiber_sample_size`：每個粗單元纖維中納入核估計之細層索引上限（固定排序截斷，確保可重現）。
- `epsilon_push_threshold`、`js_threshold`：僅供頁面標註「通過/警示」判準，不改變 `refinement.py` 內部數值計算定義。
- 主表 preset 對齊論文固定組合：`weak→medium` 與 `medium→strong`。
- 主表實驗（10.4.5/10.4.6）與 A/B 實驗（10.4.7）分為兩個按鍵與兩組 state：`res_10_4_main`、`res_10_4_ab`。
- **共用定義域 D**：實驗 B 與主表共用頁面上方之基礎域參數（`n`、`k_max`、`m_max`、`sample_limit`、`seed`、連通／度數／禁二元三角等）；僅執行按鍵與結果 state 分離。
- **A/B 參數拆分**：實驗 B 表單固定分為兩組獨立輸入——`step_a`（細化 A 參數）與 `step_b`（細化 B 參數），兩者各自包含 `kind` 與 `delta`，不共用子步驟欄位。
- **實驗 B 輸入語義（對齊 9.6-C）**：以**子步驟映射規格** `RefinementSubstepSpec`（`kind` + `delta`）指定 **R_A、R_B**，而非將 `weak`／`medium`／`strong` 當成 A/B 的「主輸入欄」。Preset 僅將主表粗/細層標籤轉成預設 `kind`（見 `preset_layer_to_substep_spec`）。
- `kind` 與內部分割簽名之固定對照（`substep_spec_to_signature_delta`）：`edge_scale_motif`→weak、`degree_split`→medium、`adjacency_motif_fine`→strong。
- A/B 可切換「手動覆寫」以直接選 `kind` 與整數 **δ（相容閾值）**。
- 進階比較表輸出鏈別推前誤差、終端熵、纖維核穩定度與跨鏈終端 JS/|ΔH|。

10.4 頁「論文 10.4.5 對齊主表」欄位對照：

- `粗層熵 H(p_Λ)` ← `single_step_Lambda_to_Lambda_prime.H_p_Lambda_bits`
- `細層熵 H(p_Λ')` ← `single_step_Lambda_to_Lambda_prime.H_p_Lambda_prime_bits`
- `推前誤差 ε_push` ← `single_step_Lambda_to_Lambda_prime.pushforward_max_error`
- `終端 JS 差異 JS_term` ← `two_path_ab_ba_10_4_7.js_divergence_bits_terminal`
- `終端熵差 ΔH_term` ← `two_path_ab_ba_10_4_7.entropy_abs_diff_terminal`
- `可交換性判定` ← 同時檢查 `ε_push <= epsilon_push_threshold` 與 `JS_term <= js_threshold`

10.4 頁「A/B 進階比較表 + 跨鏈摘要」對應論文 10.4.7：

- 鏈列（A→B、B→A）：`推前誤差`、`終端熵 bit`、`纖維核平均非零比例`、`纖維核平均熵 bit`
- 路徑識別鍵：`A_to_B.path_key="A→B"`、`B_to_A.path_key="B→A"`；各鏈結果內嵌 **`step_a`**、**`step_b`**（`{"kind": SubstepKind, "delta": int}`）供追溯映射規格（已取代舊版 `step_a_signature` 等欄位）。
- 跨鏈摘要：`js_divergence_bits_terminal_ab_ba`（`JS_term`）與 `entropy_abs_diff_terminal_ab_ba`（`ΔH_term`）

## 批次表 §10.2

`pages/10_02_Domain_10_2.py` 批次列除各層指標外，會帶回該列之掃描參數，欄名經上述映射顯示為中文。

## 測試

`tests/test_ch10_experiments.py` 含表 10-2 欄位與 `build_ch10_column_name_map` 之煙霧／單元測試；變更欄名邏輯時請一併更新測試與本文件。
