# enrichplot `nsea` / `mnsea` 可视化与新方法设计 Spec

## 1. 目标

这份文档讨论的是：

**`enrichplot` 应该如何完整消费 `nseaResult` / `mnseaResult`，并在复用现有 `gseaResult` 图形生态之外，进一步提供面向网络机制比较与解读的新方法。**

这份 spec 有两个并行目标：

1. 先把 `nseaResult` / `mnseaResult` 作为 `gseaResult` 子类可直接复用的图形方法尽量接齐。
2. 再围绕网络结构、leading-edge 重排、跨层一致性与重布线，设计只对 `nsea` / `mnsea` 才真正有意义的新方法。

## 2. 范围

这份 spec 关注：

- `nseaResult`
- `mnseaResult`
- `enrichplot` 中与 `gseaResult` 家族有关的现有图形方法
- 面向网络机制解读的新图形与 plot-facing helper

这份 spec 暂不关注：

- `enrichit` 中算法本身的修改
- `clusterProfiler` 的知识库接口
- 与网络机制无直接关系的泛型图形扩展

## 3. 设计前提

### 3.1 职责边界

`enrichplot` 的职责不是重新计算 `nsea` / `mnsea`，而是做：

- 可视化
- fortify / tidy conversion
- plot-facing helper
- explanation-ready 数据结构到图形语义的映射

换句话说：

- **计算在 `enrichit`**
- **对象语义入口在 `clusterProfiler` / `DOSE` 兼容层**
- **图形表达在 `enrichplot`**

### 3.2 基本判断

`nseaResult` 和 `mnseaResult` 虽然继承自 `gseaResult`，因此现有图形方法在技术上“很好接”，但两者并不只是“换了一个结果对象名字”。

关键差异在于：

- `nseaResult` 增加了网络扩散语义
- `mnseaResult` 增加了多层网络、collapsed score、layer score、pathway contribution、feature contribution、coupling 等 explanation-ready 结构

因此本 spec 采用两层策略：

1. **先尽量复用并接齐现有图形方法**
2. **再基于网络结构语义设计新方法**

## 4. 当前对象语义

### 4.1 `nseaResult`

当前 `nseaResult` 在 `gseaResult` 基础上额外包含：

- `network`
- `diffusion_scores`
- `mode`
- `iterations`
- `restart_prob`

它本质上仍然是单网络下的 enrichment result，因此大部分 `gseaResult` 图形语义可以直接复用。

### 4.2 `mnseaResult`

当前 `mnseaResult` 在 `gseaResult` 基础上额外包含：

- `multilayer_network`
- `layer_scores`
- `collapsed_scores`
- `layer_weights`
- `coupling_table`
- `mode`
- `iterations`
- `restart_prob`
- `collapse_method`
- `target_layer`
- `output_space`
- `pathway_contribution`
- `feature_contribution`

这意味着 `mnseaResult` 不只是“多几个注释列的 `gseaResult`”，而是已经具备：

- 层级分解语义
- 路径驱动来源语义
- feature-level explanation 语义
- 多层耦合语义

因此 `mnseaResult` 需要专门的 helper 和新图形语义。

## 5. 总体策略

### 5.1 第一阶段：先把继承自 `gseaResult` 的图接齐

原则上，凡是 `gseaResult` 的图形方法仍然成立，且 `nseaResult` / `mnseaResult` 的新增结构不会破坏原语义，都应该优先直接接入。

优先目标不是“发明更多名字”，而是：

- 让用户先能无缝进入既有工作流
- 让 `nsea` / `mnsea` 在 `enrichplot` 里变成一等对象

### 5.2 第二阶段：再做网络机制专属方法

当基础图形工作流已经接齐后，`enrichplot` 应提供新的图形家族，回答以下问题：

- 同一个 pathway 在不同网络 / 不同 layer 下，到底是机制保守还是重布线？
- pathway 名字不变时，leading-edge feature 是否已经换了？
- 变化来自富集强度，还是来自网络 wiring 改组？
- 多层网络中，到底是哪一层在驱动 term，哪一层只是跟随？

## 6. 当前 enrichplot 支持面

### 6.1 已明确支持 `mnseaResult` 的图

当前已经补上的 `mnseaResult` 图包括：

- `dotplot()`
- `heatplot()`
- `cnetplot()`
- `emapplot()`
- `gseaplot()`
- `ridgeplot()`
- `upsetplot()`
- `ssplot()`

同时已经有相应的 plot-facing helper：

- `fortify.mnseaResult()`
- `fortify_mnsea_contribution()`
- `fortify_mnsea_subnetwork()`

### 6.2 仍待补齐的继承型图

当前仍建议补齐或确认的图包括：

- `gseaplot2()`
- `gsearank()`
- `hplot()`
- `treeplot()`
- `barplot()` 的 `nsea` / `mnsea` 兼容性确认与测试

其中优先级建议为：

1. `gseaplot2()`
2. `gsearank()` / `hplot()`
3. `treeplot()`
4. `barplot()` 兼容性确认

这些“待补”图的准确状态是：

- `gseaplot2()` / `gsearank()` / `hplot()` 目前都是普通函数，内部调用普通函数 `gsInfo()`，
  **不是 S3 泛型**。它们对 `nseaResult` / `mnseaResult` 不会自动进入
  `gsInfo.mnseaResult()`；因此补齐时首先要把 `gsInfo()` 改造成 dispatch，
  否则默认路径只会取 `geneList` slot。注意一个缓解性事实：`mnseaResult@geneList`
  通常就是 collapsed scores，所以 `layer = NULL` 的默认语义不会立刻画错，
  但单 layer 视图、`geneSetID` 解析、`Description` 列等仍需显式补上。
- `treeplot()` 已有 `setMethod("treeplot", signature(x = "gseaResult"))`，
  S4 dispatch 会命中 `nseaResult` / `mnseaResult` 子类，因此它不是“未接入”。
  对 `nseaResult` 主要是补测试确认；对 `mnseaResult` 的缺口是 layer-aware
  相似性语义，根因是 `pairwise_termsim` 还没有 `mnseaResult` 方法。
- `barplot()` 只有 `barplot.enrichResult` 与 `barplot.compareClusterResult`
  两个 S3 方法，没有 `barplot.gseaResult` / 子类方法；`nseaResult` /
  `mnseaResult` 继承的是 `gseaResult` 而不是 `enrichResult`，因此会一路落到
  `graphics::barplot.default`。补兼容性的本质是增加 gseaResult 家族方法或
  等价别名，而不是仅仅“跑一下确认”。

## 7. 设计原则

### 7.1 先复用，再专门化

如果一个图本来就是针对 `gseaResult` 的单一路径或多路径可视化，且 `nseaResult` / `mnseaResult` 在不引入额外歧义的前提下可以沿用，就先接入。

### 7.2 新方法必须带来新信息

新方法不应该只是“老图换个参数名”，而应该提供现有 `gseaResult` 图看不到的信息，例如：

- rewiring 程度
- layer-specific 驱动
- feature-level 重排
- 保守机制 vs 机制转向

### 7.3 helper 优先于绘图细节

所有新图形都应建立在稳定的 plot-facing helper 之上，而不是让每个图自己直接拆 slot。

也就是说，应优先定义：

- term-level summary helper
- pathway-level contribution helper
- feature-level rewiring helper
- mechanism classification helper

再让上层图形函数消费它们。

### 7.4 不把 plotting concern 带回计算层

不建议把以下内容塞回 `enrichit`：

- node color
- node shape
- panel splitting
- legend title
- label density

这些都属于 `enrichplot` 的职责。

## 8. 数据契约

### 8.1 复用型图的数据契约

对于直接复用 `gseaResult` 语义的图，数据契约应尽量沿用：

- `result`
- `geneSets`
- `geneList`
- `permScores`
- `termsim`

只是对 `mnseaResult` 增加：

- `layer = NULL` 时使用 `collapsed_scores`
- `layer = "<single-layer>"` 时使用对应 `layer_scores[[layer]]`

关于 `geneList` 的语义需要写清楚：

- `nseaResult@geneList` 与 `mnseaResult@geneList` 分别来自各自 pipeline 的
  enrichment 输入：`nsea` 为扩散后的 ranked scores、`mnsea` 为 collapsed scores。
  因此 `layer = NULL` 的默认复用路径在语义上就是正确的 collapsed 视图；
  它们不是“原始表达排序”，不要把 `geneList` 当成与扩散无关的输入去二次处理。
- `layer = "<single-layer>"` 必须显式走 `get_mnsea_ranked_scores()` 一类入口，
  不能只读 `geneList` slot。

### 8.2 `mnsea` explanation helper 契约

`mnseaResult` 现有 helper 应被视为公共 plot-facing 契约的一部分：

#### `fortify.mnseaResult(level = "result")`

- 输出常规 term-level 结果
- 用于 `dotplot()`、`barplot()` 等总览图

#### `fortify.mnseaResult(level = "pathway")`

- 输出 pathway × layer 级别结果
- 用于 layer-aware summary 图

#### `fortify_mnsea_contribution(level = "pathway")`

- 标准列至少包括：
  - `ID`
  - `Description`
  - `layer`
  - `contribution`
  - `share`
  - `n_feature`

#### `fortify_mnsea_contribution(level = "feature")`

- 标准列至少包括：
  - `ID`
  - `Description`
  - `Feature`
  - `layer`
  - `score`
  - `abs_score`
  - `sign`
  - `is_core`

#### `fortify_mnsea_subnetwork()`

- 输出 pathway-specific node / edge tables
- 至少区分：
  - `pathway`
  - `feature`
  - `membership`
  - `intra`
  - `coupling`

## 9. 方法矩阵

| 方法 | `nseaResult` | `mnseaResult` | 策略 | 说明 |
|---|---|---|---|---|
| `fortify()` | 直接复用 | 扩展 | mixed | `mnsea` 需要 `level` 语义 |
| `dotplot()` | 直接复用 | 已扩展 | mixed | `mnsea` 支持 pathway contribution |
| `barplot()` | 待补 | 待补 | reuse | 缺 `barplot.gseaResult`，先加方法再补测试 |
| `heatplot()` | 直接复用 | 已专门化 | extend | pathway / feature 两种视图 |
| `cnetplot()` | 直接复用 | 已专门化 | extend | pathway-specific subnetwork |
| `emapplot()` | 直接复用 | 已扩展 | extend | layer-aware similarity |
| `ssplot()` | 直接复用 | 已扩展 | extend | similarity-space overview |
| `gseaplot()` | 直接复用 | 已扩展 | extend | collapsed / single-layer |
| `gseaplot2()` | 待补语义 | 待补语义 | reuse | 需先让 `gsInfo()` dispatch |
| `gsearank()` | 待补语义 | 待补语义 | reuse | 需先让 `gsInfo()` dispatch |
| `hplot()` | 待补语义 | 待补语义 | reuse | 需先让 `gsInfo()` dispatch |
| `ridgeplot()` | 直接复用 | 已扩展 | extend | collapsed / single-layer |
| `upsetplot()` | 待确认 | 已扩展 | extend | `mnsea` 已采用 pathway overlap 语义 |
| `treeplot()` | 待确认 | 待补语义 | extend | 已有 `gseaResult` 方法可选；缺 `pairwise_termsim` 的 `mnsea` 语义 |

## 10. 拟新增的新方法家族

下面这些方法不是为了“补齐已有函数”，而是为了给 `nsea` / `mnsea` 用户提供真正新的机制信息。

### 10.1 `phaseplot()`

#### 核心问题

一个 term 的变化，到底来自：

- 富集强度变化
- 还是网络结构重布线

#### 图形定义

把每个 pathway / module 放到二维平面：

- `x` 轴：enrichment shift
  - 例如 `NES` 或 `delta_NES`
- `y` 轴：rewiring score
  - 例如 edge retention 的补值、中心性漂移、leading-edge 重排程度
- 点大小：leading-edge overlap 或 `n_feature`
- 点颜色：显著性、保守性类别或 layer dominance

#### 价值

这张图可以把结果分成：

- signal changed, mechanism conserved
- signal stable, mechanism rewired
- signal strengthened and rewired
- weak and unstable

#### 适用对象

- `nseaResult`
- `mnseaResult`

但 `mnseaResult` 的价值更高。

#### 数据来源

`phaseplot()` 的 x 轴（enrichment shift）来自结果表的 `NES` 或跨 context 的
`delta_NES`；y 轴（rewiring score）必须由 `compute_rewiring_score()` 统一返回，
不允许各图层自行拼装。对于单个 `nseaResult` / `mnseaResult`，如果不存在跨
network / layer / condition 的比较上下文，y 轴应退化为 layer-to-collapsed 重布线
得分或显式声明该轴不可用，而不是静默画成无意义的值。

### 10.2 `consensusmap()`

#### 核心问题

在多个网络 / 多个 layer / 多个条件下，同一个 pathway 到底是：

- 保守
- 重布线
- 条件特异
- 矛盾

#### 图形定义

行是 pathway，列是 network / layer / condition。

每个格子至少表达两类信息：

- enrichment strength
- topology consistency

旁边再给 pathway 一个机制标签：

- `conserved`
- `rewired`
- `context_specific`
- `contradictory`

#### 价值

这是最适合作为 `mnsea` 总览主图的新方法之一。

#### 输入契约（必须先定义）

单个 `nseaResult` 只有一个 network、单个 `mnseaResult` 只有若干 layer，都没有
“多个 network / 多个 condition”的天然结构。因此 `consensusmap()` 的输入必须
显式接受一个**命名的结果列表**（例如
`list(networkA = res1, networkB = res2, ...)` 或
`list(conditionA = res1, conditionB = res2, ...)`），或对 `mnseaResult` 退化为
layer 列。没有比较上下文时不得产生不透明行为，建议直接报错并要求传入列表。

### 10.3 `rewireplot()`

#### 核心问题

同一个 pathway 名字不变时，底层驱动 feature 和边是否已经换线？

#### 图形定义

给定一个 pathway：

- 节点是 leading-edge feature
- 边是 pathway-specific subnetwork
- 用分面或颜色区分不同 network / layer / condition
- 突出：
  - 保守 feature
  - 新增 feature
  - 消失 feature
  - 桥接节点
  - coupling-mediated edges

#### 价值

它解决的是“同名 pathway 是否同机制”这个问题。

#### 适用对象与参考基准

`shared` / `gained` / `lost` / `shifted` 都需要一个明确的**参考基准**。第一版
只支持 `mnseaResult`（至少可以在 layer 之间、或 layer 与 collapsed 之间比较），
并强制要求显式传入 `reference_layer`；对 `nseaResult` 只有在提供第二个结果对象
或 reference 输入时才支持，不能在单网络结果上静默产生伪 rewiring 状态。

### 10.4 `mechanismflow()`

#### 核心问题

term 在多个网络 / 多个 layer 中的状态是如何迁移的？

#### 图形定义

把 pathway 状态定义为若干机制相位：

- conserved
- activated
- rewired
- decoupled
- collapsed

再用 sankey / river 方式展示 term 在不同条件下如何迁移。

#### 价值

它比单纯比较多个 `NES` 更接近机制演化叙事。

#### 输入契约

状态迁移至少需要两个 context，因此 `mechanismflow()` 同样需要命名的结果列表
或 `mnseaResult` 的 layer 序列。第一版允许只聚焦 `mnseaResult`（见开放问题 4），
但必须写清楚迁移序列来自哪里；不要假设单个对象隐含了“时间顺序”。

## 11. 新方法所需 helper

为了支持上述新方法，建议新增或显式稳定以下 helper。

### 11.1 `summarize_nsea_mechanism()`

输入：

- `nseaResult` 或 `mnseaResult`；如果用于跨 network / condition 比较，可接受命名列表并返回带 context 列的汇总。

输出 term-level summary table，至少包含：

- `ID`
- `Description`
- `NES`
- `p.adjust`
- `leading_edge_size`
- `leading_edge_overlap`
- `rewiring_score`
- `centrality_shift`
- `mechanism_class`

数据来源与退化行为必须显式说明：

- `leading_edge_overlap` / `rewiring_score` / `centrality_shift` 在对象里没有
  现成 slot，必须由 `compute_rewiring_score()` 等 helper 显式计算。
- 对 `nseaResult`，若没有第二个 network 或 reference，`leading_edge_overlap` /
  `centrality_shift` / `rewiring_score` 应返回 `NA` 或显式报参数错误，不能
  用恒等值（如 0）伪装成“未重布线”。
- 对 `mnseaResult`，`leading_edge_overlap` 优先来自 layer vs collapsed 比较，
  而不是重新发明一套既不同层也不同对象间的指标。

### 11.2 `extract_rewiring_features()`

输入：

- `x`
- `pathway_id`
- `layer` / `reference_layer`

输出 feature-level comparison table，至少包含：

- `Feature`
- `score`
- `abs_score`
- `sign`
- `status`
  - `shared`
  - `gained`
  - `lost`
  - `shifted`

`reference_layer` 是必需输入，不得用未文档化的“第一个 layer”作为默认值；
没有参考基准时应报错，而不是返回全部 `shared`。对 `nseaResult` 仅在提供第二个
结果对象或 reference 输入时可用。

### 11.3 `classify_mechanism_state()`

根据 enrichment shift 与 rewiring score，把 pathway 归类为：

- `conserved`
- `rewired`
- `context_specific`
- `contradictory`

### 11.4 `compute_rewiring_score()`

这是新方法家族里最核心的数值 helper，**必须先于其他 helper / 图层落地**。

第一版采用统一的透明定义：`rewiring_score = 1 - leading-edge Jaccard overlap`
（即 leading-edge overlap 的补值），范围 `[0, 1]`，`0` 表示完全保守、`1` 表示完全重布线；
`subnetwork edge overlap`、中心性漂移、coupling 占比变化等作为后续可选复合项，
由同一函数统一返回，不允许各图层各自计算。

候选定义可以来自：

- leading-edge overlap 的补值
- subnetwork edge overlap 的补值
- 节点中心性变化
- coupling edge 占比变化

最终应由一个统一函数负责返回。

## 12. API 设计建议

### 12.1 先保证新方法名清晰，而不是和旧图挤在一起

建议直接新增独立函数，而不是把所有新语义都塞进现有图的参数里。

推荐函数名：

- `phaseplot()`
- `consensusmap()`
- `rewireplot()`
- `mechanismflow()`

### 12.2 对 `mnseaResult` 的 layer 参数保持统一

建议所有涉及 `mnseaResult` 的新函数，尽量统一支持：

- `layer = NULL`
- `layer = "<single-layer>"`
- 必要时 `layer = c(...)`

但默认语义必须写清楚：

- 总览图默认使用所有 layer 的 collapsed 语义
- pathway-specific 机制图默认允许聚焦单 layer

### 12.3 pathway 选择规则保持稳定

所有 pathway-specific 图都不应让 `pathway_id = NULL` 产生不透明行为。

建议沿用当前 helper 的稳定默认选择规则：

- 数值型 `showCategory`
- 字符型 `ID`
- 字符型 `Description`

并避免仅按 `Description` 模糊匹配。

## 13. 向后兼容

### 13.1 继承型图尽量不破坏旧接口

对于 `gseaplot2()`、`gsearank()`、`hplot()`、`treeplot()` 等待补图，原则上：

- 旧参数名不变
- 对 `nseaResult` / `mnseaResult` 只是新增兼容性
- `mnseaResult` 仅在必要时增加 `layer` 等新参数

### 13.2 新方法不污染现有泛型

新方法家族建议以新函数名引入，而不是把现有 `gseaplot()` 或 `emapplot()` 继续塞进更多复杂参数。

## 14. 验证计划

### 14.1 继承型图补齐的验证

对于待补的 `gseaplot2()`、`gsearank()`、`hplot()`、`treeplot()`：

- 补 smoke tests
- 补 `layer` 语义测试
- 补稳定 pathway 选择测试
- 补无 `termsim` / 单 pathway / 双 pathway 边界测试

### 14.2 新方法家族的验证

对于 `phaseplot()`、`consensusmap()`、`rewireplot()`、`mechanismflow()`：

- 先验证 helper 输出结构
- 再验证分类或坐标逻辑
- 最后验证图形对象本身

不建议一开始就只写“图能画出来”的松散测试。

## 15. 分阶段实施顺序

### Phase 0: 建立可复现的验证基线

1. 使用 `pkgload::load_all(".")` 或重新安装当前源码作为验证入口
2. 安装建议依赖（至少 `ggupset`，`upsetplot` 测试依赖它）
3. 在干净环境复跑 `tests/testthat/test-mnsea-helpers.R` 并归档基线结果

### Phase 1: 补齐继承型图

按优先级建议：

1. 把 `gsInfo()` 改成可 dispatch（前置改造）
2. `gseaplot2()`
3. `gsearank()`
4. `hplot()`
5. 公共化 `pairwise_termsim` 的 `mnsea` 相似性语义
6. `treeplot()`
7. `barplot()` 的 gseaResult 家族方法 + 兼容性确认

### Phase 2: 先做一个最有方法学辨识度的新图

建议优先：

1. `phaseplot()`
2. `rewireplot()`

原因是：

- `phaseplot()` 定义整体框架
- `rewireplot()` 提供 pathway-specific 证据层

两者搭起来，最容易形成 `nsea` / `mnsea` 的独特方法标签。

### Phase 3: 再做总览型与演化型图

继续推进：

1. `consensusmap()`
2. `mechanismflow()`

## 16. 非目标

这份 spec 当前不打算：

- 重新设计 `enrichit` 对象结构
- 在 `enrichit` 中加入 plot metadata
- 为每一个旧图都强行发明 `nsea` / `mnsea` 特化语义
- 在第一轮就实现完整的统计重采样稳健性框架

## 17. 开放问题

当前仍需在实现前确认的问题包括：

1. **（已定第一版）**`rewiring_score` 优先采用 `1 - leading-edge Jaccard overlap`；
   是否在第二版加入 edge overlap / 中心性漂移的复合项仍需评估。
2. **（已定第一版）**`phaseplot()` 的 `y` 轴 = `compute_rewiring_score()` 的统一返回值；
   后续可与 edge overlap 并列展示。
3. `consensusmap()` 的默认分类阈值应写死还是交给用户指定？建议第一版写死，
   同时提供 `thresholds` 参数，之后再开放自动选阈值。
4. `mechanismflow()` 是否只服务于 `mnseaResult`，而不强求 `nseaResult` 支持？
   建议第一版只支持 `mnseaResult`（其 layer 序列天然提供状态迁移），
   `nseaResult` 需传入命名列表才能使用。
5. 是否要先把 `pairwise_termsim` 的 `mnsea` 语义进一步公共化，再推进 `treeplot()`？
   **建议升级为必做项**：`treeplot()`、`emapplot()`、`ssplot()` 都依赖相似性语义，
   不应各自维护一套 layer-aware 定义。
6. **（新增）**`consensusmap()` / `mechanismflow()` 的命名列表输入形态、以及
   `rewireplot()` 的 `reference_layer` 是否为必填，需要与 `enrichit` 侧确认没有
   对象结构依赖后再冻结参数签名。

## 18. 当前推荐结论

短期建议非常明确：

1. **先把 `gsInfo()` 改造成可 dispatch，再接齐 `gseaplot2()`、`gsearank()`、`hplot()`**
2. **公共化 `pairwise_termsim` 的 `mnsea` 语义，再补 `treeplot()` 与 `barplot()` 家族方法**
3. **再以 `phaseplot()` 和 `rewireplot()` 作为第一批新方法立项**

这样做的好处是：

- 用户先得到完整工作流
- `nsea` / `mnsea` 再逐步形成自己的方法门脸
- 不会把“兼容旧图”和“发明新图”两件事搅成一锅粥
