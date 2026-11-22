# 三角面壳定向 v1 与材料结算扩展 v1

输入vertices(3..1000)为唯一ASCII id和3个有界十进制字符串坐标；faces(1..2000)为唯一id和3个既有顶点id。坐标单位由调用方统一给定，体积为该单位立方。全部几何使用Fraction；输出坐标/体积为整数或分数字符串。未知顶点、三点共线/重合的退化三角形、同一组三顶点的重复面（即使朝向相反）拒绝。不同id的重合坐标不会自动焊接。

按共享无向边求面连通分量，component_id为该分量字典序最小输入face id。单面边是boundary，>2面边是nonmanifold。恰两面共享边施加方向异或约束：最终两面沿公共边必须反向。按face id确定种子，BFS求翻面；冲突边给出见证。存在非流形边或方向冲突的整个分量不发布定向面，不能把局部修补假装可用壳。

可定向的开曲面输出一致朝向及边界，但不计算体积；方向只是确定性种子选择，不声明面向外。闭合且可定向分量按Σ dot(a,cross(b,c))/6求有符号体积；负值则整壳翻面，得到正向体积。零体积明确标问题且体积null。保留全部输入面身份，记录flipped_face_ids和oriented_faces。没有新增顶点、补洞、曲面简化或几何变形。

不检测面自交、嵌套壳、非流形顶点（边流形仍可能存在）、壳间相交或同点不同id；closed仅表示所有边有两面。sum_absolute_shell_volume只在每壳均有体积时给出各独立壳绝对体积之和，不能当作材料并集、带空腔净体积或保证可打印。诊断是正常业务结果退出0；输入错误退出2。样例为单位轴四面体，初始base朝内，修正后体积1/6。

全部结果先计算后原子新建单个JSON，不覆盖已有目标。没有读取任意CAD/STL格式、自动拓扑焊接、外部服务或设备调用。

# 材料结算模式（--material）

不带 `--material` 时，`process(data)`、CLI 行为、JSON 字段及其含义与 v1 完全一致。带 `--material` 时输入格式、坐标、数量上限与身份（原顶点 ID、原面 ID、原 component_id）不变，旧报告全部字段原样保留，另加顶层 `material` 对象。材料诊断是正常业务结果（退出0）；非法输入、原子新建冲突等仍按 v1 失败语义（退出2），且失败时不产生输出文件。处理与测试离线、仅标准库、无新增运行时依赖；不焊接、不补洞、不变换坐标。

## 几何语义

所有判定使用精确坐标（`Fraction`），不使用浮点容差、近似体积或包围盒去猜包含。三处几何判定：

- `triangle_intersection_dimension(A,B)` 返回 `-1/0/1/2`，分别表示不相交、仅一点、仅一段、有面积重叠。异面三角形通过两个平面截段在交线上的投影区间判定；共面三角形做精确 Sutherland–Hodgman 裁剪后按不同点数与共线性定维。
- 点是否在闭壳内：从该点沿一条射线做奇偶交叉计数；若射线恰好穿过顶点、边或整个面所在平面则换下一条确定性方向重试，避免任何平局规则。
- 包围盒只用于「可能相交」的快速排除，绝不用于推断包含。

## 壳资格

先要求旧分析已表明壳 `closed=true`、`orientable=true` 且 `absolute_signed_volume` 非 `null`。再要求：

1. 壳内每个顶点周围入射面构成单个闭合扇环。判定按旧报告 `oriented_faces` 的绕序：每个入射面对该顶点恰有一条流入边 `prev→v` 和一条流出边 `v→next`；这些面必须通过流出顶点连到下一面的流入顶点，串成恰好访问每个入射面一次的单一闭环。
2. 任意两面实际交集维度等于共用输入顶点 ID 所决定的部分：无共同 ID 必须为 `-1`；一个共同 ID 必须为 `0`；两个共同 ID 必须为 `1`。

不满足的壳是**不合格壳**。资格的三种原因互斥、按以下顺序给出：

| `reason` | 含义 | `witness` |
| --- | --- | --- |
| `legacy_analysis_failed` | 旧分析未通过 | `{"legacy_issues":[...]}`（旧报告 `issues`，升序原样） |
| `vertex_fan_invalid` | 顶点扇环不合格 | `{"vertex_id":"..."}`（字典序最小的反例顶点） |
| `face_intersection_invalid` | 面交集不合格 | `{"face_ids":["...","..."]}`（按面顺序扫描到的第一对，两端各自是原面 ID，内部按字典序） |

不合格壳不出现在放行组中；其原旧报告 `issues` 等仍可在同一 JSON 的旧 `components` 中查到。

## 壳间关系与结算组

对每一对壳（覆盖全部输入分量，包括不合格壳；按 `component_id` 字典序扫描）：

- 两壳都合格时，先测边界接触：任一对面片实际交集维度 `>=0`（一点、沿边或共面重叠都算）即登记 `contact`，见证为该对原面 ID（升序二元素）。接触关系不推断包含。
- 无接触时按真实有界内部判定严格包含：以被包含壳字典序最小的顶点作探针，用奇偶射线法判断它是否在对方严格内部。因合格壳是面连通闭曲面、两壳边界已确证无交，一个严格内部探针即可确定整壳严格包含；探针落在对方边界上则不出包含。被包含壳的该探针顶点 ID 作为见证。
- 接触、包含或以下不确定关系都把两端并入同一结算组：
- 涉及至少一个不合格壳的壳对不判内部；只要两壳旧报告闭包围盒相交（含边界贴合）就登记 `uncertain`，带 `reason:"unqualified_shell"` 和排序后的原 `component_ids`。
- 严格包含但无法明确唯一方向时（防御分支）登记 `uncertain`，带 `reason:"containment_ambiguous"`。

关系沿多跳传递（并查集），无关系的壳单独成组。组内有任一不合格壳、任一 `contact` 或任一 `uncertain` 关系即**整组扣留**，不能发布组内任何壳的材料角色、材料朝向或净体积。互不相干的合格组仍放行。

## `material` 字段

```json
"material": {
  "released_groups": [ ... ],
  "withheld_groups": [ ... ],
  "released_net_volume": "整数或分数字符串",
  "mesh_net_volume": "整数或分数字符串或 null"
}
```

排序约定：

- `released_groups` 按 `root_component_id` 字典序；`withheld_groups` 按 `member_component_ids` 字典序。
- 放行组 `shells` 按 `(depth, component_id)` 升序；`relations` 与扣留组一致地排序。
- `relations` 排序键：`containment` 用 `(kind, container_id, contained_id, witness.vertex_id)`，其余用 `(kind, component_ids[0], component_ids[1], face_ids)`。
- 放行组内 `material_oriented_faces` 按原 `face_ids` 升序；顶点循环移位为最小原顶点 ID 起始（保持绕序方向不变）。

放行组：

```json
{
  "root_component_id": "...",
  "net_volume": "整数或分数字符串",
  "shells": [{
    "component_id": "...",
    "depth": 0,
    "parent_component_id": null,
    "role": "outer_wall",
    "material_sign": 1,
    "absolute_signed_volume": "...",
    "net_volume_contribution": "...",
    "material_oriented_faces": [{"face_id": "...", "vertices": ["...","...","..."]}],
    "flipped_face_ids_relative_to_legacy": ["..."]
  }],
  "relations": [ ... ]
}
```

以最外层壳为根，`depth` 从0起；直接容器是严格包含它且最深的壳（同深按 `component_id` 字典序最大者），根为 `null`。偶数深度 `role` 为 `outer_wall`（材料外壁），奇数深度为 `cavity_wall`（空腔壁）。净体积按 `Σ (-1)^depth × 绝对体积` 精确计算，`material_sign` 即 `(-1)^depth`。`material_oriented_faces` 按材料朝外定向：外壁保留旧报告的向外绕序，空腔壁整体反向（朝壳体有界内部），并把每个面的顶点循环移位为最小 ID 起始；`flipped_face_ids_relative_to_legacy` 记录相对旧绕序发生反转的面 ID。旧报告自身的 `oriented_faces`/`flipped_face_ids` 不被改写。

扣留组：

```json
{
  "member_component_ids": ["..."],
  "withheld_reasons": ["unqualified_shell", "contact", "uncertain_relation"],
  "qualified_member_component_ids": ["..."],
  "unqualified_shells": [{"component_id": "...", "reason": "...", "witness": { ... }}],
  "relations": [ ... ]
}
```

`withheld_reasons` 按固定顺序出现（`unqualified_shell`、`contact`、`uncertain_relation`），只含实际成立者；`member_component_ids` 与 `qualified_member_component_ids` 升序；`unqualified_shells` 按 `component_id` 升序。完整成员与见证让下游能解释为什么一个未接触问题面的内壳也一起扣留（多跳传递）。

关系对象：

```json
{"kind": "containment", "container_id": "...", "contained_id": "...", "witness": {"vertex_id": "..."}}
{"kind": "contact", "component_ids": ["...","..."], "face_ids": ["...","..."]}
{"kind": "uncertain", "component_ids": ["...","..."], "reason": "unqualified_shell"}
```

## 汇总字段

- `released_net_volume`：全部放行组 `net_volume` 的精确和；没有放行组时为字符串 `"0"`。
- `mesh_net_volume`：只在没有任何扣留组时等于 `released_net_volume`；只要存在扣留组即为 `null`（不明部分不当作空腔，也不从输入中删除）。

## 稳定性

相同几何与 ID 的输入，即使数组乱序、面三顶点循环移位或原绕序反转，新增 `material` 区块（其余旧报告字段同样）保持稳定：所有遍历顺序基于原 ID 字典序，绕序判定基于旧报告的规范化朝向，体积与长方体维数均与输入顺序无关。

# 质检交接模式（--previous + --inspection-ledger）

`--previous OLD.json --inspection-ledger LEDGER.json` 必须成对出现，只给一个是输入错误（退出 2）。两份网格各自沿用 v1 的格式、坐标边界与数量上限，面 ID 只在本版有效。不带这两个参数时，CLI、`process(data)` / `process(data, material=True)` 与报告字节不变；带参数时旧报告全部字段原样保留，另加顶层 `inspection_handoff` 对象。该模式可与 `--material` 同跑，材料结算与复检责任互不替代；开曲面也可以持有有效的旧质检凭证。所有计算离线、仅标准库，全部几何判定使用精确 `Fraction`，无浮点容差，不修改任何输入文件或工单。

## 台账与摘要

台账是严格 JSON，顶层恰有四个字段：`previous_mesh_sha256`、`current_mesh_sha256`、`previous_faces`、`current_faces`，不允许重复 JSON 键。

- `previous_faces` 每行含 `face_id`、`surface_id`、`verdict`；`verdict` 只能是 `pass`、`reject`、`uninspected`。`reject` 行必须另有全表唯一的 `ticket_id`；其他结论带工单号、工单号重复都是输入错误。
- `current_faces` 每行只含 `face_id`、`surface_id`。
- 两个面数组都必须把该版本网格的**所有**面各列一次：缺项、未知 ID、重复 ID 均为输入错误。`surface_id` 与网格 ID 同语法；同一标识可对应多个三角面，也允许只在一个版本出现。
- 两份网格先各自通过 v1 校验（未知顶点、退化面、重复三角面等同样拒绝），再校验摘要。

摘要绑定实际网格内容，不绑定面 ID 列表。规范对象保留原 `vertices`、`faces` 结构：两数组各按 ID 排序；每个三角面的三个顶点 ID 排序；坐标取解析后的精确值并写成最简整数或分数字符串。对此对象按 JSON 键排序、无分隔空白（`separators=(',', ':')`）、UTF-8 编码后取 SHA-256 小写十六进制。因此导出行顺序、等值小数写法（如 `1.0` 与 `1`）、面绕序都不影响摘要；坐标真的变动会改变摘要。摘要不符即输入错误。

## 共面与面积归属

面积归属只在**同一 `surface_id` 且严格共面**的面之间发生：不同检验面标识即使几何完全重合也互不可见；平面以法向量（首个非零分量统一为正后归一化）和平面常数为精确键，离开原平面哪怕 `0.000001` 也不共面。点接触、边接触、非共面相交都不提供面积；不使用顶点匹配、包围盒、原 ID 相等或任何容差猜测。共面后取法向量最大分量对应的坐标轴做二维投影，二维面积比等于三维面积比。

对每个新面，把同组新版面与同组旧版面所有相关三角形的支撑线切割该面，得到内部覆盖关系恒定的凸多边形分区（同一区域不会被重复计数）。每个分区按互斥优先级归入六类之一：

1. 被任意另一**新版**面覆盖的区域：`target_overlap`，不论旧证据。
2. 其余被两个或更多**旧**面覆盖的区域：`source_overlap`，即使旧结论一致也不能挑一张当凭证。
3. 其余被恰好一个旧面覆盖的区域：按该旧面的 `verdict` 记 `pass`、`reject` 或 `uninspected`。
4. 完全没有旧面覆盖的区域：`no_previous`。

六类比例的分母都是该新面自身总面积，六类同时输出（含为零者），合计恰为 1，全部为最简整数或分数字符串。

## `inspection_handoff` 字段

```json
"inspection_handoff": {
  "current_faces": [{
    "face_id": "...",
    "surface_id": "...",
    "area_fractions": {
      "target_overlap": "...", "source_overlap": "...",
      "pass": "...", "reject": "...",
      "uninspected": "...", "no_previous": "..."
    },
    "inherited": true,
    "intersecting_previous_face_ids": ["..."],
    "pass_credentials": [{"face_id": "...", "fraction": "..."}],
    "transferred_ticket_ids": ["..."]
  }],
  "reinspection_face_ids": ["..."],
  "tickets": [{
    "ticket_id": "...",
    "previous_face_id": "...",
    "transfers": [{"face_id": "...", "fraction": "..."}],
    "missing": "...",
    "ambiguous": "...",
    "fully_transferred": true
  }]
}
```

约定：

- `current_faces` 按 `face_id` 字典序；每个新面都出现。`area_fractions` 固定包含上述六个键。
- `inherited` 仅当 `pass` 比例恰为 1 时为真——整张新面是最小可签单元；缺证区域不会被自动签发部分合格。其余所有面（含 `reject`、`uninspected`、任何重叠或无旧覆盖）都进 `reinspection_face_ids`（字典序）。
- `intersecting_previous_face_ids`：同检验面、同平面上与该新面有正面积相交的全部旧面 ID，字典序；只有点/边接触的旧面不列入。
- `pass_credentials`：对 `pass` 类真正贡献凭证的旧面 ID 及各自占该新面总面积的比例，按旧面 ID 字典序。同一张新面可以合法地从多个互不重叠的旧 `pass` 面继承（例如对角线切换后各承一半）。
- `transferred_ticket_ids`：按后文规则转派进入该面的旧工单 ID，字典序。

## 旧返修工单对账

对账从**每一张**旧 `reject` 面出发，不依赖新版是否能匹配；没有任何新版承接的工单也必须列出。在同检验面、同平面上，对旧返修面做同样的平面分区，逐区统计承接它的新版面数量：

- 恰有一个新版面承接的区域，其面积转派给该面，并进入该面复检记录的 `transferred_ticket_ids`。
- 零个新版面承接的区域记 `missing`；多个新版面争抢的区域记 `ambiguous`。两部分留在旧工单侧人工追查，绝不权宜地分给其中任一面。
- 即使同处还有旧 `pass` 或其他凭证，已明确记录的缺陷不抹去；转派只看新版承接是否唯一。

`tickets` 按 `ticket_id` 字典序，每条含旧面 ID、按新面 ID 字典序排列的 `transfers`（比例分母统一为该工单的原旧面面积）、`missing` 与 `ambiguous` 比例，三者合计为 1。`fully_transferred` 仅当两项留置比例都为 0 时为真；它只表示责任全数转派，不等于返修完成。

## 失败与稳定性

输入错误（台账缺项/未知/重复 ID、摘要不符、无效网格、只给一个新参数、台账结构或工单号非法）与发布失败一律退出 2，沿用原子新建语义：不覆盖既有输出、不留半成品报告。正常得出“需复检”或“须人工追查”（`missing`/`ambiguous` 非零）仍是成功，退出 0。等价的网格重排、绕序变更、台账行乱序不改变 `inspection_handoff`：平面键、直线方向与全部遍历顺序都与输入顺序和朝向无关。

# 铺层工单模式（--layer-plan）

不带 `--layer-plan` 时，CLI、`process(data)` / `process(data, material=...)` / `process(data, handoff=...)` 与报告字节不变。带 `--layer-plan PLAN.json` 时输入网格的格式、坐标、数量上限与身份不变，旧报告全部字段原样保留，另加顶层 `layer_plan` 对象。不要求同时传 `--material`；显式传入时 `material` 区块完整保留并与之并列。可与 `--previous` / `--inspection-ledger` 同跑，但交接凭证不是铺层依据，两区块互不影响。本模式只审查工单指定的离散高度：不推断层间连续成形、不生成机床指令、不自动补支撑。所有计算离线、仅标准库，全部几何判定使用精确 `Fraction`，无浮点容差。

## 工单格式

工单是严格 JSON 对象（重复键拒绝），恰有两个字段：

- `z_levels`：1 至 32 个高度，按铺层顺序给出，必须严格递增。每项使用与顶点坐标相同的有界十进制字符串语法。
- `min_support_ratio`：同一语法的十进制字符串，解析值须大于 0、至多为 1。

工单结构非法（缺字段、多字段、层数越界、不递增、比例越界、语法非法）是输入错误，退出 2 且不产生输出文件。

## 作业层与台面

`z_levels` 的第一项代表与台面相接的作业层：台面覆盖该层的全部可结算区域，该层所有待铺岛都以台面为支撑。即使第一层没有任何岛，后续层也不重新起算台面——只有第一项是台面层。相邻两项就是相邻作业层，不要求等间距；支撑关系只在相邻作业层之间成立，跳过的层不存在。报告逐层显式保留每一作业层（可以无岛），按工单顺序排列。

## 组、壳角色与扣留

组、壳角色（`outer_wall`/`cavity_wall`）、材料朝向与整组扣留完全按材料结算规则确定（内部总是先完成结算，无论是否发布 `material` 区块）。只有放行组生成区域与铺层记录。扣留组按原分量 ID 列出并给出原因，不提供面积或支撑，也不会消失为“空层”；互不相干的放行组仍有独立计划。

## 截面区域

对每个放行组，按壳的嵌套角色解释三维材料：外壁围成的有界内部减去空腔，加回空腔中的实体岛（即按 `material_sign` 交替求和恰为 1 的点集）。在指定高度取严格处在三维材料内部的点投影到 XY 平面；铺层区域是该集合在平面内取闭包后再取内部的结果。因此：

- 只有点或线的孤立接触、零宽缝被忽略；真正的孔保留。
- 完全落在切平面上的水平外壁或空腔壁本身不算材料内部，不能因为在相邻高度的某一侧有材料就填入它覆盖的区域（例如切平面落在空腔顶/底壁上时，孔不消失；切平面只落在外盒底或顶外壁时不产生面积）。

实现不漂移切平面、不用容差合并几何、不猜测包围盒内部，也不把切到的三角形面积直接相加：用全部交线段与水平三角形边的支撑线把组的 XY 包围框切成凸格，逐格在重心处做精确三维点包含判定（先排除落在水平壁上的格），再取材料格的闭包并集的内部。包围框只用于限制排列范围，分类完全由精确判定给出。

## 岛与边界

铺层区域的连通分量（按正长度公共边相连）是待铺岛；两个区域只有边界点相接不能并岛，同一壳组同一高度可以有多个岛。每岛给出：

```json
{
  "island_id": "isl_7",
  "root_component_id": "...",
  "boundary": {
    "outer_rings": [[["0","0"],["6","0"],["6","6"],["0","6"]]],
    "hole_rings": [[["1","1"],["1","5"],["5","5"],["5","1"]]]
  },
  "net_area": "20",
  "dispatch": { ... }
}
```

- `island_id`：报告内稳定且唯一的岛 ID（`isl_` 前缀加序号），供后续层的支撑记录引用。赋值顺序：作业层按工单顺序，层内先按 `root_component_id` 字典序，同组内按规范化边界键字典序。层内 `islands` 数组同序。
- `boundary`：能还原该岛区域的精确 XY 边界。每个环是一个简单闭合多边形，顶点不重复列出末点；外环逆时针（正有向面积）、孔环顺时针，起点都是字典序最小顶点，共线点已化简。岛的区域恰是其全部环的奇偶（even-odd）填充内部；净面积等于外环面积减孔环面积。边界相接（共点）时各岛的环分别列在各自记录下，归属明确。坐标用最简整数或分数字符串。
- `net_area`：XY 平面净面积（不是原三角面在三维空间的面积），最简整数或分数字符串。

## 派发判定

第一个作业层的全部待铺岛以台面为支撑：计入面积等于岛净面积、比例为 `"1"`、无上一层候选，恒可派发。以后每岛只可使用紧邻上一作业层、同放行组且已可派发的岛作支撑；其他组、跳过的层和上一层不能派发的岛一律不算。支撑面积是本岛与这些可用区域交集的面积（同层各岛互不重叠，不重复计入），除以本岛净面积后大于或等于 `min_support_ratio` 才可派发；只碰边或点的交集面积为零。不能派发的岛仍保留完整几何记录，但下一层不得把它当成已铺。

```json
"dispatch": {
  "dispatchable": false,
  "support_source": "previous_layer",
  "counted_area": "20",
  "support_ratio": "5/9",
  "candidates": [{"island_id": "isl_3", "intersection_area": "20", "counted": true}]
}
```

- `support_source`：`build_plate`（仅第一层）或 `previous_layer`。
- `counted_area` / `support_ratio`：实际计入的精确面积及其占本岛净面积的比例，最简整数或分数字符串。
- `candidates`：仅非首层出现非空；逐一列出与上一层同组各岛的正面积交集、对方岛 ID 和是否计入（`counted` 即对方可派发）。由此可区分“没有几何承接”（候选为空）与“承接区域不能派发”（候选 `counted` 为 `false`）。候选按对方岛 ID 的赋值顺序排列。

## `layer_plan` 字段与扣留组

```json
"layer_plan": {
  "z_levels": ["1", "..."],
  "min_support_ratio": "3/4",
  "layers": [{"z": "1", "islands": [ ... ]}],
  "withheld_groups": [{"member_component_ids": ["..."], "withheld_reasons": ["..."]}]
}
```

`z_levels` 与 `min_support_ratio` 按解析后的精确值回显为最简整数或分数字符串。`layers` 与工单逐项对应。`withheld_groups` 按 `member_component_ids` 字典序，与 `material` 区块的扣留组一一对应（完整见证见 `material` 区块或旧报告 `components`）。

## 失败与稳定性

工单无效与报告发布失败沿用退出 2、原子新建、不覆盖既有目标的契约，失败时不产生输出文件；扣留或不能派发属于正常成功结果（退出 0）。仅改变顶点/面输入顺序、等值小数写法（网格或工单）或面绕序时，`layer_plan` 区块相同：所有遍历顺序基于原 ID 字典序与规范化几何键，坐标与面积都是精确值。换三角剖分或重编 ID 后不保证岛 ID 相同。

