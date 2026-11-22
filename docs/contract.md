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
