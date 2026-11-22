# data-qd-189 · 三角网格定向、闭壳体积、材料结算与质检交接

从共享边关系建立面朝向约束，修复可定向分量并用闭壳有符号体积确定方向；保留开边、非流形边和冲突见证。可选的 `--material` 模式在同一份报告上增加离线材料结算：可放行零件、嵌套空腔的净体积，以及其余壳整组扣留的可追查理由。成对使用的 `--previous` / `--inspection-ledger` 模式则在重三角化后的新版网格上生成复检清单，并逐张交代旧返修工单的责任去向。`--layer-plan` 模式按工单指定的离散高度审查铺层任务：逐层精确截面、岛级派发判定与可核对的支撑来源。

topology.py 建立连通/异或约束，shells.py 将定向结果用于体积与几何证据，geometry.py 提供精确三角形相交与点包含判定，material.py 做壳资格、接触/包含关系与结算分组，handoff.py 做规范摘要、台账校验、共面面积归属与工单对账，layerplan.py 做精确平面截面、岛提取与支撑结算，service.py 保留输入顶点身份并挂接新增区块。完整字段、单位、身份、见证与失败约定见 [数据契约](docs/contract.md)。样例是作者自制合成数据，许可见 [来源说明](fixtures/SOURCES.md)。

## 安装、测试与演示

Python >=3.11，本批实际使用本机 CPython3.12.14。运行与测试仅用标准库；构建工具精确锁定。每仓独立 .venv，公共包缓存可共享；不共享其他项目的site-packages。

下面命令已在本机已有uv、Python与公共缓存条件下离线执行通过。其他机器如未准备缓存，需先准备构建包；不保证首次安装可以离线完成。

```sh
uv venv --offline --python 3.12 .venv
uv pip install --offline --python .venv/bin/python -r requirements-build.lock
uv pip install --offline --python .venv/bin/python --no-build-isolation --no-deps .
uv pip check --python .venv/bin/python
.venv/bin/python -I -m unittest discover -s tests -v
.venv/bin/mesh-shell --source fixtures/input.json --output run/report-initial.json
.venv/bin/mesh-shell --source fixtures/input.json --output run/report-material.json --material
.venv/bin/mesh-shell --source fixtures/handoff/current.json \
  --output run/report-handoff.json \
  --previous fixtures/handoff/previous.json \
  --inspection-ledger fixtures/handoff/ledger.json
.venv/bin/mesh-shell --source fixtures/layerplan/mesh.json \
  --output run/report-layerplan.json \
  --layer-plan fixtures/layerplan/plan.json
```

实际74项测试通过（原7项 + 材料结算回归14项 + 质检交接回归29项 + 铺层工单回归24项）；CLI 演示生成JSON。旧版摘要为：

```json
{
  "all_shells_measurable": true,
  "closed_oriented_shells": 1,
  "components": 1,
  "faces": 4,
  "sum_absolute_shell_volume": "1/6"
}
```

带 `--material` 时同一 `summary` 不变，另加 `material` 区块；单壳体样例为：

```json
{
  "released_groups": [{"root_component_id": "base", "net_volume": "1/6", "...": "..."}],
  "withheld_groups": [],
  "released_net_volume": "1/6",
  "mesh_net_volume": "1/6"
}
```

输出是原子新建；演示已生成过的同名路径不能重用，复跑请换一个新文件名。非法输入与文件写入错误退出2；业务诊断（包括整组扣留、须复检、工单 `missing`/`ambiguous` 留置）退出0并以 JSON 表达，见契约。已经实际核验重复JSON键拒绝、旧结果字节不变、已有符号链接不覆盖。

## 质检交接模式

`--previous OLD.json --inspection-ledger LEDGER.json` 成对启用，可与 `--material` 同跑。台账必须逐面列出两个版本的所有面并携带两份网格的规范 SHA-256 摘要；摘要按解析后的精确坐标、排序后的顶点/面 ID 与排序后的三角面顶点计算，不受行顺序、小数等值写法与绕序影响。报告新增的 `inspection_handoff` 区块给出：

- 每个新面六类互斥区域（`target_overlap`、`source_overlap`、`pass`、`reject`、`uninspected`、`no_previous`）占该面总面积的精确比例；只有 `pass` 比例恰为 1 的整张面标记 `inherited`，其余全部进入 `reinspection_face_ids`，不做部分合格签发。
- 每个新面正面积相交的同检验面旧 ID、真正贡献 `pass` 凭证的旧面及比例，以及转派进入该面的工单 ID。
- 从每张旧 `reject` 工单出发的对账：按新面 ID 列出的转派比例，以及 `missing`（无新版承接）与 `ambiguous`（多新版争抢）留置比例；分母统一是原旧面面积，合计为 1。

面积归属只在同一 `surface_id` 且严格共面的面之间用精确半平面裁剪计算，点/边接触、跨检验面重合、偏离平面 `0.000001` 都不产生凭证。`fixtures/handoff/` 是斜面正方形换对角线划分的可跑示例：新版两张面各继承两张旧凭证、各占 `1/2`。

## 铺层工单模式

`--layer-plan PLAN.json` 按工单审查离散作业层，可与 `--material`、`--previous`/`--inspection-ledger` 同跑；交接凭证不是铺层依据。工单是恰含 `z_levels`（1 至 32 个严格递增高度）与 `min_support_ratio`（大于 0、至多为 1）的严格 JSON 对象。报告新增 `layer_plan` 区块：

- 按材料结算确定放行组与整组扣留；只给放行组生成区域与铺层记录，扣留组按原分量 ID 与原因列出，不提供面积或支撑，也不会消失为“空层”。
- 每层截面是“严格处于三维材料内部的点投影到 XY 后取闭包再取内部”的精确区域：水平壁落在切平面上不填孔，点/线接触与零宽缝被忽略，真正的孔保留。区域按连通分量分成待铺岛，只在一点相接的区域不并岛；每岛给出报告内稳定的岛 ID、可还原区域的精确 XY 环边界（含孔）与净面积。
- 第一层各岛以台面为支撑（计入面积=岛面积、比例 1）；后续每岛只累计紧邻上一层、同放行组且已可派发的岛的交集面积，比例达到门槛才可派发。非首层逐岛列出与上一层同组各岛的正面积交集、对方岛 ID 与是否计入，区分“没有几何承接”和“承接区域不能派发”。

`fixtures/layerplan/` 是 `[0,6]³` 外盒套 `[1,5]×[1,5]×[1,3]` 空腔、层高 1–5、门槛 0.75 的可跑示例：前三层各是净面积 20 的带孔岛并可派发；高度 4 的岛面积 36、支撑 20、比例 `5/9`，不能派发；高度 5 几何交集 36 但上一岛不可派发、实际支撑 0，仍不能派发。只审查工单指定的离散层，不推断层间连续成形、不生成机床指令、不自动补支撑。

## 当前边界

不自动焊接、补洞或改变坐标；旧版 `sum_absolute_shell_volume` 仍是各独立闭壳绝对体积之和，不能用于材料结算。材料模式用精确几何区分嵌套、接触、包围盒重叠与不合格壳：接触、不确定关系和任何不合格成员都会使整个结算组扣留，任何扣留都会让 `mesh_net_volume` 为 `null`。除上述精确自交判定外，不读取任意 CAD/STL 格式，也不调用外部服务或设备。
