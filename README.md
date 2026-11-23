# data-qd-189 · 三角网格定向、闭壳体积、材料结算与质检交接

从共享边关系建立面朝向约束，修复可定向分量并用闭壳有符号体积确定方向；保留开边、非流形边和冲突见证。可选的 `--material` 模式在同一份报告上增加离线材料结算：可放行零件、嵌套空腔的净体积，以及其余壳整组扣留的可追查理由。成对使用的 `--previous` / `--inspection-ledger` 模式则在重三角化后的新版网格上生成复检清单，并逐张交代旧返修工单的责任去向。`--layer-plan` 模式按工单指定的离散高度审查铺层任务：逐层精确截面、岛级派发判定与可核对的支撑来源。`--load-plan` 模式审查可调平称重台上的整批静力配载：能放行时交付一套可直接复核的目标反力与各梁实载，不能放行时区分装夹前提、单件平衡与共梁超限，并给出可独立核验的共梁超限证书。`--scan-plan` 模式做工单驱动的离线取证排程：在给定视位与拍摄预算内判断能否为全部目标面取证并逐件移出，可行时交付拍摄次数最少、逐步字典序最优的视位与移件顺序；它只审查不实拍，不产生设备指令，也不改变任何在场零件。`--feature-fit` 模式做夹具偏移下的扫描特征验收：在工单给定的角度集合与连续偏移区间内求每个角度的可注册域，用有限个轴对齐偏移矩形完整覆盖，并为每片交付可复核的统一证据（合格矩形给整片成立的完整不复用对应，不合格矩形给整片成立的 Hall 缺额），批次结论区分 `pass`/`fail`/`conditional`/`unregistered`。

topology.py 建立连通/异或约束，shells.py 将定向结果用于体积与几何证据，geometry.py 提供精确三角形相交与点包含判定，material.py 做壳资格、接触/包含关系与结算分组，handoff.py 做规范摘要、台账校验、共面面积归属与工单对账，layerplan.py 做精确平面截面、岛提取与支撑结算，loadplan.py 做装夹前提检查、一阶体积矩重心与精确线性规划配载，featurefit.py 做基准可注册域求交、检出-验收项可行盒排列与整片统一证据覆盖，service.py 保留输入顶点身份并挂接新增区块。完整字段、单位、身份、见证与失败约定见 [数据契约](docs/contract.md)。样例是作者自制合成数据，许可见 [来源说明](fixtures/SOURCES.md)。

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
.venv/bin/mesh-shell --source fixtures/loadplan/mesh.json \
  --output run/report-loadplan.json \
  --load-plan fixtures/loadplan/plan.json
.venv/bin/mesh-shell --source fixtures/scanplan/mesh.json \
  --output run/report-scanplan.json \
  --scan-plan fixtures/scanplan/plan.json
.venv/bin/mesh-shell --source fixtures/featurefit/mesh.json \
  --output run/report-featurefit.json \
  --feature-fit fixtures/featurefit/spec.json
```

实际167项测试通过（原7项 + 材料结算回归14项 + 质检交接回归29项 + 铺层工单回归24项 + 静力配载回归25项 + 离线取证回归38项 + 夹具偏移验收回归30项）；CLI 演示生成JSON。旧版摘要为：

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

## 静力配载审查模式

`--load-plan PLAN.json` 审查一批已完成几何检查的零件能否同时放上可调平称重台：每件零件落在四个垫（一个非退化轴对齐矩形的四角）上，垫的支承力按工单指定传给承载梁。工单恰含 `deck_z`（所有垫的承托高度）、`beams`（1–4 根，各有非负限额）与 `loads`（1–16 件，各有正的单位体积重量与恰四个垫）。装载项必须引用材料结算放行组的根；报告新增 `load_plan` 区块，按优先级给出四种业务结论之一：

- `clamping_failed`：所列 ID 落在被扣留组（列出对应组及原扣留原因）、零件外壳最低 Z 不恰为 `deck_z`、或垫点未实际落在根壳边界上；不给整批反力或梁超限证书。
- `part_infeasible`：忽略梁限额后单件平衡式加垫限额仍不成立，列出所有此类装载项。
- `beam_infeasible`：每件单独可行但没有满足所有梁限额的组合，交付共梁超限证书（和为 1 的非负梁权重、逐件加权载荷的精确最小值、以及严格大于加权限额的最小值之和），工艺员可用报告、原工单与单件平衡式独立重算。
- `ok`：交付同一套完整、按零件与垫 ID 可索引的目标反力，以及按梁 ID 计的全批实载；多份正确方案并存时按零件 ID 字典序逐件最小化自由参数，保证结果稳定。

每件报告净体积、重量与由一阶体积矩精确求得的重心 XY（空腔扣除体积与矩、内岛加回）。不模拟刚度、动态受力或调平过程，反力是静力可行目标而非传感器读数。`fixtures/loadplan/` 是三个 `2×2×1` 独立盒、两条对角线分接 A/B 梁、第三件全接 C 梁的可跑示例：限额 A=3、B=5、C=10 时可配载，全批实载精确为 A=3、B=5、C=4。

## 离线取证计划模式

`--scan-plan PLAN.json` 按工单审查本班次能否对全部目标面取证并逐件移出，可与 `--material` 等同跑；材料结算在内部复用，显式 `--material` 时原材料块按拍摄前完整输入原样保留。工单恰含 `targets`（1–24 个带覆盖门槛的原面 ID，门槛为大于 0、至多 1 的十进制字符串）、`poses`（1–4 个唯一视位：六向轴方向之一加平行投影窗口，去掉视位轴后其余两轴按原坐标顺序作 u、v）与 `max_passes`（1–5 的整数拍摄预算）。报告新增 `scan_plan` 区块：

- 目标面属于原扣留组时结论为 `targets_withheld`，列出扣留组、原扣留原因与触发面，不发布任何工步，也不能靠缩小目标集报完成；预算内不存在完整计划时为 `no_plan_in_budget`。两者都是退出 0 的正常业务结论，未知面等工单错误仍退出 2。
- 可行时状态为 `ok`，先按拍摄次数最短、再按每步 `(pose_id, 移出根 ID/null)` 逐步字典序（`null` 最后）择优。每步给出视位、该步开始仍在场各目标面的本次可见/新增/累计比例、拍完后全部有资格移件的根 ID 与实际移出的根（不移为 `null`）。
- 可见性同时要求材料朝外法向严格朝向相机、平行投影在窗口内、射线在严格正深度处不遇任何在场三角面。遮挡判定把投影目标面按窗口、遮挡面投影边与精确等深线切成凸格后逐格做严格深度比较，不接受投影重合、包围盒、顶点采样或容差结论；共面等深不遮挡。覆盖区域映射回原三角面的重心参数面取并集，复拍同一窗口不重复入账，比例全部以最简整数或分数字符串报告。
- 未被吊走的达标零件继续遮挡但保留证据与资格；每次吊具至多移一整件（整组消失，不止目标面），工步外不移件。计划不是实拍：不生成控制指令、不修改输入、不移动零件。

`fixtures/scanplan/` 是一个闭盒、一个带对角剖分顶面的闭盒与一个直三棱柱叠放的可跑示例：三视位 K/L/R、预算三拍时最短计划为 `(K,B),(L,A),(R,C)`，其中首拍后 A、B 同时合格但必须先吊 B，三棱柱靠 L、R 两窗拼证（`7/8` 加新增 `1/8`）才能移出。

## 夹具偏移特征验收模式

`--feature-fit SPEC.json` 在夹具偏移包络内验收扫描检出能否签收全部设计特征，可与 `--material` 等同跑；各区块互不影响。工单顶层恰含 `reference`（设计网格，同格式同上限，顶点 ID 与扫描网格互不相认）、`clamp`（0/90/180/270 中非空不重复的角度数组，以及 `dx`、`dy` 两个连续闭区间）、`datums`（1–8 条已确认基准对应，两侧各自不重复引用，三轴非负公差）、`detections`（0–12 个互不重复的扫描检出点）与 `features`（1–8 条验收项：唯一 ID、设计顶点不共用、三轴公差、候选检出名单）。扫描点先绕 +Z 按右手方向旋转再加 `(dx,dy,0)`；偏移在区间内连续可取。报告新增 `feature_fit` 区块：

- 对每个允许角度，可注册域是夹具矩形与全部基准轴对齐约束的交集（闭矩形、线、点或空），空域也逐角度列出。批次结论：全空为 `unregistered`；否则全部可注册姿态合格为 `pass`、全不合格为 `fail`、两类都有为 `conditional`，空域不参与合计。
- 每个非空域用有限个轴对齐偏移矩形无缺口覆盖（本实现给出互不相交的划分），逐轴记录精确端点与是否包含，允许线、点。合格矩形携带一套在整片每个偏移都成立的完整、不复用对应（按验收项 ID 索引）；不合格矩形携带一个非空验收项子集及它在整片上可能用到的全部检出点，点数严格少于项数（Hall 缺额），证明整片都签不满。缺额与对应都是对整片统一成立的证据，不是抽查。
- 全部坐标与判定用精确 `Fraction`：不采样整数点、不用浮点容差、不丢弃零宽区域；端点与公差用最简整数或分数字符串。下游可按契约的重算方法用原坐标与工单独立验收整片范围，不需要重新做最优匹配。

`fixtures/featurefit/` 是只允 0 度、夹具 `dx∈[-1,1]`、`dy∈[0,1]` 的可跑示例：基准都在原点、公差 `[1,1,0]`，设计 A/B 在 X=0/3（A 与原点基准共用设计顶点），检出 P/Q 在 X=1/2，两项候选都是 P、Q。整个夹具范围可注册，但仅 `dx=0` 的线合格（A→P、B→Q），批次为 `conditional`。

## 当前边界

不自动焊接、补洞或改变坐标；旧版 `sum_absolute_shell_volume` 仍是各独立闭壳绝对体积之和，不能用于材料结算。材料模式用精确几何区分嵌套、接触、包围盒重叠与不合格壳：接触、不确定关系和任何不合格成员都会使整个结算组扣留，任何扣留都会让 `mesh_net_volume` 为 `null`。除上述精确自交判定外，不读取任意 CAD/STL 格式，也不调用外部服务或设备。
