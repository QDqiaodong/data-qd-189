# data-qd-189 · 三角网格定向、闭壳体积与材料结算

从共享边关系建立面朝向约束，修复可定向分量并用闭壳有符号体积确定方向；保留开边、非流形边和冲突见证。可选的 `--material` 模式在同一份报告上增加离线材料结算：可放行零件、嵌套空腔的净体积，以及其余壳整组扣留的可追查理由。

topology.py 建立连通/异或约束，shells.py 将定向结果用于体积与几何证据，geometry.py 提供精确三角形相交与点包含判定，material.py 做壳资格、接触/包含关系与结算分组，service.py 保留输入顶点身份并挂接材料区块。完整字段、单位、身份、见证与失败约定见 [数据契约](docs/contract.md)。样例是作者自制合成数据，许可见 [来源说明](fixtures/SOURCES.md)。

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
```

实际21项测试通过（原7项 + 材料结算回归14项）；CLI 演示生成JSON。旧版摘要为：

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

输出是原子新建；演示已生成过的同名路径不能重用，复跑请换一个新文件名。非法输入与文件写入错误退出2；业务诊断（包括整组扣留）退出0并以 JSON 表达，见契约。已经实际核验重复JSON键拒绝、旧结果字节不变、已有符号链接不覆盖。

## 当前边界

不自动焊接、补洞或改变坐标；旧版 `sum_absolute_shell_volume` 仍是各独立闭壳绝对体积之和，不能用于材料结算。材料模式用精确几何区分嵌套、接触、包围盒重叠与不合格壳：接触、不确定关系和任何不合格成员都会使整个结算组扣留，任何扣留都会让 `mesh_net_volume` 为 `null`。除上述精确自交判定外，不读取任意 CAD/STL 格式，也不调用外部服务或设备。
