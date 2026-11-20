# data-qd-189 · 三角网格定向与闭壳体积

从共享边关系建立面朝向约束，修复可定向分量并用闭壳有符号体积确定方向；保留开边、非流形边和冲突见证。

topology.py 建立连通/异或约束，shells.py 将定向结果用于体积与几何证据，service.py 保留输入顶点身份。 完整字段、单位、身份与失败约定见 [数据契约](docs/contract.md)。样例是作者自制合成数据，许可见 [来源说明](fixtures/SOURCES.md)。

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
```

实际7项测试通过；CLI演示生成JSON，摘要为：

```json
{
  "all_shells_measurable": true,
  "closed_oriented_shells": 1,
  "components": 1,
  "faces": 4,
  "sum_absolute_shell_volume": "1/6"
}
```

输出是原子新建；演示已生成过的同名路径不能重用，复跑请换一个新文件名。非法输入与文件写入错误退出2；业务诊断的表达方式见契约。已经实际核验重复JSON键拒绝、旧结果字节不变、已有符号链接不覆盖。

## 当前边界

不自动焊接、补洞或检测自交；各壳绝对体积之和不等于材料并集体积。 实际后续需求可从嵌套材料壳、面片选择与剖面作业沿既有身份和结果集成；当前未实现这些扩展，也不保证固定数量的未来任务。
