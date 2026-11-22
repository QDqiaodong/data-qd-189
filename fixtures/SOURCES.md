# 样例来源

本仓全部样例与预期结果为作者自制合成数据，采用本仓 MIT 许可，不来自真实客户或实际作业。规则是 docs/contract.md 中明确的有界模型。

`handoff/` 下的三份质检交接样例（`previous.json`、`current.json`、`ledger.json`）同为自制合成数据：斜面 `(x,y,z)=(x,y,x)` 上边長为 4 的参数正方形，旧版与新版沿不同对角线划分，两个旧面均为 `pass`。台账中的两个 SHA-256 摘要是按 docs/contract.md 的规范方法对这两份网格实际计算得到的。

