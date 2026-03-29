# CASPToolkit

CASPToolkit 是一个面向蛋白质结构处理和 CASP 流程的 Python 工具包。

本项目当前采用脚本直跑模式：你可以按任务选择脚本运行，不需要统一 CLI 入口。

## 安装

```bash
git clone <repo>
cd CASPToolkit
pip install -e .
```

## 外部工具配置

部分脚本依赖外部程序（USalign、Phenix）。

推荐通过环境变量配置，不需要手动修改源码：

```bash
export CASPTOOLKIT_USALIGN_PATH=/path/to/USalign
export CASPTOOLKIT_PHENIX_CLASHSCORE_PATH=/path/to/phenix.clashscore
```

配置优先级：环境变量 > 代码默认路径。

兼容说明：历史变量 `PDBTOOLKIT_USALIGN_PATH` 和 `PDBTOOLKIT_PHENIX_CLASHSCORE_PATH` 仍可使用。

## 快速导航：按任务选择脚本

### 1) CIF 转 PDB

脚本：casptoolkit/PDBOps/cif2pdb.py

```bash
python casptoolkit/PDBOps/cif2pdb.py input.cif output.pdb
python casptoolkit/PDBOps/cif2pdb.py input_dir/ output_dir/ --num_workers 8
```

输入：单个 .cif 文件或目录。
输出：单个 .pdb 文件或输出目录。

### 2) 原子序号重编号

脚本：casptoolkit/PDBOps/renumber_atoms.py

说明：该模块主要作为函数被其他脚本调用，用于处理大结构中原子序号超过 PDB 上限的问题。

### 3) 合并多个结构文件

脚本：casptoolkit/PDBOps/merge_structure.py

```bash
python casptoolkit/PDBOps/merge_structure.py input_dir/ merged.pdb
```

输入：包含多个 PDB 文件的目录。
输出：合并后的 PDB 文件。

### 4) 链名重映射

脚本：casptoolkit/PDBOps/reassign_chain_id.py

```bash
python casptoolkit/PDBOps/reassign_chain_id.py input.pdb output.pdb ABC:DEF
```

输入：PDB 文件 + 链映射字符串（如 ABC:DEF，表示 A->D, B->E, C->F）。
输出：链名重排后的 PDB 文件。

### 5) 计算 clashscore

脚本：casptoolkit/CASP/phenix_clashscore.py

```bash
python casptoolkit/CASP/phenix_clashscore.py -d model_dir/ results.json --n_cpu 8
```

输入：单文件、目录或文件列表。
输出：JSON 结果文件。

### 6) AF3 结果 QA 评分与排序

脚本：casptoolkit/CASP/qa_af3.py

```bash
python casptoolkit/CASP/qa_af3.py input_dir/ output_dir/ --n_cpu 8
```

主要流程：
1. 解压 zip 中的 CIF 和 summary JSON。
2. CIF 转 PDB。
3. 计算 QA（默认公式：iptm*0.8 + ptm*0.2）。
4. 计算平均 pLDDT。
5. 输出 rank_*.pdb 和 qa.csv。

### 7) 模板叠合与 TM-score

脚本：casptoolkit/CASP/sup_template.py

```bash
python casptoolkit/CASP/sup_template.py model_dir/ reference.pdb --output_file tmscore.csv --sup_dir sup/ --n_cpu 8
```

输入：模型目录 + 参考结构。
输出：TM-score 表格，可选叠合结构文件。

### 8) 逐目标叠合并组装

脚本：casptoolkit/CASP/sup_assemble.py

```bash
python casptoolkit/CASP/sup_assemble.py source.pdb target_dir/ output.pdb
```

输入：source 结构 + target 子结构目录。
输出：组装后的 PDB 文件。

### 9) 同聚体叠合（An -> Am）

脚本：casptoolkit/CASP/sup_homooligo.py

```bash
python casptoolkit/CASP/sup_homooligo.py source_An.pdb target_Am.pdb output_dir/
```

输入：两个同聚体结构。
输出：按 source 各链生成的叠合结果目录。

## 常见问题

### 1) 提示找不到 USalign 或 phenix.clashscore

请先确认路径有效，并设置环境变量：

```bash
export CASPTOOLKIT_USALIGN_PATH=/path/to/USalign
export CASPTOOLKIT_PHENIX_CLASHSCORE_PATH=/path/to/phenix.clashscore
```

如果你之前已经使用旧变量名（`PDBTOOLKIT_*`），当前版本仍兼容。

### 2) 如何选择并行参数

PDBOps 脚本使用 `--num_workers`，CASP 脚本使用 `--n_cpu`。

优先从较小值开始测试（如 `--num_workers 2` 或 `--n_cpu 2`），再根据机器资源增大。

### 3) 为什么有些脚本支持目录，有些只支持文件

不同脚本面向的任务粒度不同。建议先看脚本帮助：

```bash
python casptoolkit/CASP/qa_af3.py -h
python casptoolkit/PDBOps/cif2pdb.py -h
```
