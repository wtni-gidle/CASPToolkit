# CASPToolkit

CASPToolkit 是实验室内部使用的蛋白质结构处理脚本集合，按脚本直接运行。

## 安装

```bash
git clone <repo>
cd CASPToolkit
pip install -e .
```

## 外部工具配置

部分脚本依赖 USalign 与 phenix.clashscore。请先设置：

```bash
export CASPTOOLKIT_USALIGN_PATH=/path/to/USalign
export CASPTOOLKIT_PHENIX_CLASHSCORE_PATH=/path/to/phenix.clashscore
```

## 脚本速览

### PDBOps

1) CIF 转 PDB

```bash
python casptoolkit/PDBOps/cif2pdb.py input.cif output.pdb
python casptoolkit/PDBOps/cif2pdb.py input_dir/ output_dir/ --num_workers 8
```

2) 合并 PDB

```bash
python casptoolkit/PDBOps/merge_structure.py input_dir/ merged.pdb --renumber_atoms
```

3) 链 ID 重映射

```bash
python casptoolkit/PDBOps/reassign_chain_id.py input.pdb output.pdb ABC:DEF --renumber_atoms
```

### CASP

1) clashscore

```bash
python casptoolkit/CASP/phenix_clashscore.py model_dir/ --output_path results.json --num_workers 8
```

说明：`--output_path` 仅在 `input_path` 为目录时生效。

2) 模板叠合与 TM-score

```bash
python casptoolkit/CASP/sup_template.py model_dir/ reference.pdb --output_file tmscore.csv --sup_dir sup/ --num_workers 8
```

说明：`--output_file` 与 `--sup_dir` 至少提供一个。

3) 逐目标叠合并组装

```bash
python casptoolkit/CASP/sup_assemble.py source.pdb output.pdb --target_dir target_dir/ --renumber_atoms
python casptoolkit/CASP/sup_assemble.py source.pdb output.pdb --target_file target_complex.pdb --renumber_atoms
```

4) 同聚体叠合（An -> Am）

```bash
python casptoolkit/CASP/sup_homooligo.py source_An.pdb target_Am.pdb output_dir/ --output_prefix sup_
```

5) AF3 QA 流程

```bash
python casptoolkit/CASP/qa_af3.py -h
```

## 建议

- 优先使用 `-h` 查看脚本参数。
- 并行参数从小值开始（如 `--num_workers 2`），根据机器资源调整。
