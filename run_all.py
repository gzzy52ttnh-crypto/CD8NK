#!/usr/bin/env python3
"""
一键运行全部分析流程
目录结构：
  data/
    adata/      - 数据文件（h5ad, csv, txt等）
    code/       - 分析代码（7个脚本）
    result/     - 输出结果（图表、CSV）
    run_all.py  - 本文件
    台账.md     - 分析台账

运行顺序：
  1. figure1.py           - Fig1: 单细胞图谱与NK-like细胞特征
  2. figure2.py           - Fig2: 克隆命运锁定预测响应
  3. spatial_validation.py - Spatial: 空间转录组验证
  4. figure3.py           - Fig3: 髓系微环境与SPP1+ TAM
  5. figure4.py           - Fig4: 机制链条
  6. figure5.py           - Fig5: 外部队列验证
  7. figure_supplement.py - FigS1-S3: 补充图（数据集概览、敏感性分析、额外验证）
  8. t05_bootstrap_sensitivity.py - Bootstrap内部验证 + 阈值敏感性分析
  9. table1_baseline.py   - Table 1: 患者基线特征表
 10. recompute_5_metrics.py - 早期5项指标重新计算（T04 SEM/LUSC OR/CD8表达/T08 MWU/化疗n）

使用：python3 run_all.py
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, 'code')
RESULT = os.path.join(HERE, 'result')
ADATA = os.path.join(HERE, 'adata')

PYTHON = sys.executable

os.makedirs(RESULT, exist_ok=True)
os.makedirs(ADATA, exist_ok=True)

steps = [
    {'name': 'Fig1: 单细胞图谱与NK-like细胞特征', 'script': 'figure1.py', 'depends': []},
    {'name': 'Fig2: 克隆命运锁定预测响应', 'script': 'figure2.py', 'depends': ['per_patient_metrics.csv']},
    {'name': 'Spatial: 空间转录组验证', 'script': 'spatial_validation.py', 'depends': []},
    {'name': 'Fig3: 髓系微环境与SPP1+ TAM', 'script': 'figure3.py', 'depends': ['per_patient_metrics.csv', 'spatial_patient_scores.csv']},
    {'name': 'Fig4: 机制链条', 'script': 'figure4.py', 'depends': ['myeloid_per_patient.csv', 'per_patient_metrics.csv']},
    {'name': 'Fig5: 外部队列验证', 'script': 'figure5.py', 'depends': []},
    {'name': 'FigS1-S3: 补充图（数据集概览/敏感性/额外验证）', 'script': 'figure_supplement.py',
     'depends': ['per_patient_metrics.csv', 'sig_GSE135222.csv', 'fig5_gene_match.csv']},
    {'name': 'T05 Step 2+3: Bootstrap 1000 + 阈值敏感性', 'script': 't05_bootstrap_sensitivity.py',
     'depends': []},
    {'name': 'Step 4: Table 1 患者基线特征表', 'script': 'table1_baseline.py',
     'depends': []},
    {'name': 'Recompute: 早期5项指标重新计算', 'script': 'recompute_5_metrics.py',
     'depends': []},
]

def check_dependencies(deps):
    for dep in deps:
        path = os.path.join(RESULT, dep)
        if not os.path.exists(path):
            return False, dep
    return True, None

def run_script(name, script):
    print(f"\n{'='*60}", flush=True)
    print(f"正在运行: {name}", flush=True)
    print(f"脚本: {script}", flush=True)
    print('='*60, flush=True)
    
    script_path = os.path.join(CODE, script)
    log_path = os.path.join(RESULT, f'{script.replace(".py", "")}.log')
    
    env = os.environ.copy()
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    
    with open(log_path, 'w') as log:
        proc = subprocess.run(
            [PYTHON, script_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
            timeout=None
        )
        log.write(proc.stdout)
        print(proc.stdout, flush=True)
    
    print(f"\n{'='*60}", flush=True)
    print(f"退出码: {proc.returncode}", flush=True)
    print('='*60, flush=True)
    
    return proc.returncode == 0

def main():
    print(f"数据分析目录: {HERE}", flush=True)
    print(f"Python路径: {PYTHON}", flush=True)
    print(f"步骤总数: {len(steps)}", flush=True)
    
    failed_steps = []
    
    for i, step in enumerate(steps, 1):
        ok, missing = check_dependencies(step['depends'])
        if not ok:
            print(f"\n⚠️ 步骤{i}跳过（缺少依赖: {missing}）", flush=True)
            failed_steps.append(f"步骤{i}: {step['name']}")
            continue
        
        success = run_script(f"步骤{i}: {step['name']}", step['script'])
        if not success:
            print(f"\n❌ 步骤{i}失败", flush=True)
            failed_steps.append(f"步骤{i}: {step['name']}")
    
    print(f"\n{'='*60}", flush=True)
    if failed_steps:
        print("失败的步骤:")
        for f in failed_steps:
            print(f"  - {f}")
        sys.exit(1)
    else:
        print("✅ 所有步骤运行成功!")
        print(f"结果保存在: {RESULT}")
        sys.exit(0)

if __name__ == '__main__':
    main()