#!/usr/bin/env python3
"""
一键运行 NK-like CD8+ T 细胞 NSCLC 抗 PD-1 全部分析流程
所有 18 个步骤按依赖顺序串行执行，固定使用 /opt/anaconda3/bin/python3 环境。
"""
import subprocess
import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CODE_DIR = ROOT / 'code'
RESULT_DIR = ROOT / 'result'
LOG_DIR = RESULT_DIR

PYTHON = '/opt/anaconda3/bin/python3'

STEPS = [
    {'name': 'Step 1 现象发现', 'script': 'step1_phenomenon_discovery.py', 'depends': []},
    {'name': 'Step 2 核心指标构建', 'script': 'step2_core_index_construction.py', 'depends': []},
    {'name': 'Step 2 Firth 惩罚回归验证', 'script': 'step2_firth_validation.py', 'depends': ['step2_core_index_construction.py']},
    {'name': 'Step 2.5 空间转录组验证', 'script': 'step2_5_spatial_validation.py', 'depends': []},
    {'name': 'Step 3 机制探索', 'script': 'step3_mechanism_exploration.py', 'depends': ['step2_core_index_construction.py']},
    {'name': 'Step 3b TCF7 杀伤功能失调', 'script': 'step3b_tcf7_dysfunction.py', 'depends': ['step3_mechanism_exploration.py']},
    {'name': 'Step 4 证据拓展', 'script': 'step4_evidence_generalization.py', 'depends': ['step2_core_index_construction.py', 'step2_5_spatial_validation.py']},
    {'name': 'Step 4-5 外部队列深度验证', 'script': 'step4_5_external_validation.py', 'depends': ['step2_core_index_construction.py']},
    {'name': 'Step 4b 化疗动力学 (GSE179994)', 'script': 'step4b_chemo_dynamics.py', 'depends': []},
    {'name': 'Step 4c 主队列化疗方案分层', 'script': 'step4c_chemo_stratified.py', 'depends': ['step2_core_index_construction.py']},
    {'name': 'Step 4d 化疗机制深度', 'script': 'step4d_chemo_mechanism.py', 'depends': ['step2_core_index_construction.py', 'step3_mechanism_exploration.py']},
    {'name': 'Step 4e cDC2 抗原呈递', 'script': 'step4e_cdc2_mechanism.py', 'depends': ['step3_mechanism_exploration.py']},
    {'name': 'Step 5 临床转化模型', 'script': 'step5_clinical_translation.py', 'depends': ['step2_core_index_construction.py']},
    {'name': 'Step 5.5 IRS 评分构建', 'script': 'step5_5_IRS_construction.py', 'depends': ['step2_core_index_construction.py', 'step3_mechanism_exploration.py']},
    {'name': 'Step 6 GSE241934 外部验证', 'script': 'step6_gse241934_validation.py', 'depends': []},
    {'name': 'Step 6c TCR 克隆分析', 'script': 'step6c_tcr_clonality.py', 'depends': ['step6_gse241934_validation.py']},
    {'name': 'Step 6d 机制验证 (GSE241934)', 'script': 'step6d_mechanism_validation.py', 'depends': ['step6_gse241934_validation.py']},
    {'name': 'Step 6e 克隆多样性 (GSE241934+GSE179994)', 'script': 'step6e_tcr_diversity.py', 'depends': ['step6_gse241934_validation.py', 'step4b_chemo_dynamics.py']},
]


def run_step(idx, step):
    """执行单个步骤"""
    script = CODE_DIR / step['script']
    log_file = LOG_DIR / f"{step['script'].replace('.py', '')}.log"
    print(f"\n[{idx+1}/{len(STEPS)}] {step['name']}")
    print(f"   Script: {step['script']}")
    print(f"   Log: {log_file.name}")
    try:
        result = subprocess.run(
            [PYTHON, str(script)],
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            timeout=3600
        )
        with open(log_file, 'w', encoding='utf-8') as f:
            f.write(f"=== STDOUT ===\n{result.stdout}\n=== STDERR ===\n{result.stderr}\n=== EXIT CODE: {result.returncode} ===\n")
        if result.returncode == 0:
            print(f"   OK (exit 0)")
            return True
        else:
            print(f"   FAILED (exit {result.returncode})")
            print(f"   STDERR (last 500 chars): {result.stderr[-500:]}")
            return False
    except subprocess.TimeoutExpired:
        print(f"   TIMEOUT (3600s)")
        return False
    except Exception as e:
        print(f"   ERROR: {e}")
        return False


def main():
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"ROOT: {ROOT}")
    print(f"PYTHON: {PYTHON}")
    print(f"Total steps: {len(STEPS)}")

    success_count = 0
    fail_count = 0
    failed_steps = []

    for idx, step in enumerate(STEPS):
        if not (CODE_DIR / step['script']).exists():
            print(f"\n[{idx+1}/{len(STEPS)}] SKIP {step['name']} - 脚本不存在: {step['script']}")
            fail_count += 1
            failed_steps.append(step['name'])
            continue
        ok = run_step(idx, step)
        if ok:
            success_count += 1
        else:
            fail_count += 1
            failed_steps.append(step['name'])

    print(f"\n{'='*60}")
    print(f"全部完成: {success_count}/{len(STEPS)} 成功")
    if failed_steps:
        print(f"失败步骤: {', '.join(failed_steps)}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
