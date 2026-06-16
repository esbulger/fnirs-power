#!/usr/bin/env python
"""
Thin Slurm-facing entry script for cluster-based HbO power simulations.

Small test version:
- 3 x 5 = 15 block/subject grid cells
- reduced iterations

Choose the scenario by setting SCENARIO below.
Supported:
    - "null_vs_null"
    - "signal_vs_null"
    - "signal_vs_signal"
"""

import sys
from pathlib import Path

import numpy as np

# Add the project root so the sibling fnirspower_py package is importable.
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from fnirspower_py.sim_data import SimulationConfig
from fnirspower_py.sim_runner import (
    get_slurm_array_pair,
    get_slurm_cpus,
    run_cluster_power_job,
)

if __name__ == "__main__":
    # -----------------------------------------------------------------
    # User inputs
    # -----------------------------------------------------------------
    date_tag = "9999-99-99"

    data_dir = Path("/path/to/python/data/")
    save_dir = Path("/path/to/python/results/")

    SCENARIO = "signal_vs_null"

    # ROI settings. ROI_2 is only used for signal_vs_signal contrasts.
    ROI_1 = 0
    ROI_2 = 1

    # Small test grid: 3 x 5 = 15 total Slurm jobs.
    blocks = np.array([5, 20, 40], dtype=float)
    subject_grid = np.array([5, 15, 30, 45, 60], dtype=int)

    config = SimulationConfig(
        data_mat_file=data_dir / f"{date_tag}_measurement_prediction.mat",
        data_var="M.chan.HbO",
        channel_names_file=data_dir / "dense_channel_names.mat",
        channel_positions_file=data_dir / "dense_channel_positions.mat",
        save_dir=save_dir,
        date_tag=date_tag,
        scenario=SCENARIO,
        blocks=blocks,
        subject_grid=subject_grid,
        noise_amplitude=0.2,
        n_iter=1000,
        n_permutations=256,
    )

    # Specify scenario: null, signal, or contrast.
    if SCENARIO == "null_vs_null":
        config.roi_idx_1 = ROI_1
        config.roi_idx_2 = ROI_1
        config.title = f"null_hbo_roi{config.roi_idx_1:02d}_power_test"

    elif SCENARIO == "signal_vs_null":
        config.roi_idx_1 = ROI_1
        config.roi_idx_2 = ROI_1
        config.signal_scale_1 = 1.0
        config.title = (
            f"signal_vs_null_hbo_roi{config.roi_idx_1:02d}_power_test"
        )

    elif SCENARIO == "signal_vs_signal":
        config.roi_idx_1 = ROI_1
        config.roi_idx_2 = ROI_2
        config.signal_scale_1 = 1.0
        config.signal_scale_2 = 1.0
        config.title = (
            f"signal_vs_signal_hbo_roi{config.roi_idx_1:02d}"
            f"_vs_roi{config.roi_idx_2:02d}_power_test"
        )

    else:
        raise ValueError(
            f"Unsupported SCENARIO: {SCENARIO}. Use 'null_vs_null', "
            "'signal_vs_null', or 'signal_vs_signal'."
        )

    task_id, block_idx, subj_idx = get_slurm_array_pair(
        len(config.blocks),
        len(config.subject_grid),
    )
    cpus = get_slurm_cpus(1)

    print(
        f"[Task {task_id}] scenario={config.scenario}, title={config.title}, "
        f"block_idx={block_idx} (block={config.blocks[block_idx]}), "
        f"subj_idx={subj_idx} (N={config.subject_grid[subj_idx]}), "
        f"cpus={cpus}"
    )

    out = run_cluster_power_job(
        config,
        block_idx=block_idx,
        subj_idx=subj_idx,
        cpus=cpus,
    )

    print(
        f"Saved {out['output_file']} "
        f"in {out['elapsed_sec']:.1f}s"
    )
