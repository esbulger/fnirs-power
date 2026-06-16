"""Reusable simulation runner for one Slurm array task / one grid cell."""

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import pickle
import time
import numpy as np

from .sim_data import load_prediction_dataset, build_nirs_geometry, build_nirs_adjacency
from .cluster_stats import run_cluster_permutation_test

def _make_seed(base_seed, stream_id, *indices):
    """
    Create a reproducible uint32 seed from a base seed and simulation indices.

    stream_id separates random-number streams used for different purposes,
    such as noise generation, subject sampling, and MNE permutations.
    """
    seed_sequence = np.random.SeedSequence(
        [int(base_seed), int(stream_id), *map(int, indices)]
    )

    return int(
        seed_sequence.generate_state(
            1,
            dtype=np.uint32,
        )[0]
    )

def make_noise_matrix(blocks, n_subject_pool, n_rois, n_channels, amplitude, rng):
    """Create block-scaled Gaussian noise realizations for all subjects and ROIs."""
    intra_noise = amplitude / np.sqrt(blocks)
    noise_shape = (len(blocks), n_subject_pool, n_rois, n_channels)
    return rng.randn(*noise_shape) * intra_noise[:, None, None, None]


def generate_groups(
    *,
    mode,
    subj_data,
    noise_matrix_1,
    noise_matrix_2,
    block_idx,
    subj_count,
    n_subject_pool,
    rng,
    roi_idx_1=0,
    roi_idx_2=None,
    signal_scale_1=1.0,
    signal_scale_2=1.0,
):
    """
    Generate one pair of groups for a single Monte Carlo iteration.

    subj_data is expected to have shape:
        subjects x rois x channels
    """
    if noise_matrix_2 is None:
        raise ValueError("noise_matrix_2 must be provided as an independent noise realization.")

    if roi_idx_2 is None:
        roi_idx_2 = roi_idx_1

    if subj_count > n_subject_pool:
        raise ValueError(
            f"subj_count ({subj_count}) exceeds n_subject_pool ({n_subject_pool})"
        )

    sel_null_1 = rng.permutation(n_subject_pool)[:subj_count]
    sel_null_2 = rng.permutation(n_subject_pool)[:subj_count]

    null_1 = noise_matrix_1[block_idx, sel_null_1, roi_idx_1, :]
    null_2 = noise_matrix_2[block_idx, sel_null_2, roi_idx_2, :]

    if mode == "null_vs_null":
        return null_1, null_2

    if mode == "signal_vs_null":
        sel_signal_1 = rng.permutation(n_subject_pool)[:subj_count]
        signal_1 = signal_scale_1 * subj_data[sel_signal_1, roi_idx_1, :]
        return signal_1 + null_1, null_2

    if mode == "signal_vs_signal":
        sel_signal_1 = rng.permutation(n_subject_pool)[:subj_count]
        sel_signal_2 = rng.permutation(n_subject_pool)[:subj_count]

        signal_1 = signal_scale_1 * subj_data[sel_signal_1, roi_idx_1, :]
        signal_2 = signal_scale_2 * subj_data[sel_signal_2, roi_idx_2, :]

        return signal_1 + null_1, signal_2 + null_2

    raise ValueError(f"Unknown scenario: {mode}")


def _one_iteration_worker(args):
    """
    Top-level worker so it can be pickled by ProcessPoolExecutor.
    """
    (
        nn,
        base_seed,
        block_idx,
        subj_idx,
        subj_count,
        scenario,
        subj_data,
        noise_matrix_1,
        noise_matrix_2,
        n_subject_pool,
        adjacency_sparse,
        n_permutations,
        alpha,
        roi_idx_1,
        roi_idx_2,
        signal_scale_1,
        signal_scale_2,
    ) = args

    iter_rng = np.random.RandomState(
        _make_seed(
            base_seed,
            2,
            block_idx,
            subj_idx,
            nn,
        )
    )

    permutation_seed = _make_seed(
        base_seed,
        3,
        block_idx,
        subj_idx,
        nn,
    )

    group1, group2 = generate_groups(
        mode=scenario,
        subj_data=subj_data,
        noise_matrix_1=noise_matrix_1,
        noise_matrix_2=noise_matrix_2,
        block_idx=block_idx,
        subj_count=subj_count,
        n_subject_pool=n_subject_pool,
        rng=iter_rng,
        roi_idx_1=roi_idx_1,
        roi_idx_2=roi_idx_2,
        signal_scale_1=signal_scale_1,
        signal_scale_2=signal_scale_2,
    )

    _, _, pvals, _, _ = run_cluster_permutation_test(
        X=[group1, group2],
        adjacency_sparse=adjacency_sparse,
        n_permutations=n_permutations,
        alpha=alpha,
        n_jobs=1,
        seed=permutation_seed,
    )

    hit = 0
    if pvals.size:
        flat = np.concatenate(pvals) if getattr(pvals, "dtype", None) == "O" else pvals
        hit = int(np.any(flat < alpha))

    return nn, hit


def get_slurm_array_pair(n_blocks, n_subjects):
    """Decode SLURM_ARRAY_TASK_ID into block and subject-grid indices."""
    import os

    total_jobs = n_blocks * n_subjects
    try:
        task_id = int(os.environ["SLURM_ARRAY_TASK_ID"])
    except (KeyError, ValueError):
        raise RuntimeError(
            "SLURM_ARRAY_TASK_ID not set/invalid. Use --array=0-{0}.".format(total_jobs - 1)
        )

    if not (0 <= task_id < total_jobs):
        raise RuntimeError(
            "SLURM_ARRAY_TASK_ID={0} out of range (0–{1}).".format(task_id, total_jobs - 1)
        )

    block_idx = task_id // n_subjects
    subj_idx = task_id % n_subjects
    return task_id, block_idx, subj_idx


def get_slurm_cpus(default=1):
    """Read SLURM_CPUS_PER_TASK with a simple fallback."""
    import os
    return int(os.environ.get("SLURM_CPUS_PER_TASK", default))


def run_cluster_power_job(config, block_idx, subj_idx, cpus=1):
    """
    Run one simulation job for a single (block_idx, subj_idx) pair.
    """
    tic = time.time()

    dataset = load_prediction_dataset(
        data_mat_file=config.data_mat_file,
        data_var=config.data_var,
        channel_names_file=config.channel_names_file,
        channel_names_var=config.channel_names_var,
        channel_positions_file=config.channel_positions_file,
        channel_positions_var=config.channel_positions_var,
    )

    # validation of dataset
    if dataset.subj_data.ndim != 3:
        raise ValueError(
            f"Expected loaded data to be 3D: subjects x rois x channels; got {dataset.subj_data.shape}"
        )

    if not (0 <= config.roi_idx_1 < dataset.n_rois):
        raise ValueError(
            f"config.roi_idx_1={config.roi_idx_1} is out of range for n_rois={dataset.n_rois}"
        )

    roi_idx_2_eff = config.roi_idx_2
    if roi_idx_2_eff is None:
        roi_idx_2_eff = config.roi_idx_1

    if not (0 <= roi_idx_2_eff < dataset.n_rois):
        raise ValueError(
            f"config.roi_idx_2={config.roi_idx_2} is out of range for n_rois={dataset.n_rois}"
        )

    subj_count = int(config.subject_grid[subj_idx])
    if subj_count > dataset.n_subject_pool:
        raise ValueError(
            f"Requested subj_count={subj_count}, but only {dataset.n_subject_pool} subjects are available"
        )

    # set up geometry
    info, pos, outlines = build_nirs_geometry(
        dataset.ch_names,
        dataset.ch_pos_m,
        sfreq=config.sfreq,
    )
    adjacency_sparse = build_nirs_adjacency(
        info,
        pos,
        threshold_mm=config.adjacency_threshold_mm,
    )

    # set simulation parameters and compute within-subject variance
    block = config.blocks[block_idx]
    subj_count = int(config.subject_grid[subj_idx])

    base_seed = int(config.random_seed)

    rng_noise_1 = np.random.RandomState(
        _make_seed(
            base_seed,
            0,
            block_idx,
            subj_idx,
        )
    )

    rng_noise_2 = np.random.RandomState(
        _make_seed(
            base_seed,
            1,
            block_idx,
            subj_idx,
        )
    )

    noise_matrix_1 = make_noise_matrix(
        config.blocks,
        dataset.n_subject_pool,
        dataset.n_rois,
        dataset.n_channels,
        config.noise_amplitude,
        rng_noise_1,
    )
    noise_matrix_2 = make_noise_matrix(
        config.blocks,
        dataset.n_subject_pool,
        dataset.n_rois,
        dataset.n_channels,
        config.noise_amplitude,
        rng_noise_2,
    )

    results = np.zeros((1, config.n_iter), dtype=int)

    worker_args = [
        (
            nn,
            base_seed,
            block_idx,
            subj_idx,
            subj_count,
            config.scenario,
            dataset.subj_data,
            noise_matrix_1,
            noise_matrix_2,
            dataset.n_subject_pool,
            adjacency_sparse,
            config.n_permutations,
            config.alpha,
            config.roi_idx_1,
            config.roi_idx_2,
            config.signal_scale_1,
            config.signal_scale_2,
        )
        for nn in range(config.n_iter)
    ]

    if cpus == 1:
        for args in worker_args:
            nn, hit = _one_iteration_worker(args)
            results[0, nn] = hit
    else:
        with ProcessPoolExecutor(max_workers=cpus) as exe:
            futures = {exe.submit(_one_iteration_worker, args): args[0] for args in worker_args}
            for fut in as_completed(futures):
                nn, hit = fut.result()
                results[0, nn] = hit

    save_dir = Path(config.save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    out_f = config.make_output_filename(block_idx=block_idx, subj_idx=subj_idx)
    with open(out_f, "wb") as f:
        pickle.dump(results, f)

    return {
        "results": results,
        "output_file": str(out_f),
        "elapsed_sec": time.time() - tic,
        "block": float(block),
        "subj_count": subj_count,
    }