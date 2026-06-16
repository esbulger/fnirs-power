"""
Utilities for loading per-(block, subject) result grids and computing
summary values for cluster-based power simulations.
"""

from pathlib import Path
import pickle
import re
import string
import warnings

import numpy as np


def load_pickle_array(file_path):
    """Load one pickle file and return it as a NumPy array."""
    with open(file_path, "rb") as f:
        return np.asarray(pickle.load(f))


def empirical_power_summary(arr):
    """
    Compute empirical power from one result array.

    Works for single-ROI results shaped (1, n_iter) and multi-ROI results
    shaped (n_rois, n_iter). In both cases, this returns the mean hit rate.
    """
    arr = np.asarray(arr)
    if arr.size == 0:
        raise ValueError("Cannot summarize an empty result array")
    return float(arr.mean())


def _pattern_to_regex(file_pattern):
    """
    Convert a Python format-style result filename pattern into a regex.

    Supported fields are {bb} and {ss}, with or without format specs such as
    {bb:02d} and {ss:02d}.
    """
    formatter = string.Formatter()
    regex_parts = []

    for literal_text, field_name, format_spec, conversion in formatter.parse(file_pattern):
        regex_parts.append(re.escape(literal_text))

        if field_name is None:
            continue

        if field_name == "bb":
            regex_parts.append(r"(?P<bb>\d+)")
        elif field_name == "ss":
            regex_parts.append(r"(?P<ss>\d+)")
        else:
            raise ValueError(
                f"Unsupported filename pattern field {{{field_name}}}. "
                "Only {bb} and {ss} are supported."
            )

    return re.compile("^" + "".join(regex_parts) + "$")


def find_result_pattern(data_dir, candidate_patterns):
    """
    Return the first candidate filename pattern that matches at least one file.

    Unlike a simple block00/subj00 check, this scans the directory, so it still
    works if only part of a grid exists or if indices are not contiguous.
    """
    data_dir = Path(data_dir)

    for pattern in candidate_patterns:
        regex = _pattern_to_regex(pattern)
        for f in data_dir.iterdir():
            if f.is_file() and regex.match(f.name):
                return pattern

    raise FileNotFoundError(
        "Could not find any matching result files. Tried:\n  "
        + "\n  ".join(candidate_patterns)
    )


def infer_result_grid(data_dir, file_pattern):
    """
    Infer the block and subject indices represented by result filenames.

    The block and subject indices are collected independently. Their
    Cartesian product defines the rectangular result grid. Some individual
    block-subject combinations may be absent; load_summary_grid represents
    those missing cells with NaN.

    Returns
    -------
    bb_found : list of int
        Sorted block file indices.
    ss_found : list of int
        Sorted subject-grid file indices.
    """
    data_dir = Path(data_dir)
    regex = _pattern_to_regex(file_pattern)

    bb_found = set()
    ss_found = set()

    for f in data_dir.iterdir():
        if not f.is_file():
            continue
        match = regex.match(f.name)
        if match:
            bb_found.add(int(match.group("bb")))
            ss_found.add(int(match.group("ss")))

    if not bb_found or not ss_found:
        raise FileNotFoundError(
            "Could not infer block/subject grid from files matching pattern:\n"
            f"  {file_pattern}"
        )

    return sorted(bb_found), sorted(ss_found)


def load_summary_grid(
    data_dir,
    file_pattern,
    bb_found,
    ss_found,
    summary_fn=empirical_power_summary,
):
    """
    Load result files and compute one summary value per block/subject cell.

    Missing block-subject combinations are represented by NaN and reported
    with a warning. Existing files must all have the same array shape.

    Parameters
    ----------
    data_dir : str or Path
        Directory containing result files.
    file_pattern : str
        Filename pattern using {bb} and {ss} fields.
    bb_found : sequence of int
        Block file indices defining the rows of the output grid.
    ss_found : sequence of int
        Subject-grid file indices defining the columns of the output grid.
    summary_fn : callable
        Function mapping one loaded result array to one scalar.

    Returns
    -------
    summary_grid : ndarray
        Float array with shape
        ``(len(bb_found), len(ss_found))``.

        Missing block-subject combinations contain NaN.
    """
    data_dir = Path(data_dir)
    bb_found = list(bb_found)
    ss_found = list(ss_found)

    if not bb_found:
        raise ValueError("bb_found must be non-empty")
    if not ss_found:
        raise ValueError("ss_found must be non-empty")

    # Find the first existing file rather than assuming the first inferred
    # block-subject combination exists.
    sample_file = None

    for bb in bb_found:
        for ss in ss_found:
            candidate = data_dir / file_pattern.format(bb=bb, ss=ss)

            if candidate.exists():
                sample_file = candidate
                break

        if sample_file is not None:
            break

    if sample_file is None:
        raise FileNotFoundError(
            "No result files were found for the inferred block and "
            "subject indices."
        )

    expected_shape = load_pickle_array(sample_file).shape

    # NaN explicitly represents a missing simulation result.
    summary_grid = np.full(
        (len(bb_found), len(ss_found)),
        np.nan,
        dtype=float,
    )

    missing_files = []

    for row_idx, bb in enumerate(bb_found):
        for col_idx, ss in enumerate(ss_found):
            fname = data_dir / file_pattern.format(bb=bb, ss=ss)

            if not fname.exists():
                missing_files.append(fname.name)
                continue

            arr = load_pickle_array(fname)

            if arr.shape != expected_shape:
                raise ValueError(
                    f"Unexpected shape in {fname}: {arr.shape}, "
                    f"expected {expected_shape}"
                )

            summary_grid[row_idx, col_idx] = summary_fn(arr)

    if missing_files:
        preview_count = 10
        preview = ", ".join(missing_files[:preview_count])

        if len(missing_files) > preview_count:
            preview += (
                f", ... and "
                f"{len(missing_files) - preview_count} more"
            )

        warnings.warn(
            f"{len(missing_files)} result-grid cell(s) are missing and "
            f"were set to NaN: {preview}",
            RuntimeWarning,
            stacklevel=2,
        )

    return summary_grid


def _resolve_axis_labels(labels, found_indices, axis_name):
    """
    Resolve plotting labels for inferred filename indices.

    Labels are used only when they can be aligned safely. Otherwise, inferred
    file indices are used and a warning is emitted.
    """
    found_indices = list(found_indices)

    if labels is None:
        warnings.warn(
            f"No {axis_name} labels provided. "
            f"Using inferred file indices instead: {found_indices}",
            stacklevel=3,
        )
        return np.asarray(found_indices)

    labels = np.asarray(labels)

    # Common case: labels correspond one-to-one to the inferred rows/columns.
    if len(labels) == len(found_indices):
        return labels

    # Lookup-table case: labels are indexed by the original file indices.
    if len(labels) > max(found_indices):
        warnings.warn(
            f"{axis_name} label count ({len(labels)}) does not match the number "
            f"of inferred indices ({len(found_indices)}), but labels can be "
            "indexed by inferred file indices. Using labels[found_indices].",
            stacklevel=3,
        )
        return labels[found_indices]

    warnings.warn(
        f"{axis_name} labels do not match inferred file indices. "
        f"Got {len(labels)} labels for inferred indices {found_indices}. "
        "Using inferred file indices instead.",
        stacklevel=3,
    )
    return np.asarray(found_indices)


def resolve_plot_inputs(blocks, subjects, plot_bbs, bb_found, ss_found):
    """
    Resolve plotting labels and selected block rows.

    Parameters
    ----------
    blocks : array-like or None
        Real block-count labels. If incompatible, file indices are used.
    subjects : array-like or None
        Real subject-count labels. If incompatible, file indices are used.
    plot_bbs : sequence of int or None
        Block file indices to plot. If None, all inferred block indices are plotted.
    bb_found : sequence of int
        Inferred block file indices.
    ss_found : sequence of int
        Inferred subject file indices.

    Returns
    -------
    block_labels : ndarray
        Labels aligned to rows of the summary grid.
    subject_labels : ndarray
        Labels aligned to columns of the summary grid.
    plot_rows : list of int
        Row positions in the summary grid to plot.
    """
    bb_found = list(bb_found)
    ss_found = list(ss_found)

    block_labels = _resolve_axis_labels(blocks, bb_found, "block")
    subject_labels = _resolve_axis_labels(subjects, ss_found, "subject")

    if plot_bbs is None:
        plot_rows = list(range(len(bb_found)))
    else:
        bb_to_row = {bb: row for row, bb in enumerate(bb_found)}
        missing = [bb for bb in plot_bbs if bb not in bb_to_row]
        if missing:
            raise ValueError(
                f"plot_bbs contains block indices not found in files: {missing}. "
                f"Found block indices are {bb_found}."
            )
        plot_rows = [bb_to_row[bb] for bb in plot_bbs]

    return block_labels, subject_labels, plot_rows
