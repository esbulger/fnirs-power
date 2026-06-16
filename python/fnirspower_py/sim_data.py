"""
Utilities for loading MATLAB-exported simulation data and building
MNE-compatible channel geometry / adjacency for NIRS cluster testing.

Supports both standard MATLAB .mat files readable by scipy.io.loadmat and
MATLAB v7.3 files, which are HDF5-backed and require h5py.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.io import loadmat
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
from mne import create_info
import mne
from mne._fiff.constants import FIFF
from mne.viz.utils import apply_trans, _check_sphere
from mne.viz.topomap import _get_pos_outlines


@dataclass
class SimulationConfig:
    """
    Configuration for one family of cluster-based power simulations.
    """

    # Data inputs
    data_mat_file: Path
    data_var: str
    channel_names_file: Path
    channel_names_var: str = "channel_names"
    channel_positions_file: Optional[Path] = None
    channel_positions_var: str = "channel_positions"

    # Output
    save_dir: Path = Path("utils")
    date_tag: str = ""
    title: str = ""

    # Scenario
    scenario: str = "signal_vs_null"  # "null_vs_null", "signal_vs_null", "signal_vs_signal"

    # Group 1 signal definition
    roi_idx_1: int = 0
    signal_scale_1: float = 1.0

    # Group 2 signal definition
    roi_idx_2: Optional[int] = None
    signal_scale_2: float = 1.0

    # Simulation grid
    blocks: np.ndarray = field(
        default_factory=lambda: np.linspace(5, 40, 8)
    )
    subject_grid: np.ndarray = field(
        default_factory=lambda: np.array([5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100], dtype=int)
    )

    # Simulation / stats
    noise_amplitude: float = 0.2
    alpha: float = 0.05
    n_iter: int = 1000
    n_permutations: int = 1024
    adjacency_threshold_mm: float = 20.0
    sfreq: float = 5.0
    random_seed: int = 42

    def __post_init__(self):
        self.data_mat_file = Path(self.data_mat_file)
        self.channel_names_file = Path(self.channel_names_file)
        if self.channel_positions_file is not None:
            self.channel_positions_file = Path(self.channel_positions_file)
        self.save_dir = Path(self.save_dir)

        valid_scenarios = {"null_vs_null", "signal_vs_null", "signal_vs_signal"}
        if self.scenario not in valid_scenarios:
            raise ValueError(
                "scenario must be one of {0}; got {1}".format(
                    sorted(valid_scenarios), self.scenario
                )
            )

        self.blocks = np.asarray(self.blocks, dtype=float)
        self.subject_grid = np.asarray(self.subject_grid, dtype=int)

        if self.blocks.ndim != 1 or self.blocks.size == 0:
            raise ValueError("blocks must be a non-empty 1D array")

        if self.subject_grid.ndim != 1 or self.subject_grid.size == 0:
            raise ValueError("subject_grid must be a non-empty 1D array")

        if self.n_iter <= 0:
            raise ValueError("n_iter must be positive")

        if self.n_permutations <= 0:
            raise ValueError("n_permutations must be positive")

    @property
    def output_name_core(self):
        title = self.title.strip().replace(" ", "_")
        if title:
            return title
        return "results_{0}".format(self.scenario)

    def make_output_filename(self, block_idx, subj_idx, suffix=".pkl"):
        return self.save_dir / "{0}_{1}_block{2:02d}_subj{3:02d}{4}".format(
            self.date_tag,
            self.output_name_core,
            block_idx,
            subj_idx,
            suffix,
        )


@dataclass
class PredictionDataset:
    """Container for MATLAB-exported channel-level HbO prediction data."""
    subj_data: np.ndarray  # subjects x rois x channels
    ch_names: list
    ch_pos_m: np.ndarray  # channels x xyz, in meters
    n_subject_pool: int
    n_rois: int
    n_channels: int


def _matlab_str_array_to_list(arr):
    """Convert MATLAB-loaded char/cell/string arrays into a Python list."""
    if isinstance(arr, str):
        return [arr.strip()]

    if isinstance(arr, (list, tuple)):
        out = []
        for x in arr:
            if isinstance(x, str):
                out.append(x.strip())
            else:
                out.extend(_matlab_str_array_to_list(x))
        return out

    arr = np.asarray(arr)

    if arr.dtype.kind in {"U", "S"}:
        if arr.ndim == 0:
            return [str(arr.item()).strip()]
        if arr.ndim == 1:
            return [str(x).strip() for x in arr]
        return ["".join(row).strip() for row in arr]

    if arr.dtype == object:
        out = []
        for x in arr.ravel():
            if isinstance(x, str):
                out.append(x.strip())
            else:
                x = np.asarray(x)
                if x.dtype.kind in {"U", "S"}:
                    out.append("".join(x.ravel()).strip())
                else:
                    out.append(str(x).strip())
        return out

    raise TypeError(f"Could not convert channel names with dtype {arr.dtype}")


def _get_nested_mat_field(mat, field_path):
    """
    Extract a nested MATLAB field from a scipy.loadmat dictionary.

    Example:
        field_path = "M.chan.HbO"
    """
    obj = mat

    for part in field_path.split("."):
        if isinstance(obj, dict):
            if part not in obj:
                raise KeyError(f'Missing field "{part}" while resolving "{field_path}"')
            obj = obj[part]
        elif hasattr(obj, part):
            obj = getattr(obj, part)
        else:
            raise KeyError(
                f'Cannot access field "{part}" while resolving "{field_path}". '
                f"Current object type is {type(obj)}"
            )

    return obj


def _decode_hdf5_attr(value):
    """Decode MATLAB/HDF5 attributes that may be bytes or NumPy scalars."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.ndarray) and value.size == 1:
        return _decode_hdf5_attr(value.item())
    return value


def _decode_matlab_char_codes(arr):
    """Decode a MATLAB v7.3 char dataset into Python string(s)."""
    arr = np.asarray(arr)

    # MATLAB v7.3 stores arrays with reversed dimensions relative to MATLAB.
    if arr.ndim > 1:
        arr = np.transpose(arr, axes=tuple(range(arr.ndim - 1, -1, -1)))

    def decode_one(vec):
        chars = []
        for code in np.asarray(vec).ravel():
            code = int(code)
            if code != 0:
                chars.append(chr(code))
        return "".join(chars).strip()

    if arr.ndim <= 1:
        return decode_one(arr)

    return [decode_one(row) for row in arr]


def _is_hdf5_reference_array(h5py, arr):
    """Return True if an HDF5 dataset array stores object references."""
    return h5py.check_dtype(ref=np.asarray(arr).dtype) is not None


def _hdf5_dataset_to_python(h5_file, obj):
    """Convert one h5py Dataset/Group into a Python object or NumPy array."""
    import h5py

    if isinstance(obj, h5py.Group):
        return {key: _hdf5_dataset_to_python(h5_file, obj[key]) for key in obj.keys()}

    arr = obj[()]
    matlab_class = _decode_hdf5_attr(obj.attrs.get("MATLAB_class", None))

    if _is_hdf5_reference_array(h5py, arr):
        refs = np.asarray(arr)
        if refs.ndim > 1:
            refs = np.transpose(refs, axes=tuple(range(refs.ndim - 1, -1, -1)))
        out = []
        for ref in refs.ravel():
            if ref:
                out.append(_hdf5_dataset_to_python(h5_file, h5_file[ref]))
        return out

    if matlab_class == "char":
        return _decode_matlab_char_codes(arr)

    arr = np.asarray(arr)

    # MATLAB v7.3 stores numeric arrays with dimensions reversed relative to MATLAB.
    if arr.ndim > 1:
        arr = np.transpose(arr, axes=tuple(range(arr.ndim - 1, -1, -1)))

    return arr


def _get_nested_hdf5_field(file_path, field_path):
    """
    Extract a nested MATLAB v7.3/HDF5 field.

    Example:
        field_path = "M.chan.HbO"
    """
    try:
        import h5py
    except ImportError as err:
        raise ImportError(
            "This .mat file is MATLAB v7.3/HDF5. Install h5py to read it: "
            "pip install h5py"
        ) from err

    file_path = Path(file_path)

    with h5py.File(file_path, "r") as h5_file:
        obj = h5_file
        for part in field_path.split("."):
            if isinstance(obj, h5py.Dataset):
                arr = obj[()]
                if _is_hdf5_reference_array(h5py, arr):
                    refs = np.asarray(arr).ravel()
                    if refs.size != 1:
                        raise KeyError(
                            f'Cannot descend into non-scalar reference field "{part}" '
                            f'while resolving "{field_path}"'
                        )
                    obj = h5_file[refs[0]]
                else:
                    raise KeyError(
                        f'Cannot descend into numeric dataset while resolving "{field_path}"'
                    )

            if not isinstance(obj, h5py.Group):
                raise KeyError(
                    f'Cannot access field "{part}" while resolving "{field_path}". '
                    f"Current object type is {type(obj)}"
                )

            if part not in obj:
                available = list(obj.keys())
                raise KeyError(
                    f'Missing field "{part}" while resolving "{field_path}" in {file_path}. '
                    f"Available fields here are: {available}"
                )
            obj = obj[part]

        return _hdf5_dataset_to_python(h5_file, obj)


def _load_mat_field_any_version(file_path, field_path):
    """
    Load one variable/field from either a standard MAT file or MATLAB v7.3 file.

    scipy.io.loadmat cannot read v7.3 files, so this function falls back to
    h5py when SciPy raises the v7.3 NotImplementedError.
    """
    file_path = Path(file_path)

    try:
        mat = loadmat(file_path, simplify_cells=True)
        return _get_nested_mat_field(mat, field_path)

    except (NotImplementedError, ValueError) as err:
        # MATLAB v7.3 files are HDF5-backed. SciPy usually raises
        # NotImplementedError for true MATLAB v7.3 files, but some HDF5
        # files raise ValueError before SciPy recognizes the MAT version.
        try:
            import h5py
            is_hdf5 = h5py.is_hdf5(file_path)
        except ImportError:
            is_hdf5 = False

        if not is_hdf5:
            raise

        return _get_nested_hdf5_field(file_path, field_path)


def load_prediction_dataset(
        data_mat_file,
        data_var,
        channel_names_file,
        channel_names_var,
        channel_positions_file,
        channel_positions_var,
):
    """
    Load MATLAB-exported simulation data and channel metadata.

    Expects channel-level HbO data with shape:
        (n_subject_pool, n_rois, n_channels)

    Example:
        data_var = "M.chan.HbO"

    Supports both standard MAT files and MATLAB v7.3/HDF5 MAT files.
    """
    data_mat_file = Path(data_mat_file)
    channel_names_file = Path(channel_names_file)
    channel_positions_file = Path(channel_positions_file)

    try:
        subj_data = _load_mat_field_any_version(data_mat_file, data_var)
    except KeyError as err:
        raise KeyError(
            f'Could not resolve data_var="{data_var}" in {data_mat_file}'
        ) from err

    subj_data = np.asarray(subj_data, dtype=float)

    if subj_data.ndim != 3:
        raise ValueError(
            "Expected subj_data to be 3D with shape "
            f"(subjects, rois, channels); got shape {subj_data.shape}. "
            "For MATLAB v7.3 files, check whether the HDF5 dimension order was "
            "converted correctly."
        )

    n_subject_pool, n_rois, n_channels = subj_data.shape

    try:
        ch_names_raw = _load_mat_field_any_version(channel_names_file, channel_names_var)
    except KeyError as err:
        raise KeyError(
            f'Could not resolve channel_names_var="{channel_names_var}" in '
            f"{channel_names_file}"
        ) from err

    ch_names = _matlab_str_array_to_list(ch_names_raw)

    try:
        ch_pos_raw = _load_mat_field_any_version(channel_positions_file, channel_positions_var)
    except KeyError as err:
        raise KeyError(
            f'Could not resolve channel_positions_var="{channel_positions_var}" in '
            f"{channel_positions_file}"
        ) from err

    ch_pos_m = np.asarray(ch_pos_raw, dtype=float)

    if ch_pos_m.ndim != 2:
        raise ValueError(f"Channel positions must be 2D; got shape {ch_pos_m.shape}")

    if ch_pos_m.shape[0] == 3 and ch_pos_m.shape[1] != 3:
        ch_pos_m = ch_pos_m.T

    ch_pos_m = ch_pos_m / 1000.0

    if len(ch_names) != n_channels:
        raise ValueError(
            f"Channel name count ({len(ch_names)}) does not match data channels ({n_channels}). "
            f"Data shape is {subj_data.shape}."
        )

    if ch_pos_m.shape != (n_channels, 3):
        raise ValueError(
            f"Channel positions must have shape ({n_channels}, 3); got {ch_pos_m.shape}"
        )

    return PredictionDataset(
        subj_data=subj_data,
        ch_names=ch_names,
        ch_pos_m=ch_pos_m,
        n_subject_pool=n_subject_pool,
        n_rois=n_rois,
        n_channels=n_channels,
    )


def create_nirs_info(ch_names, ch_pos_m, sfreq=5.0):
    """Create an MNE Info object containing HbO channel positions."""
    if len(ch_names) != len(ch_pos_m):
        raise ValueError(
            "ch_names and ch_pos_m must contain the same number of entries: "
            f"received {len(ch_names)} channel names and "
            f"{len(ch_pos_m)} positions."
        )

    info = create_info(
        ch_names=ch_names,
        ch_types="hbo",
        sfreq=sfreq,
    )

    for i, loc in enumerate(ch_pos_m):
        try:
            position = tuple(loc)

            if len(position) != 3:
                raise ValueError(
                    f"expected three coordinates, received {len(position)}"
                )

            info["chs"][i]["loc"][0:3] = position

        except (TypeError, ValueError, IndexError) as exc:
            raise ValueError(
                f"Could not assign coordinates for channel "
                f"{ch_names[i]!r} at index {i}: {exc}"
            ) from exc

    return info


def jitter_duplicate_channel_locations(info, tol=6, eps=1e-8, seed=42):
    """Jitter overlapping channel locations slightly to avoid numerical issues."""
    rng = np.random.RandomState(seed)
    locs = np.array([ch["loc"][0:3] for ch in info["chs"]])
    rounded = np.round(locs, tol)
    _, inv_idx, counts = np.unique(
        rounded, axis=0, return_inverse=True, return_counts=True
    )

    for grp in np.where(counts > 1)[0]:
        dup = np.where(inv_idx == grp)[0]
        for idx in dup[1:]:
            info["chs"][idx]["loc"][0:3] = locs[idx] + rng.randn(3) * eps

    return info


def compute_channel_positions(info, n_channels):
    """Compute 2D sensor positions used by MNE topographic utilities."""
    dev_head_t = info["dev_head_t"]
    chs = info["chs"][:]
    pos = np.empty((len(chs), 3))

    for ci, ch in enumerate(chs):
        pos[ci] = ch["loc"][:3]
        if ch["coord_frame"] == FIFF.FIFFV_COORD_DEVICE:
            if dev_head_t is None:
                dev_head_t = np.eye(4)
            pos[ci] = apply_trans(dev_head_t, pos[ci])

    sphere = _check_sphere(None, info)
    picks = np.arange(n_channels)
    pos, outlines = _get_pos_outlines(info, picks, sphere, to_sphere=True)
    return pos, outlines


def build_nirs_geometry(ch_names, ch_pos_m, sfreq=5.0, jitter_duplicates=True):
    """Build MNE Info plus projected sensor positions for NIRS channels."""
    info = create_nirs_info(ch_names, ch_pos_m, sfreq=sfreq)
    if jitter_duplicates:
        info = jitter_duplicate_channel_locations(info)
    pos, outlines = compute_channel_positions(info, len(ch_names))
    return info, pos, outlines


def build_nirs_adjacency(info, pos, threshold_mm=20.0):
    """Build a distance-thresholded channel adjacency matrix for clustering."""
    adjacency, _ = mne.channels.find_ch_adjacency(info, ch_type=None)
    adjacency_array = adjacency.toarray()
    distances_mm = squareform(pdist(pos)) * 1000.0
    adjacency_array[distances_mm > threshold_mm] = 0
    return csr_matrix(adjacency_array)
