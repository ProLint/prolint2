"""Color configuration and utilities.

This module provides color schemes, amino acid colors, and color
interpolation functions for ProLint visualizations.
"""

import json
from pathlib import Path
from typing import Dict, List, Tuple


def load_theme() -> Dict:
    """Load the ProLint color theme from JSON file.

    Returns
    -------
    dict
        Theme configuration containing colors, gradients, and amino acid mappings.
    """
    theme_path = Path(__file__).parent / "theme.json"
    with open(theme_path, "r") as f:
        return json.load(f)


# Load theme at module import
_THEME = load_theme()

# Export commonly used theme sections
COLORS: Dict = _THEME["colors"]
COLOR_SCALES: Dict[str, List[str]] = _THEME["colorScales"]
GRADIENTS: Dict = _THEME["gradients"]
AMINO_ACID_COLORS: Dict[str, str] = _THEME["aminoAcidColors"]
AMINO_ACID_ONE_LETTER: Dict[str, str] = _THEME["aminoAcidOneLetter"]
UNIT_LABELS: Dict[str, str] = _THEME["unitLabels"]


def hex_to_rgb(hex_color: str) -> Tuple[int, int, int]:
    """Convert hex color string to RGB tuple.

    Parameters
    ----------
    hex_color : str
        Hex color string (e.g., "#FF5733" or "FF5733").

    Returns
    -------
    tuple of int
        RGB values as (red, green, blue), each 0-255.
    """
    hex_color = hex_color.lstrip("#")
    return (
        int(hex_color[0:2], 16),
        int(hex_color[2:4], 16),
        int(hex_color[4:6], 16),
    )


def color_to_tuple(
    color: Dict, normalize: bool = True
) -> Tuple[float, float, float, float]:
    """Convert color dict to RGBA tuple.

    Parameters
    ----------
    color : dict
        Color dict with "r", "g", "b", "a" keys (RGB in 0-255, alpha in 0-1).
    normalize : bool, default=True
        If True, normalize RGB values to 0-1 range. Alpha is unchanged.

    Returns
    -------
    tuple of float
        RGBA values as (red, green, blue, alpha).
    """
    r, g, b, a = color["r"], color["g"], color["b"], color["a"]
    if normalize:
        return (r / 255.0, g / 255.0, b / 255.0, a)
    return (float(r), float(g), float(b), a)


def interpolate_gradient(
    value: float,
    min_val: float,
    max_val: float,
    gradient_name: str = "sharedContacts",
) -> Tuple[float, float, float, float]:
    """Interpolate a color from a gradient based on value.

    Parameters
    ----------
    value : float
        Value to map to a color.
    min_val : float
        Minimum value of the range.
    max_val : float
        Maximum value of the range.
    gradient_name : str, default="sharedContacts"
        Name of the gradient to use from GRADIENTS.

    Returns
    -------
    tuple of float
        Interpolated RGBA color values (0-1 range).
    """
    gradient = GRADIENTS.get(gradient_name, GRADIENTS["sharedContacts"])

    # Normalize value to 0-1 range
    if max_val > min_val:
        normalized = (value - min_val) / (max_val - min_val)
    else:
        normalized = 0.5
    normalized = max(0.0, min(1.0, normalized))

    # Find the two stops to interpolate between
    for i in range(len(gradient) - 1):
        pos1 = gradient[i]["position"]
        pos2 = gradient[i + 1]["position"]

        if pos1 <= normalized <= pos2:
            t = (normalized - pos1) / (pos2 - pos1) if pos2 > pos1 else 0

            c1 = gradient[i]["color"]
            c2 = gradient[i + 1]["color"]

            # Interpolate each channel
            r = c1["r"] + (c2["r"] - c1["r"]) * t
            g = c1["g"] + (c2["g"] - c1["g"]) * t
            b = c1["b"] + (c2["b"] - c1["b"]) * t
            a = c1["a"] + (c2["a"] - c1["a"]) * t

            return (r / 255.0, g / 255.0, b / 255.0, a)

    # Fallback to last color
    return color_to_tuple(gradient[-1]["color"])


def get_color_for_value(
    value: float, scale_name: str = "prolint", vmin: float = 0, vmax: float = 1
) -> str:
    """Get a hex color for a value from a discrete color scale.

    Parameters
    ----------
    value : float
        Value to map to a color.
    scale_name : str, default="prolint"
        Name of the color scale from COLOR_SCALES.
    vmin : float, default=0
        Minimum value of the range.
    vmax : float, default=1
        Maximum value of the range.

    Returns
    -------
    str
        Hex color string.
    """
    colors = COLOR_SCALES.get(scale_name, COLOR_SCALES["prolint"])
    normalized = (value - vmin) / (vmax - vmin) if vmax > vmin else 0.5
    normalized = max(0, min(1, normalized))
    idx = int(normalized * (len(colors) - 1))
    return colors[idx]


def get_unit_label(unit: str) -> str:
    """Get display label for a time unit.

    Parameters
    ----------
    unit : str
        Time unit code (e.g., "us", "ns", "ps").

    Returns
    -------
    str
        Human-readable unit label (e.g., "μs", "ns", "ps").
    """
    return UNIT_LABELS.get(unit, unit)
