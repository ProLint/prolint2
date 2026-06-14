"""ProLint theme and color utilities.

This module provides color schemes, gradients, and styling functions
for consistent visualization across ProLint plots.
"""

from typing import Tuple

# Import colors and utilities from the single source of truth
from prolint.config.colors import (
    COLORS,
    COLOR_SCALES,
    GRADIENTS,
    AMINO_ACID_COLORS,
    AMINO_ACID_ONE_LETTER,
    UNIT_LABELS,
    hex_to_rgb,
    interpolate_gradient,
    get_color_for_value,
    get_unit_label,
)


# =============================================================================
# Matplotlib Color Map Utilities
# =============================================================================


def get_colormap(name: str = "viridis"):
    """Get a matplotlib colormap by name.

    Parameters
    ----------
    name : str, default="viridis"
        Color scale name from COLOR_SCALES.

    Returns
    -------
    LinearSegmentedColormap
        Matplotlib colormap object.

    Raises
    ------
    ValueError
        If name is not a valid color scale.
    """
    from matplotlib.colors import LinearSegmentedColormap

    if name not in COLOR_SCALES:
        raise ValueError(
            f"Unknown color scale: {name}. Available: {list(COLOR_SCALES.keys())}"
        )

    colors = COLOR_SCALES[name]
    return LinearSegmentedColormap.from_list(f"prolint_{name}", colors)


def interpolate_color(
    value: float,
    min_val: float,
    max_val: float,
    gradient_name: str = "sharedContacts",
) -> Tuple[float, float, float, float]:
    """Interpolate a color from a gradient based on value.

    Parameters
    ----------
    value : float
        Value to interpolate.
    min_val : float
        Minimum value in range.
    max_val : float
        Maximum value in range.
    gradient_name : str, default="sharedContacts"
        Gradient name from GRADIENTS.

    Returns
    -------
    tuple of (float, float, float, float)
        RGBA color tuple with values in 0-1 range.
    """
    return interpolate_gradient(value, min_val, max_val, gradient_name)


# =============================================================================
# Plot Styling Defaults
# =============================================================================

PLOT_STYLE = {
    "figure.facecolor": COLORS["background"]["default"],
    "axes.facecolor": COLORS["background"]["subtle"],
    "axes.edgecolor": COLORS["border"]["default"],
    "axes.labelcolor": COLORS["text"]["primary"],
    "axes.titlecolor": COLORS["text"]["primary"],
    "xtick.color": COLORS["text"]["secondary"],
    "ytick.color": COLORS["text"]["secondary"],
    "grid.color": COLORS["border"]["light"],
    "text.color": COLORS["text"]["primary"],
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
}


def apply_prolint_style():
    """Apply ProLint plotting style to matplotlib.

    Updates matplotlib rcParams with ProLint theme colors and fonts.
    Called automatically by plotter classes.
    """
    import matplotlib.pyplot as plt

    plt.rcParams.update(PLOT_STYLE)


# =============================================================================
# Re-export everything for backwards compatibility
# =============================================================================

__all__ = [
    # Core colors
    "COLORS",
    "COLOR_SCALES",
    "GRADIENTS",
    "AMINO_ACID_COLORS",
    "AMINO_ACID_ONE_LETTER",
    "UNIT_LABELS",
    # Color utilities
    "hex_to_rgb",
    "interpolate_gradient",
    "interpolate_color",
    "get_color_for_value",
    "get_unit_label",
    # Matplotlib utilities
    "get_colormap",
    "apply_prolint_style",
    "PLOT_STYLE",
]
