"""Network plotters for shared contact visualization.

This module provides network graph visualizations for residue
contact correlations.
"""

from typing import List, Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from prolint.analysis.base import AnalysisResult
from prolint.plotting.base import BasePlotter, PlottingRegistry
from prolint.plotting.theme import COLORS, interpolate_color, apply_prolint_style


MAX_NETWORK_RESIDUES = 100


class NetworkPlotter(BasePlotter):
    """Plotter for residue contact correlation networks.

    Visualizes query residues as nodes connected by edges when they
    share contacts with the same database molecule.

    See Also
    --------
    SharedContactsAnalysis : Generates correlation matrix data
    HeatmapPlotter : Alternative matrix visualization
    """

    name = "network"
    required_analysis = "shared_contacts"
    description = "Network graph of query residues sharing database contacts"

    @classmethod
    def validate_result(cls, result: AnalysisResult) -> None:
        """Validate that result contains required network data."""
        required_keys = ["matrix", "labels"]
        missing = [k for k in required_keys if k not in result.data]
        if missing:
            raise ValueError(
                f"AnalysisResult missing required keys for {cls.name}: {missing}. "
                f"Expected result from '{cls.required_analysis}' analysis."
            )

    @classmethod
    def plot(
        cls,
        result: AnalysisResult,
        threshold: int = 0,
        selected_residues: Optional[List[int]] = None,
        max_residues: int = MAX_NETWORK_RESIDUES,
        figsize: Tuple[float, float] = (10, 10),
        node_size: int = 500,
        font_size: int = 8,
        show_edge_labels: bool = False,
        title: str = "Shared Contacts Network",
        highlight_nodes: Optional[List[int]] = None,
    ) -> Tuple[Figure, Axes]:
        """Create network graph visualization.

        Parameters
        ----------
        result : AnalysisResult
            Result from shared_contacts analysis.
        threshold : int, default=0
            Minimum shared frames to draw an edge.
        selected_residues : list of int, optional
            Subset of residues to display.
        max_residues : int, default=100
            Maximum residues before raising error.
        highlight_nodes : list of int, optional
            Residue IDs to highlight.
        show_edge_labels : bool, default=False
            Whether to show edge weights.

        Returns
        -------
        tuple of (Figure, Axes)
            Matplotlib figure and axes objects.

        Raises
        ------
        ImportError
            If networkx is not installed.
        ValueError
            If too many residues to display.
        """
        try:
            import networkx as nx
        except ImportError:
            raise ImportError(
                "NetworkX is required for network visualization. "
                "Install it with: pip install networkx"
            )

        cls.validate_result(result)
        apply_prolint_style()

        # Extract data from result
        matrix = np.asarray(result.data["matrix"])
        labels = result.data["labels"]

        # Filter to selected residues if specified
        if selected_residues is not None and len(selected_residues) > 0:
            selected_set = set(selected_residues)
            selected_indices = [i for i, l in enumerate(labels) if l in selected_set]
            labels = [labels[i] for i in selected_indices]
            matrix = matrix[np.ix_(selected_indices, selected_indices)]

        n = len(labels)

        # Enforce maximum residue limit to avoid overcrowding
        if n > max_residues:
            raise ValueError(
                f"Number of residues ({n}) exceeds maximum allowed ({max_residues}). "
                f"Please use 'selected_residues' to select a subset of residues to display. "
                f"You can increase 'max_residues' if needed, but this may result in an overcrowded plot."
            )

        # Handle empty case
        if n == 0:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(
                0.5,
                0.5,
                "No nodes to display",
                ha="center",
                va="center",
                fontsize=12,
                color=COLORS["text"]["secondary"],
            )
            ax.axis("off")
            return fig, ax

        # Build nodes
        nodes = [
            {"id": f"res-{labels[i]}", "label": str(labels[i]), "residue_id": labels[i]}
            for i in range(n)
        ]

        # Build edges (only upper triangle, non-zero values above threshold)
        edges = []
        for i in range(n):
            for j in range(i + 1, n):
                value = matrix[i][j]
                if value > threshold:
                    edges.append(
                        {
                            "source": f"res-{labels[i]}",
                            "target": f"res-{labels[j]}",
                            "value": int(value),
                        }
                    )

        # Create NetworkX graph
        G = nx.Graph()

        # Add nodes
        for node in nodes:
            G.add_node(
                node["id"],
                label=node.get("label", node["id"]),
                residue_id=node.get("residue_id"),
            )

        # Calculate edge statistics
        edge_values = [e.get("value", 1) for e in edges]
        max_value = max(edge_values) if edge_values else 1
        min_value = min(edge_values) if edge_values else 0

        # Add edges with attributes
        for edge in edges:
            value = edge.get("value", 1)
            width = 1 + (value / max_value) * 4 if max_value > 0 else 2
            color = interpolate_color(value, min_value, max_value)
            G.add_edge(
                edge["source"], edge["target"], weight=value, width=width, color=color
            )

        # Create figure
        fig, ax = plt.subplots(figsize=figsize)

        # Use circular layout
        pos = nx.circular_layout(G)

        # Draw edges first
        edge_widths = [G[u][v].get("width", 2) for u, v in G.edges()]
        edge_colors = [
            G[u][v].get("color", COLORS["neutral"]["400"]) for u, v in G.edges()
        ]

        nx.draw_networkx_edges(
            G,
            pos,
            ax=ax,
            width=edge_widths,
            edge_color=edge_colors,
            alpha=0.7,
        )

        # Prepare node colors
        highlight_set = set(highlight_nodes or [])
        node_colors = []
        node_sizes = []
        for node_id in G.nodes():
            residue_id = G.nodes[node_id].get("residue_id")
            if residue_id in highlight_set:
                node_colors.append(COLORS["data"]["highlight"])
                node_sizes.append(node_size * 1.5)
            else:
                node_colors.append(COLORS["data"]["query"])
                node_sizes.append(node_size)

        # Draw nodes
        nx.draw_networkx_nodes(
            G,
            pos,
            ax=ax,
            node_color=node_colors,
            node_size=node_sizes,
            edgecolors=COLORS["data"]["query_dark"],
            linewidths=2,
        )

        # Draw labels
        node_labels = {node: G.nodes[node].get("label", node) for node in G.nodes()}
        nx.draw_networkx_labels(
            G,
            pos,
            ax=ax,
            labels=node_labels,
            font_size=font_size,
            font_weight="bold",
            font_color=COLORS["text"]["primary"],
        )

        # Draw edge labels if requested
        if show_edge_labels:
            edge_labels = {(u, v): str(G[u][v].get("weight", "")) for u, v in G.edges()}
            nx.draw_networkx_edge_labels(
                G,
                pos,
                ax=ax,
                edge_labels=edge_labels,
                font_size=font_size - 2,
                font_color=COLORS["text"]["secondary"],
            )

        ax.set_title(title, fontsize=12, fontweight="semibold")
        ax.axis("off")

        # Add stats
        stats_text = f"{len(nodes)} nodes, {len(edges)} edges"
        if edge_values:
            stats_text += f"\nMax shared: {max_value}, Avg: {np.mean(edge_values):.1f}"
        ax.text(
            0.02,
            0.02,
            stats_text,
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment="bottom",
            color=COLORS["text"]["secondary"],
        )

        plt.tight_layout()
        return fig, ax


# Register plotter
PlottingRegistry.register("network", NetworkPlotter)
