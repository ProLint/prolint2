import os
from collections import defaultdict
import plotly.graph_objects as go
import numpy as np

class LipidProteinPlotter:
    def __init__(self, residue_index, interactions, gap=200, ylabel=None, title=None):
        self.residue_index = residue_index
        self.interactions = interactions
        self.gap = gap
        self.ylabel = ylabel if ylabel else "Interactions"
        self.title = title

    def plot_residue_data(self, fn=None):
        bar_color = "#176BA0"

        if fn is None:
            fn = os.path.join(os.getcwd(), "Figure_interactions.html")

        # check for chain breaks
        gray_areas = defaultdict(list)  # show grey area to indicate missing residues
        chain_starts = [0]  # plot in separate figures if the gap between two adjacent residues is larger than 50
        for idx in np.arange(1, len(self.residue_index)):
            if self.residue_index[idx] - self.residue_index[idx - 1] < 0:
                chain_starts.append(idx)
            elif self.residue_index[idx] - self.residue_index[idx - 1] > self.gap:
                chain_starts.append(idx)
            elif 1 < self.residue_index[idx] - self.residue_index[idx - 1] <= self.gap:
                gray_areas[chain_starts[-1]].append([self.residue_index[idx - 1] + 1, self.residue_index[idx] - 1])
        chain_starts.append(len(self.residue_index))

        # plot
        for chain_idx in np.arange(len(chain_starts[:-1])):
            df = self.interactions[chain_starts[chain_idx]:chain_starts[chain_idx + 1]]
            resi_selected = self.residue_index[chain_starts[chain_idx]:chain_starts[chain_idx + 1]]

            fig = go.Figure()
            fig.update_layout(
                yaxis=dict(title=self.ylabel),
                xaxis=dict(title="Residue Index"),
                title=self.title,
                template="plotly_white"
            )

            fig.add_trace(go.Bar(x=resi_selected, y=df, marker_color=bar_color))

            # plot missing residue area
            if chain_starts[chain_idx] in gray_areas.keys():
                for gray_area in gray_areas[chain_starts[chain_idx]]:
                    fig.add_shape(
                        type="rect",
                        xref="x",
                        yref="paper",
                        x0=gray_area[0],
                        y0=0,
                        x1=gray_area[1],
                        y1=1,
                        fillcolor="#c0c0c0",
                        opacity=0.3,
                        layer="below"
                    )

            fig.update_xaxes(
                tickvals=resi_selected,
                tickmode="array",
                ticktext=resi_selected,
                tickangle=45,
                tickfont=dict(size=8)
            )
            fig.update_yaxes(tickfont=dict(size=8))

            if len(chain_starts) == 2:
                fig.write_html(fn)
            else:
                name, ext = os.path.splitext(fn)
                fig.write_html("{}_{}{}".format(name, chain_idx, ext))

