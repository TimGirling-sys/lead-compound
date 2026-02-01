"""
Visualization module for compound evolution pathways.

Creates publication-quality visualizations showing:
- Compound similarity networks
- Lead compound clusters
- Optimization pathways with chemical structures
"""

from typing import List, Dict, Tuple, Optional, Any
from pathlib import Path
import io
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec

import networkx as nx
from PIL import Image

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from .fingerprints import Compound
from .clustering import Cluster
from .pathway import EvolutionPathway, LeadCompound, LeadType


class EvolutionVisualizer:
    """
    Creates visualizations of compound evolution pathways.

    Generates publication-quality figures showing:
    - Network view of compound similarities
    - Cluster membership with lead compounds highlighted
    - Pathway arrows from initial to final leads
    - Chemical structure insets

    Attributes:
        figsize: Default figure size (width, height) in inches
        dpi: Resolution for output images
        style: Visual style ('light' or 'dark')
    """

    # Color schemes for different lead types
    COLORS = {
        'initial': '#2ECC71',      # Green
        'intermediate': '#3498DB',  # Blue
        'final': '#E74C3C',         # Red
        'analog': '#95A5A6',        # Gray
        'edge': '#BDC3C7',          # Light gray
        'pathway': '#9B59B6',       # Purple
        'background': '#FFFFFF',    # White
        'text': '#2C3E50',          # Dark blue-gray
    }

    # Generation colors for visual distinction
    GENERATION_CMAP = LinearSegmentedColormap.from_list(
        'generations',
        ['#2ECC71', '#F1C40F', '#E67E22', '#E74C3C'],  # Green to Red
        N=10
    )

    def __init__(self,
                 figsize: Tuple[int, int] = (16, 12),
                 dpi: int = 150,
                 style: str = 'light'):
        """
        Initialize the visualizer.

        Args:
            figsize: Figure dimensions (width, height) in inches
            dpi: Output resolution
            style: 'light' or 'dark' color scheme
        """
        self.figsize = figsize
        self.dpi = dpi
        self.style = style

        if style == 'dark':
            self.COLORS['background'] = '#1A1A2E'
            self.COLORS['text'] = '#EAEAEA'
            self.COLORS['edge'] = '#4A4A6A'

    def visualize_pathway(self,
                         pathway: EvolutionPathway,
                         compounds: List[Compound],
                         similarity_matrix: np.ndarray,
                         clusters: List[Cluster],
                         output_path: Optional[str] = None,
                         show_structures: bool = True,
                         show_all_compounds: bool = True,
                         title: str = "Compound Evolution Pathway") -> plt.Figure:
        """
        Create comprehensive pathway visualization.

        Generates a multi-panel figure showing:
        1. Main network view with pathway highlighted
        2. Chemical structures of key leads
        3. Legend and annotations

        Args:
            pathway: Evolution pathway to visualize
            compounds: All compounds
            similarity_matrix: Pairwise similarities
            clusters: Compound clusters
            output_path: Optional path to save the figure
            show_structures: Whether to show chemical structures
            show_all_compounds: Show all compounds or just leads
            title: Figure title

        Returns:
            matplotlib Figure object
        """
        # Create figure with custom layout
        fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
        fig.patch.set_facecolor(self.COLORS['background'])

        if show_structures:
            # Layout: main network (left), structures (right)
            gs = gridspec.GridSpec(
                3, 4,
                width_ratios=[3, 1, 1, 1],
                height_ratios=[1, 1, 1],
                hspace=0.3,
                wspace=0.2
            )
            ax_main = fig.add_subplot(gs[:, 0])
            ax_structures = [
                fig.add_subplot(gs[0, 1:]),
                fig.add_subplot(gs[1, 1:]),
                fig.add_subplot(gs[2, 1:])
            ]
        else:
            ax_main = fig.add_subplot(111)
            ax_structures = []

        # Draw main network
        self._draw_network(
            ax_main, pathway, compounds, similarity_matrix,
            clusters, show_all_compounds
        )

        # Draw chemical structures
        if show_structures and ax_structures:
            self._draw_structure_panel(ax_structures, pathway)

        # Add title and legend
        fig.suptitle(
            title,
            fontsize=16,
            fontweight='bold',
            color=self.COLORS['text'],
            y=0.98
        )

        # Add legend
        self._add_legend(fig)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        if output_path:
            fig.savefig(
                output_path,
                dpi=self.dpi,
                bbox_inches='tight',
                facecolor=self.COLORS['background'],
                edgecolor='none'
            )

        return fig

    def _draw_network(self,
                     ax: plt.Axes,
                     pathway: EvolutionPathway,
                     compounds: List[Compound],
                     similarity_matrix: np.ndarray,
                     clusters: List[Cluster],
                     show_all_compounds: bool) -> None:
        """
        Draw the compound similarity network.

        Args:
            ax: Matplotlib axes
            pathway: Evolution pathway
            compounds: All compounds
            similarity_matrix: Pairwise similarities
            clusters: Compound clusters
            show_all_compounds: Whether to show all compounds
        """
        ax.set_facecolor(self.COLORS['background'])

        # Build the visualization graph
        G = nx.Graph()
        compound_to_idx = {c.smiles: i for i, c in enumerate(compounds)}
        lead_smiles = {l.compound.smiles for l in pathway.all_leads}

        # Determine which compounds to include
        if show_all_compounds:
            compounds_to_show = compounds
        else:
            # Only show leads and their direct analogs
            compounds_to_show = []
            for lead in pathway.all_leads:
                compounds_to_show.append(lead.compound)
                lead_idx = compound_to_idx[lead.compound.smiles]
                for i, c in enumerate(compounds):
                    if similarity_matrix[lead_idx, i] >= 0.7:
                        if c not in compounds_to_show:
                            compounds_to_show.append(c)

        # Add nodes
        for compound in compounds_to_show:
            is_lead = compound.smiles in lead_smiles
            G.add_node(
                compound.smiles,
                compound=compound,
                is_lead=is_lead
            )

        # Add edges (only between similar compounds)
        shown_smiles = {c.smiles for c in compounds_to_show}
        for i, c1 in enumerate(compounds_to_show):
            idx1 = compound_to_idx[c1.smiles]
            for j, c2 in enumerate(compounds_to_show[i + 1:], i + 1):
                idx2 = compound_to_idx[c2.smiles]
                sim = similarity_matrix[idx1, idx2]
                if sim >= 0.5:  # Lower threshold for visualization
                    G.add_edge(c1.smiles, c2.smiles, weight=sim)

        if len(G.nodes) == 0:
            ax.text(0.5, 0.5, "No compounds to display",
                   ha='center', va='center', transform=ax.transAxes)
            return

        # Calculate layout - use spring layout with custom parameters
        try:
            pos = nx.spring_layout(
                G,
                k=2.0 / np.sqrt(len(G.nodes)),
                iterations=100,
                weight='weight',
                seed=42
            )
        except Exception:
            pos = nx.kamada_kawai_layout(G)

        # Assign colors and sizes to nodes
        node_colors = []
        node_sizes = []
        node_edge_colors = []
        node_edge_widths = []

        for node in G.nodes():
            is_lead = G.nodes[node].get('is_lead', False)

            if is_lead:
                # Find the lead type
                lead_type = LeadType.INTERMEDIATE
                for lead in pathway.all_leads:
                    if lead.compound.smiles == node:
                        lead_type = lead.lead_type
                        break

                if lead_type == LeadType.INITIAL:
                    node_colors.append(self.COLORS['initial'])
                    node_sizes.append(800)
                elif lead_type == LeadType.FINAL:
                    node_colors.append(self.COLORS['final'])
                    node_sizes.append(1000)
                else:
                    node_colors.append(self.COLORS['intermediate'])
                    node_sizes.append(600)

                node_edge_colors.append(self.COLORS['text'])
                node_edge_widths.append(3)
            else:
                node_colors.append(self.COLORS['analog'])
                node_sizes.append(100)
                node_edge_colors.append(self.COLORS['analog'])
                node_edge_widths.append(0.5)

        # Draw edges first (below nodes)
        edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
        edge_alphas = [0.1 + 0.5 * w for w in edge_weights]

        nx.draw_networkx_edges(
            G, pos, ax=ax,
            edge_color=self.COLORS['edge'],
            alpha=0.3,
            width=[0.5 + 2 * w for w in edge_weights]
        )

        # Draw pathway edges with special styling
        if pathway.pathway_graph:
            pathway_edges = []
            for u, v in pathway.pathway_graph.edges():
                if u in G.nodes() and v in G.nodes():
                    pathway_edges.append((u, v))

            if pathway_edges:
                nx.draw_networkx_edges(
                    G, pos, ax=ax,
                    edgelist=pathway_edges,
                    edge_color=self.COLORS['pathway'],
                    width=4,
                    alpha=0.8,
                    arrows=True,
                    arrowsize=20,
                    arrowstyle='-|>',
                    connectionstyle='arc3,rad=0.1'
                )

        # Draw nodes
        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            node_color=node_colors,
            node_size=node_sizes,
            edgecolors=node_edge_colors,
            linewidths=node_edge_widths,
            alpha=0.9
        )

        # Add labels for lead compounds
        labels = {}
        for lead in pathway.all_leads:
            if lead.compound.smiles in G.nodes():
                labels[lead.compound.smiles] = lead.compound.name

        nx.draw_networkx_labels(
            G, pos, ax=ax,
            labels=labels,
            font_size=8,
            font_color=self.COLORS['text'],
            font_weight='bold'
        )

        # Add cluster hulls
        self._draw_cluster_hulls(ax, clusters, compounds, pos, compound_to_idx)

        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.3, 1.3)
        ax.axis('off')
        ax.set_title(
            "Compound Similarity Network",
            fontsize=12,
            color=self.COLORS['text'],
            pad=10
        )

    def _draw_cluster_hulls(self,
                           ax: plt.Axes,
                           clusters: List[Cluster],
                           compounds: List[Compound],
                           pos: Dict[str, Tuple[float, float]],
                           compound_to_idx: Dict[str, int]) -> None:
        """
        Draw convex hulls around compound clusters.

        Args:
            ax: Matplotlib axes
            clusters: List of clusters
            compounds: All compounds
            pos: Node positions
            compound_to_idx: Compound index mapping
        """
        from scipy.spatial import ConvexHull
        from matplotlib.patches import Polygon

        colors = plt.cm.Set3(np.linspace(0, 1, len(clusters)))

        for i, cluster in enumerate(clusters):
            if len(cluster.compounds) < 3:
                continue

            # Get positions of cluster compounds
            points = []
            for compound in cluster.compounds:
                if compound.smiles in pos:
                    points.append(pos[compound.smiles])

            if len(points) < 3:
                continue

            points = np.array(points)

            try:
                hull = ConvexHull(points)
                hull_points = points[hull.vertices]

                # Expand hull slightly for visual clarity
                center = np.mean(hull_points, axis=0)
                expanded = center + 1.2 * (hull_points - center)

                polygon = Polygon(
                    expanded,
                    alpha=0.1,
                    facecolor=colors[i % len(colors)],
                    edgecolor=colors[i % len(colors)],
                    linewidth=2,
                    linestyle='--',
                    zorder=0
                )
                ax.add_patch(polygon)
            except Exception:
                # Skip if hull calculation fails
                pass

    def _draw_structure_panel(self,
                             axes: List[plt.Axes],
                             pathway: EvolutionPathway) -> None:
        """
        Draw chemical structure panels for lead compounds.

        Args:
            axes: List of axes for structure panels
            pathway: Evolution pathway
        """
        # Collect leads to show
        leads_to_show = []

        if pathway.initial_leads:
            leads_to_show.append(
                ("Initial Lead", pathway.initial_leads[0], self.COLORS['initial'])
            )

        if pathway.intermediate_leads:
            # Show first intermediate
            leads_to_show.append(
                ("Intermediate", pathway.intermediate_leads[0],
                 self.COLORS['intermediate'])
            )

        if pathway.final_lead:
            leads_to_show.append(
                ("Final Lead", pathway.final_lead, self.COLORS['final'])
            )

        for i, ax in enumerate(axes):
            ax.set_facecolor(self.COLORS['background'])
            ax.axis('off')

            if i < len(leads_to_show):
                title, lead, color = leads_to_show[i]

                # Render molecule image
                mol_img = self._render_molecule(lead.compound.mol, size=(400, 300))

                if mol_img is not None:
                    ax.imshow(mol_img)

                # Add title with colored background
                ax.set_title(
                    f"{title}\n{lead.compound.name}",
                    fontsize=10,
                    fontweight='bold',
                    color=self.COLORS['text'],
                    bbox=dict(
                        facecolor=color,
                        alpha=0.3,
                        edgecolor=color,
                        boxstyle='round,pad=0.3'
                    )
                )

                # Add metrics
                metrics_text = (
                    f"Analogs: {lead.analog_count}\n"
                    f"Generation: {lead.generation}\n"
                    f"Cluster size: {len(lead.cluster.compounds)}"
                )
                ax.text(
                    0.02, 0.02, metrics_text,
                    transform=ax.transAxes,
                    fontsize=8,
                    verticalalignment='bottom',
                    color=self.COLORS['text'],
                    bbox=dict(
                        facecolor=self.COLORS['background'],
                        alpha=0.8,
                        edgecolor=self.COLORS['edge'],
                        boxstyle='round,pad=0.2'
                    )
                )

    def _render_molecule(self,
                        mol: Chem.Mol,
                        size: Tuple[int, int] = (300, 200),
                        highlight_atoms: Optional[List[int]] = None) -> Optional[np.ndarray]:
        """
        Render a molecule as an image.

        Args:
            mol: RDKit molecule
            size: Image dimensions (width, height)
            highlight_atoms: Optional atoms to highlight

        Returns:
            Numpy array of the image or None if rendering fails
        """
        try:
            # Generate 2D coordinates if needed
            if mol.GetNumConformers() == 0:
                AllChem.Compute2DCoords(mol)

            # Create drawer
            drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])

            # Configure drawing options
            opts = drawer.drawOptions()
            opts.addStereoAnnotation = True
            opts.addAtomIndices = False
            opts.legendFontSize = 12

            if self.style == 'dark':
                opts.setBackgroundColour((0.1, 0.1, 0.18, 1))

            # Draw molecule
            if highlight_atoms:
                drawer.DrawMolecule(
                    mol,
                    highlightAtoms=highlight_atoms
                )
            else:
                drawer.DrawMolecule(mol)

            drawer.FinishDrawing()

            # Convert to image
            img_data = drawer.GetDrawingText()
            img = Image.open(io.BytesIO(img_data))

            return np.array(img)

        except Exception as e:
            print(f"Warning: Failed to render molecule: {e}")
            return None

    def _add_legend(self, fig: plt.Figure) -> None:
        """
        Add legend to the figure.

        Args:
            fig: Matplotlib figure
        """
        legend_elements = [
            mpatches.Patch(
                facecolor=self.COLORS['initial'],
                edgecolor=self.COLORS['text'],
                label='Initial Lead',
                linewidth=2
            ),
            mpatches.Patch(
                facecolor=self.COLORS['intermediate'],
                edgecolor=self.COLORS['text'],
                label='Intermediate Lead',
                linewidth=2
            ),
            mpatches.Patch(
                facecolor=self.COLORS['final'],
                edgecolor=self.COLORS['text'],
                label='Final Lead',
                linewidth=2
            ),
            mpatches.Patch(
                facecolor=self.COLORS['analog'],
                edgecolor=self.COLORS['analog'],
                label='Analog',
                linewidth=1
            ),
            plt.Line2D(
                [0], [0],
                color=self.COLORS['pathway'],
                linewidth=3,
                label='Evolution Path'
            )
        ]

        fig.legend(
            handles=legend_elements,
            loc='lower center',
            ncol=5,
            fontsize=9,
            frameon=True,
            fancybox=True,
            shadow=True,
            facecolor=self.COLORS['background'],
            edgecolor=self.COLORS['edge']
        )

    def visualize_pathway_linear(self,
                                pathway: EvolutionPathway,
                                output_path: Optional[str] = None,
                                title: str = "Lead Optimization Pathway") -> plt.Figure:
        """
        Create a linear pathway diagram showing progression.

        Shows the pathway as a horizontal flow from initial to final lead,
        with chemical structures and connecting arrows.

        Args:
            pathway: Evolution pathway to visualize
            output_path: Optional path to save the figure
            title: Figure title

        Returns:
            matplotlib Figure object
        """
        all_leads = pathway.all_leads
        n_leads = len(all_leads)

        if n_leads == 0:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, "No pathway to display",
                   ha='center', va='center')
            return fig

        # Calculate figure size based on number of leads
        fig_width = max(12, n_leads * 4)
        fig, ax = plt.subplots(figsize=(fig_width, 8), dpi=self.dpi)
        fig.patch.set_facecolor(self.COLORS['background'])
        ax.set_facecolor(self.COLORS['background'])

        # Sort leads by generation
        sorted_leads = sorted(all_leads, key=lambda l: l.generation)

        # Calculate positions
        x_positions = np.linspace(0.1, 0.9, n_leads)
        y_position = 0.5

        # Draw each lead
        for i, (x, lead) in enumerate(zip(x_positions, sorted_leads)):
            # Determine color
            if lead.lead_type == LeadType.INITIAL:
                color = self.COLORS['initial']
            elif lead.lead_type == LeadType.FINAL:
                color = self.COLORS['final']
            else:
                color = self.COLORS['intermediate']

            # Render molecule
            mol_img = self._render_molecule(lead.compound.mol, size=(250, 200))

            if mol_img is not None:
                # Create image annotation
                imagebox = OffsetImage(mol_img, zoom=0.6)
                ab = AnnotationBbox(
                    imagebox,
                    (x, y_position),
                    frameon=True,
                    bboxprops=dict(
                        facecolor=self.COLORS['background'],
                        edgecolor=color,
                        linewidth=3,
                        boxstyle='round,pad=0.1'
                    )
                )
                ax.add_artist(ab)

            # Add label
            label = f"{lead.compound.name}\n({lead.lead_type.value})"
            ax.text(
                x, y_position - 0.35,
                label,
                ha='center',
                va='top',
                fontsize=10,
                fontweight='bold',
                color=self.COLORS['text'],
                bbox=dict(
                    facecolor=color,
                    alpha=0.3,
                    edgecolor=color,
                    boxstyle='round,pad=0.3'
                )
            )

            # Add generation indicator
            ax.text(
                x, y_position + 0.35,
                f"Gen {lead.generation}",
                ha='center',
                va='bottom',
                fontsize=9,
                color=self.COLORS['text']
            )

            # Draw arrow to next lead
            if i < n_leads - 1:
                arrow_start = x + 0.08
                arrow_end = x_positions[i + 1] - 0.08

                ax.annotate(
                    '',
                    xy=(arrow_end, y_position),
                    xytext=(arrow_start, y_position),
                    arrowprops=dict(
                        arrowstyle='-|>',
                        color=self.COLORS['pathway'],
                        lw=3,
                        mutation_scale=20
                    )
                )

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        ax.set_title(
            title,
            fontsize=14,
            fontweight='bold',
            color=self.COLORS['text'],
            pad=20
        )

        plt.tight_layout()

        if output_path:
            fig.savefig(
                output_path,
                dpi=self.dpi,
                bbox_inches='tight',
                facecolor=self.COLORS['background']
            )

        return fig

    def create_summary_report(self,
                             pathway: EvolutionPathway,
                             compounds: List[Compound],
                             clusters: List[Cluster],
                             output_path: str) -> None:
        """
        Generate a comprehensive HTML summary report.

        Args:
            pathway: Evolution pathway
            compounds: All compounds
            clusters: Compound clusters
            output_path: Path for the HTML report
        """
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Compound Evolution Analysis Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Arial, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background: {self.COLORS['background']};
            color: {self.COLORS['text']};
        }}
        h1, h2, h3 {{
            color: {self.COLORS['text']};
        }}
        .lead-card {{
            border: 2px solid;
            border-radius: 10px;
            padding: 15px;
            margin: 10px 0;
            background: white;
        }}
        .initial {{ border-color: {self.COLORS['initial']}; }}
        .intermediate {{ border-color: {self.COLORS['intermediate']}; }}
        .final {{ border-color: {self.COLORS['final']}; }}
        .stats {{
            display: grid;
            grid-template-columns: repeat(4, 1fr);
            gap: 15px;
            margin: 20px 0;
        }}
        .stat-box {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
        }}
        .stat-value {{
            font-size: 24px;
            font-weight: bold;
            color: {self.COLORS['intermediate']};
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            padding: 10px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background: {self.COLORS['intermediate']};
            color: white;
        }}
        .pathway-arrow {{
            text-align: center;
            font-size: 24px;
            color: {self.COLORS['pathway']};
            margin: 10px 0;
        }}
    </style>
</head>
<body>
    <h1>Compound Evolution Analysis Report</h1>

    <div class="stats">
        <div class="stat-box">
            <div class="stat-value">{len(compounds)}</div>
            <div>Total Compounds</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{len(clusters)}</div>
            <div>Clusters</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{len(pathway.all_leads)}</div>
            <div>Lead Compounds</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{max(l.generation for l in pathway.all_leads) if pathway.all_leads else 0}</div>
            <div>Generations</div>
        </div>
    </div>

    <h2>Evolution Pathway</h2>
"""

        # Add leads in order
        for lead in sorted(pathway.all_leads, key=lambda l: l.generation):
            lead_class = lead.lead_type.value
            html_content += f"""
    <div class="lead-card {lead_class}">
        <h3>{lead.lead_type.value.title()} Lead - {lead.compound.name}</h3>
        <p><strong>SMILES:</strong> {lead.compound.smiles}</p>
        <p><strong>Generation:</strong> {lead.generation}</p>
        <p><strong>Cluster Size:</strong> {len(lead.cluster.compounds)}</p>
        <p><strong>Close Analogs:</strong> {lead.analog_count}</p>
        <p><strong>Centrality Score:</strong> {lead.centrality_score:.3f}</p>
    </div>
    <div class="pathway-arrow">▼</div>
"""

        html_content += """
    <h2>Cluster Summary</h2>
    <table>
        <tr>
            <th>Cluster ID</th>
            <th>Size</th>
            <th>Centroid</th>
            <th>Avg Similarity</th>
            <th>Density</th>
        </tr>
"""

        for cluster in clusters:
            html_content += f"""
        <tr>
            <td>{cluster.id}</td>
            <td>{len(cluster.compounds)}</td>
            <td>{cluster.centroid.name if cluster.centroid else 'N/A'}</td>
            <td>{cluster.avg_internal_similarity:.3f}</td>
            <td>{cluster.density:.3f}</td>
        </tr>
"""

        html_content += """
    </table>
</body>
</html>
"""

        with open(output_path, 'w') as f:
            f.write(html_content)
