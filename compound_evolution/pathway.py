"""
Pathway analysis module for identifying lead compound evolution.

This module analyzes compound clusters to identify:
- Initial lead compounds (starting points)
- Intermediate leads (optimization steps)
- Final lead compounds (optimized results)
- The progression pathway connecting them
"""

from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass, field
from enum import Enum
import numpy as np

import networkx as nx

from .fingerprints import Compound, FingerprintCalculator
from .clustering import Cluster, CompoundClusterer


class LeadType(Enum):
    """Classification of lead compounds in the optimization pathway."""
    INITIAL = "initial"       # Starting point - first/parent scaffold
    INTERMEDIATE = "intermediate"  # Optimization step
    FINAL = "final"           # Optimized lead compound


@dataclass
class LeadCompound:
    """
    A lead compound in the optimization pathway.

    Attributes:
        compound: The Compound object
        lead_type: Classification (initial, intermediate, final)
        cluster: The cluster this lead represents
        analog_count: Number of close analogs around this lead
        centrality_score: How central this compound is in similarity network
        generation: Estimated generation in the optimization (0 = initial)
    """
    compound: Compound
    lead_type: LeadType
    cluster: Cluster
    analog_count: int = 0
    centrality_score: float = 0.0
    generation: int = 0

    def __hash__(self):
        return hash(self.compound.smiles)


@dataclass
class EvolutionPathway:
    """
    The complete optimization pathway from initial to final lead.

    Attributes:
        initial_leads: Starting lead compound(s)
        intermediate_leads: Optimization steps
        final_lead: The optimized lead compound
        pathway_graph: NetworkX graph showing compound relationships
        generation_map: Maps compounds to their generation number
    """
    initial_leads: List[LeadCompound] = field(default_factory=list)
    intermediate_leads: List[LeadCompound] = field(default_factory=list)
    final_lead: Optional[LeadCompound] = None
    pathway_graph: Optional[nx.DiGraph] = None
    generation_map: Dict[str, int] = field(default_factory=dict)

    @property
    def all_leads(self) -> List[LeadCompound]:
        """Get all leads in order: initial -> intermediate -> final."""
        leads = self.initial_leads.copy()
        leads.extend(self.intermediate_leads)
        if self.final_lead:
            leads.append(self.final_lead)
        return leads


class PathwayAnalyzer:
    """
    Analyzes compound clusters to identify the optimization pathway.

    The analyzer uses multiple heuristics to identify lead compounds:
    1. Cluster size - leads have many analogs around them
    2. Centrality - leads are structurally central
    3. Network position - initial leads connect to many clusters
    4. Density - final leads have tight clusters of optimized variants

    Attributes:
        similarity_threshold: Threshold for "close analog" (default 0.7)
        fp_calculator: Fingerprint calculator for similarity
    """

    def __init__(self, similarity_threshold: float = 0.7):
        """
        Initialize the pathway analyzer.

        Args:
            similarity_threshold: Similarity cutoff for close analogs
        """
        self.similarity_threshold = similarity_threshold
        self.fp_calculator = FingerprintCalculator()

    def analyze(self,
               compounds: List[Compound],
               clusters: List[Cluster],
               similarity_matrix: np.ndarray) -> EvolutionPathway:
        """
        Analyze compounds to identify the optimization pathway.

        Args:
            compounds: All compounds from the patent
            clusters: Clustered compounds
            similarity_matrix: Pairwise similarity matrix

        Returns:
            EvolutionPathway describing the lead compound evolution
        """
        # Build compound similarity network
        compound_graph = self._build_compound_graph(
            compounds, similarity_matrix
        )

        # Calculate centrality scores for all compounds
        centrality_scores = self._calculate_centrality(
            compound_graph, compounds
        )

        # Identify lead compounds from each cluster
        leads = self._identify_leads(
            clusters, centrality_scores, similarity_matrix, compounds
        )

        # Determine the pathway progression
        pathway = self._determine_pathway(
            leads, compounds, similarity_matrix, compound_graph
        )

        return pathway

    def _build_compound_graph(self,
                             compounds: List[Compound],
                             similarity_matrix: np.ndarray) -> nx.Graph:
        """
        Build a similarity network connecting compounds.

        Args:
            compounds: List of compounds
            similarity_matrix: Pairwise similarities

        Returns:
            NetworkX graph with similarity-weighted edges
        """
        G = nx.Graph()

        # Add all compounds as nodes
        for i, compound in enumerate(compounds):
            G.add_node(
                i,
                compound=compound,
                smiles=compound.smiles,
                name=compound.name
            )

        # Add edges for similar compounds
        n = len(compounds)
        for i in range(n):
            for j in range(i + 1, n):
                sim = similarity_matrix[i, j]
                if sim >= self.similarity_threshold:
                    G.add_edge(i, j, weight=sim, similarity=sim)

        return G

    def _calculate_centrality(self,
                             G: nx.Graph,
                             compounds: List[Compound]) -> Dict[str, float]:
        """
        Calculate multiple centrality measures for compounds.

        Uses a combination of:
        - Degree centrality (number of connections)
        - Betweenness centrality (bridge position)
        - Eigenvector centrality (connection to important nodes)

        Args:
            G: Compound similarity graph
            compounds: List of compounds

        Returns:
            Dictionary mapping SMILES to combined centrality score
        """
        if len(G.nodes) == 0:
            return {}

        # Calculate various centrality measures
        try:
            degree_cent = nx.degree_centrality(G)
        except Exception:
            degree_cent = {i: 0 for i in range(len(compounds))}

        try:
            betweenness_cent = nx.betweenness_centrality(G, weight='weight')
        except Exception:
            betweenness_cent = {i: 0 for i in range(len(compounds))}

        try:
            eigenvector_cent = nx.eigenvector_centrality_numpy(
                G, weight='weight'
            )
        except Exception:
            eigenvector_cent = {i: 0 for i in range(len(compounds))}

        # Normalize and combine
        def normalize(d):
            values = list(d.values())
            if not values:
                return d
            min_v, max_v = min(values), max(values)
            if max_v == min_v:
                return {k: 0.5 for k in d}
            return {k: (v - min_v) / (max_v - min_v) for k, v in d.items()}

        degree_norm = normalize(degree_cent)
        between_norm = normalize(betweenness_cent)
        eigen_norm = normalize(eigenvector_cent)

        # Combined score with weights
        centrality = {}
        for i, compound in enumerate(compounds):
            score = (
                0.4 * degree_norm.get(i, 0) +
                0.3 * between_norm.get(i, 0) +
                0.3 * eigen_norm.get(i, 0)
            )
            centrality[compound.smiles] = score

        return centrality

    def _identify_leads(self,
                       clusters: List[Cluster],
                       centrality_scores: Dict[str, float],
                       similarity_matrix: np.ndarray,
                       all_compounds: List[Compound]) -> List[LeadCompound]:
        """
        Identify potential lead compounds from clusters.

        Each cluster's centroid is a potential lead. The leads are
        scored based on cluster properties and centrality.

        Args:
            clusters: List of compound clusters
            centrality_scores: Centrality score for each compound
            similarity_matrix: Pairwise similarities
            all_compounds: All compounds

        Returns:
            List of LeadCompound objects
        """
        leads = []
        compound_to_idx = {c.smiles: i for i, c in enumerate(all_compounds)}

        for cluster in clusters:
            if cluster.centroid is None:
                continue

            # Count close analogs
            centroid_idx = compound_to_idx.get(cluster.centroid.smiles)
            if centroid_idx is None:
                continue

            analog_count = 0
            for compound in all_compounds:
                idx = compound_to_idx.get(compound.smiles)
                if idx is not None and idx != centroid_idx:
                    if similarity_matrix[centroid_idx, idx] >= self.similarity_threshold:
                        analog_count += 1

            lead = LeadCompound(
                compound=cluster.centroid,
                lead_type=LeadType.INTERMEDIATE,  # Will be refined later
                cluster=cluster,
                analog_count=analog_count,
                centrality_score=centrality_scores.get(
                    cluster.centroid.smiles, 0
                )
            )
            leads.append(lead)

        return leads

    def _determine_pathway(self,
                          leads: List[LeadCompound],
                          compounds: List[Compound],
                          similarity_matrix: np.ndarray,
                          compound_graph: nx.Graph) -> EvolutionPathway:
        """
        Determine the optimization pathway from leads.

        Uses heuristics to identify:
        - Initial lead: High centrality, connects to many clusters
        - Final lead: Dense cluster of recent optimizations
        - Intermediates: Connected leads between initial and final

        Args:
            leads: Potential lead compounds
            compounds: All compounds
            similarity_matrix: Pairwise similarities
            compound_graph: Compound similarity network

        Returns:
            EvolutionPathway with classified leads
        """
        pathway = EvolutionPathway()

        if not leads:
            return pathway

        compound_to_idx = {c.smiles: i for i, c in enumerate(compounds)}

        # Score each lead for being initial vs final
        # Initial leads: high connectivity to diverse structures
        # Final leads: tight cluster of very similar compounds
        initial_scores = []
        final_scores = []

        for lead in leads:
            # Initial lead characteristics:
            # - High centrality (hub position)
            # - Moderate cluster size (spawned many directions)
            # - Lower cluster density (more diverse analogs)
            initial_score = (
                0.5 * lead.centrality_score +
                0.3 * min(lead.analog_count / 20, 1.0) +
                0.2 * (1 - lead.cluster.density)  # Lower density = more diverse
            )
            initial_scores.append(initial_score)

            # Final lead characteristics:
            # - High cluster density (tight optimization)
            # - Many close analogs (intensive exploration)
            # - Compounds might be at "edge" of network (endpoint)
            final_score = (
                0.4 * lead.cluster.density +
                0.3 * lead.cluster.avg_internal_similarity +
                0.3 * min(lead.analog_count / 10, 1.0)
            )
            final_scores.append(final_score)

        # Identify initial lead(s)
        initial_idx = np.argmax(initial_scores)
        leads[initial_idx].lead_type = LeadType.INITIAL
        leads[initial_idx].generation = 0
        pathway.initial_leads.append(leads[initial_idx])

        # Identify final lead
        # Exclude initial lead from consideration
        final_scores_adj = final_scores.copy()
        final_scores_adj[initial_idx] = -1
        final_idx = np.argmax(final_scores_adj)

        if final_idx != initial_idx and len(leads) > 1:
            leads[final_idx].lead_type = LeadType.FINAL
            pathway.final_lead = leads[final_idx]

        # Remaining leads are intermediates
        for i, lead in enumerate(leads):
            if i != initial_idx and i != final_idx:
                lead.lead_type = LeadType.INTERMEDIATE
                pathway.intermediate_leads.append(lead)

        # Build pathway graph and assign generations
        pathway.pathway_graph = self._build_pathway_graph(
            pathway, compounds, similarity_matrix
        )

        # Determine generations using BFS from initial lead
        self._assign_generations(pathway, similarity_matrix, compound_to_idx)

        return pathway

    def _build_pathway_graph(self,
                            pathway: EvolutionPathway,
                            compounds: List[Compound],
                            similarity_matrix: np.ndarray) -> nx.DiGraph:
        """
        Build a directed graph showing the progression pathway.

        Args:
            pathway: The evolution pathway
            compounds: All compounds
            similarity_matrix: Pairwise similarities

        Returns:
            Directed graph with pathway edges
        """
        G = nx.DiGraph()
        compound_to_idx = {c.smiles: i for i, c in enumerate(compounds)}

        all_leads = pathway.all_leads
        if not all_leads:
            return G

        # Add lead nodes
        for lead in all_leads:
            G.add_node(
                lead.compound.smiles,
                lead=lead,
                lead_type=lead.lead_type.value,
                cluster_size=len(lead.cluster.compounds)
            )

        # Connect leads based on similarity (potential evolution path)
        for i, lead1 in enumerate(all_leads):
            idx1 = compound_to_idx.get(lead1.compound.smiles)
            if idx1 is None:
                continue

            for lead2 in all_leads[i + 1:]:
                idx2 = compound_to_idx.get(lead2.compound.smiles)
                if idx2 is None:
                    continue

                sim = similarity_matrix[idx1, idx2]
                if sim >= self.similarity_threshold * 0.7:  # Lower threshold for leads
                    # Determine direction based on lead types
                    if lead1.lead_type == LeadType.INITIAL:
                        G.add_edge(
                            lead1.compound.smiles,
                            lead2.compound.smiles,
                            weight=sim
                        )
                    elif lead2.lead_type == LeadType.INITIAL:
                        G.add_edge(
                            lead2.compound.smiles,
                            lead1.compound.smiles,
                            weight=sim
                        )
                    elif lead1.lead_type == LeadType.FINAL:
                        G.add_edge(
                            lead2.compound.smiles,
                            lead1.compound.smiles,
                            weight=sim
                        )
                    elif lead2.lead_type == LeadType.FINAL:
                        G.add_edge(
                            lead1.compound.smiles,
                            lead2.compound.smiles,
                            weight=sim
                        )
                    else:
                        # Both intermediate - use higher centrality as source
                        if lead1.centrality_score >= lead2.centrality_score:
                            G.add_edge(
                                lead1.compound.smiles,
                                lead2.compound.smiles,
                                weight=sim
                            )
                        else:
                            G.add_edge(
                                lead2.compound.smiles,
                                lead1.compound.smiles,
                                weight=sim
                            )

        return G

    def _assign_generations(self,
                           pathway: EvolutionPathway,
                           similarity_matrix: np.ndarray,
                           compound_to_idx: Dict[str, int]) -> None:
        """
        Assign generation numbers to leads using BFS from initial.

        Args:
            pathway: Evolution pathway to update
            similarity_matrix: Pairwise similarities
            compound_to_idx: Mapping from SMILES to index
        """
        if not pathway.initial_leads:
            return

        G = pathway.pathway_graph
        if G is None or len(G.nodes) == 0:
            return

        # BFS from initial lead
        visited = set()
        queue = [(pathway.initial_leads[0].compound.smiles, 0)]
        pathway.generation_map[pathway.initial_leads[0].compound.smiles] = 0

        while queue:
            current_smiles, generation = queue.pop(0)
            if current_smiles in visited:
                continue
            visited.add(current_smiles)

            # Find the lead for this compound
            for lead in pathway.all_leads:
                if lead.compound.smiles == current_smiles:
                    lead.generation = generation
                    pathway.generation_map[current_smiles] = generation
                    break

            # Add neighbors
            if current_smiles in G:
                for neighbor in G.neighbors(current_smiles):
                    if neighbor not in visited:
                        queue.append((neighbor, generation + 1))

        # Ensure final lead has highest generation
        if pathway.final_lead:
            max_gen = max(
                lead.generation for lead in pathway.all_leads
            )
            pathway.final_lead.generation = max_gen

    def get_compound_generation(self,
                               compound: Compound,
                               pathway: EvolutionPathway,
                               similarity_matrix: np.ndarray,
                               all_compounds: List[Compound]) -> int:
        """
        Determine which generation a compound belongs to.

        A compound belongs to the generation of its most similar lead.

        Args:
            compound: Compound to classify
            pathway: Evolution pathway
            similarity_matrix: Pairwise similarities
            all_compounds: All compounds

        Returns:
            Generation number (0 = initial, higher = later optimization)
        """
        compound_to_idx = {c.smiles: i for i, c in enumerate(all_compounds)}
        compound_idx = compound_to_idx.get(compound.smiles)

        if compound_idx is None:
            return -1

        # Find most similar lead
        best_sim = -1
        best_generation = 0

        for lead in pathway.all_leads:
            lead_idx = compound_to_idx.get(lead.compound.smiles)
            if lead_idx is not None:
                sim = similarity_matrix[compound_idx, lead_idx]
                if sim > best_sim:
                    best_sim = sim
                    best_generation = lead.generation

        return best_generation
