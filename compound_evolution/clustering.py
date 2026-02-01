"""
Clustering module for identifying structural groups of compounds.

This module implements various clustering algorithms to group similar
compounds together, which helps identify lead compound series.
"""

from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass, field
import numpy as np

from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
import networkx as nx

from .fingerprints import Compound, FingerprintCalculator


@dataclass
class Cluster:
    """
    Represents a cluster of structurally similar compounds.

    Attributes:
        id: Unique cluster identifier
        compounds: List of compounds in this cluster
        centroid: The compound most representative of this cluster
        centroid_index: Index of centroid in compounds list
        avg_internal_similarity: Average pairwise similarity within cluster
        density: How tightly grouped the compounds are
    """
    id: int
    compounds: List[Compound] = field(default_factory=list)
    centroid: Optional[Compound] = None
    centroid_index: int = 0
    avg_internal_similarity: float = 0.0
    density: float = 0.0

    def __len__(self):
        return len(self.compounds)

    @property
    def size(self) -> int:
        return len(self.compounds)


class CompoundClusterer:
    """
    Clusters compounds based on structural similarity.

    Supports multiple clustering approaches:
    - Hierarchical clustering with distance threshold
    - Network-based clustering using community detection
    - Butina clustering (sphere exclusion)

    Attributes:
        similarity_threshold: Minimum similarity for compounds to be
                            considered "related" (default 0.7)
        fp_calculator: FingerprintCalculator instance for similarity
    """

    def __init__(self, similarity_threshold: float = 0.7):
        """
        Initialize the clusterer.

        Args:
            similarity_threshold: Tanimoto similarity cutoff (0.0-1.0)
                                Higher values create smaller, tighter clusters
        """
        self.similarity_threshold = similarity_threshold
        self.fp_calculator = FingerprintCalculator()

    def cluster_hierarchical(self,
                            compounds: List[Compound],
                            similarity_matrix: np.ndarray,
                            method: str = 'average') -> List[Cluster]:
        """
        Perform hierarchical clustering on compounds.

        Uses agglomerative clustering with the specified linkage method.
        Distance is defined as (1 - similarity).

        Args:
            compounds: List of compounds to cluster
            similarity_matrix: Precomputed similarity matrix
            method: Linkage method ('single', 'complete', 'average', 'ward')

        Returns:
            List of Cluster objects
        """
        # Convert similarity to distance
        distance_matrix = 1.0 - similarity_matrix

        # Ensure diagonal is exactly 0
        np.fill_diagonal(distance_matrix, 0)

        # Handle any numerical issues
        distance_matrix = np.clip(distance_matrix, 0, 1)

        # Convert to condensed form for scipy
        condensed = squareform(distance_matrix, checks=False)

        # Perform hierarchical clustering
        linkage_matrix = linkage(condensed, method=method)

        # Cut tree at threshold distance
        distance_threshold = 1.0 - self.similarity_threshold
        cluster_labels = fcluster(
            linkage_matrix,
            t=distance_threshold,
            criterion='distance'
        )

        return self._build_clusters(
            compounds, cluster_labels, similarity_matrix
        )

    def cluster_network(self,
                       compounds: List[Compound],
                       similarity_matrix: np.ndarray) -> List[Cluster]:
        """
        Perform network-based clustering using community detection.

        Creates a similarity network where edges connect compounds
        with similarity above threshold, then uses the Louvain
        algorithm to detect communities.

        Args:
            compounds: List of compounds to cluster
            similarity_matrix: Precomputed similarity matrix

        Returns:
            List of Cluster objects
        """
        # Build similarity graph
        G = nx.Graph()

        # Add nodes
        for i, compound in enumerate(compounds):
            G.add_node(i, compound=compound)

        # Add edges for similar compounds
        n = len(compounds)
        for i in range(n):
            for j in range(i + 1, n):
                sim = similarity_matrix[i, j]
                if sim >= self.similarity_threshold:
                    G.add_edge(i, j, weight=sim)

        # Detect communities using Louvain algorithm
        try:
            communities = nx.community.louvain_communities(
                G, weight='weight', resolution=1.0
            )
        except Exception:
            # Fallback to greedy modularity if Louvain fails
            communities = list(nx.community.greedy_modularity_communities(G))

        # Convert communities to cluster labels
        cluster_labels = np.zeros(n, dtype=int)
        for cluster_id, community in enumerate(communities):
            for node in community:
                cluster_labels[node] = cluster_id + 1

        # Handle isolated nodes (no edges above threshold)
        max_label = max(cluster_labels) if len(cluster_labels) > 0 else 0
        for i in range(n):
            if cluster_labels[i] == 0:
                max_label += 1
                cluster_labels[i] = max_label

        return self._build_clusters(
            compounds, cluster_labels, similarity_matrix
        )

    def cluster_butina(self,
                      compounds: List[Compound],
                      similarity_matrix: np.ndarray) -> List[Cluster]:
        """
        Perform Butina clustering (sphere exclusion).

        The Butina algorithm:
        1. Calculate number of neighbors for each compound
        2. Pick compound with most neighbors as cluster center
        3. Remove center and its neighbors from pool
        4. Repeat until all compounds assigned

        This is particularly good for identifying lead series
        as cluster centers are naturally the most "central" compounds.

        Args:
            compounds: List of compounds to cluster
            similarity_matrix: Precomputed similarity matrix

        Returns:
            List of Cluster objects
        """
        n = len(compounds)
        cluster_labels = np.zeros(n, dtype=int)
        assigned = set()
        current_cluster = 0

        # Calculate neighbor counts
        def get_neighbors(idx: int) -> Set[int]:
            neighbors = set()
            for j in range(n):
                if j != idx and j not in assigned:
                    if similarity_matrix[idx, j] >= self.similarity_threshold:
                        neighbors.add(j)
            return neighbors

        while len(assigned) < n:
            # Find unassigned compound with most neighbors
            best_idx = -1
            best_neighbor_count = -1

            for i in range(n):
                if i not in assigned:
                    neighbors = get_neighbors(i)
                    if len(neighbors) > best_neighbor_count:
                        best_neighbor_count = len(neighbors)
                        best_idx = i

            if best_idx == -1:
                break

            # Create new cluster
            current_cluster += 1
            neighbors = get_neighbors(best_idx)

            # Assign center and neighbors to cluster
            cluster_labels[best_idx] = current_cluster
            assigned.add(best_idx)

            for neighbor in neighbors:
                cluster_labels[neighbor] = current_cluster
                assigned.add(neighbor)

        # Handle any remaining unassigned (shouldn't happen but safety check)
        for i in range(n):
            if cluster_labels[i] == 0:
                current_cluster += 1
                cluster_labels[i] = current_cluster

        return self._build_clusters(
            compounds, cluster_labels, similarity_matrix
        )

    def _build_clusters(self,
                       compounds: List[Compound],
                       labels: np.ndarray,
                       similarity_matrix: np.ndarray) -> List[Cluster]:
        """
        Build Cluster objects from cluster labels.

        Args:
            compounds: List of all compounds
            labels: Cluster label for each compound
            similarity_matrix: Pairwise similarity matrix

        Returns:
            List of Cluster objects with centroids identified
        """
        unique_labels = np.unique(labels)
        clusters = []

        for cluster_id in unique_labels:
            # Get compounds in this cluster
            indices = np.where(labels == cluster_id)[0]
            cluster_compounds = [compounds[i] for i in indices]

            if len(cluster_compounds) == 0:
                continue

            cluster = Cluster(
                id=int(cluster_id),
                compounds=cluster_compounds
            )

            # Find centroid (compound most similar to all others in cluster)
            if len(indices) == 1:
                cluster.centroid = cluster_compounds[0]
                cluster.centroid_index = 0
                cluster.avg_internal_similarity = 1.0
                cluster.density = 1.0
            else:
                # Calculate average similarity to cluster for each member
                avg_sims = []
                for i, idx_i in enumerate(indices):
                    sims = []
                    for idx_j in indices:
                        if idx_i != idx_j:
                            sims.append(similarity_matrix[idx_i, idx_j])
                    avg_sims.append(np.mean(sims) if sims else 0)

                centroid_local_idx = np.argmax(avg_sims)
                cluster.centroid = cluster_compounds[centroid_local_idx]
                cluster.centroid_index = centroid_local_idx

                # Calculate cluster metrics
                all_sims = []
                for i, idx_i in enumerate(indices):
                    for j, idx_j in enumerate(indices):
                        if i < j:
                            all_sims.append(similarity_matrix[idx_i, idx_j])

                cluster.avg_internal_similarity = np.mean(all_sims) if all_sims else 1.0
                cluster.density = max(avg_sims) if avg_sims else 1.0

            clusters.append(cluster)

        # Sort clusters by size (largest first)
        clusters.sort(key=lambda c: len(c.compounds), reverse=True)

        # Reassign IDs based on new order
        for new_id, cluster in enumerate(clusters):
            cluster.id = new_id

        return clusters

    def get_inter_cluster_similarity(self,
                                     cluster1: Cluster,
                                     cluster2: Cluster,
                                     similarity_matrix: np.ndarray,
                                     all_compounds: List[Compound]) -> float:
        """
        Calculate similarity between two clusters.

        Uses the maximum similarity between cluster centroids or
        the average of top-k similarities between cluster members.

        Args:
            cluster1: First cluster
            cluster2: Second cluster
            similarity_matrix: Full pairwise similarity matrix
            all_compounds: Complete list of compounds (for index mapping)

        Returns:
            Inter-cluster similarity score
        """
        # Build index mapping
        compound_to_idx = {c.smiles: i for i, c in enumerate(all_compounds)}

        # Get indices for each cluster
        indices1 = [compound_to_idx[c.smiles] for c in cluster1.compounds]
        indices2 = [compound_to_idx[c.smiles] for c in cluster2.compounds]

        # Calculate all pairwise similarities between clusters
        similarities = []
        for i in indices1:
            for j in indices2:
                similarities.append(similarity_matrix[i, j])

        if not similarities:
            return 0.0

        # Return max similarity (closest pair)
        return max(similarities)

    def build_cluster_graph(self,
                           clusters: List[Cluster],
                           similarity_matrix: np.ndarray,
                           all_compounds: List[Compound],
                           edge_threshold: float = 0.5) -> nx.Graph:
        """
        Build a graph connecting related clusters.

        Args:
            clusters: List of clusters
            similarity_matrix: Pairwise similarity matrix
            all_compounds: All compounds
            edge_threshold: Minimum inter-cluster similarity for edge

        Returns:
            NetworkX graph with clusters as nodes
        """
        G = nx.Graph()

        # Add cluster nodes
        for cluster in clusters:
            G.add_node(
                cluster.id,
                size=len(cluster.compounds),
                centroid=cluster.centroid,
                density=cluster.density
            )

        # Add edges between related clusters
        for i, c1 in enumerate(clusters):
            for c2 in clusters[i + 1:]:
                sim = self.get_inter_cluster_similarity(
                    c1, c2, similarity_matrix, all_compounds
                )
                if sim >= edge_threshold:
                    G.add_edge(c1.id, c2.id, weight=sim)

        return G
