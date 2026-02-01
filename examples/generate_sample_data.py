"""
Generate sample SDF data for demonstrating compound evolution analysis.

This script creates a realistic set of compounds that simulate a lead
optimization campaign, such as a kinase inhibitor SAR (Structure-Activity
Relationship) study.

The generated data demonstrates:
1. An initial lead compound (hit from screening)
2. First-generation analogs exploring different positions
3. Second-generation compounds with combined optimizations
4. Final optimized lead with refined properties
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from typing import List, Tuple


def generate_kinase_inhibitor_series() -> List[Tuple[str, str]]:
    """
    Generate a kinase inhibitor lead optimization series.

    Returns a list of (SMILES, name) tuples representing a realistic
    SAR campaign starting from a pyrimidine-based kinase inhibitor.
    """
    compounds = []

    # ================================================================
    # GENERATION 0: Initial Hit/Lead
    # ================================================================
    # Simple pyrimidine scaffold - the starting point
    initial_lead = "c1ccc(Nc2ncnc3ccccc23)cc1"  # 4-anilinoquinazoline
    compounds.append((initial_lead, "Initial_Lead_001"))

    # ================================================================
    # GENERATION 1: First-round analogs - exploring the scaffold
    # ================================================================

    # Analogs varying the aniline ring substituents
    gen1_analogs = [
        # Halogen substitutions on aniline
        ("Fc1ccc(Nc2ncnc3ccccc23)cc1", "Gen1_Fluoro_para"),
        ("Clc1ccc(Nc2ncnc3ccccc23)cc1", "Gen1_Chloro_para"),
        ("Fc1ccccc1Nc1ncnc2ccccc12", "Gen1_Fluoro_ortho"),
        ("Fc1cccc(Nc2ncnc3ccccc23)c1", "Gen1_Fluoro_meta"),

        # Methyl substitutions
        ("Cc1ccc(Nc2ncnc3ccccc23)cc1", "Gen1_Methyl_para"),
        ("Cc1ccccc1Nc1ncnc2ccccc12", "Gen1_Methyl_ortho"),

        # Electron-donating groups
        ("COc1ccc(Nc2ncnc3ccccc23)cc1", "Gen1_Methoxy_para"),
        ("Nc1ccc(Nc2ncnc3ccccc23)cc1", "Gen1_Amino_para"),

        # Disubstituted anilines
        ("Fc1ccc(Nc2ncnc3ccccc23)cc1F", "Gen1_DiFluoro"),
        ("Clc1cc(Nc2ncnc3ccccc23)ccc1Cl", "Gen1_DiChloro"),
    ]
    compounds.extend(gen1_analogs)

    # Analogs varying the quinazoline core
    gen1_core_analogs = [
        # Substituents on quinazoline
        ("c1ccc(Nc2ncnc3cc(Cl)ccc23)cc1", "Gen1_Core_Chloro"),
        ("c1ccc(Nc2ncnc3cc(F)ccc23)cc1", "Gen1_Core_Fluoro"),
        ("c1ccc(Nc2ncnc3cc(C)ccc23)cc1", "Gen1_Core_Methyl"),
        ("c1ccc(Nc2ncnc3ccc(OC)cc23)cc1", "Gen1_Core_Methoxy"),

        # Different nitrogen positions
        ("c1ccc(Nc2ncc3ccccc3n2)cc1", "Gen1_Quinoxaline"),
    ]
    compounds.extend(gen1_core_analogs)

    # ================================================================
    # GENERATION 2: Second-round - combining favorable modifications
    # ================================================================

    # Best aniline + core modifications combined
    gen2_combined = [
        # 3-Fluoro aniline emerged as favorable
        ("Fc1cccc(Nc2ncnc3cc(Cl)ccc23)c1", "Gen2_3F_6Cl"),
        ("Fc1cccc(Nc2ncnc3cc(F)ccc23)c1", "Gen2_3F_6F"),
        ("Fc1cccc(Nc2ncnc3ccc(OC)cc23)c1", "Gen2_3F_7OMe"),

        # 4-Chloro aniline combinations
        ("Clc1ccc(Nc2ncnc3cc(Cl)ccc23)cc1", "Gen2_4Cl_6Cl"),
        ("Clc1ccc(Nc2ncnc3cc(F)ccc23)cc1", "Gen2_4Cl_6F"),

        # 3,4-Difluoro combinations (potency driver)
        ("Fc1ccc(Nc2ncnc3cc(Cl)ccc23)cc1F", "Gen2_34F_6Cl"),
        ("Fc1ccc(Nc2ncnc3cc(OC)ccc23)cc1F", "Gen2_34F_6OMe"),
    ]
    compounds.extend(gen2_combined)

    # Adding solubilizing groups (addressing ADMET)
    gen2_solubility = [
        ("Fc1cccc(Nc2ncnc3cc(CN4CCNCC4)ccc23)c1", "Gen2_Piperazine"),
        ("Fc1cccc(Nc2ncnc3cc(CNCCN)ccc23)c1", "Gen2_Ethylenediamine"),
        ("Fc1cccc(Nc2ncnc3cc(C(=O)N)ccc23)c1", "Gen2_Carboxamide"),
    ]
    compounds.extend(gen2_solubility)

    # ================================================================
    # GENERATION 3: Final optimization - refined analogs
    # ================================================================

    # Optimized lead series with best combination
    gen3_optimized = [
        # The emerging optimized scaffold: 3-F-aniline + 6-position modification
        ("Fc1cccc(Nc2ncnc3cc(CN4CCN(C)CC4)ccc23)c1", "Gen3_Lead_NMePiperazine"),
        ("Fc1cccc(Nc2ncnc3cc(CN4CCOCC4)ccc23)c1", "Gen3_Lead_Morpholine"),
        ("Fc1cccc(Nc2ncnc3cc(CNCC(C)C)ccc23)c1", "Gen3_Lead_Isobutylamine"),

        # Fine-tuning the aniline
        ("Fc1cc(F)cc(Nc2ncnc3cc(CN4CCN(C)CC4)ccc23)c1", "Gen3_35F_NMePip"),
        ("Fc1ccc(F)c(Nc2ncnc3cc(CN4CCN(C)CC4)ccc23)c1", "Gen3_23F_NMePip"),

        # Different linkers
        ("Fc1cccc(Nc2ncnc3cc(CCN4CCN(C)CC4)ccc23)c1", "Gen3_EthylLinker"),
        ("Fc1cccc(Nc2ncnc3cc(OCN4CCN(C)CC4)ccc23)c1", "Gen3_OxyLinker"),
    ]
    compounds.extend(gen3_optimized)

    # ================================================================
    # GENERATION 4: Final Lead and close analogs (dense cluster)
    # ================================================================

    # The final optimized lead - extensive analog coverage
    final_lead = "Fc1cccc(Nc2ncnc3cc(CN4CCN(C)CC4)ccc23)c1"  # Same as Gen3_Lead_NMePiperazine
    # Note: Already added above, these are the final close analogs

    gen4_final_analogs = [
        # Very close analogs of the final lead (fine-tuning)
        ("Fc1cccc(Nc2ncnc3cc(CN4CCN(CC)CC4)ccc23)c1", "Final_NEtPiperazine"),
        ("Fc1cccc(Nc2ncnc3cc(CN4CCN(C(C)C)CC4)ccc23)c1", "Final_NiPrPiperazine"),
        ("Fc1cccc(Nc2ncnc3cc(CN4CCN(CC(F)(F)F)CC4)ccc23)c1", "Final_NTFM_Pip"),
        ("Fc1cccc(Nc2ncnc3cc(CN4CCN(Cc5ccccc5)CC4)ccc23)c1", "Final_NBnPiperazine"),

        # Stereoisomers and constrained analogs
        ("Fc1cccc(Nc2ncnc3cc(C[C@H]4CNCCN4C)ccc23)c1", "Final_Chiralpip_R"),
        ("Fc1cccc(Nc2ncnc3cc(C[C@@H]4CNCCN4C)ccc23)c1", "Final_Chiralpip_S"),

        # Final backup compounds
        ("Fc1cccc(Nc2ncnc3cc(CN4CCC(N(C)C)CC4)ccc23)c1", "Final_4NMe2_pip"),
        ("Fc1cccc(Nc2ncnc3cc(CN4CCN(C)C4)ccc23)c1", "Final_NMe_Pyr"),
    ]
    compounds.extend(gen4_final_analogs)

    # ================================================================
    # Add some "dead-end" analogs that were explored but not pursued
    # ================================================================
    dead_ends = [
        # Large lipophilic groups (poor solubility)
        ("c1ccc(Nc2ncnc3cc(Cc4ccccc4)ccc23)cc1", "DeadEnd_Benzyl"),
        ("c1ccc(Nc2ncnc3cc(CCc4ccccc4)ccc23)cc1", "DeadEnd_Phenethyl"),

        # Too polar
        ("O=C(O)c1ccc(Nc2ncnc3ccccc23)cc1", "DeadEnd_Carboxylic"),

        # Different core (explored but abandoned)
        ("c1ccc(Nc2ccc3nccnc3n2)cc1", "DeadEnd_Pyrido"),
    ]
    compounds.extend(dead_ends)

    return compounds


def create_sdf_file(compounds: List[Tuple[str, str]], output_path: str) -> None:
    """
    Create an SDF file from a list of SMILES and names.

    Args:
        compounds: List of (SMILES, name) tuples
        output_path: Path for output SDF file
    """
    writer = Chem.SDWriter(output_path)

    for smiles, name in compounds:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Warning: Could not parse SMILES for {name}: {smiles}")
            continue

        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)

        # Set molecule name
        mol.SetProp("_Name", name)

        # Add some properties
        mol.SetProp("SMILES", smiles)
        mol.SetProp("MolWeight", f"{Descriptors.MolWt(mol):.2f}")
        mol.SetProp("LogP", f"{Descriptors.MolLogP(mol):.2f}")
        mol.SetProp("HBD", str(Descriptors.NumHDonors(mol)))
        mol.SetProp("HBA", str(Descriptors.NumHAcceptors(mol)))
        mol.SetProp("TPSA", f"{Descriptors.TPSA(mol):.2f}")

        writer.write(mol)

    writer.close()
    print(f"Created SDF file with {len(compounds)} compounds: {output_path}")


def main():
    """Generate the sample SDF file."""
    import os

    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, "sample_kinase_inhibitors.sdf")

    print("Generating sample kinase inhibitor lead optimization series...")
    compounds = generate_kinase_inhibitor_series()
    print(f"Generated {len(compounds)} compounds")

    create_sdf_file(compounds, output_path)
    print(f"\nSample data saved to: {output_path}")
    print("\nYou can analyze this file using:")
    print(f"  compound-evolution analyze {output_path} -o results/")


if __name__ == "__main__":
    main()
