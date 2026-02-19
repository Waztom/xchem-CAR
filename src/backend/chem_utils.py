"""Chemistry utilities using RDKit for molecular manipulations.

This module contains functions for working with SMILES, SMARTS, molecular
properties, structure comparison, and chemical reactions.
"""

from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import MolStandardize
from rdkit.Chem import inchi

import itertools
import logging
import pandas as pd
import datetime

from .exceptions import (
    ChemistryError,
    MolecularPropertyError,
    SMILESParsingError,
    SMARTSReactionError,
    SVGGenerationError,
)

logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Molecular properties
# -----------------------------------------------------------------------------


def get_mws(smiles: list[str]) -> list[float]:
    """Gets the molecular weights of a list of compounds SMILES

    Parameters
    ----------
    smiles: list[str]
        The SMILES to calculate molecular weights for
    Returns
    -------
    MWs: list[float]
        The list of molecular weights
    """
    try:
        MWs = [
            Descriptors.MolWt(Chem.MolFromSmiles(smi)) for smi in smiles if smi != ""
        ]
        return MWs
    except Exception as e:
        raise MolecularPropertyError(
            f"Failed to compute molecular weights for {smiles}"
        ) from e


def get_inchi_key(smiles: str) -> str:
    """Gets the inchikey for a compound SMILES

    Parameters
    ----------
    smiles: str
        The SMILESs to convert to an inchikey

    Returns
    -------
    inchikeys: str
        The inchikeys
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Could not create molecule from SMILES: {smiles}")
            return None

        # Use the correct function from the inchi module
        inchikey = Chem.inchi.MolToInchiKey(mol)
        return inchikey
    except Exception as e:
        logger.error(f"get_inchi_key yielded error for {smiles}: {e}")
        return None


def get_molecular_formula(smiles: list) -> list:
    """Gets the molecular formula of a list of compounds SMILES
    Parameters
    ----------
    smiles: list[str]
        The SMILES to calculate molecular formula for
    Returns
    -------
    formula: list[str]
    """
    try:
        formula = [
            rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smi)) for smi in smiles
        ]
        return formula
    except Exception as e:
        raise MolecularPropertyError(
            f"Failed to compute molecular formula for {smiles}"
        ) from e


# -----------------------------------------------------------------------------
# SMILES manipulation
# -----------------------------------------------------------------------------


def canon_smiles(smiles: str) -> str:
    """Function to canonicalise SMILES

    Parameters
    ----------
    smiles: str
        The SMILES to be canonicalised

    Returns
    -------
    canon_smiles: str
        The canonicalised SMILES
    status: bool
        Returns False if the input smiles could not be canonicalised
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            canon_smiles = Chem.MolToSmiles(mol)
            return canon_smiles
        else:
            return False
    except Exception as e:
        raise SMILESParsingError(
            f"Failed to canonicalize SMILES '{smiles}'"
        ) from e


def strip_salts(smiles: str, return_details: bool = False):
    """Strips salts from a SMILES string by returning the largest molecular fragment.

    Parameters
    ----------
    smiles : str
        The input SMILES string potentially containing salts
    return_details : bool, optional
        If True, returns additional information about salt stripping

    Returns
    -------
    str or tuple
        If return_details is False: The SMILES string with salts removed
        If return_details is True: A tuple (desalted_smiles, salts_removed, salt_fragments)
        Returns the original SMILES if processing fails

    Examples
    --------
    >>> strip_salts("CC(=O)O.Na")  # Sodium acetate
    'CC(=O)O'
    >>> strip_salts("CC(=O)O.[Na+]", True)  # Sodium acetate with details
    ('CC(=O)O', True, ['[Na+]'])
    """
    try:
        # Canonicalize input SMILES first
        canonical_smiles = canon_smiles(smiles)
        if not canonical_smiles:
            logger.warning(f"Could not canonicalize SMILES: {smiles}")
            return (smiles, False, []) if return_details else smiles

        # Convert to RDKit molecule
        mol = Chem.MolFromSmiles(canonical_smiles)
        if mol is None:
            logger.warning(f"Could not parse SMILES: {canonical_smiles}")
            return (canonical_smiles, False, []) if return_details else canonical_smiles

        # Get fragments
        fragments = Chem.GetMolFrags(mol, asMols=True)
        if len(fragments) <= 1:
            # No salts present
            return (canonical_smiles, False, []) if return_details else canonical_smiles

        # Find the largest fragment by molecular weight
        fragment_weights = [Descriptors.MolWt(frag) for frag in fragments]
        largest_idx = fragment_weights.index(max(fragment_weights))
        main_fragment = fragments[largest_idx]

        # Get salt fragments
        salt_fragments = []
        for i, frag in enumerate(fragments):
            if i != largest_idx:
                salt_fragments.append(Chem.MolToSmiles(frag))

        # Convert to canonical SMILES
        desalted_smiles = Chem.MolToSmiles(main_fragment)

        logger.info(f"Removed salts from {canonical_smiles} -> {desalted_smiles}")
        logger.info(f"Salt fragments: {salt_fragments}")

        if return_details:
            return desalted_smiles, True, salt_fragments
        else:
            return desalted_smiles

    except Exception as e:
        logger.error(f"Error stripping salts from {smiles}: {e}")
        return (smiles, False, []) if return_details else smiles


# -----------------------------------------------------------------------------
# SMARTS matching
# -----------------------------------------------------------------------------


def match_smarts(smiles: str, smarts: str) -> bool:
    """Matches a SMILES pattern to a SMARTS pattern

    Parameters
    ----------
    smiles: str
        The SMILES pattern to match
    smarts: str
        The SMARTS pattern to match

    Returns
    -------
    status: bool
        Returns True if the SMILES pattern matches the SMARTS pattern
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            return True
        else:
            return False
    except Exception as e:
        raise SMARTSReactionError(
            f"Failed to match SMILES '{smiles}' against SMARTS '{smarts}'"
        ) from e


def check_reactant_smarts(reactant_SMILES: tuple, reaction_SMARTS: list) -> list:
    """Checks if reactant pair can produce a product

    Parameters
    ----------
    reactant_SMILES_pair: tuple
        The pair of reactant smiles to check
    reaction_SMARTS: str
        The reaction SMARTS pattern used to check the reactant SMILES

    Returns
    -------
    product_mols: list
        The list of product mols formed between reactant SMILES from SMARTS pattern
    status: None
        Returns None if no product mols are formed
    """
    reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactant_SMILES if smi != ""]
    if len(reaction_SMARTS) == 1:
        SMARTS_pattern = reaction_SMARTS[0]
        rxn = AllChem.ReactionFromSmarts(SMARTS_pattern)
        product_mols = None
        if len(reactant_mols) == 1:
            try:
                products = rxn.RunReactants(reactant_mols)
                product_mols = [product[0] for product in products]
                if product_mols:
                    product_smiles = set(
                        [Chem.MolToSmiles(mol) for mol in product_mols]
                    )
                    product_mols = [Chem.MolFromSmiles(smi) for smi in product_smiles]
            except Exception as e:
                logger.warning(f"check_reactant_smarts yielded error: {e}")

        if len(reactant_mols) > 1:
            logger.debug("The number of reactant mols is > 1")
            for reactant_permutation in list(itertools.permutations(reactant_mols)):
                try:
                    products = rxn.RunReactants(reactant_permutation)
                    product_mols = [product[0] for product in products]
                    if product_mols:
                        logger.debug(f"The product mols are {product_mols}")
                        product_smiles = set(
                            [Chem.MolToSmiles(mol) for mol in product_mols]
                        )
                        product_mols = [
                            Chem.MolFromSmiles(smi) for smi in product_smiles
                        ]
                        break
                    if not product_mols:
                        continue  # reactants were in wrong order so no product
                except Exception as e:
                    logger.warning(f"check_reactant_smarts yielded error: {e}")
                    logger.debug(f"Permutation: {reactant_permutation}")
                    continue
        if product_mols:
            return [product_mols[0]]
        else:
            logger.debug(f"No products for SMARTS={SMARTS_pattern}, reactants={reactant_SMILES}")
            return None

    if len(reaction_SMARTS) > 1:
        product_mols = []
        if len(reactant_mols) == 1:
            logger.debug("The number of reactant mols is 1")
            try:
                product_mol = None
                for SMARTS_pattern in reaction_SMARTS:
                    rxn = AllChem.ReactionFromSmarts(SMARTS_pattern)
                    if product_mol is None:
                        products = rxn.RunReactants(reactant_mols)
                    else:
                        products = rxn.RunReactants(product_mol)
                    product_mol = [products[0][0]]
                    logger.debug(f"The product mol is {product_mol}")
                    logger.debug(
                        f"The product smiles is {Chem.MolToSmiles(product_mol[0])}"
                    )
                    if product_mol:
                        product_mols.append(product_mol[-1])
            except Exception as e:
                logger.warning(
                    f"check_reactant_smarts yielded error: {e} "
                    f"(reactants={reactant_SMILES}, SMARTS={SMARTS_pattern})"
                )
        if len(reactant_mols) > 1:
            logger.debug("The number of reactant mols is > 1")
            try:
                for reactant_permutation in list(itertools.permutations(reactant_mols)):
                    product_mol = None
                    for SMARTS_pattern in reaction_SMARTS:
                        rxn = AllChem.ReactionFromSmarts(SMARTS_pattern)
                        if product_mol is None:
                            products = rxn.RunReactants(reactant_mols)
                        else:
                            products = rxn.RunReactants(product_mol)
                        product_mol = [products[0][0]]
                        logger.debug(f"The product mol is {product_mol}")
                        logger.debug(
                            f"The product smiles is {Chem.MolToSmiles(product_mol[0])}"
                        )
                        if product_mol:
                            product_mols.extend(product_mol[-1])
                        if not product_mol:
                            continue
            except Exception as e:
                logger.warning(f"check_reactant_smarts yielded error: {e}")
                logger.debug(f"Permutation: {reactant_permutation}")
        if len(product_mols) != 0:
            return product_mols
        else:
            logger.debug(f"No products for SMARTS={SMARTS_pattern}, reactants={reactant_SMILES}")
            return None


# -----------------------------------------------------------------------------
# Reaction utilities
# -----------------------------------------------------------------------------


def get_addition_order(
    product_smi: str, reactant_SMILES: tuple, reaction_SMARTS: list
) -> list:
    """Gets reactant pair addition order as SMILES that yields the expected
       product via the reaction SMARTS pattern

    Parameters
    ----------
    product_smi: str
        The product SMILES
    reactant_SMILES_pair: tuple
        The reactant SMILES pair for a reaction
    reaction_SMARTS: lists
        The reaction SMARTS pattern. List for multiple SMARTS patterns/reaction transformations.

    Returns
    -------
    reactant_SMILES_pair: list
        The list of ordered reactant smiles
    status: None
        None if no order can create the input product
    """
    try:
        reaction_SMARTS = reaction_SMARTS[0]
        reaction_SMARTS_reactants = reaction_SMARTS.split(">>")[0].split(".")
        rxn = AllChem.ReactionFromSmarts(reaction_SMARTS)
        reactant_mols = [
            Chem.MolFromSmiles(smi) for smi in reactant_SMILES if smi != ""
        ]
        ordered_smis = None

        if reaction_SMARTS_reactants == "":
            ordered_smis = [canon_smiles(smi) for smi in reactant_SMILES if smi != ""]

        if len(reactant_mols) == 1:
            ordered_smis = [canon_smiles(smi) for smi in reactant_SMILES if smi != ""]

        if len(reactant_mols) > 1:

            if Chem.MolFromSmarts(reaction_SMARTS_reactants[0]) == Chem.MolFromSmarts(
                reaction_SMARTS_reactants[1]
            ):
                ordered_smis = [
                    canon_smiles(smi) for smi in reactant_SMILES if smi != ""
                ]
            else:
                for reactant_permutation in list(itertools.permutations(reactant_mols)):
                    try:
                        products = rxn.RunReactants(reactant_permutation)
                        product_mols = [product[0] for product in products]
                        if not product_mols:
                            continue  # reactants were in wrong order so no product
                    except Exception as e:
                        logger.warning(f"get_addition_order yielded error: {e}")
                        logger.debug(f"Permutation: {reactant_permutation}")
                        continue
                    product_smis = [
                        Chem.MolToSmiles(m) for m in product_mols if m is not None
                    ]
                    if product_smi in product_smis:
                        ordered_smis = [
                            Chem.MolToSmiles(m) for m in reactant_permutation
                        ]
        if ordered_smis is not None:
            return ordered_smis
        else:
            logger.debug(
                f"Addition order not found for product SMILES {product_smi}, "
                f"SMARTS={reaction_SMARTS}, reactants={reactant_SMILES}"
            )
            return None
    except Exception as e:
        raise ChemistryError(
            f"Failed to determine addition order for product '{product_smi}'"
        ) from e



# -----------------------------------------------------------------------------
# Molecular manipulation
# -----------------------------------------------------------------------------


def atom_remover(mol, rxn):
    """Remove atoms from a molecule using a reaction pattern.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input molecule to modify
    rxn : rdkit.Chem.rdChemReactions.ChemicalReaction
        The reaction pattern to apply to the molecule

    Returns
    -------
    rdkit.Chem.rdchem.Mol or None
        The modified molecule with atoms removed, or None if the reaction fails

    Examples
    --------
    >>> from rdkit import Chem
    >>> from rdkit.Chem import AllChem
    >>> mol = Chem.MolFromSmiles('CC(=O)O')
    >>> rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])[OH]>>[C:1](=[O:2])')
    >>> result = atom_remover(mol, rxn)
    """
    try:
        ps = rxn.RunReactants((mol,))

        logger.debug(f"Attempting to run reaction on molecule: {Chem.MolToSmiles(mol)}")

        if not ps:
            logger.warning("Could not run the reaction, returning original molecule")
            return Chem.Mol(mol)

        for p in ps:
            res = Chem.RemoveHs(p[0])
            logger.info(f"Successfully removed atoms, result: {Chem.MolToSmiles(res)}")
            return res

    except Exception as e:
        logger.error(f"Error in atom_remover: {str(e)}")
        return None


def get_frags(mols: list, smarts: str) -> list:
    """Get the fragments of a list of molecules"
    Parameters
    ----------
    frags: list[rdkit.Chem.rdchem.Mol]
        The molecules to fragment

    Returns
    -------
    frags: list[rdkit.Chem.rdchem.Mol]
        The fragments of the input molecules
    """
    frag_mols = []
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
    except Exception as e:
        raise SMARTSReactionError(
            f"Failed to parse SMARTS pattern '{smarts}'"
        ) from e

    for mol in mols:
        try:
            ps = rxn.RunReactants((mol,))
            if not ps:
                frag_mols.append(None)
                continue
            for p in ps:
                res = Chem.RemoveHs(p[0])
                frag_mols.append(res)
        except Exception as e:
            logger.error(f"Error in get_frags for molecule: {e}")
            frag_mols.append(None)
            continue
    return frag_mols


def remove_radicals(mol):
    """Remove radicals from a molecule by adding hydrogens.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input molecule that may contain radicals

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        The molecule with radicals removed, or None if the operation fails

    Examples
    --------
    >>> from rdkit import Chem
    >>> mol = Chem.MolFromSmiles('[CH2]CC')
    >>> result = remove_radicals(mol)
    >>> Chem.MolToSmiles(result)
    'CCC'
    """
    try:
        for atom in mol.GetAtoms():
            if atom.GetNumRadicalElectrons() > 0:
                atom.SetNumRadicalElectrons(0)
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)

        Chem.SanitizeMol(mol)

        return mol

    except Exception as e:
        logger.error(f"Error removing radicals: {str(e)}")
        return None


# -----------------------------------------------------------------------------
# Structure comparison
# -----------------------------------------------------------------------------


def are_equivalent_structures(smiles1, smiles2, match_tautomers=True, similarity_threshold=0.9):
    """
    Determines if two chemical structures are equivalent using multiple approaches.
    Handles tautomers and different structural representations of the same compound.

    Parameters
    ----------
    smiles1 : str
        First SMILES string
    smiles2 : str
        Second SMILES string
    match_tautomers : bool
        Whether to attempt tautomer matching (more computationally expensive)
    similarity_threshold : float
        Threshold for fingerprint similarity (0-1)

    Returns
    -------
    bool
        True if structures are considered equivalent, False otherwise
    """
    try:
        # Handle empty or None inputs
        if not smiles1 or not smiles2:
            return False
            
        # Method 1: Direct canonicalized SMILES comparison after stripping salts
        clean1 = strip_salts(canon_smiles(smiles1))
        clean2 = strip_salts(canon_smiles(smiles2))
        
        if clean1 == clean2:
            logger.debug(f"Structures matched by direct SMILES comparison")
            return True
        
        mol1 = Chem.MolFromSmiles(clean1)
        mol2 = Chem.MolFromSmiles(clean2)
        
        if not mol1 or not mol2:
            logger.warning(f"Could not parse molecules from SMILES")
            return False
            
        # Method 3: Tautomer standardization (if enabled)
        if match_tautomers:
            try:
                enumerator = MolStandardize.tautomerEnumerator.TautomerEnumerator()
                
                # Get canonical tautomers
                canon_taut1 = enumerator.Canonicalize(mol1)
                canon_taut2 = enumerator.Canonicalize(mol2)
                
                if canon_taut1 and canon_taut2:
                    smiles_taut1 = Chem.MolToSmiles(canon_taut1)
                    smiles_taut2 = Chem.MolToSmiles(canon_taut2)
                    
                    if smiles_taut1 == smiles_taut2:
                        logger.debug(f"Structures matched by tautomer canonicalization")
                        return True
                        
                    # Try full enumeration as a last resort
                    tautomers1 = {Chem.MolToSmiles(x) for x in enumerator.Enumerate(mol1)}
                    tautomers2 = {Chem.MolToSmiles(x) for x in enumerator.Enumerate(mol2)}
                    
                    if tautomers1.intersection(tautomers2):
                        logger.debug(f"Structures matched by tautomer enumeration")
                        return True
            except Exception as e:
                logger.debug(f"Tautomer comparison failed: {str(e)}")
        
        # If we get here, no match was found
        return False
                
    except Exception as e:
        logger.error(f"Error comparing structures: {str(e)}")
        return False  # Better to return False than None on error


# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------


def create_svg_string(smiles: str) -> str:
    """Function that creates a SVG image string from smiles string

    Parameters
    ----------
    smiles: string
        The SMILES to create an SVG image string from

    Returns
    -------
    svg_string: string
        The SVG image string
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(100, 50)
        drawer.SetFontSize(8)
        drawer.SetLineWidth(1)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg_string = drawer.GetDrawingText()
        return svg_string
    except Exception as e:
        raise SVGGenerationError(
            f"Failed to generate SVG for SMILES '{smiles}'"
        ) from e


def create_reaction_svg_string(smarts: str) -> str:
    """Function that creates a SVG image string from smarts string

    Parameters
    ----------
    smarts: string
        The SMARTS reaction pattern to create an SVG image
        string from

    Returns
    -------
    svg_string: string
        The SVG image string of the SMARTS pattern
    """
    try:
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(900, 200)
        drawer.DrawReaction(smarts)
        drawer.FinishDrawing()
        svg_string = drawer.GetDrawingText()
        return svg_string
    except Exception as e:
        raise SVGGenerationError(
            f"Failed to generate reaction SVG for SMARTS '{smarts}'"
        ) from e


# -----------------------------------------------------------------------------
# Combinatorial chemistry
# -----------------------------------------------------------------------------


def combi_chem(
    reactant_1_SMILES: list, reactant_2_SMILES: list, are_product_SMILES: bool = False
) -> list:
    """Gets all possible combinations between two uneven lists of
       reactants

    Parameters
    ----------
    reactant_1_SMILES: list
        The list of reactant one smiles
    reactant_2_SMILES: list
        The second list of reactant two smiles
    are_product_SMILES: boolean
        Set to True if reactant_2_SMILES is list of products from
        previous reaction step

    Returns
    -------
    all_possible_combinations: list
        All possible reactant combinations possible
        between reactat 1 and reactant two lists
        as a list of tuples
    """
    try:
        if len(reactant_1_SMILES) == 0:
            reactant_2_SMILES_canon = [canon_smiles(smi) for smi in reactant_2_SMILES]
            if not are_product_SMILES:
                reactant_2_SMILES_canon = list(dict.fromkeys(reactant_2_SMILES_canon))
            all_possible_combinations = list(
                itertools.product([""], reactant_2_SMILES_canon)
            )
        if len(reactant_2_SMILES) == 0:
            reactant_1_SMILES_canon = [canon_smiles(smi) for smi in reactant_1_SMILES]
            if not are_product_SMILES:
                reactant_1_SMILES_canon = list(dict.fromkeys(reactant_1_SMILES_canon))
            all_possible_combinations = list(
                itertools.product([""], reactant_1_SMILES_canon)
            )
        if len(reactant_1_SMILES) != 0 and len(reactant_2_SMILES) != 0:
            reactant_1_SMILES_canon = [canon_smiles(smi) for smi in reactant_1_SMILES]
            reactant_2_SMILES_canon = [canon_smiles(smi) for smi in reactant_2_SMILES]
            if not are_product_SMILES:
                reactant_2_SMILES_canon = list(dict.fromkeys(reactant_2_SMILES_canon))
            all_possible_combinations = list(
                itertools.product(
                    list(dict.fromkeys(reactant_1_SMILES_canon)), reactant_2_SMILES_canon
                )
            )
        return all_possible_combinations
    except Exception as e:
        raise ChemistryError(
            f"Failed to compute combinatorial chemistry combinations for "
            f"reactant lists '{reactant_1_SMILES}' and '{reactant_2_SMILES}'"
        ) from e    


def create_combi_chem_csv(csv_input_file: str, out_dir: str):
    """Creates a .csv file for all the combinations possible for a given input
        of reactant SMILES pairs

    Parameters
    ----------
    csv_input_file: str
        The path to .csv file to read the reactant SMILES pairs and reactant classes from
    out_dir: str
        The directory to write the csv to
    """
    from .recipe_utils import get_recipe_smarts

    try:
        output_list = []
        input_df = pd.read_csv(csv_input_file)
        grouped_df = input_df.groupby(["reactant_class", "reaction_recipe"])
        for group in grouped_df:
            reaction_classes = group[1]["reactant_class"].tolist()
            reaction_recipes = group[1]["reaction_recipe"].tolist()
            reactant_1_SMILES = group[1]["reactant_1"].tolist()
            reactant_2_SMILES = group[1]["reactant_2"].tolist()
            reaction_SMARTS = get_recipe_smarts(
                reaction_classes[0], reaction_recipes[0]
            )

            all_possible_combinations = combi_chem(reactant_1_SMILES, reactant_2_SMILES)
            product_smiles = []
            for reactant_pair in all_possible_combinations:
                product_mols = check_reactant_smarts(
                    reactant_SMILES=reactant_pair, reaction_SMARTS=reaction_SMARTS
                )
                product_smiles.append(Chem.MolToSmiles(product_mols[0]))

            all_possible_reactant_1_smiles = [x[0] for x in all_possible_combinations]
            all_possible_reactant_2_smiles = [x[1] for x in all_possible_combinations]

            output_list.append(
                [
                    all_possible_reactant_1_smiles,
                    all_possible_reactant_2_smiles,
                    product_smiles,
                    reaction_classes,
                    reaction_recipes,
                ]
            )
        out_df = pd.DataFrame(
            output_list,
            columns=[
                "reactant_1_SMILES",
                "reactant_2_SMILES",
                "product_SMILES",
                "reaction_class",
                "reaction_recipe",
            ],
        )
        out_df.to_csv(
            out_dir
            + "{}-combi-chem.csv".format(
                datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
            )
        )

    except Exception as e:
        raise ChemistryError(
            f"Failed to create combinatorial chemistry CSV from '{csv_input_file}'"
        ) from e
