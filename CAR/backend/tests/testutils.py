from unittest import TestCase
from rdkit import Chem
from rdkit.Chem import AllChem

from car.utils import (
    getMWs,
    getInchiKey,
    calculateMolsFromConc,
    calculateMassFromMols,
    canonSmiles,
    combiChem,
    createSVGString,
    createReactionSVGString,
    getAddtionOrder,
    checkReactantSMARTS,
    getChemicalName,
    getPubChemCAS,
    getPubChemCompound,
)

from .testdata.indata.testutils import (
    snar_reactant_smiles_one,
    snar_reactant_smiles_two,
)
from .testdata.outdata.testutils import (
    snar_combo_equal,
    snar_combo_unequal,
    svg_str,
    svg_reaction_str,
)


class ChemistryFunctionsTestCase(TestCase):
    def setUp(self) -> None:
        self.smiles = "C1COC2=C(C3=C(C(=C21)CCN)OCC3)Br"
        self.snar_reactant_smiles_one = snar_reactant_smiles_one
        self.snar_reactant_smiles_two = snar_reactant_smiles_two
        self.snar_reactant_smiles_tuple = (
            snar_reactant_smiles_one[0],
            snar_reactant_smiles_two[0],
        )
        self.snar_encoded_smarts = [
            "[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]"
        ]
        self.snar_product_smiles = "O=C(O)Cc1ccc(Nc2ccccc2)cc1F"
        self.snar_product_mols = Chem.MolFromSmiles(self.snar_product_smiles)
        self.svg_str = self.strip_white_space(str=svg_str)
        self.reaction_smarts = AllChem.ReactionFromSmarts(
            "{}.{}>>{}".format(
                self.snar_reactant_smiles_one[0],
                self.snar_reactant_smiles_two[0],
                self.snar_product_smiles,
            )
        )
        self.reaction_svg_str = self.strip_white_space(svg_reaction_str)
        self.boc_tbs_reaction_smarts = [
            "[#7:2]-[#6](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])>>[#7:2]",
            "[#6:1]-[#8:2]-[#14]([#6H3])([#6H3])-[#6]([#6H3])([#6H3])([#6H3])>>[#6:1]-[#8:2]-[#1]",
        ]
        self.boc_tbs_reactant_smiles = (
            "CC(C)(C)OC(=O)NC(CCN)C(=O)Nc1ccccc1Cc1cccc2c(cc(C)nc12)CO[Si](C)(C)C(C)(C)C",
            "",
        )
        self.boc_tbs_product_smiles = [
            "Cc1cc(CO[Si](C)(C)C(C)(C)C)c2cccc(Cc3ccccc3NC(=O)C(N)CCN)c2n1",
            "[H]OCc1cc(C)nc2c(Cc3ccccc3NC(=O)C(N)CCN)cccc12",
        ]
        self.ester_double_reaction_smarts = [
            "[#6:1](=[#8:2])-[#8:3][#6H1,#6H2,#6H3]>>[#6:1](=[#8:2])-[#8H:3]",
            "[#6:1](=[#8:2])-[#8:3][#6H1,#6H2,#6H3]>>[#6:1](=[#8:2])-[#8H:3]",
        ]
        self.ester_double_reactant_smiles = (
            "O=C(N[C@H](Cc1cccc2ccc(C)nc21)C(=O)OC)c1cc([NH]n1)C(=O)OC",
            "",
        )
        self.ester_double_product_smiles = [
            "COC(=O)c1cc(C(=O)N[C@H](Cc2cccc3ccc(C)nc23)C(=O)O)n[nH]1",
            "Cc1ccc2cccc(C[C@@H](NC(=O)c3cc(C(=O)O)[nH]n3)C(=O)O)c2n1",
        ]

    def strip_white_space(self, str):
        return str.replace(" ", "").replace("\t", "").replace("\n", "")

    def test_get_mws(self):
        test_mws = getMWs(smiles=self.snar_reactant_smiles_tuple)
        self.assertAlmostEqual(
            first=test_mws[0],
            second=172.13,
            places=2,
            msg="incorrect molecular weight calculated",
        )
        self.assertAlmostEqual(
            first=test_mws[1],
            second=93.13,
            places=2,
            msg="incorrect molecular weight calculated",
        )

    def test_get_mws_fail(self):
        test_mws = getMWs(smiles=["OT Chemistry is possible"])
        self.assertEqual(
            test_mws,
            None,
            "incorrect capture of bad SMILES input",
        )

    def test_get_inchi_key(self):
        test_inchi_key = getInchiKey(smiles=self.smiles)
        self.assertEqual(
            test_inchi_key,
            "YZDFADGMVOSVIX-UHFFFAOYSA-N",
            "incorrect InChI key calculated",
        )

    def test_get_inchi_key_fail(self):
        test_inchi_key = getInchiKey(smiles="OT Chemistry is possible")
        self.assertEqual(
            test_inchi_key,
            None,
            "incorrect capture of bad SMILES input",
        )

    def test_calculate_mols_from_conc(self):
        test_product_mols = calculateMolsFromConc(
            target_concentration=0.5,
            target_volume=1,
        )

        self.assertAlmostEqual(
            first=test_product_mols,
            second=5e-10,
            places=4,
            msg="incorrect product mols calculated",
        )

    def test_calculate_mols_fail(self):
        test_product_mols = calculateMolsFromConc(
            target_concentration="090", target_volume="OT Chemistry is possible"
        )

        self.assertEqual(
            test_product_mols,
            None,
            "incorrect capture of bad target mass and SMILES input",
        )

    def test_calculate_mass_from_mols(self):
        mass = calculateMassFromMols(mols=0.5, SMILES=self.snar_product_smiles)
        self.assertAlmostEqual(
            mass,
            122626.5,
            places=1,
            msg="Incorrect mass calculated from mols and SMILES",
        )

    def test_calculate_mass_from_mols_fail(self):
        mass = calculateMassFromMols(mols=0.5, SMILES="OT Chemistry is possible")
        self.assertEqual(
            mass,
            None,
            "Incorrect capture of bad SMILES input",
        )

    def test_canon_smiles(self):
        test_canon_smiles = canonSmiles(smiles=self.smiles)
        self.assertEqual(
            test_canon_smiles,
            "NCCc1c2c(c(Br)c3c1OCC3)OCC2",
            "incorrect canonicalisation of SMILES",
        )

    def test_canon_smiles_fail(self):
        test_canon_smiles = canonSmiles(smiles="OT Chemistry is possible")
        self.assertEqual(
            test_canon_smiles,
            False,
            "incorrect capture of bad SMILES input",
        )

    def test_combi_chem(self):
        test_all_possible_combinations = combiChem(
            reactant_1_SMILES=self.snar_reactant_smiles_one,
            reactant_2_SMILES=self.snar_reactant_smiles_two,
        )
        self.assertEqual(
            test_all_possible_combinations,
            snar_combo_equal,
            "incorrect combinatorial product of two equal length lists of smiles",
        )

    def test_combi_chem_fail(self):
        test_all_possible_combinations = combiChem(
            reactant_1_SMILES=self.snar_reactant_smiles_one[0:1],
            reactant_2_SMILES=self.snar_reactant_smiles_two,
        )
        self.assertEqual(
            test_all_possible_combinations,
            snar_combo_unequal,
            "incorrect combinatorial product of two non-equal length lists of smiles",
        )

    def test_create_svg_string(self):
        test_svg_str = createSVGString(smiles=self.smiles)
        test_svg_str = self.strip_white_space(test_svg_str)
        self.assertEqual(
            test_svg_str,
            self.svg_str,
            "incorrect creation of a svg string from SMILES",
        )

    def test_create_reaction_svg_string(self):
        test_svg_str = createReactionSVGString(smarts=self.reaction_smarts)
        test_svg_str = self.strip_white_space(test_svg_str)
        self.assertEqual(
            test_svg_str,
            self.reaction_svg_str,
            "incorrect creation of a reaction svg string from reaction SMARTS",
        )

    def test_get_addition_order_success(self):
        test_ordered_smis = getAddtionOrder(
            product_smi=self.snar_product_smiles,
            reactant_SMILES=self.snar_reactant_smiles_tuple,
            reaction_SMARTS=self.snar_encoded_smarts,
        )
        self.assertEqual(
            test_ordered_smis,
            ["Nc1ccccc1", "O=C(O)Cc1ccc(F)cc1F"],
            "incorrect addtion order for a SNAr encoded recipe SMARTS",
        )

    def test_get_addition_order_fail(self):
        first_reactant_smiles = "OT Chemistry is possible"
        second_reactant_smiles = self.snar_reactant_smiles_one[0]
        reactant_SMILES = (first_reactant_smiles, second_reactant_smiles)
        test_ordered_smis = getAddtionOrder(
            product_smi=self.snar_product_smiles,
            reactant_SMILES=reactant_SMILES,
            reaction_SMARTS=self.snar_encoded_smarts,
        )
        self.assertEqual(
            test_ordered_smis,
            None,
            "incorrect SMILES input should return None for addtion order",
        )

    def test_check_reactant_smarts_success(self):
        test_product_mols = checkReactantSMARTS(
            reactant_SMILES=self.snar_reactant_smiles_tuple,
            reaction_SMARTS=self.snar_encoded_smarts,
        )
        test_product_smiles = [Chem.MolToSmiles(mol) for mol in test_product_mols]

        self.assertEqual(
            len(test_product_mols),
            2,
            "incorrect length of product mols for testing reaction SMARTS",
        )
        self.assertEqual(
            set(test_product_smiles),
            set(["O=C(O)Cc1ccc(Nc2ccccc2)cc1F", "O=C(O)Cc1ccc(F)cc1Nc1ccccc1"]),
            "incorrect product SMILES match for testing reaction SMARTS",
        )

    def test_more_one_reaction_smarts_success_boc_tbs(self):
        test_product_mols = checkReactantSMARTS(
            reactant_SMILES=self.boc_tbs_reactant_smiles,
            reaction_SMARTS=self.boc_tbs_reaction_smarts,
        )
        test_product_smiles = [Chem.MolToSmiles(mol) for mol in test_product_mols]
        self.assertEqual(
            len(test_product_mols),
            2,
            "incorrect length of product mols for testing reaction SMARTS",
        )
        self.assertEqual(
            set(test_product_smiles),
            set(self.boc_tbs_product_smiles),
            "incorrect product SMILES match for testing reaction SMARTS",
        )

    def test_more_one_reaction_smarts_success_double_ester(self):
        test_product_mols = checkReactantSMARTS(
            reactant_SMILES=self.ester_double_reactant_smiles,
            reaction_SMARTS=self.ester_double_reaction_smarts,
        )
        test_product_smiles = [Chem.MolToSmiles(mol) for mol in test_product_mols]
        self.assertEqual(
            len(test_product_mols),
            2,
            "incorrect length of product mols for testing reaction SMARTS",
        )
        self.assertEqual(
            set(test_product_smiles),
            set(self.ester_double_product_smiles),
            "incorrect product SMILES match for testing reaction SMARTS",
        )


class PubChemFunctionsTestCase(TestCase):
    def setUp(self) -> None:
        self.smiles = "C1COC2=C(C3=C(C(=C21)CCN)OCC3)Br"
        self.inchikey = "YZDFADGMVOSVIX-UHFFFAOYSA-N"
        self.compound = getPubChemCompound(inchikey=self.inchikey)

    def test_get_pubchem_compound(self):
        compound = getPubChemCompound(inchikey=self.inchikey)
        self.assertEqual(
            compound.cid,
            10265873,
            "incorrect PubChem id for getting a compound from PubChem",
        )

    def test_get_pubchem_compound_fail(self):
        test_compound = getPubChemCompound(inchikey="OT chemistry is possible")
        self.assertEqual(
            test_compound,
            None,
            "PubChem search should yield None response",
        )

    def test_get_pubchem_cas(self):
        test_cas = getPubChemCAS(compound=self.compound)
        self.assertEqual(test_cas, "733720-95-1", "incorrect CAS number returned")

    def test_get_chemical_name(self):
        test_name = getChemicalName(inchikey=self.inchikey)
        self.assertEqual(
            test_name,
            "2-(4-bromo-2,3,6,7-tetrahydrofuro[2,3-f][1]benzofuran-8-yl)ethanamine",
            "incorrect PubChem IUPAC name for compound",
        )

    def test_get_chemical_name_fail(self):
        test_name = getChemicalName(inchikey="OT chemistry is possible")
        self.assertEqual(
            test_name,
            None,
            "PubChem name search should fail and return None",
        )

    def test_get_pubchem_cas(self):
        test_cas = getPubChemCAS(compound=self.compound)
        self.assertEqual(test_cas, "178557-21-6", "incorrect CAS number returned")

    def get_pubchem_cas_fail(self):
        test_cas = getPubChemCAS(compound=None)
        self.assertEqual(test_cas, None, "CAS should be None")
