from unittest import TestCase

from backend.conversions import calculateMolsFromConc, calculateMassFromMols


class ConversionsFunctionsTestCase(TestCase):
    def setUp(self) -> None:
        self.snar_product_smiles = "O=C(O)Cc1ccc(Nc2ccccc2)cc1F"

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
