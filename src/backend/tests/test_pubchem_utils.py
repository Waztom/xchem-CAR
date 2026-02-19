from unittest import TestCase

from backend.pubchem_utils import get_chemical_name, get_pubchem_cas, get_pubchem_compound


class PubChemFunctionsTestCase(TestCase):
    def setUp(self) -> None:
        self.smiles = "C1COC2=C(C3=C(C(=C21)CCN)OCC3)Br"
        self.inchikey = "YZDFADGMVOSVIX-UHFFFAOYSA-N"
        self.compound = get_pubchem_compound(inchikey=self.inchikey)

    def test_get_pubchem_compound(self):
        compound = get_pubchem_compound(inchikey=self.inchikey)
        self.assertEqual(
            compound.cid,
            10265873,
            "incorrect PubChem id for getting a compound from PubChem",
        )

    def test_get_pubchem_compound_fail(self):
        test_compound = get_pubchem_compound(inchikey="OT chemistry is possible")
        self.assertEqual(
            test_compound,
            None,
            "PubChem search should yield None response",
        )

    def test_get_pubchem_cas(self):
        test_cas = get_pubchem_cas(compound=self.compound)
        self.assertEqual(test_cas, "178557-21-6", "incorrect CAS number returned")

    def test_get_chemical_name(self):
        test_name = get_chemical_name(inchikey=self.inchikey)
        self.assertEqual(
            test_name,
            "2-(4-bromo-2,3,6,7-tetrahydrofuro[2,3-f][1]benzofuran-8-yl)ethanamine",
            "incorrect PubChem IUPAC name for compound",
        )

    def test_get_chemical_name_fail(self):
        test_name = get_chemical_name(inchikey="OT chemistry is possible")
        self.assertEqual(
            test_name,
            None,
            "PubChem name search should fail and return None",
        )

    def test_get_pubchem_cas_fail(self):
        test_cas = get_pubchem_cas(compound=None)
        self.assertEqual(test_cas, None, "CAS should be None")
