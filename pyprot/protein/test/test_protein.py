from unittest import TestCase
from pyprot.protein.protein import Protein
from pyprot.protein.aminoacid import AminoAcid


class TestProtein(TestCase):
    def test_getDescription_staysAfterCopy(self):
        description = "this is a description!"
        prot = Protein("A", description)
        prot2 = Protein(prot)
        self.assertEquals(prot2.getDescription(), description)

    def test_setNameMode(self):
        prot = Protein("A")
        self.assertTrue(len(str(prot)) == 1)
        prot.setNameMode("medium")
        self.assertTrue(len(str(prot)) == 3)
        prot.setNameMode("long")
        self.assertTrue(len(str(prot)) > 3)

    def test_setSeparator(self):
        prot = Protein("AE")
        prot.setSeparator("~")
        self.assertEquals(len(str(prot).split("~")), 2)

    def test_insert(self):
        prot = Protein("KL")
        prot.insert(0, "M")
        self.assertEquals(prot, Protein("MKL"))

    def test_extend(self):
        prot = Protein("OK")
        prot.extend(Protein("L"))
        prot.extend("M")
        self.assertEquals(prot, Protein("OKLM"))

    def test_remove(self):
        prot = Protein("TRANKIL")
        prot.remove(AminoAcid("A"))
        self.assertEquals(prot, Protein("TRNKIL"))

    def test_delete(self):
        prot = Protein("YEP")
        del prot[1]
        self.assertEquals(prot, Protein("YP"))

    def test_count(self):
        prot = Protein("ABCAX")
        self.assertEquals(prot.count(AminoAcid("A")), 2)

    def test_copyEqualsOriginal(self):
        prot1 = Protein("methionine")
        prot2 = Protein(prot1)
        self.assertEquals(prot1, prot2)

    def test_differentNamesEqual(self):
        prot1 = Protein(["Methionine", "pyrrolysine", "CYSTEINE", "K"])
        prot2 = Protein("MOCK")
        self.assertEquals(prot1, prot2)
