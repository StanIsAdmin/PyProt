from unittest import TestCase
from pyprot.base.sequence import Sequence
from pyprot.base.aminoacid import AminoAcid


class TestSequence(TestCase):
    def test_getDescription_staysAfterCopy(self):
        description = "this is a description!"
        prot = Sequence("A", description)
        prot2 = Sequence(prot)
        self.assertEquals(prot2.getDescription(), description)

    def test_setNameMode(self):
        prot = Sequence("A")
        self.assertTrue(len(str(prot)) == 1)
        prot.setNameMode("medium")
        self.assertTrue(len(str(prot)) == 3)
        prot.setNameMode("long")
        self.assertTrue(len(str(prot)) > 3)

    def test_setSeparator(self):
        prot = Sequence("AE")
        prot.setSeparator("~")
        self.assertEquals(len(str(prot).split("~")), 2)

    def test_insert(self):
        prot = Sequence("KL")
        prot.insert(0, "M")
        self.assertEquals(prot, Sequence("MKL"))

    def test_extend(self):
        prot = Sequence("OK")
        prot.extend(Sequence("L"))
        prot.extend("M")
        self.assertEquals(prot, Sequence("OKLM"))

    def test_remove(self):
        prot = Sequence("TRANKIL")
        prot.remove(AminoAcid("A"))
        self.assertEquals(prot, Sequence("TRNKIL"))

    def test_delete(self):
        prot = Sequence("YEP")
        del prot[1]
        self.assertEquals(prot, Sequence("YP"))

    def test_count(self):
        prot = Sequence("ABCAX")
        self.assertEquals(prot.count(AminoAcid("A")), 2)

    def test_copyEqualsOriginal(self):
        prot1 = Sequence("methionine")
        prot2 = Sequence(prot1)
        self.assertEquals(prot1, prot2)

    def test_differentNamesEqual(self):
        prot1 = Sequence(["Methionine", "pyrrolysine", "CYSTEINE", "K"])
        prot2 = Sequence("MOCK")
        self.assertEquals(prot1, prot2)
