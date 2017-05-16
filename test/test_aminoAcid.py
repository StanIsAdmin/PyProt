from unittest import TestCase

from pyprot.protein import AminoAcid


class TestAminoAcid(TestCase):
    def test_getAllNames(self):
        for name in AminoAcid.getAllNames():
            try:
                aa = AminoAcid(name)
            except ValueError:
                self.fail()

    def test_getName(self):
        aa = AminoAcid("ala")
        self.assertEquals("A", aa.getName())

    def test_isGap(self):
        aa = AminoAcid("gap")
        self.assertTrue(aa.isGap())

if __name__ == "__main__":
    TestAminoAcid.main()