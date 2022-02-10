import unittest, os, random, time, copy
from PHANTASM.taxonomy.TaxRank import TaxRank


class TestTaxRank(unittest.TestCase):
    # constants
    SPE_STRING = 'species'
    SPE_PLURAL = 'species'
    GEN_STRING = 'genus'
    GEN_PLURAL = 'genera'
    FAM_STRING = 'family'
    FAM_PLURAL = 'families'
    ORD_STRING = 'order'
    ORD_PLURAL = 'orders'
    CLS_STRING = 'class'
    CLS_PLURAL = 'classes'
    PHY_STRING = 'phylum'
    PHY_PLURAL = 'phyla'
    DOM_STRING = 'domain'
    DOM_PLURAL = 'domains'
    SUP_STRING = 'superkingdom'


    # tests
    def testA_isValidRank(self) -> None:
        """ Tests for expected functionality of isValidRank()
        """
        # expected strings should be valid
        self.assertTrue(TaxRank.isValidRank(TestTaxRank.SPE_STRING))
        self.assertTrue(TaxRank.isValidRank(TestTaxRank.GEN_STRING))
        self.assertTrue(TaxRank.isValidRank(TestTaxRank.FAM_STRING))
        self.assertTrue(TaxRank.isValidRank(TestTaxRank.ORD_STRING))
        self.assertTrue(TaxRank.isValidRank(TestTaxRank.CLS_STRING))
        self.assertTrue(TaxRank.isValidRank(TestTaxRank.PHY_STRING))
        self.assertTrue(TaxRank.isValidRank(TestTaxRank.DOM_STRING))
        self.assertTrue(TaxRank.isValidRank(TestTaxRank.SUP_STRING))

        # unexpected strings should be invalid (species: singular == plural)
        self.assertFalse(TaxRank.isValidRank('fake'))
        self.assertFalse(TaxRank.isValidRank('subspecies'))
        self.assertFalse(TaxRank.isValidRank('subphylum'))
        self.assertFalse(TaxRank.isValidRank('kingdom'))
        self.assertFalse(TaxRank.isValidRank(TestTaxRank.GEN_PLURAL))
        self.assertFalse(TaxRank.isValidRank(TestTaxRank.FAM_PLURAL))
        self.assertFalse(TaxRank.isValidRank(TestTaxRank.ORD_PLURAL))
        self.assertFalse(TaxRank.isValidRank(TestTaxRank.CLS_PLURAL))
        self.assertFalse(TaxRank.isValidRank(TestTaxRank.PHY_PLURAL))

    def testB_rankCeiling(self) -> None:
        """ Tests that there is no rank above 'domain'
        """
        rank = TaxRank(TestTaxRank.DOM_STRING)
        self.assertEqual(rank.getRankAbove(), None)
    
    def testC_rankFloor(self) -> None:
        """ Tests that there is no rank below 'species'
        """
        rank = TaxRank(TestTaxRank.SPE_STRING)
        self.assertEqual(rank.getRankBelow(), None)
    
    def testD_comparisonOperators(self) -> None:
        """ Tests __eq__, __ne__, __gt__, __ge__, __lt__, and __le__
        """
        # make all possible ranks
        spe = TaxRank(TestTaxRank.SPE_STRING)
        gen = TaxRank(TestTaxRank.GEN_STRING)
        fam = TaxRank(TestTaxRank.FAM_STRING)
        ord = TaxRank(TestTaxRank.ORD_STRING)
        cls = TaxRank(TestTaxRank.CLS_STRING)
        phy = TaxRank(TestTaxRank.PHY_STRING)
        dom = TaxRank(TestTaxRank.DOM_STRING)

        # test equality
        self.assertEqual(spe, spe)
        self.assertEqual(gen, gen)
        self.assertEqual(fam, fam)
        self.assertEqual(ord, ord)
        self.assertEqual(cls, cls)
        self.assertEqual(phy, phy)
        self.assertEqual(dom, dom)

        # test inequality (all possible, non-redundant combinations)
        self.assertNotEqual(spe, gen)
        self.assertNotEqual(spe, fam)
        self.assertNotEqual(spe, ord)
        self.assertNotEqual(spe, cls)
        self.assertNotEqual(spe, phy)
        self.assertNotEqual(spe, dom)
        self.assertNotEqual(gen, fam)
        self.assertNotEqual(gen, ord)
        self.assertNotEqual(gen, cls)
        self.assertNotEqual(gen, phy)
        self.assertNotEqual(gen, dom)
        self.assertNotEqual(fam, ord)
        self.assertNotEqual(fam, cls)
        self.assertNotEqual(fam, phy)
        self.assertNotEqual(fam, dom)
        self.assertNotEqual(ord, cls)
        self.assertNotEqual(ord, phy)
        self.assertNotEqual(ord, dom)
        self.assertNotEqual(phy, dom)

        # test __lt__ and __gt__ for expected hierarchy
        # species
        self.assertLess(spe, gen)
        self.assertLess(spe, fam)
        self.assertLess(spe, ord)
        self.assertLess(spe, cls)
        self.assertLess(spe, phy)
        self.assertLess(spe, dom)
        # genus
        self.assertGreater(gen, spe)
        self.assertLess(gen, fam)
        self.assertLess(gen, ord)
        self.assertLess(gen, cls)
        self.assertLess(gen, phy)
        self.assertLess(gen, dom)
        # family
        self.assertGreater(fam, spe)
        self.assertGreater(fam, gen)
        self.assertLess(fam, ord)
        self.assertLess(fam, cls)
        self.assertLess(fam, phy)
        self.assertLess(fam, dom)
        # order
        self.assertGreater(ord, spe)
        self.assertGreater(ord, gen)
        self.assertGreater(ord, fam)
        self.assertLess(ord, cls)
        self.assertLess(ord, phy)
        self.assertLess(ord, dom)
        # class
        self.assertGreater(cls, spe)
        self.assertGreater(cls, gen)
        self.assertGreater(cls, fam)
        self.assertGreater(cls, ord)
        self.assertLess(cls, phy)
        self.assertLess(cls, dom)
        # phylum
        self.assertGreater(phy, spe)
        self.assertGreater(phy, gen)
        self.assertGreater(phy, fam)
        self.assertGreater(phy, ord)
        self.assertGreater(phy, cls)
        self.assertLess(phy, dom)
        # domain
        self.assertGreater(dom, spe)
        self.assertGreater(dom, gen)
        self.assertGreater(dom, fam)
        self.assertGreater(dom, ord)
        self.assertGreater(dom, cls)
        self.assertGreater(dom, phy)

        # test __ge__ and __le__ for the ceiling and floor, respectively
        allRanks = [spe, gen, fam, ord, cls, phy, dom]
        for rank in allRanks:
            self.assertGreaterEqual(dom, rank)
            self.assertLessEqual(spe, rank)
    
    def testE_superkingdomIsDomain(self) -> None:
        """ Tests that 'superkingdom' behaves appropriately as input
        """
        domain = TaxRank(TestTaxRank.DOM_STRING)
        superkingdom = TaxRank(TestTaxRank.SUP_STRING)
        self.assertEqual(domain, superkingdom)

    def testF_getRankAboveBelow(self) -> None:
        """ Tests getRankAbove and getRankBelow for expepcted behavior
        """
        # make all possible ranks
        SPE = TaxRank(TestTaxRank.SPE_STRING)
        GEN = TaxRank(TestTaxRank.GEN_STRING)
        FAM = TaxRank(TestTaxRank.FAM_STRING)
        ORD = TaxRank(TestTaxRank.ORD_STRING)
        CLS = TaxRank(TestTaxRank.CLS_STRING)
        PHY = TaxRank(TestTaxRank.PHY_STRING)
        DOM = TaxRank(TestTaxRank.DOM_STRING)

        # species to genus
        newRank = SPE.getRankAbove()
        self.assertGreater(newRank, SPE)
        self.assertEqual(newRank, GEN)
        self.assertEqual(SPE, TaxRank(TestTaxRank.SPE_STRING))

        # genus to family
        newRank = newRank.getRankAbove()
        self.assertGreater(newRank, GEN)
        self.assertEqual(newRank, FAM)

        # family to order
        newRank = newRank.getRankAbove()
        self.assertGreater(newRank, FAM)
        self.assertEqual(newRank, ORD)

        # order to class
        newRank = newRank.getRankAbove()
        self.assertGreater(newRank, ORD)
        self.assertEqual(newRank, CLS)

        # class to phylum
        newRank = newRank.getRankAbove()
        self.assertGreater(newRank, CLS)
        self.assertEqual(newRank, PHY)

        # phylum to domain
        newRank = newRank.getRankAbove()
        self.assertGreater(newRank, PHY)
        self.assertEqual(newRank, DOM)

        # domain to phylum
        newRank = DOM.getRankBelow()
        self.assertLess(newRank, DOM)
        self.assertEqual(newRank, PHY)
        self.assertEqual(DOM, TaxRank(TestTaxRank.DOM_STRING))

        # phylum to class
        newRank = newRank.getRankBelow()
        self.assertLess(newRank, PHY)
        self.assertEqual(newRank, CLS)

        # class to order
        newRank = newRank.getRankBelow()
        self.assertLess(newRank, CLS)
        self.assertEqual(newRank, ORD)

        # order to family
        newRank = newRank.getRankBelow()
        self.assertLess(newRank, ORD)
        self.assertEqual(newRank, FAM)

        # family to genus
        newRank = newRank.getRankBelow()
        self.assertLess(newRank, FAM)
        self.assertEqual(newRank, GEN)

        # genus to species
        newRank = newRank.getRankBelow()
        self.assertLess(newRank, GEN)
        self.assertEqual(newRank, SPE)

    def testG_addition(self) -> None:
        """ Tests that adding an int to a TaxRank  works as expected
        """
        SPE = TaxRank(TestTaxRank.SPE_STRING)
        GEN = TaxRank(TestTaxRank.GEN_STRING)
        FAM = TaxRank(TestTaxRank.FAM_STRING)
        ORD = TaxRank(TestTaxRank.ORD_STRING)
        CLS = TaxRank(TestTaxRank.CLS_STRING)
        PHY = TaxRank(TestTaxRank.PHY_STRING)
        DOM = TaxRank(TestTaxRank.DOM_STRING)

        self.assertEqual((SPE + 0), SPE)
        self.assertEqual((SPE + 1), GEN)
        self.assertEqual((SPE + 2), FAM)
        self.assertEqual((SPE + 3), ORD)
        self.assertEqual((SPE + 4), CLS)
        self.assertEqual((SPE + 5), PHY)
        self.assertEqual((SPE + 6), DOM)
        self.assertEqual((SPE + 7), None)
        self.assertEqual((SPE + -1), None)
    
    def testH_subtraction1(self):
        """ Tests that subtracting an int from a TaxRank works as expected
        """
        # make all ranks
        SPE = TaxRank(TestTaxRank.SPE_STRING)
        GEN = TaxRank(TestTaxRank.GEN_STRING)
        FAM = TaxRank(TestTaxRank.FAM_STRING)
        ORD = TaxRank(TestTaxRank.ORD_STRING)
        CLS = TaxRank(TestTaxRank.CLS_STRING)
        PHY = TaxRank(TestTaxRank.PHY_STRING)
        DOM = TaxRank(TestTaxRank.DOM_STRING)
    
        # check subtraction
        self.assertEqual((DOM - 0), DOM)
        self.assertEqual((DOM - 1), PHY)
        self.assertEqual((DOM - 2), CLS)
        self.assertEqual((DOM - 3), ORD)
        self.assertEqual((DOM - 4), FAM)
        self.assertEqual((DOM - 5), GEN)
        self.assertEqual((DOM - 6), SPE)
        self.assertEqual((DOM - 7), None)
        self.assertEqual((DOM - -1), None)

    def testI_subtraction2(self):
        """ Tests that subtracting a TaxRank from a TaxRank works as expected
        """
        # make all ranks
        SPE = TaxRank(TestTaxRank.SPE_STRING)
        GEN = TaxRank(TestTaxRank.GEN_STRING)
        FAM = TaxRank(TestTaxRank.FAM_STRING)
        ORD = TaxRank(TestTaxRank.ORD_STRING)
        CLS = TaxRank(TestTaxRank.CLS_STRING)
        PHY = TaxRank(TestTaxRank.PHY_STRING)
        DOM = TaxRank(TestTaxRank.DOM_STRING)

        # subtract forward (larger minus smaller)
        self.assertEqual((DOM - DOM), 0)
        self.assertEqual((DOM - PHY), 1)
        self.assertEqual((DOM - CLS), 2)
        self.assertEqual((DOM - ORD), 3)
        self.assertEqual((DOM - FAM), 4)
        self.assertEqual((DOM - GEN), 5)
        self.assertEqual((DOM - SPE), 6)

        # subtract reverse (smaller minus larger)
        self.assertEqual((SPE - DOM), -6)
        self.assertEqual((GEN - DOM), -5)
        self.assertEqual((FAM - DOM), -4)
        self.assertEqual((ORD - DOM), -3)
        self.assertEqual((CLS - DOM), -2)
        self.assertEqual((PHY - DOM), -1)
        self.assertEqual((DOM - DOM), 0)

    def testJ_increment(self):
        """ Tests for the expected functionality of increment()
        """
        # make all ranks
        SPE = TaxRank(TestTaxRank.SPE_STRING)
        GEN = TaxRank(TestTaxRank.GEN_STRING)
        FAM = TaxRank(TestTaxRank.FAM_STRING)
        ORD = TaxRank(TestTaxRank.ORD_STRING)
        CLS = TaxRank(TestTaxRank.CLS_STRING)
        PHY = TaxRank(TestTaxRank.PHY_STRING)
        DOM = TaxRank(TestTaxRank.DOM_STRING)

        # increment should modify calling object
        rank = copy.deepcopy(SPE)
        self.assertEqual(rank, SPE)

        rank.increment()
        self.assertEqual(rank, GEN)

        rank.increment()
        self.assertEqual(rank, FAM)

        rank.increment()
        self.assertEqual(rank, ORD)

        rank.increment()
        self.assertEqual(rank, CLS)

        rank.increment()
        self.assertEqual(rank, PHY)

        rank.increment()
        self.assertEqual(rank, DOM)

        # increment should raise an error when trying to go above the ceiling
        try: 
            rank.increment()
            failed = False
        except:
            failed = True
        self.assertTrue(failed)

    def testK_decrement(self):
        """ Tests for the expected functionality of increment()
        """
        # make all ranks
        SPE = TaxRank(TestTaxRank.SPE_STRING)
        GEN = TaxRank(TestTaxRank.GEN_STRING)
        FAM = TaxRank(TestTaxRank.FAM_STRING)
        ORD = TaxRank(TestTaxRank.ORD_STRING)
        CLS = TaxRank(TestTaxRank.CLS_STRING)
        PHY = TaxRank(TestTaxRank.PHY_STRING)
        DOM = TaxRank(TestTaxRank.DOM_STRING)

        # decrement should modify calling object
        rank = copy.deepcopy(DOM)
        self.assertEqual(rank, DOM)

        rank.decrement()
        self.assertEqual(rank, PHY)

        rank.decrement()
        self.assertEqual(rank, CLS)

        rank.decrement()
        self.assertEqual(rank, ORD)

        rank.decrement()
        self.assertEqual(rank, FAM)

        rank.decrement()
        self.assertEqual(rank, GEN)

        rank.decrement()
        self.assertEqual(rank, SPE)

        # decrement should raise an error when trying to go below the floor
        try: 
            rank.decrement()
            failed = False
        except:
            failed = True
        self.assertTrue(failed)

    def testL_plural(self):
        """ Tests that plural() produces the expected strings
        """
        # make all ranks
        SPE = TaxRank(TestTaxRank.SPE_STRING)
        GEN = TaxRank(TestTaxRank.GEN_STRING)
        FAM = TaxRank(TestTaxRank.FAM_STRING)
        ORD = TaxRank(TestTaxRank.ORD_STRING)
        CLS = TaxRank(TestTaxRank.CLS_STRING)
        PHY = TaxRank(TestTaxRank.PHY_STRING)
        DOM = TaxRank(TestTaxRank.DOM_STRING)

        # check that each rank produces the expected pluralized string
        self.assertEqual(SPE.plural(), TestTaxRank.SPE_PLURAL)
        self.assertEqual(GEN.plural(), TestTaxRank.GEN_PLURAL)
        self.assertEqual(FAM.plural(), TestTaxRank.FAM_PLURAL)
        self.assertEqual(ORD.plural(), TestTaxRank.ORD_PLURAL)
        self.assertEqual(CLS.plural(), TestTaxRank.CLS_PLURAL)
        self.assertEqual(PHY.plural(), TestTaxRank.PHY_PLURAL)
        self.assertEqual(DOM.plural(), TestTaxRank.DOM_PLURAL)



