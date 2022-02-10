import unittest
from PHANTASM.taxonomy.processLpsn import makeLpsnD
from PHANTASM.taxonomy.TaxRank import TaxRank
from param import CSV_1, CSV_2, CSV_3, CSV_4

class TestLpsnData(unittest.TestCase):
    # constants
    RANK_KEYS = {'species', 'genus', 'family', 'order', 'class', 'phylum'}
    SYN_KEYS = {'species_syn', 'genus_syn', 'family_syn', 'order_syn', 'class_syn', 'phylum_syn'}
    ENTRY_KEYS = {'type', 'status', 'parent'}

    def setUp(self):
        self.lpsnD = makeLpsnD(CSV_1, CSV_2, CSV_3, CSV_4)


    def testA_expectedKeys(self):
        """ Tests if the top-level keys (RANK_KEYS + SYN_KEYS) are present in
            lpsnD.keys()
        """
        # check that the top-level keys match the expected keys
        expected = set()
        expected.update(TestLpsnData.RANK_KEYS)
        expected.update(TestLpsnData.SYN_KEYS)
        observed = set(self.lpsnD.keys())
        self.assertEqual(expected, observed)


    def testB_synonymsFormattedCorrectly(self):
        """ Tests for expected formatting in the synonyms dictionary
        """
        # for each synonym dictionary
        for sk in TestLpsnData.SYN_KEYS:
            # for each synonym
            for synName in self.lpsnD[sk].keys():
                # affirm synName is a string and its value is a string
                self.assertIsInstance(synName, str)
                self.assertIsInstance(self.lpsnD[sk][synName], str)


    def testC_ranksFormattedCorrectly(self):
        """ Tests for expected formatting in the ranks dictionary
        """
        # constants
        TYPE   = 'type'
        STATUS = 'status'
        PARENT = 'parent'

        # for each rank dictionary
        for rk in TestLpsnData.RANK_KEYS:
            # for each entry in the rank
            for sciName in self.lpsnD[rk].keys():
                # affirm sciName is a string and the entry's length is 3
                self.assertIsInstance(sciName, str)

                entry:dict = self.lpsnD[rk][sciName]

                self.assertEqual(len(entry),len(TestLpsnData.ENTRY_KEYS))
                
                # affirm the keys are sound
                self.assertEqual(entry.keys(), TestLpsnData.ENTRY_KEYS)
                
                self.assertIsInstance(entry[STATUS], str)
                self.assertIsInstance(entry[PARENT], str)

                # species type material must be a list            
                if rk == 'species':
                    self.assertIsInstance(entry[TYPE], list)
                
                # genus type material must be a string
                elif rk == 'genus':
                    self.assertIsInstance(entry[TYPE], str)
                
                # ranks greater than genus are allowed to have 'None' type material
                else:
                    if entry[TYPE] is not None:
                        self.assertIsInstance(entry[TYPE], str)


    def testD_synonymsAbsentFromRankDictionaries(self):
        """ Looks for keys that are shared between the synonym and rank dictio-
            naries (there shouldn't be any shared keys).
        """
        # for each synonym dictionary
        for sk in TestLpsnData.SYN_KEYS:
            rk = sk[:-4]  # drop '_syn' to get rank key
            
            # make sure all synonyms are missing from its respective rank dictionary
            for synName in self.lpsnD[sk].keys():
               self.assertNotIn(synName, self.lpsnD[rk].keys())


    def testE_synonymsPreferredNamesPresentInRankDictionaries(self):
        # for each synonym dictionary
        for sk in TestLpsnData.SYN_KEYS:
            rk = sk[:-4]  # drop '_syn' to get rank key
            
            # make sure all preferred names are present in the rank dictionary
            for synName in self.lpsnD[sk].keys():
                sciName = self.lpsnD[sk][synName]
                self.assertIn(sciName, self.lpsnD[rk].keys())


    def testF_parentsExistInHigherRankDict(self):
        # constants
        MISSING_PAR = ""

        # for ranks of class and below (domains have no LPSN dictionary)
        for rk in TestLpsnData.RANK_KEYS:
            # use the TaxRank class to get the parental rank
            parentRank = TaxRank(rk)
            parentRank.increment()

            # if the parent is a domain, then there are only two valid names
            if parentRank == TaxRank('domain'):
                parentKeys = ['Bacteria', 'Archaea']
            
            # get all of the possible parent names
            else:
                parentRank = str(parentRank)
                parentKeys = self.lpsnD[parentRank].keys()

            # for each entry
            for sciName in self.lpsnD[rk].keys():
                # get the parent name
                parentName = self.lpsnD[rk][sciName]['parent']

                # make sure there are no missing parents
                errorMsg = "no parent for " + sciName
                self.assertNotEqual(parentName, MISSING_PAR, errorMsg)

                # assert that the parent could be looked up
                errorMsg = "parent of " + sciName + " (" + parentName + \
                                                     ") could not be looked up"
                self.assertIn(parentName, parentKeys, errorMsg)

# this test is important for my homebrewed dataset. modify when API goes live
# only check phylum-class-order-familiy-gen relationships; all others are based on LPSN data
    def testG_typeExistsInLowerRankDict(self):
        ranks = ['phylum', 'class', 'order', 'family', 'genus']
        for i in range(len(ranks)-1):
            parentD:dict = self.lpsnD[ranks[i]]
            childD:dict  = self.lpsnD[ranks[i+1]]

        # for each entry in the parent dictionary
            for parentName in parentD.keys():
                # confirm that type is None or it exists in the child dictionary
                typeMat = parentD[parentName]['type']
                if typeMat is not None:
                    self.assertIn(typeMat, childD.keys())

# this test is important for my homebrewed dataset. modify when API goes live
# only check phylum-class-order-familiy-gen relationships; all others are based on LPSN data
    def testH_typeMaterialListsExpectedParent(self):
        ranks = ['phylum', 'class', 'order', 'family', 'genus']
        for i in range(len(ranks)-1):
            parentD:dict = self.lpsnD[ranks[i]]
            childD:dict  = self.lpsnD[ranks[i+1]]

            # for each entry in parent dictionary
            for parentName in parentD.keys():
                # confirm that type is None or that its parent matches
                typeMat = parentD[parentName]['type']
                if typeMat is not None:
                    self.assertEqual(parentName, childD[typeMat]['parent'])
                

    def testI_speciesTypeMaterial(self):
        # for each species
        speD = self.lpsnD['species']
        for species in speD.keys():
            typeMat = speD[species]['type']
            self.assertIsInstance(typeMat, list)

            for strain in typeMat:
                self.assertIsInstance(strain, str)


