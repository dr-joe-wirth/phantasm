from __future__ import annotations

class TaxRank:
    """ TaxRank:
            This class is designed to store the taxonomic rank of a Taxonomy 
            object. It facilitates the direct comparison of two taxonomic ranks 
            to one another. The constructor accepts a string indicating a taxo-
            nomic rank. The strings allowed as input to the TaxRank constructor
            are as follows (listed from lowest to highest):
                'species'
                'genus'
                'family'
                'order'
                'class'
                'phylum'
                'domain'
    """
    # constants used to determine how different ranks relate to one another, 
    # and how each rank's string changes when pluralized.
    __RANK_TO_NUM:dict   = { 'species':      0,
                             'genus':        1,
                             'family':       2,
                             'order':        3,
                             'class':        4,
                             'phylum':       5,
                             'domain':       6,
                             'superkingdom': 6 }  # easier interface with NCBI
    
    __NUM_TO_RANK:dict   = { 0: 'species',
                             1: 'genus',
                             2: 'family',
                             3: 'order',
                             4: 'class',
                             5: 'phylum',
                             6: 'domain' }
    
    __NUM_TO_PLURAL:dict = { 0: 'species',
                             1: 'genera',
                             2: 'families',
                             3: 'orders',
                             4: 'classes',
                             5: 'phyla',
                             6: 'domains' }
    
    # overloads
    def __init__(self, rank:str) -> None:
        """ __init__:
                Accepts a string as input. Constructs a TaxRank object if the 
                string is an allowed taxonomic rank. Does not return
        """
        # convert 'superkingdom' to 'domain' for interfacing with ncbi taxonomy
        # if rank == 'superkingdom':
        #     rank = 'domain'

        # initialize member variables if the input is valid
        if TaxRank.isValidRank(rank):
            self.__rankNum:int = TaxRank.__RANK_TO_NUM[rank]
            self.__rankStr:str = TaxRank.__NUM_TO_RANK[self.__rankNum]
            
        else:
            raise ValueError("Invalid input.")   
    

    def __str__(self) -> str:
        """ __str__:
                Accepts no inputs. Returns a string indicating the taxonomic
                rank.
        """
        return self.__rankStr
    

    def __repr__(self) -> str:
        """ __repr__:
                Accepts no inputs. Returns a string indicating the taxonomic
                rank.
        """
        return str(self)
    
        
    def __gt__(self, other:TaxRank) -> bool:
        """ __gt__:
                Directly compares two TaxRank objects. Returns a boolean indi-
                cating whether or not one rank is greater (higher rank) than 
                another.
        """
        # use the values of the rankNum member variables to compare ranks
        return self.__rankNum > other.__rankNum
    

    def __ge__(self, other:TaxRank) -> bool:
        """ __ge__:
                Directly compares two TaxRank objects. Returns a boolean indi-
                cating whether or not one rank is greater (higher rank) than or
                equal to another.
        """
        # use the values of the rankNum member variables to compare ranks
        return self.__rankNum >= other.__rankNum


    def __lt__(self, other:TaxRank) -> bool:
        """ __lt__:
                Directly compares two TaxRank objects. Returns a boolean indi-
                cating whether or not one rank is less (lower rank) than 
                another.
        """
        # use the values of the rankNum member variables to compare ranks
        return self.__rankNum < other.__rankNum
    

    def __le__(self, other:TaxRank) -> bool:
        """ __le__:
                Directly compares two TaxRank objects. Returns a boolean indi-
                cating whether or not one rank is less (lower rank) than or 
                equal to another.
        """
        # use the values of the rankNum member variables to compare ranks
        return self.__rankNum <= other.__rankNum


    def __eq__(self, other:TaxRank) -> bool:
        """ __eq__:
                Directly compares two TaxRank objects. Returns a boolean indi-
                cating whether or not the ranks are equivalent to one another.
        """
        # use the values of the rankNum member variables to compare ranks
        return self.__rankNum == other.__rankNum
    

    def __ne__(self, other:TaxRank) -> bool:
        """ __ne__:
                Directly compares two TaxRank objects. Returns a boolean indi-
                cating whether or not the ranks are not equivalent to one 
                another.
        """
        return not self == other


    def __sub__(self, other:TaxRank) -> TaxRank:
        """ __sub__:
                Allows one rank to be subtracted from another. This is mainly 
                provided as an easy way to determine how many taxonomic ranks
                separate the two ranks being compared. It also allows for an 
                int to be subtracted from a rank in order to obtain the ranks
                below it. Returns an integer if substracting two TaxRank 
                objects. Returns a TaxRank object if subtracting an integer 
                from a TaxRank object.
        """
        # return the difference between the rankNum member variables
        if type(other) == TaxRank:
            return self.__rankNum - other.__rankNum
        
        # return the rank 
        elif type(other) == int:
            newNum = self.__rankNum - other

            # make a new TaxRank at the new level unless it is not possible.
            if newNum in self.__NUM_TO_RANK.keys():
                return TaxRank(self.__NUM_TO_RANK[newNum])
            else:
                return None
        else:
            raise TypeError("Can only subtract an int or " + 
                            "a TaxRank object from a TaxRank object.")


    def __add__(self, number:int) -> TaxRank:
        """ __add__:
                Allows an integer to be added to a rank object. This is mainly
                provided as an easy way to obtain the ranks above the object. 
                Unlike subtraction, there is no reason to provide 'rank + rank'
                functionality. Returns a TaxRank object.
        """
        # perform the addition
        newNum = self.__rankNum + number

        # return a TaxRank object at the new level unless it is not possible.
        if newNum in self.__NUM_TO_RANK.keys():
            return TaxRank(self.__NUM_TO_RANK[newNum])
        else:
            return None


    def __hash__(self) -> int:
        """ __hash__:
                Accepts no inputs. Allows TaxRank objects to be stored in sets
                and/or used as dict.keys(). Returns an int.
        """
        return self.__rankNum



    # member functions
    def isValidRank(rankString:str) -> bool:
        """ isAllowed:
                Accepts a string as input. Returns a boolean indicating whether
                or not the string could be successfully converted to a TaxRank
                object.
        """
        # the string is only valid if it is a key in the RANK_TO_NUM dictionary
        return rankString in TaxRank.__RANK_TO_NUM.keys()

    
    def getRankBelow(self) -> TaxRank:
        """ getRankBelow:
                Accepts no inputs. Returns a TaxRank object that is exactly one
                taxonomic rank below the current taxonomic rank. If the current
                taxonomic rank is 'species', then returns None. DOES NOT MODIFY
                THE CALLING OBJECT!
        """
        # no ranks are below species, so return None
        if self.__rankStr == 'species': return None
        
        # subtract one from the rank and return it
        else: return self - 1
    

    def getRankAbove(self) -> TaxRank:
        """ getRankBelow:
                Accepts no inputs. Returns a TaxRank object that is exactly one
                taxonomic rank above the current taxonomic rank. If the current
                taxonomic rank is 'phylum', then returns None. DOES NOT MODIFY
                THE CALLING OBJECT!
        """
        # no ranks are above domain, so return None
        if self.__rankStr == 'domain': return None

        # add one to the rank and return it
        else: return self + 1


    def increment(self) -> None:
        """ increment:
                Accepts no inputs. Increments the calling object by one taxono-
                mic rank. Does not return.
        """
        # get the rank above
        temp = self.getRankAbove()
        
        # update member variables
        self.__rankStr = temp.__rankStr
        self.__rankNum = temp.__rankNum


    def decrement(self) -> None:
        """ decrement:
                Accepts no inputs. Decrements the calling object by one taxono-
                mic rank. Does not return.
        """
        # get the rank below
        temp = self.getRankBelow()

        # update member variables
        self.__rankStr = temp.__rankStr
        self.__rankNum = temp.__rankNum


    def plural(self) -> str:
        """ plural:
                Accepts no inputs. Returns a string with the pluralized spell-
                ing of the rank.
        """
        # use self.rankNum as a key for the NUM_TO_PLURAL dictionary
        return TaxRank.__NUM_TO_PLURAL[self.__rankNum]

