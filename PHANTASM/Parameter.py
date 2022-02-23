from __future__ import annotations
import os

class Parameters:
    """ This class was created to facilitate interfacing with xenoGI. xenoGI
        relies on a 'parameter dictionary' as a way to provide a shared set of
        filenames, paths, and other parameters. This class stores all of the 
        values necessary to build the parameter dictionary so long as the user
        provides a few input values.

            Use the 'toDict' method to convert a Parameters object to the xeno-
            GI parameter dictionary.

            Use the 'fromDict' method to convert a xenoGI parameter dictionary
            to a Parameters object.
    """

    # private member variables for automated object creation

    # these paths need to be nested within the working directory
    __genbankFilePath = os.path.join('wgs', '*')
    __fileNameMapFN = 'wgsHumanMap.txt'
    __fastaFilePath = os.path.join('fasta', '*.fa')
    __strainInfoFN = 'strainInfo.txt'
    __geneInfoFN = 'geneInfo.txt'
    __geneOrderFN = 'geneOrder.txt'
    __problemGenbankFN = 'problemGenbankFiles.txt'
    __blastFilePath = os.path.join('blast', '*.out')
    __aabrhFN = 'aabrhHardCore.out'
    __speciesTreeFN = 'speciesTree.nwk'
    __makeSpeciesTreeWorkingDir = 'makeSpeciesTreeWorkDir'
    __concatenatedAlignmentFN = os.path.join(__makeSpeciesTreeWorkingDir, 'aabrhHardCore_concatenated.afa')
    __famToGeneKeyFN = os.path.join(__makeSpeciesTreeWorkingDir, 'aabrhHardCoreFamToGeneKey.txt')
    __phyloMarkersFN = 'putativePhylogeneticMarkers.txt'
    __taxonomyObjectFilePath = '*.tax'
    __phyloMarkerFaaFN = 'phylogeneticMarker.faa'
    __blastpResultFN = 'phylogeneticMarker.blastp'
    __aaiFN = 'aai_matrix.txt'
    __aaiHeatmapFN = 'aai_heatmap.pdf'
    __aniFN = 'ani_matrix.txt'
    __aniWorkDir = 'aniWorkdir'
    __aniBlastFilePath = os.path.join(__aniWorkDir, 'blast', '*.blastn')

    # these objects will never change
    __dnaBasedGeneTrees = False
    __blastFileJoinStr = '_-VS-_'
    __blastCLine = 'blastp -matrix BLOSUM62 -gapopen 11 -gapextend 1 -seg yes -outfmt "6 qseqid sseqid evalue qlen qstart qend slen sstart send pident score" -evalue '
    __evalueThresh = 1e-8
    __alignCoverThresh = 0.65
    __percIdentThresh = 0.35
    __aabrhHardCoreGeneTreeFileStem = 'aabrhHardCoreFam'

    # stuff needed for calcScores (not necessary for phantasm)
    __gapOpen = 12
    __gapExtend = 1
    __matrix = 'parasail.blosum62'
    

    def __init__(self, email:str, workdir:str, blastExecutDirPath:str, \
                 musclePath:str, fastTreePath:str, numProcesses:int=1, \
                                        maxNumTreeLeaves:int=50) -> Parameters:
        """ __init__:
                Accepts an email address, a working directory, a directory con-
                taining blast+ executables, a path to a MUSCLE executable, a
                path to a FastTree executable, the number of processors, and
                the maximum number of leaves for species trees as inputs. Popu-
                lates its member variables using the inputs as several private
                static class variables. Returns the newly constructed object.
        """
        # import the user specified variables; make sure all paths are absolute
        self.email = email
        self.workdir = os.path.abspath(workdir)
        self.blastExecutDirPath = os.path.abspath(blastExecutDirPath)
        self.musclePath = os.path.abspath(musclePath)
        self.fastTreePath = os.path.abspath(fastTreePath)
        self.numProcesses = numProcesses
        self.maxNumTreeLeaves = maxNumTreeLeaves

        # add the workdir to the path of selected static class members
        self.genbankFilePath = os.path.join(workdir, Parameters.__genbankFilePath)
        self.fileNameMapFN = os.path.join(workdir, Parameters.__fileNameMapFN)
        self.fastaFilePath = os.path.join(workdir, Parameters.__fastaFilePath)
        self.strainInfoFN = os.path.join(workdir, Parameters.__strainInfoFN)
        self.geneInfoFN = os.path.join(workdir, Parameters.__geneInfoFN)
        self.geneOrderFN = os.path.join(workdir, Parameters.__geneOrderFN)
        self.problemGenbankFN = os.path.join(workdir, Parameters.__problemGenbankFN)
        self.blastFilePath = os.path.join(workdir, Parameters.__blastFilePath)
        self.aabrhFN = os.path.join(workdir, Parameters.__aabrhFN)
        self.speciesTreeFN = os.path.join(workdir, Parameters.__speciesTreeFN)
        self.makeSpeciesTreeWorkingDir = os.path.join(workdir, Parameters.__makeSpeciesTreeWorkingDir)
        self.concatenatedAlignmentFN = os.path.join(workdir, Parameters.__concatenatedAlignmentFN)
        self.famToGeneKeyFN = os.path.join(workdir, Parameters.__famToGeneKeyFN)
        self.phyloMarkersFN = os.path.join(workdir, Parameters.__phyloMarkersFN)
        self.taxonomyObjectFilePath = os.path.join(workdir, Parameters.__taxonomyObjectFilePath)
        self.phyloMarkerFaaFN = os.path.join(workdir, Parameters.__phyloMarkerFaaFN)
        self.blastpResultFN = os.path.join(workdir, Parameters.__blastpResultFN)
        self.aaiFN = os.path.join(workdir, Parameters.__aaiFN)
        self.aaiHeatmapFN = os.path.join(workdir, Parameters.__aaiHeatmapFN)
        self.aniFN = os.path.join(workdir, Parameters.__aniFN)
        self.aniWorkDir = os.path.join(workdir, Parameters.__aniWorkDir)
        self.aniBlastFilePath = os.path.join(workdir, Parameters.__aniBlastFilePath)

        # import remaining variables from static class members
        self.dnaBasedGeneTrees = Parameters.__dnaBasedGeneTrees
        self.blastFileJoinStr = Parameters.__blastFileJoinStr
        self.blastCLine = Parameters.__blastCLine
        self.evalueThresh = Parameters.__evalueThresh
        self.alignCoverThresh = Parameters.__alignCoverThresh
        self.percIdentThresh = Parameters.__percIdentThresh
        self.aabrhHardCoreGeneTreeFileStem = Parameters.__aabrhHardCoreGeneTreeFileStem

        # stuff required for making calcscores work
        # self.gapOpen = Parameters.__gapOpen
        # self.gapExtend = Parameters.__gapExtend
        # self.matrix = Parameters.__matrix


    def __repr__(self) -> str:
        # returns the string representation
        return str(self)


    def __str__(self) -> str:
        # returns a string of the object in dict format
        return str(self.toDict())
    

    def toDict(self) -> dict:
        # converts a Parameters object to the paramD dictionary format
        return vars(self)
    
    
    def fromDict(paramD:dict) -> Parameters:
        # converts a paramD dictionary to a Parameters object
        return Parameters(paramD['email'],
                          paramD['workdir'],
                          paramD['blastExecutDirPath'], 
                          paramD['musclePath'],
                          paramD['fastTreePath'],
                          paramD['numProcesses'],
                          paramD['maxNumTreeLeaves'])


