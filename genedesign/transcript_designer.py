import random

from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript

from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker
from genedesign.checkers.codon_checker import CodonChecker

class TranscriptDesigner:
    """
    Previous: Reverse translates a protein sequence into a DNA sequence and chooses an RBS using the highest CAI codon for each amino acid.
    Improved version: 
    - Uses multiple synonymous codons (guided random)
    - Rejects candidates that fail cehckers (forbidden, hairpin, promoter, internal RBS)
    """

    def __init__(self):
        self.aminoAcidToCodon = {}
        self.rbsChooser = None
        self.checkers = []
        self.codonChecker = None

    def initiate(self) -> None:
        """
        Initializes the codon table and the RBS chooser.
        """
        self.rbsChooser = RBSChooser()
        self.rbsChooser.initiate()

        # Initialize all checkers
        forbidden = ForbiddenSequenceChecker()
        forbidden.initiate()

        promoter = PromoterChecker()
        promoter.initiate()

        internal_rbs = InternalRBSChecker()
        internal_rbs.initiate()

        self.codonChecker = CodonChecker()
        self.codonChecker.initiate()

        self.checkers = [forbidden, promoter, internal_rbs]        

        # Multiple codon options (simple "guided random" approach)
        self.aminoAcidToCodon = {
            'A': ["GCG", "GCC", "GCT"],
            'C': ["TGC", "TGT"],
            'D': ["GAT", "GAC"],
            'E': ["GAA", "GAG"],
            'F': ["TTC", "TTT"],
            'G': ["GGT", "GGC", "GGA"],
            'H': ["CAC", "CAT"],
            'I': ["ATC", "ATT", "ATA"],
            'K': ["AAA", "AAG"],
            'L': ["CTG", "CTT", "CTC", "TTA"],
            'M': ["ATG"],
            'N': ["AAC", "AAT"],
            'P': ["CCG", "CCA", "CCC", "CCT"],
            'Q': ["CAG", "CAA"],
            'R': ["CGT", "CGC", "AGA"],
            'S': ["TCT", "TCC", "AGC", "AGT"],
            'T': ["ACC", "ACG", "ACT", "ACA"],
            'V': ["GTT", "GTG", "GTC", "GTA"],
            'W': ["TGG"],
            'Y': ["TAC", "TAT"]
        }

    def _passes_checks(self, dna: str, codons: list[str]) -> bool:
        passed_hairpin, _ = hairpin_checker(dna)
        if not passed_hairpin:
            return False

        for checker in self.checkers:
            ok, _ = checker.run(dna)
            if not ok:
                return False

        codons_above_board, _, _, _ = self.codonChecker.run(codons)
        if not codons_above_board:
            return False

        return True

    def _choose_codon(self, aa: str, used_codons: dict[str, int]) -> str:
        options = self.aminoAcidToCodon[aa]

        # sort codons by how often they were used
        ranked = sorted(options, key=lambda c: used_codons.get(c, 0))

        # prefer the most common codon but allow alternatives
        if random.random() < 0.7:
            return ranked[0]
        else:
            return random.choice(ranked[:min(2, len(ranked))])

    def run(self, peptide: str, ignores: set) -> Transcript:
        best_codons = None

        """
        Try multiple random designs; accept first that passes all checks. 
        If none pass after a reasonable number of attempts, 
        fallback to a design that only uses the most common codon for each amino acid (previous approach).
        """

        for _ in range(200):  # Try up to 200 random designs
            # Translate peptide to codons
            used_codons = {}
            codons = []

            for aa in peptide:
                codon = self._choose_codon(aa, used_codons)
                codons.append(codon)
                used_codons[codon] = used_codons.get(codon, 0) + 1

            # Append the stop codon (TAA in this case)
            codons.append("TAA")

            # Build the CDS from the codons
            cds = ''.join(codons)

            selectedRBS = self.rbsChooser.run(cds, ignores)
            transcript_dna = selectedRBS.utr.upper() + cds

            if self._passes_checks(transcript_dna, codons):
                best_codons = codons
                break

        #Fallback if no design passed checks
        if best_codons is None:
            best_codons = [self.aminoAcidToCodon[aa][0] for aa in peptide]  # Use the most common codon for each amino acid
            best_codons.append("TAA")  # Append stop codon

        cds = ''.join(best_codons)
        # Choose an RBS
        selectedRBS = self.rbsChooser.run(cds, ignores)

        # Return the Transcript object
        return Transcript(selectedRBS, peptide, best_codons)

if __name__ == "__main__":
    # Example usage of TranscriptDesigner
    peptide = "MYPFIRTARMTV"
    
    designer = TranscriptDesigner()
    designer.initiate()

    ignores = set()
    transcript = designer.run(peptide, ignores)
    
    # Print out the transcript information
    print(transcript)
