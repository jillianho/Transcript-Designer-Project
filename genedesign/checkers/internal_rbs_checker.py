class InternalRBSChecker: 
    def __init__(self):
        self.motifs = []
        self.ignore_prefix = 9
    
    def initiate(self):
        self.motifs = [
            "AGGAGG",
            "GGAGG",
            "AAGGAG",
        ]
    
    def run(self, dnaseq):
        
        seq = dnaseq.upper()
        
        for motif in self.motifs:
            found_index = -1

            for i in range(self.ignore_prefix, len(seq) - len(motif) + 1):
                if seq[i:i+len(motif)] == motif:
                    found_index = i
                    break

            if found_index != -1:
                left = max(0, found_index -5)
                right = min(len(seq), found_index + len(motif) + 5)
            
                return False, seq[left:right]
    
        return True, None