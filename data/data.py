from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio.Align import AlignInfo
import sys
import math
import os
import time

IUPAC_REV = {('A',): 'A',
             ('A', 'C'): 'M',
             ('A', 'C', 'G'): 'V',
             ('A', 'C', 'T'): 'H',
             ('A', 'G'): 'R',
             ('A', 'G', 'T'): 'D',
             ('A', 'C', 'G', 'T'): 'N',
             ('A', 'T'): 'W',
             ('C',): 'C',
             ('C', 'G'): 'S',
             ('C', 'G', 'T'): 'B',
             ('C', 'T'): 'Y',
             ('G',): 'G',
             ('G', 'T'): 'K',
             ('T',): 'T',}


class Data:
    def __init__(self, fasta_file_path, muscle_path='./third_party/muscle.exe', debug=True):
        self.debug = debug
        self.unaligned_sequences = self.read_fasta(fasta_file_path)
        self.muscle_aligned_sequences = self.muscle_align(fasta_file_path, muscle_path)
        self.summary_align = self.get_summary_align(self.out_aligned_file_path)
        self.align_consensus = self.get_consensus(self.summary_align)
        self.write_consensus_to_file(fasta_file_path)
        self.total_multiple_alignment_time = self.alignment_time + self.consensus_from_alignment_time

    def get_summary_align(self, path):
        alignment = AlignIO.read(path, 'fasta')
        summary_align = AlignInfo.SummaryInfo(alignment)
        return summary_align

    def read_fasta(self, path):
        self.data_log(f'Reading fasta file at {path}')
        in_file = open(path, 'r')

        data=''
        sequences = {}
        name_list=[]
        seq_list=[]

        for line in in_file:
            line=line.strip()
            for i in line:
                if i=='>':
                    tokens = line.split()
                    name = tokens[0][1:]
                    name_list.append(name)
                    if data:
                        seq_list.append(data)
                        data=''
                    break
                else:
                    line=line.upper()
            if all([k==k.upper() for k in line]):
                data=data+line

        if data:
            seq_list.append(data)

        for i, name in enumerate(name_list):
            sequences[name] = seq_list[i]

        return sequences

    def muscle_align(self, fasta_file_path, muscle_path):
        path, filename = os.path.split(fasta_file_path)

        self.out_aligned_file_path = os.path.join(path, filename.split('.')[0] + '-aligned.fasta')
        self.data_log('Starting to align with muscle.')
        start = time.time()

        muscle_cline = MuscleCommandline(muscle_path, input=fasta_file_path, out=self.out_aligned_file_path)
        _ = muscle_cline()

        end = time.time()

        self.alignment_time = end - start
        self.data_log(f'Writing muscle aligned sequences to {self.out_aligned_file_path}.')
        # Get consensus sequence and put in separate file.
        alignment = self.read_fasta(self.out_aligned_file_path)
        return alignment

    def get_consensus(self, summary_align):
        self.data_log('Obtaining consensus sequence from alignment...')
        start = time.time()
        tabulated = self._tabulate(summary_align)
        consensus = ''.join([self._get_consensus_base(d) for d in tabulated]).upper()
        end = time.time()
        self.consensus_from_alignment_time = end - start
        return consensus
    
    def write_consensus_to_file(self, fasta_file_path):
        path, filename = os.path.split(fasta_file_path)
        out_consensus_path = os.path.join(path, filename.split('.')[0] + '-consensus.txt')
        self.data_log(f'Writing consensus to {out_consensus_path}.')
        with open(out_consensus_path, 'w+') as f:
            f.write(f'>consensus\n{self.align_consensus}\n')

    def _get_consensus_base(self, tabdict, gap='-', errorchar='X', use_ambi=True):
        """Given a dictionary representing character frequencies
        at a single position, returns the most common char at
        that position subject to the rules below.
        plu            plurality for calling a consensus character
        use_ambi - uses IUPAC ambiguity codes if possible in place of errorchar
        """

        # remove any base that has a zero count.
        # convert fractional frequencies to log2 scale
        tabdict = dict([(k,math.log(v,2)) for k,v in tabdict.items() if v > 0])

        # if our frequency table no remaining entries, error
        if len(tabdict) == 0:
            return errorchar

        # if our frequency table has only one entry, we are done.
        if len(tabdict) == 1:
            return list(tabdict.keys())[0]

        # order nucleotides by increasing frequency
        sortedkeys = sorted(tabdict.keys(), key=lambda k: tabdict[k], reverse=True)

        # Collect the nucleotides that constitute a majority of samples at a site.
        # We want to emit ambiguity codes that summarize a position
        # while avoiding truly low-frequency bases.  An alignment motif would convey
        # more information, and perhaps that's what one should use, but a consensus
        # sequence must express each site in only a single letter so some loss of
        # information is inevitable.
        #
        # Our cheesy algorithm is to include all bases whose frequency is
        # at least half that of the next highest frequency nucleotide.
        nuc = [sortedkeys[0]]

        for i,v in enumerate(sortedkeys):
            if i == 0:
                continue
                
            # calculate log fold change between successive frequencies.
            logf = tabdict[sortedkeys[i - 1]]- tabdict[v]

            if logf <= 1:
                nuc.append(v)
            else:
                # otherwise stop looking.
                break
        
        if len(nuc) > 1:
            if use_ambi:
                return IUPAC_REV.get(tuple(sorted(nuc)), errorchar)
            else:
                return errorchar
        else:
            return sortedkeys[0]

    def _tabulate(self, seqList, start = 0, end = None):
        """List of dictionaries, each holding nucleotide frequences for a single site in multiple alignment."""
        if end is None: 
            end = len(seqList.alignment._records[0].seq) 
    
        if start < 0 or end > len(seqList.alignment._records[0].seq): 
            raise ValueError("Start (%s) and end (%s) are not in the range %s to %s" 
                            % (start, end, 0, len(seqList.alignment._records[0].seq))) 

        all_letters = seqList._get_all_letters()
        chars_to_ignore = []
        dictList = []
        for residue_num in range(start, end): 
            dictList.append(seqList._get_letter_freqs(residue_num, 
                                                seqList.alignment._records, 
                                                all_letters, chars_to_ignore))

        return dictList
    
    def data_log(self, log):
        if(self.debug):
            print(f'[ DATA ]: {log}')