"""
Author - K H M Jayasinghe
Date - 21.02.2021

"""



class Sequence:


    def fasta_split(self, filename):
        """Give name of the sequence and sequence containing python dictionary

        input- FASTA file containing multiple DNA,RNA,or amino acid sequences
        output- Dictionary containing sequence as value and name of the sequence as key
        """

        self.seq_Dict = {}
        newKey = ""

        with open(filename, 'r') as sequence:  # Read sequence store in the variable
            for line in sequence:

                line = line.strip('\n')
                if '>' in line:  # differentiate keys of the dictionary using fasta head
                    newKey = line
                    self.seq_Dict[newKey] = ''
                    continue
                # add values to the key
                self.seq_Dict[newKey] = self.seq_Dict[newKey] + line

            return self.seq_Dict



    def get_Seq_Type(self):
        """Give name of the sequence and type of the sequence containing python dictionary

        input- FASTA file containing multiple DNA,RNA,or amino acid sequences
        output- Dictionary containing, sequence type as value and name of the sequence as key
        """

        self.sequence_Identifing_dict = {}  # define new key to store sequence and sequence type

        # use dictionary of fasta split method to get type of each sequence
        for key, value in self.seq_Dict.items():

            new_key = key  # use current key to assign value
            if 'M' in value:
                self.sequence_Identifing_dict[new_key] = 'amino acid sequence'

            elif 'U' in value:
                self.sequence_Identifing_dict[new_key] = 'RNA sequence'

            else:
                self.sequence_Identifing_dict[new_key] = 'DNA sequence'

        return self.sequence_Identifing_dict



    def get_global_pairwise_similarity_score(self):
        dict = self.seq_Dict
        from itertools import combinations
        from Bio import Align

        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'

        keys = tuple(dict)
        combs = list(combinations(keys, 2)) # create combination with all possible pairs
        #print(combs)
        score_result = ' '
        for combination in combs:
            seq1 = dict[combination[0]]     # get each pair to get pairwise similarity score
            seq2 = dict[combination[1]]

            score = aligner.score(seq1, seq2)
            score_result = score_result+ (str('pairwise similarity between ' + str(combination[0]) + ' and ' + str(combination[1]) + ' is ' + str(score)+'\n'))

        #print(score_result)

        save_result = open("global_pairwise_similarity_score.txt", "w")
        save_result.write(score_result)
        save_result.close()



    def pairwise_dot_matrix(self,seq_one, seq_two):
        import seaborn as sns
        import matplotlib.pyplot as plt
        sns.set_theme(style="whitegrid")

        seq_oneAr = []     # define arrays for two sequences
        seq_twoAr = []
        for na in seq_one:
            seq_oneAr.append(na)
        for nac in seq_two:             # get each single letter for x axis and y axis for draw the plot
            seq_twoAr.append(nac)
        sns.scatterplot(x=seq_oneAr, y=seq_twoAr)
        plt.show()

    def sequence_finder(self,accession_number):

        """Find sequence from NCBI using accession number and save it in FASTA file

         input- accession number of particular sequence
         output- sequence containing FASTA file
         """

        from Bio import Entrez
        from Bio import SeqIO
        Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
        # Get user given accession number directly into id
        handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")

        output = '<' + str(record.description) + '\n' + str(record.seq)  # Build FASTA file
        save_result = open("result_of_sequence_finder.fasta", "w")
        save_result.write(str(output))      # save output as string
        save_result.close()



    def blast(self):

        """Do appropriate BLAST for given DNA ,RNA or amino acid sequence

         input- FASTA file containing multiple DNA,RNA,or amino acid sequences
         output- BLAST result containing xml file
         """

        from Bio.Blast import NCBIWWW


        self.get_Seq_Type()
        self.result = ' '

        for key, value in self.sequence_Identifing_dict.items():    # Get previous seq type containing dict to identify seq
            if (value == 'RNA sequence') or (value == 'DNA sequence'):
                blast_type = 'blastn'
                blast_db = 'nt'
                                                                        # change blast input variable according to the sequences
            elif (value == 'amino acid sequence'):
                blast_type = 'blastp'
                blast_db = 'pr'

            result_handle = NCBIWWW.qblast(blast_type, blast_db, str(self.seq_Dict[key]))
            self.result = self.result + result_handle.read()
        save_result = open("result_of_blast.xml", "w")
        save_result.write(self.result)



    def Calculate_melting_temperature(self):
        """Calculate the melting temperature of a DNA and RNA from given fasta sequences

        input- FASTA file containing multiple DNA,RNA,or amino acid sequences
        output- melting temperature of DNA and RNA sequences
        """

        from Bio.SeqUtils import MeltingTemp as mt
        from Bio.Seq import Seq

        for key, value in self.sequence_Identifing_dict.items():    # Get previous seq type containing dict to identify seq

            if (value == 'RNA sequence') or (value == 'DNA sequence'): # calculate Tm only for DNA & RNA sequences
                self.melting_temperature = ('%0.2f' % mt.Tm_Wallace(str(self.seq_Dict[key])))   # use sequence containing dict to get particular seq according to the key
                print(key + ' melting temperature is ' + self.melting_temperature)




    def calculate_isoelectricPoint(self):
        """Calculate the isoelectric point.
        Uses the module IsoelectricPoint to calculate the pI of a protein.

        input- FASTA file containing multiple DNA,RNA,or amino acid sequences
        output- isoelectric point and it's charge at pH 7, only for amino acid sequences
        """

        from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP

        for key, value in self.sequence_Identifing_dict.items(): # Get previous seq type containing dict to identify seq

            if value == 'amino acid sequence': # filter out amino acid sequences
                protein = IP(str(self.seq_Dict[key])) # use sequence containing dict to get particular seq according to the key

                print("IEP of peptide {} is {:.2f}"
                      .format(protein.sequence, protein.pi()))

                print("It's charge at pH 7 is {:.2f}"
                      .format(protein.charge_at_pH(7.0)))



    def aromaticity(self):

        """Calculate the aromaticity according to Lobry, 1994.
        Calculates the aromaticity value of a protein according to Lobry, 1994.
        It is simply the relative frequency of Phe+Trp+Tyr.
        """
        aromatic_aas = ['Y', 'W', 'F']  # define array for aromatic aa
        self.aromatic_aa_count = 0 # define variable for count aromatic aa
        self.length_of_aa_seq = 0 # define variable for count length of aa seq

        for key, value in self.sequence_Identifing_dict.items():       # Get previous seq type containing dict to identify seq
            if (value == 'amino acid sequence'):        # filter out amino acid sequences
                for aa in self.seq_Dict[key]:       # use sequence containing dict to get particular seq according to the key

                    if aa in aromatic_aas:
                        self.aromatic_aa_count += 1
                self.length_of_aa_seq = len(self.seq_Dict[key])

                print('aromaticity of the ' + key + ' is ', self.aromatic_aa_count / self.length_of_aa_seq * 100)   # calculate aromaticity



s1 = Sequence()
print(s1.fasta_split('seqList.fasta'))
# print(s1.get_Seq_Type())
# print(s1.get_global_pairwise_similarity_score())
# print(s1.pairwise_dot_matrix('ATGrr','AGCrr'))
# print(s1.sequence_finder("EU490707"))
print(s1.blast())
# print(s1.Calculate_melting_temperature())
# print(s1.aromaticity())
# print(s1.Calculate_isoelectricPoint())
