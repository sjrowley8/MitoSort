'''
Class which stores and performs functions on variants from a sample

Created on Aug 30, 2016

@author: Jared Evans evans.jared@mayo.edu
'''
from math import log, floor, ceil
import re
import pysam
import RNA


rCRS_LEN = 16569
pileup_acgt_index = {'A_f':0, 'A_r':1, 'C_f':2, 'C_r':3, 'G_f':4, 'G_r':5, 'T_f':6, 'T_r':7, 'D':8}

class Variants:
    '''
    A Class for all variant info for a sample
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.snvs = {} # pos_ref_alt -> heteroplasmy
        self.indels = {} # pos_ref_alt -> heteroplasmy
        self.variants_allelic_depth = {} # pos_ref_alt -> refsupport,altsupport
        self.pileup = []
        for i in range(0,rCRS_LEN):
            acgt_counts = [0,0,0,0,0,0,0,0,0] # A_f, A_r, C_f, C_r, G_f, G_r, T_f, T_r, Del
            self.pileup.append(acgt_counts)
        self.pileup_indel = {} # pos_ref_alt -> (ref_f, ref_r, alt_f, alt_r)
        self.genotype = {} # pos -> genotype
        self.variants_phylotree = {} # pos_ref_alt -> variant in phylotree format
        self.variants_reverse_phylotree = {} # Key: phylotree format. Value: pos_ref_alt format
        self.variants_haplogrep = {} # The haplogrep label for each variant (yes, no, hotspot, localPrivateMut, globalPrivateMut)
        self.variants_genbank = {}
        self.variants_polyphen2_pred = {}
        self.variants_polyphen2_score = {}
        self.variants_sift_pred = {}
        self.variants_sift_score = {}
        self.variants_cadd_score = {}
        self.variants_provean_pred = {}
        self.variants_provean_score = {}
        self.variants_mutationassessor_pred = {}
        self.variants_mutationassessor_score = {}
        self.variants_disease = {}
        self.variants_locus = {}
        self.variants_locus_type = {}
        self.variants_effect = {}
        self.variants_aminoacid_change = {}
        self.variants_aminoacid_position = {}
        self.variants_codon_change = {}
        self.variants_mitomap = {}
        self.variants_somatic = {}
        self.variants_backfill = {} # pos_ref_alt -> heteroplasmy
        self.variants_rnafold = {}
        
        
    def pileup_increment_read_count(self, ref_position, nucleotide, forward_strand=True, increment=1):
        '''
        Given a pos and nucl, increment the number of reads supporting that nucl 
        '''
        strand = 'f'
        if not forward_strand:
            strand = 'r'
        key = nucleotide+"_"+strand
        if nucleotide == 'D':
            key = nucleotide
        self.pileup[ref_position][pileup_acgt_index[key]] += increment
        
    def pileup_indel_increment_read_count(self, ref_position, ref_allele, alt_allele, forward_strand=True, increment_alt=True, increment=1):
        '''
        Given a pos, ref, and alt allele, increment the number of reads supporting the INDEL 
        '''
        if 'N' not in ref_allele and 'N' not in alt_allele:
            key = str(ref_position)+"_"+ref_allele+"_"+alt_allele
            if not key in self.pileup_indel:
                self.pileup_indel[key] = [0,0,0,0]
            index = 0
            if increment_alt:
                index = 2
            if forward_strand:
                self.pileup_indel[key][index] += 1
            else:
                self.pileup_indel[key][index+1] += 1
        
    def get_depth(self, ref_position):
        '''
        Return the total read depth for a position
        '''
        return int(sum(self.pileup[int(ref_position)]))

    def get_nucleotide_support(self, ref_position, nucleotide):
        '''
        Return the total number of F and R reads supporting a nucleotide
        '''
        return self.pileup[int(ref_position)][pileup_acgt_index[nucleotide+"_f"]] + self.pileup[int(ref_position)][pileup_acgt_index[nucleotide+"_r"]]
        
    def strand_bias(self, ref_position, ref_allele, alt_allele):
        '''
        Calculate and return GATKs StrandOddsRatio
        '''
        ref_f = 0.01
        ref_r = 0.01
        alt_f = 0.01
        alt_r = 0.01
        
        if(len(ref_allele) == 1 and len(alt_allele) == 1):
            # SNV
            ref_f += self.pileup[ref_position][pileup_acgt_index[ref_allele+"_f"]]
            ref_r += self.pileup[ref_position][pileup_acgt_index[ref_allele+"_r"]]
            alt_f += self.pileup[ref_position][pileup_acgt_index[alt_allele+"_f"]]
            alt_r += self.pileup[ref_position][pileup_acgt_index[alt_allele+"_r"]]
        elif(len(ref_allele) > 1 and len(alt_allele) == 1):
            # DEL
            indel_support = self.pileup_indel[str(ref_position)+"_"+ref_allele+"_"+alt_allele]
            #ref_f += self.pileup[ref_position][pileup_acgt_index[alt_allele+"_f"]]
            #ref_r += self.pileup[ref_position][pileup_acgt_index[alt_allele+"_r"]]
            ref_f += indel_support[0]
            ref_r += indel_support[1]
            alt_f += indel_support[2]
            alt_r += indel_support[3]
        else:
            # INS
            indel_support = self.pileup_indel[str(ref_position)+"_"+ref_allele+"_"+alt_allele]
            #ref_f += self.pileup[ref_position][pileup_acgt_index[ref_allele+"_f"]]
            #ref_r += self.pileup[ref_position][pileup_acgt_index[ref_allele+"_r"]]
            ref_f += indel_support[0]
            ref_r += indel_support[1]
            alt_f += indel_support[2]
            alt_r += indel_support[3]
        
        R = (ref_f*alt_r)/(ref_r*alt_f)
        #refRatio = max(ref_f,ref_r)/min(ref_f,ref_r)
        #altRatio = max(alt_f,alt_r)/min(alt_f,alt_r)
        #return log((refRatio/altRatio)*(R+(1/R)))
        return log((((ref_f/ref_r)/(alt_r/alt_f))+((ref_r/ref_f)/(alt_f/alt_r)))*(R+(1/R)))
        
    def add_variant(self, ref_position, ref_allele, alt_allele, freq_support):
        '''
        Add a called variant to the variant list and store the ref and alt supporting reads
        '''
        if(len(ref_allele) == 1 and len(alt_allele) == 1):
            # SNV
            self.snvs[str(ref_position)+"_"+ref_allele+"_"+alt_allele] = freq_support
        else:
            # INDEL
            self.indels[str(ref_position)+"_"+ref_allele+"_"+alt_allele] = freq_support
        
        ad = self.get_allelic_depth(ref_position, ref_allele, alt_allele)
        self.variants_allelic_depth[str(ref_position)+"_"+ref_allele+"_"+alt_allele] = ad


    def get_allelic_depth(self, ref_position, ref_allele, alt_allele):
        '''
        Calculate the ref and alt supporting read depths 
        '''
        ref_support = 0
        alt_support = 0
        variant = str(ref_position)+"_"+ref_allele+"_"+alt_allele
        if(len(ref_allele) == 1 and len(alt_allele) == 1):
            # SNV
            ref_support = self.get_nucleotide_support(ref_position, ref_allele)
            alt_support = self.get_nucleotide_support(ref_position, alt_allele)
        elif len(ref_allele) > 1 and len(alt_allele) == 1:
            # DEL
            if variant in self.pileup_indel.keys():
                ref_f, ref_r, alt_f, alt_r = self.pileup_indel[variant]
                alt_support = int(alt_f) + int(alt_r)
                ref_support = int(ref_f) + int(ref_r)
            else:
                ref_support = self.get_nucleotide_support(ref_position, ref_allele[0])
        else:
            # INS
            if variant in self.pileup_indel.keys():
                ref_f, ref_r, alt_f, alt_r = self.pileup_indel[variant]
                ref_support = int(ref_f) + int(ref_r)
                alt_support = int(alt_f) + int(alt_r)
            else:
                ref_support = self.get_nucleotide_support(ref_position, ref_allele[0])
        return str(ref_support) + "," + str(alt_support)


    def calculate_genotypes(self, min_het):
        '''
        Calculate the genotype for each variant
        '''
        min_het = 0.01
        # SNVs
        for variant, freq in self.snvs.items():
            pos, ref, alt = variant.split("_")
            genotype = alt + "/" + alt
            ref_support, alt_support = self.variants_allelic_depth[variant].split(",")
            depth = float(self.get_depth(pos))
            if (float(ref_support)/depth) >= min_het:
                genotype = ref + "/" + alt
            self.genotype[variant] = genotype
        # INDELs
        for variant, freq in self.indels.items():
            pos, ref, alt = variant.split("_")
            genotype = alt + "/" + alt
            ref_support, alt_support = self.variants_allelic_depth[variant].split(",")
            depth = float(self.get_depth(pos))
            if (float(ref_support)/depth) >= min_het:
                genotype = ref + "/" + alt
            self.genotype[variant] = genotype
                    
    
    def get_phylotree_format(self, min_het_for_haplogrouping):
        '''
        Converts and stores all variants in phylotree format (eg 123A 123d 123.1A)
        Returns phylotree format homoplasmic variants to be used by haplogrep
        '''
        mutations = []
        tmp_insertions = {}
        
        # SNVs
        for variant in self.snvs.keys():
            pos, ref, alt = variant.split("_")
            pos = int(pos) + 1
            self.variants_phylotree[variant] = str(pos)+alt
            self.variants_reverse_phylotree[str(pos)+alt] = variant
            if self.snvs[variant] >= min_het_for_haplogrouping:
                # only return mutations for haplogrouping if above heteroplasmy threshold
                mutations.append(str(pos)+alt)
        for indel in self.indels.keys():
            pos, ref, alt = indel.split("_")
            pos = int(pos) + 1
            phylo_format = ""
            if len(ref) > 1 and len(alt) == 1:
                # DEL
                if len(ref) == 2:
                    phylo_format = str(pos)+"d"
                else:
                    stop = pos + len(ref)
                    phylo_format = str(pos)+"-"+str(stop)+"d"
            else:
                # INS
                phylo_format = str(pos)+".1"+alt[1:len(alt)]
                if self.indels[indel] >= min_het_for_haplogrouping:
                    tmp_insertions[str(pos-1) + "_" + ref + "_" + alt] = phylo_format # for temp hack (see below)
            self.variants_phylotree[indel] = phylo_format
            self.variants_reverse_phylotree[phylo_format] = indel
            if self.indels[indel] >= min_het_for_haplogrouping:
                mutations.append(phylo_format)
        # sort mutation list by position
        mutations.sort(key=self.get_position_sort_key)
        ## there is a bug in haplogrep if 2 INS have the same position so this is a temp hack to get around that. Just keep the INS with highest heteroplasmy:
        insertions_to_keep = {}
        for vari1, phylo1 in tmp_insertions.items():
            pos1, ref1, alt1 = vari1.split("_")
            keep = 1
            for vari2, phylo2 in tmp_insertions.items():
                pos2, ref2, alt2 = vari2.split("_")
                if pos1 == pos2 and vari1 != vari2:
                    if self.indels[vari1] <= self.indels[vari2]:
                        keep = 0
            if keep:
                insertions_to_keep[vari1] = 1
        for vari, phylo in tmp_insertions.items():
            if not vari in insertions_to_keep:
                # remove problem variants
                mutations.remove(phylo)
        ## end temp hack
        
        return mutations
    
    def get_position_sort_key(self, variant):
        '''
        function to return the correct Key that should be sorted on
        '''
        split_variant = filter(None,re.split(r'(\d+)',variant))
        return int(split_variant[0])
    
    def get_sorted_mutation_list(self):
        '''
        Return a combined list of SNVs/INDELs in sorted order
        '''
        tmp_mutations = []
        for key in self.snvs.keys():
            tmp_mutations.append(key)
        for key in self.indels.keys():
            tmp_mutations.append(key)
    
        return sorted(tmp_mutations, key=self.get_position_sort_key)
        
    def calculate_coding_effect(self, variant, codon_table, coding_genes, ref_fasta):
        '''
        Calculate the variant effect, amino acid pos, amino acid change for coding genes
        '''
        effect = ""
        aa_change = ""
        aa_pos = ""
        pos, ref, alt = variant.split("_")
        pos = int(pos) + 1
        if variant in self.variants_locus:
            for locus in self.variants_locus[variant].split(","):
                if locus in coding_genes:
                    start, stop, sequence = coding_genes[locus].split(",")
                    start = int(start)
                    stop = int(stop)
                    dna_distance = pos - start
                    aa_distance = ceil(dna_distance/3)
                    codon_position = dna_distance % 3
                    ref_codon = ref_nucl = ref_fasta.fetch(reference="chrM",start=(start+dna_distance)-codon_position-1,end=(start+dna_distance)+(3-codon_position)-1)
                    ref_codon.upper()
                    alt_upper = alt.upper()
                    alt_codon_list = list(ref_codon)
                    alt_codon_list[codon_position] = alt_upper
                    alt_codon = "".join(alt_codon_list)
                    ref_aa = codon_table[ref_codon]
                    alt_aa = codon_table[alt_codon]
                    aa_change = ref_aa + "-" + alt_aa
                    aa_pos = aa_distance
                    if ref_aa == alt_aa:
                        effect = "synonymous"
                    else:
                        effect = "non-synonymous"
                        
        self.variants_effect[variant] = effect
        self.variants_aminoacid_change[variant] = aa_change
        self.variants_aminoacid_position[variant] = str(aa_pos)
        
    def annotate_rna_variant(self,variant,rna_start,rna_stop,rna_sequence):
        '''
        Annotate whether the RNA secondary structure is predicted to change with this mutation
        '''
        annotation = "Preserved"
        pos, ref, alt = variant.split("_")
        # if INDEL then assume the structure will change
        if len(ref) > 1 or len(alt) > 1:
            annotation = "Modified"
        else:
            rna_variant_sequence = rna_sequence[0:((int(pos)+1)-int(rna_start))] + alt + rna_sequence[((int(pos)+1)-int(rna_start))+1:len(rna_sequence)]
            fc_ref = RNA.fold_compound(rna_sequence)
            fc_var = RNA.fold_compound(rna_variant_sequence)
            (structure_ref, mfe_ref) = fc_ref.mfe()
            (structure_var, mfe_var) = fc_var.mfe()
            if structure_ref != structure_var:
                annotation = "Modified"
        self.variants_rnafold[variant] = annotation
        
