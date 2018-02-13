'''
Class used to store and access info for each sample

Created on Aug 17, 2016

@author: Jared Evans evans.jared@mayo.edu
'''



import os, sys, csv
import logging
import pickle
from variants import Variants


class Sample:
    ''' A class for Sample data '''
    def __init__(self):
        ''' Constructor '''
        self.sname = None
        self.group = None # iPSC line, family ID, or tumor/normal subject ID depending on the analysis type
        self.clone = None # ipsc
        self.passage_number = None # ipsc
        self.sample_type = None # possible values: iPSC, FB, normal, tumor, family
        self.in_file = []
        self.haplogroup = None
        self.haplogroup_score = None
        self.haplogroup_yes_total = 0
        self.haplogroup_no_total = 0
        self.haplogroup_gpm_total = 0
        self.haplogroup_lpm_total = 0
        self.variants = Variants()
        self.coverage = []
        self.file_type = None
        self.is_aligned = False
        self.relative_copy_number = 0
        self.read_length = 0
        self.ipsc_risk_category = "-"
        self.ipsc_risk_score = 0
        self.ipsc_highest_risk_variant = "-"
        self.mitochondrial_reads = 0
        self.nuclear_reads = 0
        self.idxstats_file = None
        self.is_pedigree = False # pedigree analysis
        self.is_somatic = False # tumor/normal analysis
        self.is_ipsc = False # ipsc analysis

    def load_sample_info(self, sample_info, ipsc, somatic):
        ''' Load initial sample_info.txt variables. '''
        if ipsc:
            self.group = sample_info[0]
            self.clone = sample_info[1]
            self.passage_number = sample_info[2]
            self.in_file = sample_info[3].split(',')
            self.sname = self.group + "." + self.clone + ".p" + self.passage_number
            self.is_ipsc = True
            self.sample_type = "iPSC"
            if self.clone.lower() == "fibroblast" or self.clone.lower() == "fb":
                self.sample_type = "FB"
        elif somatic:
            self.sname = sample_info[0]
            self.group = sample_info[1]
            self.sample_type = sample_info[2]
            self.in_file = sample_info[3].split(',')
            self.is_somatic = True
            if len(sample_info) == 5:
                self.idxstats_file = sample_info[4]
        else:
            self.sname = sample_info[0]
            self.in_file = sample_info[1].split(',')
            if len(sample_info) == 3:
                self.idxstats_file = sample_info[2]
                
        # check filetype extension
        filename, ext = os.path.splitext(self.in_file[0])
        if ext.lower() == ".bam":
            self.is_aligned = True
        # check if files exist
        for input_file in self.in_file:
            if not os.path.isfile(input_file):
                logging.error("Input file from sample_info doesn't exist: " + input_file)
                sys.exit(1)


    def write_vcf(self, output_dir):
        '''
        write variants from this sample to VCF file
        '''
        f = open(output_dir + "/variants/vcf/" + self.sname + ".vcf","w")
        f.write("##fileformat=VCFv4.3\n")
        f.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n")
        f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("##FORMAT=<ID=HG,Number=1,Type=String,Description=\"Haplogrep mutation type\">\n")
        f.write("##FORMAT=<ID=HP,Number=1,Type=Float,Description=\"Heteroplasmy fraction\">\n")
        f.write("##INFO=<ID=LOC,Number=1,Type=String,Description=\"Locus of mutation\">\n")
        f.write("##INFO=<ID=GF,Number=1,Type=Float,Description=\"Percentage of GenBank Alt Allele Frequency\">\n")
        f.write("##INFO=<ID=PT,Number=1,Type=String,Description=\"Mutation in PhyloTree format\">\n")
        f.write("##META=<ID=Haplogroup,Number=1,Type=String,Description=\"The predicted Haplogroup the sample belongs to\">\n")
        f.write("##META=<ID=Haplogroup_score,Number=1,Type=Float,Description=\"The haplogrep score associated with the predicted Haplogroup\">\n")
        f.write("##SAMPLE=<ID=" + self.sname + ",Haplogroup=" + self.haplogroup + ",Haplogroup_score=" + self.haplogroup_score + ">\n")
        f.write("##contig=<ID=chrM,length=16569>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + self.sname + "\n")
        for variant in self.variants.get_sorted_mutation_list():
            pos, ref, alt = variant.split("_")
            ad = self.variants.variants_allelic_depth[variant]
            locus = "." 
            if variant in self.variants.variants_locus:
                locus = self.variants.variants_locus[variant]
            phylotree = self.variants.variants_phylotree[variant]
            haplogroup = "."
            if variant in self.variants.variants_haplogrep:
                haplogroup = self.variants.variants_haplogrep[variant]
            het = "."
            if len(ref) == 1 and len(alt) == 1:
                het = self.variants.snvs[variant]
            else:
                het = self.variants.indels[variant]
            #depth = "%.0f" % self.variants.get_depth(pos)
            depth = self.variants.get_depth(pos)
            genbank_freq = "."
            if variant in self.variants.variants_genbank:
                genbank_freq = self.variants.variants_genbank[variant]
            genotype = "0/1"
            if variant in self.variants.genotype:
                base1, base2 = self.variants.genotype[variant].split("/")
                if base1 == alt:
                    genotype = "1/1"
            if het < 1:
                het = "%.4f" % het
            else:
                het = "1.0" 
            f.write("chrM\t" + str(int(pos)+1) + "\t.\t" + ref + "\t" + alt + "\t.\t.\tLOC=" + locus + ";GF=" + genbank_freq[0:len(genbank_freq)-1]  + ";PT=" + phylotree + "\tGT:AD:DP:HG:HP\t" + genotype + ":" + ad + ":" + str(depth) + ":" + haplogroup + ":" + str(het) + "\n")
        f.close()

    def write_annotation_tables(self, output_dir):
        '''
        write variants and annotations to output tables
        '''
        pm = open(output_dir + "/variants/annotation/" + self.sname + ".private_mutations.csv","w")
        av = open(output_dir + "/variants/annotation/" + self.sname + ".all_variants.csv","w")
        pm_writer = csv.writer(pm, delimiter=",", quoting=csv.QUOTE_MINIMAL)
        av_writer = csv.writer(av, delimiter=",", quoting=csv.QUOTE_MINIMAL)
        if self.is_ipsc:
            pm_writer.writerow(["Line","Clone","Passage_Number","Haplogroup","Haplogroup_Quality","Private_Category","Pos","Ref","Alt","Ref_Supporting_Reads","Alt_Supporting_Reads","Total_Depth","Alt_Supporting_Percent","Genotype","Variant_Type","Locus","MITOMAP_Poly","Codon_Change","Amino_Acid_Change","Amino_Acid_Position","Disease","GenBank_Allele_Frequency","PolyPhen2_Prediction","PolyPhen2_Score","SIFT_Prediction","SIFT_Score","PROVEAN_Prediction","PROVEAN_Score","CADD_Score"])
            av_writer.writerow(["Line","Clone","Passage_Number","Haplogroup","Haplogroup_Quality","Haplogrep_Category","Pos","Ref","Alt","Ref_Supporting_Reads","Alt_Supporting_Reads","Total_Depth","Alt_Supporting_Percent","Genotype","Variant_Type","Locus","MITOMAP_Poly","GenBank_Allele_Frequency"])
        else:
            pm_writer.writerow(["Sample","Haplogroup","Haplogroup_Quality","Private_Category","Pos","Ref","Alt","Ref_Supporting_Reads","Alt_Supporting_Reads","Total_Depth","Alt_Supporting_Percent","Genotype","Variant_Type","Locus","MITOMAP_Poly","Codon_Change","Amino_Acid_Change","Amino_Acid_Position","Disease","GenBank_Allele_Frequency","PolyPhen2_Prediction","PolyPhen2_Score","SIFT_Prediction","SIFT_Score","PROVEAN_Prediction","PROVEAN_Score","CADD_Score"])
            av_writer.writerow(["Sample","Haplogroup","Haplogroup_Quality","Haplogrep_Category","Pos","Ref","Alt","Ref_Supporting_Reads","Alt_Supporting_Reads","Total_Depth","Alt_Supporting_Percent","Genotype","Variant_Type","Locus","MITOMAP_Poly","GenBank_Allele_Frequency"])
        
        for variant in self.variants.get_sorted_mutation_list():
            pos, ref, alt = variant.split("_")
            ref_support, alt_support = self.variants.variants_allelic_depth[variant].split(",")
            depth = self.variants.get_depth(pos)
            haplogrep_cat = "."
            if variant in self.variants.variants_haplogrep:
                haplogrep_cat = self.variants.variants_haplogrep[variant]
            het = "."
            if len(ref) == 1 and len(alt) == 1:
                het = self.variants.snvs[variant]
            else:
                het = self.variants.indels[variant]
            genotype = ""
            if variant in self.variants.genotype:
                genotype = self.variants.genotype[variant]
            #genotype = ref + "/" + alt
            #if het > 0.95:
            #    genotype = alt + "/" + alt
            alt_perc = "%.2f" % (het*100)
            variant_type = "SNV"
            if len(ref) > len(alt):
                variant_type = "DEL"
            if len(ref) < len(alt):
                variant_type = "INS"
            locus = ""
            if variant in self.variants.variants_locus:
                locus = self.variants.variants_locus[variant]
            mitomap = ""
            if variant in self.variants.variants_mitomap:
                mitomap = "YES"
            genbank_freq = ""
            if variant in self.variants.variants_genbank:
                genbank_freq = self.variants.variants_genbank[variant]
                genbank_freq = genbank_freq[0:len(genbank_freq)-1]
            #effect = ""
            #if variant in self.variants.variants_effect:
            #    effect = self.variants.variants_effect[variant]
            amino_change = ""
            if variant in self.variants.variants_aminoacid_change:
                amino_change = self.variants.variants_aminoacid_change[variant]
            amino_pos = ""
            if variant in self.variants.variants_aminoacid_position:
                amino_pos = self.variants.variants_aminoacid_position[variant]
            codon_change = ""
            if variant in self.variants.variants_codon_change:
                codon_change = self.variants.variants_codon_change[variant]
            disease = ""
            if variant in self.variants.variants_disease:
                disease = self.variants.variants_disease[variant]
            polyphen_pred = ""
            polyphen_score = ""
            if variant in self.variants.variants_polyphen2_score:
                polyphen_pred = self.variants.variants_polyphen2_pred[variant]
                polyphen_score = self.variants.variants_polyphen2_score[variant]
            sift_cat = ""
            sift_score = ""
            if variant in self.variants.variants_sift_score:
                sift_cat = self.variants.variants_sift_pred[variant]
                sift_score = self.variants.variants_sift_score[variant]
            cadd_score = ""
            if variant in self.variants.variants_cadd_score:
                cadd_score = self.variants.variants_cadd_score[variant]
            provean_pred = ""
            provean_score = ""
            if variant in self.variants.variants_provean_score:
                provean_pred = self.variants.variants_provean_pred[variant]
                provean_score = self.variants.variants_provean_score[variant]
            mutassessor_pred = ""
            mutassessor_score = ""
            if variant in self.variants.variants_mutationassessor_score:
                mutassessor_pred = self.variants.variants_mutationassessor_pred[variant]
                mutassessor_score = self.variants.variants_mutationassessor_score[variant]
            rna_fold = ""
            if variant in self.variants.variants_rnafold:
                rna_fold = self.variants.variants_rnafold[variant]
            if self.is_ipsc:
                av_writer.writerow([self.group, self.clone, self.passage_number, self.haplogroup, str(self.haplogroup_score), haplogrep_cat, str(int(pos)+1), ref, alt, ref_support, alt_support, str(depth), alt_perc, genotype, variant_type, locus, mitomap, genbank_freq])
            else:
                av_writer.writerow([self.sname, self.haplogroup, str(self.haplogroup_score), haplogrep_cat, str(int(pos)+1), ref, alt, ref_support, alt_support, str(depth), alt_perc, genotype, variant_type, locus, mitomap, genbank_freq])
            if haplogrep_cat == "globalPrivateMut" or haplogrep_cat == "localPrivateMut" or haplogrep_cat == "backfilledPrivateMut":
                if self.is_ipsc:
                    pm_writer.writerow([self.group, self.clone, self.passage_number, self.haplogroup, str(self.haplogroup_score), haplogrep_cat[:-10], str(int(pos)+1), ref, alt, ref_support, alt_support, str(depth), alt_perc, genotype, variant_type, locus, mitomap, codon_change, amino_change, amino_pos, disease, genbank_freq, polyphen_pred, polyphen_score, sift_cat, sift_score, provean_pred, provean_score, mutassessor_pred, mutassessor_score, cadd_score, rna_fold])
                else:
                    pm_writer.writerow([self.sname, self.haplogroup, str(self.haplogroup_score), haplogrep_cat[:-10], str(int(pos)+1), ref, alt, ref_support, alt_support, str(depth), alt_perc, genotype, variant_type, locus, mitomap, codon_change, amino_change, amino_pos, disease, genbank_freq, polyphen_pred, polyphen_score, sift_cat, sift_score, provean_pred, provean_score, mutassessor_pred, mutassessor_score, cadd_score, rna_fold])
        av.close()
        pm.close()            
        

    def write_pileup(self, output_dir):
        '''
        write pileup files
        '''
        f = open(output_dir + "/variants/pileup/" + self.sname + ".pileup_counts_snvs.txt","w")
        f.write("Pos\tA_f\tA_r\tC_f\tC_r\tG_f\tG_r\tT_f\tT_r\tDel\n")
        pos = 1
        for row in self.variants.pileup:
            delim = "\t"
            f.write(str(pos) + "\t" + str(row[0]) + "\t" + str(row[1]) + "\t" + str(row[2]) + "\t" + str(row[3]) + "\t" + str(row[4]) + "\t" + str(row[5]) + "\t" + str(row[6]) + "\t" + str(row[7]) + "\t" + str(row[8]) + "\n")
            pos += 1
        f.close()
        
        f = open(output_dir + "/variants/pileup/" + self.sname + ".pileup_counts_indels.txt","w")
        f.write("Pos\tRef\tAlt\tRef_f\tRef_r\tAlt_f\tAlt_r\n")
        sorted_indels = sorted(self.variants.pileup_indel.keys(),key=self.variants.get_position_sort_key)
        for indel in sorted_indels:
            pos, ref, alt = indel.split("_")
            ref_fwd, ref_rev, alt_fwd, alt_rev = self.variants.pileup_indel[indel]
            f.write(str(int(pos)+1) + "\t" + ref + "\t" + alt + "\t" + str(ref_fwd) + "\t" + str(ref_rev) + "\t" + str(alt_fwd) + "\t" + str(alt_rev) + "\n") 
        f.close()




