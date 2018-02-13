'''
Analysis module for the MitoSort package

Created on Aug 11, 2016

@author: Jared Evans evans.jared@mayo.edu 
'''

import ConfigParser
import pickle
import sys, os, csv, subprocess
import logging
import pysam
from shutil import copytree, copy
from sample import Sample
import html
import parsing
import plotting

def backfill_mutations(sample_objects, output_dir, min_het_to_backfill):
    ''' Backfill mutations across all samples '''
    logging.info("Backfilling mutations")
    # get master list of all mutations for backfilled VCF
    all_mutations = {}
    sample_metadata = ""
    sample_names = ""
    group_analysis = False
    for sobj in sample_objects:
        mutations = sobj.variants.get_sorted_mutation_list()
        for mutation in mutations:
            pos, ref, alt = mutation.split("_")
            all_mutations[mutation] = sobj.variants.variants_phylotree[mutation]
        sample_metadata += "##SAMPLE=<ID=" + sobj.sname + ",Haplogroup=" + sobj.haplogroup + ",Haplogroup_score=" + sobj.haplogroup_score + ">\n"
        if sample_names == "":
            sample_names = sobj.sname
        else:
            sample_names += "\t" + sobj.sname
        if sobj.is_ipsc or sobj.is_somatic or sobj.is_pedigree:
            group_analysis = True
    f = open(output_dir + "/variants/vcf/all_samples.vcf","w")
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
    f.write(sample_metadata)
    f.write("##contig=<ID=chrM,length=16569>\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_names + "\n")
    # backfill each mutation
    for mutation in sorted(all_mutations.keys(), key=sample_objects[0].variants.get_position_sort_key):
        pos, ref, alt = mutation.split("_")
        f.write("chrM\t" + str(int(pos)+1) + "\t.\t" + ref + "\t" + alt + "\t.\t.")
        locus = "."
        phylotree = "."
        genbank_freq = "."
        row_to_write = ""
        for sobj in sample_objects:
            # check for INFO annotations 
            if mutation in sobj.variants.variants_locus:
                locus = sobj.variants.variants_locus[mutation]
            phylotree = all_mutations[mutation]
            if mutation in sobj.variants.variants_genbank:
                genbank_freq = sobj.variants.variants_genbank[mutation]
            # calculate sample-specific values
            depth = float(sobj.variants.get_depth(pos))
            heteroplasmy = 0
            ad = sobj.variants.get_allelic_depth(pos, ref, alt)
            ref_support, alt_support = ad.split(",")
            if depth > 0:
                heteroplasmy = (int(alt_support)/depth)
            haplogrep_mut_type = "."
            if mutation in sobj.variants.variants_haplogrep:
                haplogrep_mut_type = sobj.variants.variants_haplogrep[mutation]
            if len(ref) == 1 and len(alt) == 1:
                if mutation not in sobj.variants.snvs.keys():
                    sobj.variants.variants_backfill[mutation] = heteroplasmy
            else:
                if mutation not in sobj.variants.indels.keys():
                    sobj.variants.variants_backfill[mutation] = heteroplasmy
            genotype = "./."
            if mutation in sobj.variants.genotype:
                allele_genotype = sobj.variants.genotype[mutation]
                base1, base2 = sobj.variants.genotype[mutation].split("/")
                if base1 == alt and base2 == alt:
                    genotype = "1/1"
                if base1 == ref and base2 == alt:
                    genotype = "0/1"
            elif depth > 0:
                genotype = "0/0"
            if heteroplasmy > 0 and heteroplasmy < 1:
                heteroplasmy = "%.4f" % heteroplasmy
            elif heteroplasmy == 0:
                heteroplasmy = "0"
            elif heteroplasmy == 1:
                heteroplasmy = "1.0"
            row_to_write += "\t" + genotype + ":" + ad + ":" + str(int(depth)) + ":" + haplogrep_mut_type + ":" + str(heteroplasmy)
        f.write("\tLOC=" + locus + ";GF=" + genbank_freq[0:len(genbank_freq)-1]  + ";PT=" + phylotree + "\tGT:AD:DP:HG:HP" + row_to_write + "\n")
    f.close()   
        
    # for samples in the same group, add their backfilled private mutations to their mutation lists
    if group_analysis:
        groups = {}
        # get list of all groups
        for sobj in sample_objects:
            groups[sobj.group] = 1
        for group in groups.keys():
            group_private_muts = {}
            # get master list of mutations to backfill for this group
            for sobj in sample_objects:
                if sobj.group == group:
                    mutations = sobj.variants.get_sorted_mutation_list()
                    for mutation in mutations:
                        if mutation in sobj.variants.variants_haplogrep.keys():
                            pos, ref, alt = mutation.split("_")
                            het = 0
                            if len(ref) == 1 and len(alt) == 1:
                                het = sobj.variants.snvs[mutation]
                            else:
                                het = sobj.variants.indels[mutation]
                            if het >= min_het_to_backfill and (sobj.variants.variants_haplogrep[mutation] == "globalPrivateMut" or sobj.variants.variants_haplogrep[mutation] == "localPrivateMut"):
                                group_private_muts[mutation] = sobj.variants.variants_phylotree[mutation]
            # backfill each of the mutations
            for sobj in sample_objects:
                if sobj.group == group:
                    for mut in group_private_muts.keys():
                        pos, ref, alt = mut.split("_")
                        if len(ref) == 1 and len(alt) == 1:
                            if mut not in sobj.variants.snvs.keys():
                                depth = float(sobj.variants.get_depth(pos))
                                ref_support = sobj.variants.get_nucleotide_support(pos,ref)
                                alt_support = sobj.variants.get_nucleotide_support(pos,alt)
                                heteroplasmy = 0
                                if depth > 0:
                                    heteroplasmy = (alt_support/depth)
                                sobj.variants.variants_backfill[mut] = heteroplasmy
                                # add mutation to sample's list, but give private cat as "backfilled"
                                sobj.variants.add_variant(pos, ref, alt, heteroplasmy)
                                sobj.variants.variants_haplogrep[mut] = "backfilledPrivateMut"
                        else:
                            if mut not in sobj.variants.indels.keys() and mut in sobj.variants.pileup_indel.keys():
                                depth = float(sobj.variants.get_depth(pos))
                                indel_support = sobj.variants.pileup_indel[mut]
                                alt_support = indel_support[2] + indel_support[3]
                                heteroplasmy = 0
                                if depth > 0:
                                    heteroplasmy = (alt_support/depth)
                                sobj.variants.variants_backfill[mut] = heteroplasmy
                                # add mutation to sample's list, but give private cat as "backfilled"
                                sobj.variants.add_variant(pos, ref, alt, heteroplasmy)
                                sobj.variants.variants_haplogrep[mut] = "backfilledPrivateMut"
                    
            # print backfilled private mutation tables to be used for plotting
            f = open(output_dir + "/variants/plot/" + group + ".private_mutations.txt","w")
            f.write("Sample\tSample_Type\tPosition\tVariant\tVariant_Type\tHeteroplasmy\tMut_Type\n")
            for mut in group_private_muts.keys():
                pos, ref, alt = mut.split("_")
                var_type = "SNV"
                if len(ref) > len(alt):
                    var_type = "DEL"
                if len(ref) < len(alt):
                    var_type = "INS"
                phylo = group_private_muts[mut]
                for sobj in sample_objects:
                    if sobj.group == group:
                        sample_type = sobj.sample_type
                        het = 0
                        if mut in sobj.variants.snvs.keys():
                            het = sobj.variants.snvs[mut]
                        if mut in sobj.variants.indels.keys():
                            het = sobj.variants.indels[mut]
                        if mut in sobj.variants.variants_backfill.keys():
                            het = sobj.variants.variants_backfill[mut]
                        if mut in sobj.variants.variants_haplogrep.keys():
                            mut_type = sobj.variants.variants_haplogrep[mut]
                        
                        het = "%.2f" % (het*100)
                        if sample_type == "iPSC":
                            sample_name = "s." + sobj.clone + ".p" + str(sobj.passage_number)
                        elif sample_type == "FB":
                            sample_name = "s." + sample_type + ".p" + str(sobj.passage_number)
                        else:
                            sample_name = "s." + sobj.sname
                        
                        f.write(sample_name + "\t" + sample_type + "\t" + str(int(pos)+1) + "\t" + \
                                phylo + "\t" + var_type + "\t" + str(het) + "\t" + mut_type + "\n")
            f.close()
            # rewrite annotation tables to include backfilled variants
            for sobj in sample_objects:
                sobj.write_annotation_tables(output_dir)

def calculate_coverage(sobj, output_dir, gene_file ,log_scale, coverage_resolution, rcrs_len):
    ''' Generate coverage plots and output files '''
    logging.info("Calculating coverage for sample: " + sobj.sname)
    wig_file = output_dir + "/coverage/wig/" + sobj.sname + ".wig"
    f = open(wig_file,"w")
    f.write("fixedStep chrom=chrM start=1 step=" + str(coverage_resolution) + " span=" + str(coverage_resolution) + "\n")
    coverage = []
    no_cov = 0
    counter = 0
    step_sum = 0.0
    for i in range(0,rcrs_len):
        counter += 1
        i_cov = float(sobj.variants.get_depth(i))
        step_sum += i_cov
        coverage.append(i_cov)
        if i_cov == 0:
            no_cov += 1
        if counter == coverage_resolution:
            f.write(str(int(step_sum/coverage_resolution)) + "\n")
            step_sum = 0.0
            counter = 0
    f.close()
    coverage.sort()

    s = open(output_dir + "/coverage/stats/" + sobj.sname + ".coverage_stats.csv","w")
    s_writer = csv.writer(s, delimiter=",", quoting=csv.QUOTE_MINIMAL)
    mt_copies = "-"
    if sobj.relative_copy_number != 0:
        mt_copies = str("%.2f" % sobj.relative_copy_number)
    avg_cov = "%.2f" % (sum(coverage)/float(rcrs_len)) 
    #perc_no_cov = "%.2f" % ((no_cov/float(rcrs_len))*100)
    q1_pos = int(rcrs_len*0.25)
    med_pos = int(rcrs_len*0.5)
    q2_pos = int(rcrs_len*0.75)
    s_writer.writerow(["Sample","Mitochondrial_Reads","Nuclear_Reads","Mitochondrial_Copies",
                       "Average_Coverage","Minimum_Coverage","First_Quartile_Coverage",
                       "Median_Coverage","Third_Quartile_Coverage","Maximum_Coverage","Bases_With_No_Coverage"])
    s_writer.writerow([sobj.sname, str(sobj.mitochondrial_reads), str(sobj.nuclear_reads), 
                       mt_copies, str(avg_cov), str(int(coverage[0])), str(int(coverage[q1_pos])), 
                       str(int(coverage[med_pos])), str(int(coverage[q2_pos])), 
                       str(int(coverage[rcrs_len-1])), str(no_cov)])
    s.close()
    
    # plot coverage
    plotting.plotCoverage(wig_file, output_dir + '/coverage/plot/' + sobj.sname + '.coverage.png',
                  sobj.sname, gene_file, log_scale)
    # plot_cmd = [rscript, rlib + "/plot_coverage.R", wig_file, output_dir + "/coverage/plot/" + sobj.sname + ".coverage.png", sobj.sname, gene_file, log_scale]
    # subprocess.call(plot_cmd)


def ipsc_annotation(sample_objects):
    ''' Calculate iPSC quality annotations '''
    
    for sobj in sample_objects:
        # Only score the iPSCs, not the Fibroblasts
        if sobj.sample_type == "iPSC":
            sobj.ipsc_risk_category = "low"
            # find the line's associated FB to compare with
            fb_sobj = None
            for sobj2 in sample_objects:
                if sobj.group == sobj2.group and sobj2.sample_type == "FB":
                    fb_sobj = sobj2
            # SNVs
            for variant in sobj.variants.snvs:
                variant_heteroplasmy = sobj.variants.snvs[variant]
                if variant in sobj.variants.variants_haplogrep:
                    # Only GlobalPrivateMutations are used in the scoring
                    if sobj.variants.variants_haplogrep[variant] == "globalPrivateMut":
                        fb_heteroplasmy = 0
                        if fb_sobj:
                            if variant in fb_sobj.variants.snvs:
                                fb_heteroplasmy = fb_sobj.variants.snvs[variant]
                        # check if homoplasmic in FB
                        homoplasmic_in_fb = False
                        if fb_heteroplasmy > 0.90:
                            homoplasmic_in_fb = True
                        # check if either known disease causing or predicted damaging
                        disease_mutation = False
                        predicted_deleterious = False
                        if variant in sobj.variants.variants_disease:
                            disease_mutation = True
                        # check if polymorphic in mitomap
                        mitomap_poly = False
                        if variant in sobj.variants.variants_mitomap:
                            mitomap_poly = True
                        # check if RNA variant and predicted to modify secondary structure
                        if variant in sobj.variants.variants_rnafold:
                            if sobj.variants.variants_rnafold[variant] == "Modified":
                                predicted_deleterious = True
                        # check if protein_coding and predicted damaging
                        predicted_damaging_count = 0
                        if variant in sobj.variants.variants_polyphen2_pred:
                            if sobj.variants.variants_polyphen2_pred[variant] == "probably_damaging" or sobj.variants.variants_polyphen2_pred[variant] == "possibly_damaging":
                                predicted_damaging_count += 1
                        if variant in sobj.variants.variants_sift_pred:
                            if sobj.variants.variants_sift_pred[variant] == "deleterious":
                                predicted_damaging_count += 1
                        if variant in sobj.variants.variants_provean_pred:
                            if sobj.variants.variants_provean_pred[variant] == "deleterious":
                                predicted_damaging_count += 1
                        if variant in sobj.variants.variants_mutationassessor_pred:
                            if sobj.variants.variants_mutationassessor_pred[variant] == "high_impact" or sobj.variants.variants_mutationassessor_pred[variant] == "medium_impact":
                                predicted_damaging_count += 1
                        # if protein_coding variant and predicted damaging in 3/4 of tools then consider deleterious
                        if predicted_damaging_count >= 3:
                            predicted_deleterious = True
                        # classify clone quality
                        if disease_mutation or (predicted_deleterious and not mitomap_poly and not homoplasmic_in_fb):
                            # ipsc score is just max heteroplasmy of all damaging mutations
                            if variant_heteroplasmy > float(sobj.ipsc_risk_score):
                                if variant_heteroplasmy > 0 and variant_heteroplasmy < 1:
                                    sobj.ipsc_risk_score = "%.3f" % variant_heteroplasmy
                                else:
                                    sobj.ipsc_risk_score = variant_heteroplasmy
                                if variant in sobj.variants.variants_phylotree:
                                    sobj.ipsc_highest_risk_variant = sobj.variants.variants_phylotree[variant]
                            if variant_heteroplasmy > 0.7:
                                sobj.ipsc_risk_category = "high"
                            elif variant_heteroplasmy >= 0.3 and variant_heteroplasmy <= 0.7 and sobj.ipsc_risk_category != "high":
                                sobj.ipsc_risk_category = "unknown"
            # INDELs
            for variant in sobj.variants.indels:
                variant_heteroplasmy = sobj.variants.indels[variant]
                if variant in sobj.variants.variants_haplogrep:
                    # Only GlobalPrivateMutations are used in the scoring
                    if sobj.variants.variants_haplogrep[variant] == "globalPrivateMut" and variant != "309_T_TC":
                        fb_heteroplasmy = 0
                        if fb_sobj:
                            if variant in fb_sobj.variants.indels:
                                fb_heteroplasmy = fb_sobj.variants.indels[variant]
                        # check if homoplasmic in FB
                        homoplasmic_in_fb = False
                        if fb_heteroplasmy > 0.90:
                            homoplasmic_in_fb = True
                        # check if either known disease causing or predicted damaging
                        disease_mutation = False
                        predicted_deleterious = False
                        if variant in sobj.variants.variants_disease:
                            disease_mutation = True
                        # check if polymorphic in mitomap
                        mitomap_poly = False
                        if variant in sobj.variants.variants_mitomap:
                            mitomap_poly = True
                        # check if RNA variant and predicted to modify secondary structure
                        if variant in sobj.variants.variants_rnafold:
                            if sobj.variants.variants_rnafold[variant] == "Modified":
                                predicted_deleterious = True
                        # check if INDEL is frameshift and protein_coding
                        if variant in sobj.variants.variants_locus_type:
                            if sobj.variants.variants_locus_type[variant] == "protein_coding":
                                pos, ref, alt = variant.split("_")
                                alt_len = 0
                                if len(alt) > len(ref):
                                    alt_len = len(alt)-1
                                else:
                                    alt_len = len(ref)-1
                                if alt_len%3 != 0:
                                    # frameshift
                                    predicted_deleterious = True
                        # classify clone quality
                        if disease_mutation or (predicted_deleterious and not mitomap_poly and not homoplasmic_in_fb):
                            # ipsc score is just max heteroplasmy of all damaging mutations
                            if variant_heteroplasmy > float(sobj.ipsc_risk_score):
                                if variant_heteroplasmy > 0 and variant_heteroplasmy < 1:
                                    sobj.ipsc_risk_score = "%.3f" % variant_heteroplasmy
                                else:
                                    sobj.ipsc_risk_score = variant_heteroplasmy
                                if variant in sobj.variants.variants_phylotree:
                                    sobj.ipsc_highest_risk_variant = sobj.variants.variants_phylotree[variant]
                            if variant_heteroplasmy > 0.7:
                                sobj.ipsc_risk_category = "high"
                            elif variant_heteroplasmy >= 0.3 and variant_heteroplasmy <= 0.7 and sobj.ipsc_risk_category != "high":
                                sobj.ipsc_risk_category = "unknown"
    


def annotate_variants(sobj, mito_genes_file, mitomap_counts_file, mitomap_poly_coding_file, mitomap_poly_control_file, mitomap_mut_coding_control_file, mitomap_mut_rna_file, mitomap_mut_somatic_file, mitimpact_file, mito_ref):
    ''' Annotate variants for a sample. '''
    logging.info("Annotating variants for sample: " + sobj.sname)
    # open reference FASTA
    try:
        ref_fasta = pysam.FastaFile(filename=mito_ref)
        #print ref_fasta.fetch(reference="chrM",start=0,end=1)
    except:
        logging.error("Opening reference FASTA " + mito_ref + " " + str(sys.exc_info()[0]))
        sys.exit(1)
        
    # load locus annotation 
    locus = []
    try:
        with open(mito_genes_file,'r') as f:
            reader = csv.reader(f,delimiter='\t',quoting=csv.QUOTE_NONE)
            for row in reader:
                if row[0] != "Map Locus":
                    sep = "\t"
                    locus.append(sep.join(row))
    except:
        logging.error("Reading " + mito_genes_file)
        sys.exit(1)
    # load mitomap allele frequencies 
    mitomap_freq = {}
    try:
        with open(mitomap_counts_file,'r') as f:
            reader = csv.reader(f,delimiter=',',quoting=csv.QUOTE_NONE)
            for row in reader:
                if row[0] != "rCRS":
                    mitomap_freq[row[1] + "_" + row[0] + "_" + row[2]] = row[6] # GB allele freq percent
    except:
        logging.error("Reading " + mitomap_counts_file)
        sys.exit(1)
    # load mitmap poly coding 
    mitomap_poly_coding = {}
    try:
        with open(mitomap_poly_coding_file,'r') as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                if row[0] != "Position":
                    alleles = row[2].split("-")
                    if len(alleles) > 1:
                        mitomap_poly_coding[row[0] + "_" + alleles[0] + "_" + alleles[1]] = row[3] + "_" + row[5] # CodonNumber_Effect:AminoAcidChange
    except:
        logging.error("Reading " + mitomap_poly_coding_file)
        print sys.exc_info()[0]
        raise
        sys.exit(1)
    # load mitomap poly control 
    mitomap_poly_control = {}
    try:
        with open(mitomap_poly_control_file,'r') as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                if row[0] != "Position":
                    alleles = row[2].split("-")
                    if len(alleles) > 1:
                        mitomap_poly_control[row[0] + "_" + alleles[0] + "_" + alleles[1]] = 1 # boolean to indicate that it is present in MitoMap
    except:
        logging.error("Reading " + mitomap_poly_control_file)
        print sys.exc_info()[0]
        raise
        sys.exit(1)
    # load mitomap mut coding/control
    mitomap_mut_coding_control = {}
    try:
        with open(mitomap_mut_coding_control_file,'r') as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                if row[0] != "Position":
                    alleles = row[4].split("-")
                    if len(alleles) > 1:
                        mitomap_mut_coding_control[row[0] + "_" + alleles[0] + "_" + alleles[1]] = row[2] # disease
    except:
        logging.error("Reading " + mitomap_mut_coding_control_file)
        print sys.exc_info()[0]
        raise
        sys.exit(1)
    # load mitomap mut RNA 
    mitomap_mut_rna = {}
    try:
        with open(mitomap_mut_rna_file,'r') as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                if row[0] != "Position":
                    mitomap_mut_rna[row[3]] = row[2] # disease
    except:
        logging.error("Reading " + mitomap_mut_rna_file)
        print sys.exc_info()[0]
        raise
        sys.exit(1)
    # load mitomap mut Somatic 
    mitomap_mut_somatic = {}
    try:
        with open(mitomap_mut_somatic_file,'r') as f:
            reader = csv.reader(f,delimiter=',')
            for row in reader:
                if row[0] != "Position":
                    alleles = row[2].split("-")
                    if len(alleles) > 1:
                        mitomap_mut_rna[row[0] + "_" + alleles[0] + "_" + alleles[1]] = row[6] # cell or tissue type
    except:
        logging.error("Reading " + mitomap_mut_somatic_file)
        print sys.exc_info()[0]
        raise
        sys.exit(1)
     
    # load MitImpact
    mitimpact = {}
    try:
        with open(mitimpact_file,'r') as f:
            reader = csv.reader(f,delimiter='\t')
            for row in reader:
                if row[0] != "Chr":
                    # replace all commoas with decimals
                    row = [r.replace(",",".") for r in row]
                    mitimpact[row[1]+"_"+row[2]+"_"+row[3]+"_"+row[4]] = row[5] + "," + row[6] + "," + row[7] + "," + row[8] + "," + row[9] + "," + row[10] + "," + row[11] + "," + row[12] + "," + row[13] + "," + row[14] + "," + row[15] + "," + row[16] + "," + row[17] + "," + row[18] + "," + row[19] # pos_ref_alt_gene -> AA_position,AA_ref,AA_alt,Codon_substitution,PolyPhen2_prediction,PolyPhen2_score,SIFT_prediction,SIFT_score,PROVEAN_prediction,PROVEAN_score,MutationAssessor_prediction,MutationAssessor_score,CADD_score,CADD_phred_score,CADD_prediction
    except:
        logging.error("Reading " + mitomap_mut_rna_file)
        print sys.exc_info()[0]
        raise
        sys.exit(1)
    
    
    # Annotate SNVs
    for variant in sobj.variants.snvs.keys():
        pos, ref, alt = variant.split("_")
        pos = int(pos) + 1
        # get Locus annotation
        for row in locus:
            locus_info = row.split("\t")
            if int(locus_info[1]) <= pos and int(locus_info[2]) >= pos:
                if variant in sobj.variants.variants_locus:
                    sobj.variants.variants_locus[variant] = sobj.variants.variants_locus[variant] + "," + locus_info[0]
                else:
                    sobj.variants.variants_locus[variant] = locus_info[0]
                if locus_info[5] == "tRNA":
                    sobj.variants.variants_locus_type = "tRNA"
                elif locus_info[5] == "rRNA":
                    sobj.variants.variants_locus_type = "rRNA"
                elif locus_info[5] == "protein_coding":
                    sobj.variants.variants_locus_type = "protein_coding"
                # if tRNA or rRNA then calculate RNAfold annotation
                if locus_info[5] == "tRNA" or locus_info[5] == "rRNA":
                    rna_sequence = ref_fasta.fetch(reference="chrM",start=(int(locus_info[1])-1),end=int(locus_info[2]))
                    sobj.variants.annotate_rna_variant(variant, locus_info[1], locus_info[2], rna_sequence)
        # get mitomap genbank allele freqs
        if str(pos) + "_" + ref + "_" + alt in mitomap_freq:
            sobj.variants.variants_genbank[variant] = mitomap_freq[str(pos) + "_" + ref + "_" + alt]
        # get codon number/amino acid change from mitomap poly coding
        if str(pos) + "_" + ref + "_" + alt in mitomap_poly_coding:
            #coding_annotation = mitomap_poly_coding[str(pos) + "_" + ref + "_" + alt].split("_")
            #effect = coding_annotation[1].split(":")
            #sobj.variants.variants_aminoacid_position[variant] = coding_annotation[0]
            #sobj.variants.variants_effect[variant] = effect[0]
            #if len(effect) > 1:
            #    sobj.variants.variants_aminoacid_change[variant] = effect[1]
            #else:
            #    sobj.variants.variants_aminoacid_change[variant] = "-"
            sobj.variants.variants_mitomap[variant] = 1;
        # nothing to get from mitomap poly control other than if it is in mitomap or not
        if str(pos) + "_" + ref + "_" + alt in mitomap_poly_control:
            sobj.variants.variants_mitomap[variant] = 1;
        # get disease info from mitomap mut coding/control annotation
        if str(pos) + "_" + ref + "_" + alt in mitomap_mut_coding_control:
            sobj.variants.variants_disease[variant] = mitomap_mut_coding_control[str(pos) + "_" + ref + "_" + alt]
        # get disease info from mitomap mut RNA annotation
        if ref + str(pos) + alt in mitomap_mut_rna:
            sobj.variants.variants_disease[variant] = mitomap_mut_rna[ref + str(pos) + alt]
        # get somatic cell or tissue annotation for tumor variants
        if str(pos) + "_" + ref + "_" + alt in mitomap_mut_somatic:
            sobj.variants.variants_somatic[variant] = mitomap_mut_somatic[str(pos) + "_" + ref + "_" + alt]
        # get annotations from MitImpact
        if variant in sobj.variants.variants_locus:
            for gene in sobj.variants.variants_locus[variant].split(","):
                if str(pos) + "_" + ref + "_" + alt + "_" + gene in mitimpact:
                    mitimpact_ann = mitimpact[str(pos) + "_" + ref + "_" + alt + "_" + gene].split(",")
                    # AA_position,AA_ref,AA_alt,Codon_substitution,PolyPhen2_prediction,PolyPhen2_score,SIFT_prediction,SIFT_score,PROVEAN_prediction,PROVEAN_score,MutationAssessor_prediction,MutationAssessor_score,CADD_score,CADD_phred_score,CADD_prediction
                    if variant in sobj.variants.variants_aminoacid_position:
                        sobj.variants.variants_aminoacid_position[variant] += "," + mitimpact_ann[0]
                        sobj.variants.variants_aminoacid_change[variant] += "," + mitimpact_ann[1]+"/"+mitimpact_ann[2]
                        sobj.variants.variants_codon_change[variant] += "," + mitimpact_ann[3]
                        sobj.variants.variants_polyphen2_pred[variant] += "," + mitimpact_ann[4]
                        sobj.variants.variants_polyphen2_score[variant] += "," + mitimpact_ann[5]
                        sobj.variants.variants_sift_pred[variant] += "," + mitimpact_ann[6]
                        sobj.variants.variants_sift_score[variant] += "," + mitimpact_ann[7]
                        sobj.variants.variants_provean_pred[variant] += "," + mitimpact_ann[8]
                        sobj.variants.variants_provean_score[variant] += "," + mitimpact_ann[9]
                        sobj.variants.variants_mutationassessor_pred[variant] += "," + mitimpact_ann[10]
                        sobj.variants.variants_mutationassessor_score[variant] += "," + mitimpact_ann[11]
                        sobj.variants.variants_cadd_score[variant] += "," + mitimpact_ann[13]
                    else:
                        sobj.variants.variants_aminoacid_position[variant] = mitimpact_ann[0]
                        sobj.variants.variants_aminoacid_change[variant] = mitimpact_ann[1]+"/"+mitimpact_ann[2]
                        sobj.variants.variants_codon_change[variant] = mitimpact_ann[3]
                        sobj.variants.variants_polyphen2_pred[variant] = mitimpact_ann[4]
                        sobj.variants.variants_polyphen2_score[variant] = mitimpact_ann[5]
                        sobj.variants.variants_sift_pred[variant] = mitimpact_ann[6]
                        sobj.variants.variants_sift_score[variant] = mitimpact_ann[7]
                        sobj.variants.variants_provean_pred[variant] = mitimpact_ann[8]
                        sobj.variants.variants_provean_score[variant] = mitimpact_ann[9]
                        sobj.variants.variants_mutationassessor_pred[variant] = mitimpact_ann[10]
                        sobj.variants.variants_mutationassessor_score[variant] = mitimpact_ann[11]
                        sobj.variants.variants_cadd_score[variant] = mitimpact_ann[13]
    
    # Annotate INDELs
    for variant in sobj.variants.indels.keys():
        pos, ref, alt = variant.split("_")
        pos = int(pos) + 1
        # get Locus annotation
        for row in locus:
            locus_info = row.split("\t")
            if int(locus_info[1]) <= pos and int(locus_info[2]) >= pos:
                if variant in sobj.variants.variants_locus:
                    sobj.variants.variants_locus[variant] = sobj.variants.variants_locus[variant] + "," + locus_info[0]
                else:
                    sobj.variants.variants_locus[variant] = locus_info[0]
                if locus_info[5] == "tRNA":
                    sobj.variants.variants_locus_type = "tRNA"
                elif locus_info[5] == "rRNA":
                    sobj.variants.variants_locus_type = "rRNA"
                elif locus_info[5] == "protein_coding":
                    sobj.variants.variants_locus_type = "protein_coding"
                # if tRNA or rRNA then calculate RNAfold annotation
                if locus_info[5] == "tRNA" or locus_info[5] == "rRNA":
                    rna_sequence = ref_fasta.fetch(reference="chrM",start=(int(locus_info[1])-1),end=int(locus_info[2]))
                    sobj.variants.annotate_rna_variant(variant, locus_info[1], locus_info[2], rna_sequence)
        if len(ref) > len(alt): # DEL
            # get mitomap genbank allele freqs
            if str(pos) + "_" + ref[len(alt):] + "_" + ":" in mitomap_freq:
                sobj.variants.variants_genbank[variant] = mitomap_freq[str(pos) + "_" + ref[len(alt):] + "_" + ":"]
            # get codon number/amino acid change from mitomap poly coding
            if str(pos) + "_" + ref[len(alt):] + "_" + "del" in mitomap_poly_coding:
                coding_annotation = mitomap_poly_coding[str(pos) + "_" + ref[len(alt):] + "_" + "del"].split("_")
                effect = coding_annotation[1].split(":")
                sobj.variants.variants_aminoacid_position[variant] = coding_annotation[0]
                sobj.variants.variants_effect[variant] = effect[0]
                if len(effect) > 1:
                    sobj.variants.variants_aminoacid_change[variant] = effect[1]
                else:
                    sobj.variants.variants_aminoacid_change[variant] = "-"
                sobj.variants.variants_mitomap[variant] = 1;
            # nothing to get from mitomap poly control other than if it is in mitomap or not
            if str(pos) + "_" + ref[len(alt):] + "_" + "del" in mitomap_poly_control:
                sobj.variants.variants_mitomap[variant] = 1;
            # get disease info from mitomap mut coding/control annotation
            if str(pos) + "_" + ref[len(alt):] + "_" + "del" in mitomap_mut_coding_control:
                sobj.variants.variants_disease[variant] = mitomap_mut_coding_control[str(pos) + "_" + ref[len(alt):] + "_" + "del"]
            # get disease info from mitomap mut RNA annotation
            if ref + str(pos) + ":" in mitomap_mut_rna:
                sobj.variants.variants_disease[variant] = mitomap_mut_rna[ref + str(pos) + ":"]
            # get somatic cell or tissue annotation for tumor variants
            if str(pos) + "_" + ref[len(alt):] + "_" + "del" in mitomap_mut_somatic:
                sobj.variants.variants_somatic[variant] = mitomap_mut_somatic[str(pos) + "_" + ref[len(alt):] + "_" + "del"]
                
            
        else: # INS
            # get mitomap genbank allele freqs
            if str(pos) + "_" + ref + "_" + alt in mitomap_freq:
                sobj.variants.variants_genbank[variant] = mitomap_freq[str(pos) + "_" + ref + "_" + alt]
            # get codon number/amino acid change from mitomap poly coding
            if str(pos) + "_" + ref + "_" + alt in mitomap_poly_coding:
                coding_annotation = mitomap_poly_coding[str(pos) + "_" + ref + "_" + alt].split("_")
                effect = coding_annotation[1].split(":")
                sobj.variants.variants_aminoacid_position[variant] = coding_annotation[0]
                sobj.variants.variants_effect[variant] = effect[0]
                if len(effect) > 1:
                    sobj.variants.variants_aminoacid_change[variant] = effect[1]
                else:
                    sobj.variants.variants_aminoacid_change[variant] = "-"
                sobj.variants.variants_mitomap[variant] = 1;
            # nothing to get from mitomap poly control other than if it is in mitomap or not
            if str(pos) + "_" + ref + "_" + alt in mitomap_poly_control:
                sobj.variants.variants_mitomap[variant] = 1;
            # get disease info from mitomap mut coding/control annotation
            if str(pos) + "_" + ref + "_" + alt in mitomap_mut_coding_control:
                sobj.variants.variants_disease[variant] = mitomap_mut_coding_control[str(pos) + "_" + ref + "_" + alt]
            # get disease info from mitomap mut RNA annotation
            if ref + str(pos) + alt in mitomap_mut_rna:
                sobj.variants.variants_disease[variant] = mitomap_mut_rna[ref + str(pos) + alt]
            # get somatic cell or tissue annotation for tumor variants
            if str(pos) + "_" + ref + "_" + alt in mitomap_mut_somatic:
                sobj.variants.variants_somatic[variant] = mitomap_mut_somatic[str(pos) + "_" + ref + "_" + alt]
    
    

def assign_haplogroups(sobj, output, java, haplogrep, phylotree, min_heteroplasmy, rcrs_len):
    ''' Assign haplogroups for a sample. '''
    logging.info("Assigning haplogroup for sample: " + sobj.sname)
    f = open(output+"/variants/haplogroups/"+sobj.sname+".haplogrep.in.hsd","w")
    f.write("SampleId\tRange\tHaplogroup\tPolymorphisms\n")
    variants = sobj.variants.get_phylotree_format(min_heteroplasmy)
    f.write(sobj.sname+"\t\"1-"+str(rcrs_len)+"\"\t?\t"+"\t".join(variants)+"\n")
    f.close()
    
    # run haplogrep
    try:
        haplogrep_cmd = [java,"-Xmx3g","-Xms2g","-jar",haplogrep,"--in",output+"/variants/haplogroups/"+sobj.sname+".haplogrep.in.hsd","--out",output+"/variants/haplogroups/"+sobj.sname+".haplogrep.out.txt","--phylotree","17","--format","hsd"]
        subprocess.call(haplogrep_cmd)
    except:
        logging.error("Unable to run Haplogrep for sample " + sobj.sname)
        sys.exit(1)

    # read in phylotree polymorphisms
    phylotree_polymorphisms = {}
    hotspot_mutations = {}
    try:
        with open(phylotree,'r') as f:
            reader = csv.reader(f,delimiter='\t',quoting=csv.QUOTE_NONE)
            for row in reader:
                if row[0] == "haplogroup":
                    # skip header
                    continue 
                mutations = row[2].split(",")
                if row[0] == "hotspot_mutations":
                    for mut in mutations:
                        hotspot_mutations[mut] = 1
                else:
                    for mut in mutations:
                        phylotree_polymorphisms[mut] = 1
                    
    except:
        logging.error("Reading file " + phylotree)
        logging.error(sys.exc_info()[0])
        sys.exit(1)

    # read in haplogrep output
    variants_haplogrep = {}
    found_polys = []
    not_found_polys = []
    try:
        with open(output+"/variants/haplogroups/"+sobj.sname+".haplogrep.out.txt",'r') as f:
            reader = csv.reader(f,delimiter='\t',quoting=csv.QUOTE_NONE)
            for row in reader:
                if row[0] == sobj.sname:
                    sobj.haplogroup = row[2]
                    sobj.haplogroup_score = row[3]
                    not_found_polys = row[4].split(" ")
                    found_polys = row[5].split(" ")
                    # filter any empty strings from arrays due to leading spaces
                    not_found_polys = filter(None,not_found_polys)
                    found_polys = filter(None,found_polys)
    except:
        logging.error("Reading file " + output+"/variants/haplogroups/"+sobj.sname+".haplogrep.out.txt")
        logging.error(sys.exc_info()[0])
        sys.exit(1)
        
    if len(found_polys) > 0:
        sobj.haplogroup_yes_total = len(found_polys)
    if len(not_found_polys) > 0:
        sobj.haplogroup_no_total = len(not_found_polys)
        
    # annotate all mutations
    for variant in sobj.variants.variants_phylotree.keys():
        phylotree_format = sobj.variants.variants_phylotree[variant]
        haplogrep_annotation = ""
        if phylotree_format in hotspot_mutations:
            variants_haplogrep[variant] = "hotspot"
        elif phylotree_format in found_polys:
            variants_haplogrep[variant] = "yes"
        elif phylotree_format in phylotree_polymorphisms:
            variants_haplogrep[variant] = "localPrivateMut"
            sobj.haplogroup_lpm_total += 1
        else:
            variants_haplogrep[variant] = "globalPrivateMut"
            sobj.haplogroup_gpm_total += 1
             
    sobj.variants.variants_haplogrep = variants_haplogrep
    #print sobj.variants.variants_haplogrep
    f = open(output+"/variants/haplogroups/"+sobj.sname+".haplogroup_summary.csv","w")
    f_writer = csv.writer(f, delimiter=",", quoting=csv.QUOTE_MINIMAL)
    f_writer.writerow(["Sample","Haplogroup","Quality_Score","Polymorphisms_Found","Polymorphisms_Missing","Global_Private_Mutations","Local_Private_Mutations"])
    f_writer.writerow([sobj.sname, sobj.haplogroup, str(sobj.haplogroup_score), str(sobj.haplogroup_yes_total), str(sobj.haplogroup_no_total), str(sobj.haplogroup_gpm_total), str(sobj.haplogroup_lpm_total)])
    f.close()

def find_contig_name(samfile):
    ''' Find which MT contig name is being using in the BAM '''
    possible_mt_names = ["chrM","chrMT","M","MT","rCRS","NC_012920"]
    mt_name = ""
    for name in possible_mt_names:
        if samfile.get_tid(name) >= 0:
            mt_name = name
    if mt_name == "":
        logging.error("Couldn't find mitochondrial contig name in BAM " + sobj.in_file[0] + " among commonly used names: " + possible_mt_names)
        sys.exit(1)
    return mt_name


def call_variants(sobj, output, min_het, min_mq, min_bq, max_sb, mito_ref, rcrs_len, max_pileup_depth):
    ''' Call variants using the pileup format. '''
    logging.info("Calling variants for sample: " + sobj.sname)
    # open reference FASTA
    try:
        ref_fasta = pysam.FastaFile(filename=mito_ref)
        #print ref_fasta.fetch(reference="chrM",start=0,end=1)
    except:
        logging.error("Opening reference FASTA " + mito_ref + " " + str(sys.exc_info()[0]))
        sys.exit(1)
        
    samfile = pysam.AlignmentFile(sobj.in_file[0], "rb")
    # check if index is present and if not try to index bam
    if not samfile.has_index():
        logging.warning("BAM file " + sobj.in_file[0] + " is not indexed. Will attempt to index file...")
        try:
            pysam.index(sobj.in_file[0])
        except:
            logging.error("Unable to index " + sobj.in_file[0])
            sys.exit(1)
            
    # find the MT contig name being used in the BAM
    contig = find_contig_name(samfile)

    # check if contig is expected length of rCRS
    if samfile.lengths[samfile.get_tid(contig)] != rcrs_len:
        logging.error("Mitochondrial contig in BAM " + sobj.in_file[0] + " does not match expected rCRS length of " + rcrs_len)
        sys.exit(1)

    

    # calculate coverage
#    for i in range(0,(rcrs_len/COVERAGE_RESOLUTION)+1):
#        sobj.coverage.append(samfile.count(reference=contig,start=i*COVERAGE_RESOLUTION,end=(i*COVERAGE_RESOLUTION)+1,read_callback="all"))
    
    # retrieve and store pileup
    for pileupcolumn in samfile.pileup(reference=contig,max_depth=max_pileup_depth,stepper="all"):
        ref_nucl = ref_fasta.fetch(reference="chrM",start=pileupcolumn.reference_pos,end=pileupcolumn.reference_pos+1)
        indels = {}
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_refskip:
                nucl_index = pileupread.query_position_or_next - pileupread.alignment.query_alignment_start
                if pileupread.alignment.mapping_quality >= min_mq:
                    forward_strand = True
                    if pileupread.alignment.is_reverse:
                        forward_strand = False
                    if pileupread.is_del and pileupread.alignment.query_alignment_qualities[nucl_index] >= min_bq:
                        sobj.variants.pileup_increment_read_count(pileupcolumn.reference_pos, 'D')
                    elif pileupread.alignment.query_alignment_qualities[nucl_index] >= min_bq:
                        # SNV
                        nucl = pileupread.alignment.query_alignment_sequence[nucl_index]
                        nucl.upper()
                        if nucl != 'N':
                            sobj.variants.pileup_increment_read_count(pileupcolumn.reference_pos,nucl,forward_strand)
                    if pileupread.indel > 0:
                        # INS
                        ref = ref_nucl
                        alt = pileupread.alignment.query_alignment_sequence[nucl_index:nucl_index+pileupread.indel+1]
                        # check base qualities of all bases in the INS
                        if min(pileupread.alignment.query_alignment_qualities[nucl_index:nucl_index+pileupread.indel+1]) >= min_bq:
                            sobj.variants.pileup_indel_increment_read_count(pileupcolumn.reference_pos,ref,alt,forward_strand)
                            indels[ref + "_" + alt] = 1
                    elif pileupread.indel < 0:
                        # DEL
                        ref = ref_fasta.fetch(reference="chrM",start=pileupcolumn.reference_pos,end=pileupcolumn.reference_pos-pileupread.indel+1)
                        alt = ref_nucl
                        # check base quality of leading nucleotide
                        if pileupread.alignment.query_alignment_qualities[nucl_index] >= min_bq:
                            sobj.variants.pileup_indel_increment_read_count(pileupcolumn.reference_pos,ref,alt,forward_strand)
                            indels[ref + "_" + alt] = 1
                    # store largest read length for this sample
                    if pileupread.alignment.query_length > sobj.read_length:
                        sobj.read_length = pileupread.alignment.query_length
                        
        # go back and fill in ref supporting reads for INDELs
        for indel_key in indels.keys():
            ref, alt = indel_key.split("_")
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_refskip:
                    nucl_index = pileupread.query_position_or_next - pileupread.alignment.query_alignment_start
                    if pileupread.alignment.mapping_quality >= min_mq and nucl_index+len(ref) <= pileupread.alignment.reference_end:
                        forward_strand = True
                        if pileupread.alignment.is_reverse:
                            forward_strand = False
                        if ref == pileupread.alignment.query_alignment_sequence[nucl_index:nucl_index+len(ref)] and min(pileupread.alignment.query_alignment_qualities[nucl_index:nucl_index+len(ref)]) >= min_bq and pileupread.indel == 0:
                            sobj.variants.pileup_indel_increment_read_count(pileupcolumn.reference_pos,ref,alt,forward_strand,False)

        # call SNVs
        depth = float(sobj.variants.get_depth(pileupcolumn.reference_pos))
        if depth > 0:
            for n in ["A","C","G","T"]:
                n_support = sobj.variants.get_nucleotide_support(pileupcolumn.reference_pos,n)
                #depth = float(sobj.variants.get_depth(pileupcolumn.reference_pos))
                if n != ref_nucl and (n_support/depth) >= min_het and ref_nucl != 'N':
                    sb = sobj.variants.strand_bias(pileupcolumn.reference_pos,ref_nucl,n)
                    # only worry about strand bias with low heteroplasmy
                    if (n_support/depth) < 0.10 and sb > max_sb:
                        continue
                    sobj.variants.add_variant(pileupcolumn.reference_pos, ref_nucl, n, (n_support/depth))
                    #print ref_nucl+str(pileupcolumn.reference_pos)+n+"\t"+str(n_support/depth)+"\t"+str(sb)
                    #print sobj.pileup[pileupcolumn.reference_pos]
        
    # call INDELs
    for indel, support in sobj.variants.pileup_indel.items():
        pos, ref, alt = indel.split("_")
        pos = int(pos)
        depth = float(sobj.variants.get_depth(pos))
        if depth > 0:
            sb = 0.0
            sb = sobj.variants.strand_bias(pos,ref,alt)
            #if len(ref) == 1 and len(alt) > 1:
                # INS
            #    sb = sobj.variants.strand_bias(pos,ref,alt)
            #else:
                # DEL
            #    sb = sobj.variants.strand_bias(pos,ref,alt)
            alt_support = support[2] + support[3]
            if alt_support/depth >= min_het:
                if alt_support/depth < 0.10 and sb > max_sb:
                        continue
                sobj.variants.add_variant(pos, ref, alt, (alt_support/depth))
                #print pos+"\t"+ref+"\t"+alt+"\t"+str(sum(support)/depth)+"\t"+str(sb)
                #print support
    
    # calculate genotypes fom variant calls
    sobj.variants.calculate_genotypes(min_het)

    # calculate relative copy number for WGS or exome data using method from MitoCounter publication
    idxstats_rows = []
    if sobj.idxstats_file:
        # option to pass in idxstats file in sample_info
        try:
            with open(sobj.idxstats_file,'r') as f:
                reader = csv.reader(f,delimiter='\t',quoting=csv.QUOTE_NONE)
                for row in reader:
                    idxstats_rows.append("\t".join(row))
        except:
            logging.error("Can't read file: " + sobj.idxstats_file)
            sys.exit(1)
    else:
        idxstats_rows = pysam.idxstats(sobj.in_file[0]).strip().split("\n")
    mt_mapped_reads = 0.0
    nuclear_mapped_reads = 0.0
    nuclear_genome_size = 0.0
    for row in idxstats_rows:
        read_counts = row.split("\t")
        if read_counts[0] == contig:
            mt_mapped_reads += int(read_counts[2])
        else:
            nuclear_mapped_reads += int(read_counts[2])
            nuclear_genome_size += int(read_counts[1])
    # only calculate the relative copy number if there are 10 times more nuclear reads
    if nuclear_mapped_reads > (mt_mapped_reads*10):
        sobj.relative_copy_number = (mt_mapped_reads*sobj.read_length*2*nuclear_genome_size)/((nuclear_mapped_reads*sobj.read_length)*rcrs_len)
    sobj.mitochondrial_reads = int(mt_mapped_reads)
    sobj.nuclear_reads = int(nuclear_mapped_reads)


def align_reads(sobj, output_dir, alignment_ref, bwa, samtools):
    ''' Align fastqs to reference sequence '''
    logging.info("Aligning Reads for sample: " + sobj.sname)
    rg = "@RG\\tID:"+sobj.sname+"\\tSM:"+sobj.sname+"\\tLB:"+sobj.sname+"\\tPL:Illumina\\tCN:Mayo"
    fastq1 = sobj.in_file[0]
    fastq2 = ""
    if len(sobj.in_file) > 1:
        fastq2 = sobj.in_file[1]
    aligned_bam = output_dir+"/"+sobj.sname+".unsorted.bam"
    sorted_bam = output_dir + "/" + sobj.sname + ".bam"
    # align with BWA-MEM
    try:
        f = open(aligned_bam, "w")
        bwa_cmd = [bwa,"mem","-M","-R",rg,alignment_ref,fastq1]
        if len(sobj.in_file) > 1:
            bwa_cmd = [bwa,"mem","-M","-R",rg,alignment_ref,fastq1,fastq2]
        samtools_cmd = [samtools,"view","-bS","-"]
        bwa_process = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE)
        samtools_process = subprocess.Popen(samtools_cmd, stdin=bwa_process.stdout, stdout=f)
        samtools_process.wait()
        f.close()
    except:
        logging.error("Unable to align fastqs: " + fastq1 + " " + fastq2)
        logging.error(sys.exc_info()[0])
        sys.exit(1)
    # sort with samtools
    try:
        pysam.sort("-o", sorted_bam, aligned_bam)
        pysam.index(sorted_bam)
    except:
        logging.error("Unable to sort and index " + sorted_bam)
        logging.error(sys.exc_info()[0])
        sys.exit(1)

    try:
        os.remove(aligned_bam)
    except:
        logging.error("Unable to remove " + aligned_bam)
        logging.error(sys.exc_info()[0])
        sys.exit(1)

    sobj.in_file = []
    sobj.in_file.append(sorted_bam)

def write_variants(sobj, output_dir):
    ''' write variants to vcf and annotation table files '''
    sobj.write_vcf(output_dir)
    sobj.write_annotation_tables(output_dir)
    sobj.write_pileup(output_dir)


def run_sample_analysis(sobj, config, output_dir, wgs):
    ''' Run the necessary methods for per-sample variant calling, annotation, and coverage '''
    if not sobj.is_aligned:
        alignment_ref = config.get("REF","MITO_REF")
        if wgs:
            alignment_ref = config.get("REF","GENOME_REF")
        align_reads(sobj, output_dir + "/alignment", alignment_ref, config.get("TOOLS","BWA"), config.get("TOOLS","SAMTOOLS"))
    else:
        if not wgs:
            # copy bams to output if not large WGS files
            try:
                copy(sobj.in_file[0],output_dir + "/alignment/")
                if os.path.isfile(sobj.in_file[0] + ".bai"):
                    copy(sobj.in_file[0] + ".bai",output_dir + "/alignment/")
            except:
                logging.error("Unable to copy bam to output folder " + sobj.in_file[0])
                logging.error(sys.exc_info()[0])
                sys.exit(1)
    
    call_variants(sobj,output_dir,
                  config.getfloat("PARAMS","HETEROPLASMY_FRACTION"),
                  config.getint("PARAMS","MAPPING_QUALITY"),
                  config.getint("PARAMS","BASE_QUALITY"),
                  config.getint("PARAMS","STRAND_BIAS_FILTER"),
                  config.get("REF","MITO_REF"),
                  config.getint("PARAMS", "RCRS_LEN"),
                  config.getint("PARAMS", "MAX_PILEUP_DEPTH"))
    
    assign_haplogroups(sobj,output_dir,
                   config.get("TOOLS","JAVA"),
                   config.get("TOOLS","HAPLOGREP"),
                   config.get("REF","PHYLOTREE"),
                   config.getfloat("PARAMS","MIN_HET_FOR_HAPLOGROUPING"),
                   config.getint("PARAMS", "RCRS_LEN"))
    
    annotate_variants(sobj,
                      config.get("REF","MITO_GENES"),
                      config.get("REF","MITOMAP_COUNTS"),
                      config.get("REF","MITOMAP_POLY_CODING"),
                      config.get("REF","MITOMAP_POLY_CONTROL"),
                      config.get("REF","MITOMAP_MUT_CODING_CONTROL"),
                      config.get("REF","MITOMAP_MUT_RNA"),
                      config.get("REF","MITOMAP_MUT_SOMATIC"),
                      config.get("REF","MITIMPACT"),
                      config.get("REF","MITO_REF"))
    
    write_variants(sobj,output_dir)
    
    calculate_coverage(sobj,output_dir,
                       config.get("REF","MITO_GENES_PLOT"),
                       config.getboolean("PARAMS","LOG_SCALE_COVERAGE_PLOT"),
                       config.getint("PARAMS", "COVERAGE_RESOLUTION"),
                       config.getint("PARAMS", "RCRS_LEN"))



def cluster_submit_sample_analysis(sobj, config, output_dir, wgs):
    '''
    Submit the sample analysis job and return job id
    '''
    # write a shell script
    mitosort_path = config.get("TOOLS","MITOSORT")
    python = config.get("TOOLS","PYTHON")
    matplotlib_path = config.get("TOOLS", "MATPLOTLIB")
    pysam_path = config.get("TOOLS", "PYSAM")
    viennarna_path = config.get("TOOLS", "VIENNARNA") 
    shell_script = output_dir + "/docs/scripts/mito_analysis." + sobj.sname + ".sh"
    f = open(shell_script,"w")
    f.write("#!/bin/sh\n\n")
    f.write("echo $(date)\n")
    #f.write("PYTHONPATH=" + mitosort_path + ":" + matplotlib + ":$PYTHONPATH \n")
    f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/:/usr/local/biotools/python/2.7.10/lib/ \n")
    f.write("export PYTHONPATH=" + mitosort_path + ":" + matplotlib_path + ":" + pysam_path + ":" + viennarna_path + "\n")
    f.write(python + " -c 'import analysis; analysis.cluster_run_sample_analysis(\"" + sobj.sname + "\",\"" + output_dir + "\",\"" + str(sobj.is_ipsc) + "\",\"" + str(sobj.is_somatic) + "\",\"" + str(wgs) + "\") '\n")
    #f.write(python + " -c 'import parallel; parallel.cluster_run_sample_analysis(\"" + sobj.sname + "\",\"" + output_dir + "\",\"" + str(sobj.is_ipsc) + "\",\"" + str(sobj.is_somatic) + "\",\"" + str(wgs) + "\") '\n")
    f.write("echo $(date)\n")
    f.close()

    qsub = config.get("CLUSTER","QSUB")
    queue = config.get("CLUSTER","QUEUE")
    memory = config.get("CLUSTER","MEMORY")
    email = config.get("CLUSTER","EMAIL")

    #qsub_cmd = [qsub,"-V","-wd",output_dir+"/docs/logs/","-q",queue,"-M",email,"-m","a","-l",memory,"-l","h_stack=10M",shell_script]
    qsub_cmd = [qsub, "-wd",output_dir+"/docs/logs/","-q",queue,"-M",email,"-m","a","-l",memory,"-l","h_stack=10M",shell_script]
    qsub_output = subprocess.check_output(qsub_cmd)
    print qsub_output.strip()
    job_id = qsub_output.split(" ")[2]

    return job_id


def cluster_run_sample_analysis(sample_name, output_dir, ipsc, somatic, wgs):
    '''
    Run the analysis for a sample, to be called from the SGE job
    '''
    # load necessary objects
    logging.basicConfig(format='%(levelname)s\t%(asctime)s - %(message)s', level=logging.INFO)
    config_file = output_dir + "/docs/config/config.txt"
    sample_info_file = output_dir + "/docs/config/sample_info.txt"
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)
    
    ipsc = bool(ipsc)
    somatic = bool(somatic)
    wgs = bool(wgs)

    sample_info = parsing.parse_sample_info(sample_info_file, ipsc, somatic)
    this_sobj = None
    for s in sample_info:
        sobj = Sample()
        sobj.load_sample_info(s, ipsc, somatic)
        if sobj.sname == sample_name:
            this_sobj = sobj
            break

    # run usual analysis
    run_sample_analysis(this_sobj, config, output_dir, wgs)
    # dump necessary objects to disk for use in the next job
    #this_sobj.dump_variants_to_file(output_dir)
    output_stream = open(output_dir + "/docs/scripts/" + this_sobj.sname + ".sample.pkl", 'wb')
    pickle.dump(this_sobj, output_stream)
    output_stream.close()
    
def cluster_submit_allsample_reports(sample_objects, config, output_dir, group_analysis, job_ids):
    '''
    submit the allsample reports SGE job
    ''' 
    # write a shell script
    mitosort_path = config.get("TOOLS","MITOSORT")
    python = config.get("TOOLS","PYTHON")
    matplotlib_path = config.get("TOOLS", "MATPLOTLIB")
    pysam_path = config.get("TOOLS", "PYSAM")
    viennarna_path = config.get("TOOLS", "VIENNARNA") 
    ipsc = False
    somatic = False
    pedigree = False
    if group_analysis:
        for sobj in sample_objects:
            if sobj.is_ipsc:
                ipsc = True
            if sobj.is_somatic:
                somatic = True
            if sobj.is_pedigree:
                pedigree = True
    
    shell_script = output_dir + "/docs/scripts/mito_analysis.allsample_reports.sh"
    f = open(shell_script,"w")
    f.write("#!/bin/sh\n\n")
    f.write("echo $(date)\n")
    f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/:/usr/local/biotools/python/2.7.10/lib/ \n")
    f.write("export PYTHONPATH=" + mitosort_path + ":" + matplotlib_path + ":" + pysam_path + ":" + viennarna_path + "\n")
    f.write(python + " -c 'import analysis; analysis.cluster_run_allsample_reports(\"" + output_dir + "\",\"" + str(ipsc) + "\",\"" + str(somatic) + "\",\"" + str(pedigree) + "\") '\n")
    #f.write(python + " -c 'import parallel; parallel.cluster_run_allsample_reports(\"" + output_dir + "\",\"" + str(ipsc) + "\",\"" + str(somatic) + "\",\"" + str(pedigree) + "\") '\n")
    f.write("echo $(date)\n")
    f.close()

    qsub = config.get("CLUSTER","QSUB")
    queue = config.get("CLUSTER","QUEUE")
    memory = config.get("CLUSTER","MEMORY")
    email = config.get("CLUSTER","EMAIL")

    s = ","
    hold_jobs = s.join(job_ids)

    qsub_cmd = [qsub, "-wd",output_dir+"/docs/logs/","-q",queue,"-hold_jid",hold_jobs,"-M",email,"-m","a","-l",memory,"-l","h_stack=10M",shell_script]
    qsub_output = subprocess.check_output(qsub_cmd)
    print qsub_output.strip()
    
    
def cluster_run_allsample_reports(output_dir, ipsc, somatic, pedigree):
    '''
    run the allsample reports code, called from an SGE job
    '''
    # load necessary objects
    logging.basicConfig(format='%(levelname)s\t%(asctime)s - %(message)s', level=logging.INFO)
    config_file = output_dir + "/docs/config/config.txt"
    sample_info_file = output_dir + "/docs/config/sample_info.txt"
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)
    
    group_analysis = False
    
    ipsc = bool(ipsc)
    somatic = bool(somatic)
    pedigree = bool(pedigree)

    if ipsc | somatic | pedigree:
        group_analysis = True
    
    sample_info = parsing.parse_sample_info(sample_info_file, ipsc, somatic)
    sample_objects = []
    for s in sample_info:
        sobj = Sample()
        sobj.load_sample_info(s, ipsc, somatic)
        # load object data from pickle dump
        #sobj.load_variants_from_file(output_dir)
        output_stream = open(output_dir + "/docs/scripts/" + sobj.sname + ".sample.pkl", 'rb')
        sobj = pickle.load(output_stream)
        output_stream.close()
        os.remove(output_dir + "/docs/scripts/" + sobj.sname + ".sample.pkl")
        sample_objects.append(sobj)
    
    # run usual analysis from mito module
    create_allsample_reports(sample_objects, config, output_dir, group_analysis, ipsc)


def create_allsample_reports(sample_objects, config, output_dir, group_analysis, ipsc):
    ''' Create final reports using results from all samples, including final html '''
    
    backfill_mutations(sample_objects, output_dir, config.getfloat("PARAMS", "MIN_HET_TO_BACKFILL"))
    if group_analysis:
        # plot per-group priv mutation plots  

        # rscript = config.get("TOOLS","RSCRIPT")
        # rlib = config.get("TOOLS","RLIB")
        
        gene_file = config.get("REF","MITO_GENES_PLOT")
        groups = {}

        for sobj in sample_objects:
            groups[sobj.group] = 1
        for group in groups.keys():
            # plot_cmd = [rscript, rlib + "/plot_group_mutations.R", output_dir + "/variants/plot/" + group + ".private_mutations.txt", output_dir + "/variants/plot/" + group + ".private_mutations.png", group, gene_file]
            # subprocess.call(plot_cmd)
            plotting.plotMutations(output_dir + '/variants/plot/' + group + '.private_mutations.txt', 
                                               output_dir + '/variants/plot/' + group + '.private_mutations.png',
                                               group, gene_file)
    # if iPSC then score and write iPSC qualities
    # look into making this write into its own function.
    if ipsc:
        ipsc_annotation(sample_objects)
        with open(output_dir + "/iPSC_Risk.csv","w") as ipsc_file:
            ipsc_writer = csv.writer(ipsc_file, delimiter=",", quoting=csv.QUOTE_MINIMAL)
            ipsc_writer.writerow(["Line","Clone","Passage_Number","Quality_Category","Quality_Score","Highest_Risk_Variant"])
            for sobj in sample_objects:
                ipsc_writer.writerow([sobj.group, sobj.clone, sobj.passage_number, sobj.ipsc_risk_category, sobj.ipsc_risk_score, sobj.ipsc_highest_risk_variant])
        plotting.plot_risk(output_dir + "/iPSC_Risk.csv", output_dir + "/iPSC_risk.png")

    # Copy HTML supporting files
    html_files = config.get("TOOLS","HTML_FILES")
    destination = output_dir + "/docs/html"
    try:
        copytree(html_files, destination)
    except:
        logging.error("Couldn't copy html files " + html_files + " to " + destination)
        logging.error(sys.exc_info()[0])
    help_document = config.get("TOOLS","MITOSORT") + "/docs/MitoSort_help_manual.pdf"
    try:
        copy(help_document, output_dir + "/docs/MitoSort_help_manual.pdf")
    except:
        logging.error("Couldn't copy help document " + help_document)
        logging.error(sys.exc_info()[0])
    logging.info("Generating final HTML report")
    # create main html report
    html.write_html_report(sample_objects, config, output_dir, output_dir + "/mito_report.html", output_dir + "/docs/html/")
