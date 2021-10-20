import os
import sys
import csv
import glob
import gzip 
import pysam
import numpy as np
import pandas as pd
import dask
import dask.dataframe as dd

from find_canonical import gene_to_canonical_transcript
from synonAnnotationFunction import synonAnnotate

configfile: "config/config.json"

perl_lib = config["perl_lib"]
perl_exec = config["perl_exec"] # path, not the file itself
vep_exec = config["vep_exec"]
vep_cache = config["vep_cache"]
vep_plugins = config["vep_plugins"]
ancestral = config["ancestral"]
dbNSFP = config["dbNSFP"]
dbNSFP_fields = config["dbNSFP_fields"]
revel = config["revel"]
cov_file = config["cov_file"]
cov_dir = config["cov_dir"]
ref_dir = config["ref_dir"]
GIAB_file = config["GIAB_file"]
gnomad_file = config["gnomad_file"]
data_dir = config["data_dir"]
rate_dir_v4 = config["rate_dir"]
region_file = config["region_file"]
region_type = config["region_type"]

base_set = ["A", "C", "T", "G"]
chrom_set = [str(x) for x in range(1, 23)]

GPCR_file = "/net/home/emkoch/GPCR/DATA/GPCR_gene_info.tsv"
GPCR_table = pd.read_csv(GPCR_file, sep="\t")
GPCR_gene_ids = set(GPCR_table.gene_id)

## NOTE: Put this into config file!!
carlson_rate_folder = '/net/home/vseplyarsky/GNOMAD_Tufts/data/Carlson/by_chr/'

def check_GPCR_ids(fname):
    for GPCR_gene_id in GPCR_gene_ids:
        if GPCR_gene_id in fname:
            return True
    return False

def flip_nt(nt):
    if nt.upper() == "A":
        flip_nt = "T"
    elif nt.upper() == "T":
        flip_nt = "A"
    elif nt.upper() == "C":
        flip_nt = "G"
    elif nt.upper() == "G":
        flip_nt = "C"
    else:
        print(nt, " not a real nucleotide!")
        flip_nt = "N"
    return flip_nt

def nt_to_CT(nt):
    if nt.upper() == "A":
        return "T"
    if nt.upper() == "G":
        return "C"
    return nt.upper()

# Here the reference nt is for the genome not Vova's context
# If the reference is A or G, then Vova's C/T based reference
# will need flipping
def mut_to_CT(ref_nt, mut_nt):
    if ref_nt.upper() == "A" or ref_nt.upper() == "G":
        return flip_nt(mut_nt)
    return mut_nt.upper()

def get_unique_names(region_file):
    names = []
    with open(region_file , "r") as f_in:
        region_reader = csv.reader(f_in, delimiter="\t")
        next(region_reader)
        for region in region_reader:
            name = region[0]+"_"+region[3]
            if not name in names:
                names.append(name)
    return names

region_names = get_unique_names(region_file)
print(region_names)

# split_fnames = glob.glob(os.path.join(data_dir, "gene_splits/fake_transcript_variants_v4_*_processed.tsv.gz"))
# split_fnames = [split_fname.split("_processed")[0] for split_fname in split_fnames]
# GPCR_fnames = [split_fname for split_fname in split_fnames if check_GPCR_ids(split_fname)]

# cov_fnames = glob.glob(os.path.join(data_dir, "coverage", "gnomad.exomes.coverage.summary.*"))

rule all:
    input:
        os.path.join(data_dir, "fake_{region_type}_variants_syn.tsv.gz".format(region_type=region_type))
        # [os.path.join(data_dir, "region_splits/fake_{region_type}_variants_".format(region_type=region_type) +
        #               region_name + "_syn.tsv.gz")
        #  for region_name in region_names]
        # os.path.join(data_dir, "fake_{region_type}_variants_all.tsv.gz".format(region_type=region_type))
        #os.path.join(data_dir, "fake_full_transcript_variants_v4.tsv.gz"),
        #os.path.join(data_dir, "fake_transcript_variants_v4.tsv.gz")
        # os.path.join(data_dir,"gene_splits","fake_transcript_variants_v4_ENSG00000109927_ENST00000392793.6_giab.tsv.gz")
        #os.path.join(data_dir, "fake_transcript_variants_v4_all.tsv.gz"),
        #os.path.join(data_dir, "fake_transcript_variants_v4_GPCR.tsv.gz"),
        # [split_fname + "_vep.tsv.gz" for split_fname in split_fnames]
#        os.path.join(data_dir, "fake_transcript_variants_sorted_all_syn_counts_v5.gz"),
        # os.path.join(data_dir, "fake_transcript_variants_sorted_v3_aa_syn_GIAB_v5.gz")
        # os.path.join(data_dir, "fake_transcript_variants_syn_v5.tsv.gz"),
        # os.path.join(data_dir, "fake_transcript_variants_ms_v5.tsv.gz")
        # os.path.join(data_dir, "fake_transcript_variants_sorted_all_ms_counts_v5.gz")
        #os.path.join(data_dir, "coverage", "gnomad.exomes.lifted.coverage.summary.tsv.gz")

rule lift_coverage:
    input:
        cov_file
    output:
        os.path.join(cov_dir, "gnomad.exomes.lifted.coverage.summary.tsv.gz")
    run:
        from pyliftover import LiftOver
        lo = LiftOver("hg19", "hg38")
        with gzip.open(input[0], "rt") as f_in, gzip.open(output[0], "wt", newline="") as f_out:
            in_reader = csv.reader(f_in, delimiter="\t")
            out_writer = csv.writer(f_out, delimiter="\t")
            next(in_reader)
            for row in in_reader:
                chrom = row[0]
                try:
                    chrom_num = int(chrom)
                except ValueError:
                    break
                pos = row[1]
                print(chrom, pos)
                new_coords = lo.convert_coordinate("chr" + str(chrom), int(pos), "+")
                if len(new_coords) != 0:
                    out_writer.writerow([new_coords[0][0][3:], new_coords[0][1]] + row[2:])
        # Split the lifted over file by chromosome
        # This part is still a bad hack since we don't verify the right files are created
        os.system("awk -F '\t' '{{print $0 > " + input[0].split(".tsv")[0] + "$1}}'")
                    
rule fake_variants:
    input:
        region_file
    output:
        os.path.join(data_dir, "fake_{region_type}_variants.tsv.gz".format(region_type=region_type))
    run:
        ref_seqs = dict()
        for chr_num in range(1, 23):
            ref_seqs[str(chr_num)] = pysam.FastaFile(filename=os.path.join(ref_dir,
                                                                "chr{}.fa.gz".format(chr_num)))
        with open(input[0] , "r") as f_in, gzip.open(output[0], "wt", newline="") as f_out:
            fake_writer = csv.writer(f_out, delimiter="\t", lineterminator="\n")
            region_reader = csv.reader(f_in, delimiter="\t")
            next(region_reader)
            for region in region_reader:
                chrom = region[4]
                start = int(region[5])
                end = int(region[6])
                ref_seq = ref_seqs[chrom]
                ref_region = ref_seq.fetch(reference="chr"+chrom, start=start-1, end=end)
                for ii, coord in enumerate(range(start, end+1)):
                    ref_base = ref_region[ii].upper()
                    for base in base_set:
                        if base != ref_base:
                            fake_writer.writerow([chrom, coord, coord, ref_base+"/"+base, "+",
                                                  region[0]+"_"+region[3]+"_"+ref_base+"_"+base])

rule sort_fakes:
    input:
        os.path.join(data_dir, "fake_{region_type}_variants.tsv.gz".format(region_type=region_type))
    output:
        os.path.join(data_dir, "fake_{region_type}_variants_sorted.tsv.gz".format(region_type=region_type))
    shell:
        "gzip -cd {input} | sort -T /net/scratch/ -k 1,1 -k 2,2 -n | tr -d '\r' | gzip -c > {output}"

rule split_fakes_by_region:
    input:
        os.path.join(data_dir, "fake_{region_type}_variants.tsv.gz".format(region_type=region_type))
    output:
        [os.path.join(data_dir, "region_splits", "fake_{region_type}_variants_".format(region_type=region_type) +
                      region_name + ".tsv.gz")
         for region_name in region_names]
    params:
        region_type = region_type
    run:
        fstem = os.path.join(data_dir, "region_splits", "fake_{region_type}_variants_".format(region_type=params.region_type))
        command =  "./split_script.sh " + input[0] + " " + fstem
        print(command)
        os.system(command)
        
rule vep_splits:
    input:
        os.path.join(data_dir, "region_splits/fake_{region_type}_variants_{{region_name}}.tsv.gz".format(region_type=region_type))
    output:
        os.path.join(data_dir, "region_splits/fake_{region_type}_variants_{{region_name}}_vep.tsv.gz".format(region_type=region_type))
    params:
        fname = os.path.join(data_dir, "region_splits/fake_{region_type}_variants_{{region_name}}.tsv".format(region_type=region_type))
    shell:
        # --plugin dbNSFP,{dbNSFP},{dbNSFP_fields}
        # --plugin AncestralAllele,{ancestral}
        # Might like to use this but seems to make certain jobs take a very long time
        "export PERL5LIB={perl_lib}; export PATH={perl_exec}:$PATH; gzip -cd {input} > {params.fname}; "
        "{vep_exec} --cache --dir_cache {vep_cache} --dir_plugins {vep_plugins} --plugin AncestralAllele,{ancestral} --plugin dbNSFP,{dbNSFP},{dbNSFP_fields} --everything --force_overwrite  --no_stats --verbose -i {params.fname} -o STDOUT | gzip -c > {output}; "
        "rm {params.fname}"

rule process_splits:
    input:
        os.path.join(data_dir, "region_splits/fake_{region_type}_variants_{{region_name}}_vep.tsv.gz".format(region_type=region_type))
    output:
        temp(os.path.join(data_dir,
                     "region_splits/fake_{region_type}_variants_{{region_name}}_processed.tsv.gz".format(region_type=region_type)))
    run:
        vep_out = pd.read_csv(input[0], sep="\t",
                              names=["Uploaded_variation", "Location", "Allele", "Gene",
                                     "Feature", "Feature_type", "Consequence",
                                     "cDNA_position", "CDS_position", "Protein_position",
                                     "Amino_acids", "Codons", "Existing_variation", "Extra"],
                              comment="#")
        outfile_uncomp = os.path.splitext(output[0])[0]
        with open(outfile_uncomp, "w", newline="") as f_out:
            fake_writer = csv.writer(f_out, delimiter="\t", lineterminator="\n")
            fake_writer.writerow(["Gene", "Canonical_transcript", "Chrom", "Pos",
                                  "Allele_ref", "Allele", "Consequence",
                                  "cDNA_position", "CDS_position", "Protein_position",
                                  "Noncanonical_consequences", "Amino_acids", "Codons",
                                  "Existing_variation", "Extra"])
            current_gene = None
            current_canonical = None
            current_pos = None
            current_allele = None
            noncanonical_consequences = ""
            current_row = [None]*15
            for ii in range(len(vep_out)):
                row_gene, row_canonical, row_ref, row_alt = vep_out.Uploaded_variation[ii].split("_")
                row_gene_vep = vep_out.Gene[ii]
                row_feature_vep = vep_out.Feature[ii]
                if (current_pos != vep_out.Location[ii].split(":")[1]) or (current_allele != vep_out.Allele[ii]):
                    current_row[9] = noncanonical_consequences[1:]
                    if (current_pos is not None) and (current_row[0] is not None):
                        fake_writer.writerow(current_row)
                    noncanonical_consequences = ""
                    current_row = [None]*15
                    current_gene = row_gene
                    current_canonical = row_canonical.split(".")[0]
                    current_pos = vep_out.Location[ii].split(":")[1]
                    current_allele = vep_out.Allele[ii]
                if row_feature_vep == current_canonical:
                    current_row[0] = row_gene
                    current_row[1] = row_canonical
                    current_row[2] = vep_out.Location[ii].split(":")[0]
                    current_row[3] = vep_out.Location[ii].split(":")[1]
                    current_row[4] = row_ref
                    current_row[5] = vep_out.Allele[ii]
                    current_row[6] = vep_out.Consequence[ii]
                    current_row[7] = vep_out.cDNA_position[ii]
                    current_row[8] = vep_out.CDS_position[ii]
                    current_row[10] = vep_out.Protein_position[ii]
                    current_row[11] = vep_out.Amino_acids[ii]
                    current_row[12] = vep_out.Codons[ii]
                    current_row[13] = vep_out.Existing_variation[ii]
                    current_row[14] = vep_out.Extra[ii]
                else:
                    noncanonical_consequences += (";" + vep_out.Feature[ii] + "," + vep_out.Feature_type[ii] + "," +
                                                  vep_out.Consequence[ii])
            
            current_row[9] = noncanonical_consequences[1:]
            if current_row[0] is not None:
                fake_writer.writerow(current_row)
        os.system("gzip " + outfile_uncomp)


rule add_mut_rate:
    input:
        os.path.join(data_dir, "region_splits/fake_{region_type}_variants_{{region_name}}_processed.tsv.gz".format(region_type=region_type))
    output:
        os.path.join(data_dir, "region_splits/fake_{region_type}_variants_{{region_name}}_mu.tsv.gz".format(region_type=region_type))
    run:
        with gzip.open(input[0], "rt") as f_in, gzip.open(output[0], "wt", newline="") as f_out:
            mut_fname_v4 = "1_rate_v5"
            mut_file_v4 = open(os.path.join(rate_dir_v4, mut_fname_v4), "r")
            mut_reader_v4 = csv.reader(mut_file_v4, delimiter=" ")
            mut_current_v4 = next(mut_reader_v4)
            pos_current_v4 = int(mut_current_v4[0])
            in_reader = csv.reader(f_in, delimiter="\t")
            out_writer = csv.writer(f_out, delimiter="\t", lineterminator="\n")

            header = next(in_reader)
            out_writer.writerow(header + ["mu", "mu_quality", "context"])
            
            quality = None
            context_current = None

            alt_dict_v4 = {}
            
            for row in in_reader:
                chrom = row[2]
                try:
                    pos = int(row[3])
                except ValueError:
                    continue
                ref = row[4]
                alt = row[5]
                # print(chrom, pos, ref, alt)
                
                ## Update v4 mutation rate reader if not on the right chromosome
                if chrom != mut_fname_v4.split("_")[0]:
                    mut_file_v4.close()
                    mut_fname_v4 = "{}_rate_v5".format(chrom)
                    mut_file_v4 = open(os.path.join(rate_dir_v4, mut_fname_v4), "r")
                    mut_reader_v4 = csv.reader(mut_file_v4, delimiter=" ")
                    mut_current_v4 = next(mut_reader_v4)
                    pos_current_v4 = int(mut_current_v4[0])

                # Locate current mutation in v4
                if pos_current_v4 != pos:
                    pos_current_v4 = int(mut_current_v4[0])
                    # print("locating mutation position: " + str(pos) + " from " + str(pos_current_v4))
                    while pos_current_v4 < pos:
                        mut_current_v4 = next(mut_reader_v4)
                        pos_current_v4 = int(mut_current_v4[0])
                    ref_current_v4 = mut_current_v4[1][2]
                    pos_ind_v4 = pos_current_v4
                    if ref_current_v4 != nt_to_CT(ref) and (pos_ind_v4 == pos):
                        print("reference doesn't match, {}:{}".format(chrom, pos))
                    alt_dict_v4 = {}
                    while (pos_ind_v4 == pos) and (ref_current_v4 == nt_to_CT(ref)):
                        quality = mut_current_v4[4]
                        context_current = mut_current_v4[1][0:5]
                        alt_dict_v4[mut_current_v4[1][6]] = float(mut_current_v4[3])
                        mut_current_v4 = next(mut_reader_v4)
                        pos_ind_v4 = int(mut_current_v4[0])
                        
                try:
                    out_writer.writerow(row + [alt_dict_v4[mut_to_CT(ref, alt)], quality, context_current])
                except KeyError:
                    out_writer.writerow(row + ["NA", "NA", "NA"])
        
rule add_coverage:
    input:
        os.path.join(data_dir, "region_splits/fake_{region_type}_variants_{{region_name}}_mu.tsv.gz".format(region_type=region_type))
    output:
        os.path.join(data_dir, "region_splits/fake_{region_type}_variants_{{region_name}}_cov.tsv.gz".format(region_type=region_type))
    run:
        with gzip.open(input[0], "rt") as f_in, gzip.open(output[0], "wt", newline="") as f_out:
            cov_fname = "gnomad.exomes.lifted.coverage.summary.1"
            cov_file = open(os.path.join(cov_dir, cov_fname), "r")
            cov_reader = csv.reader(cov_file, delimiter="\t")
            cov_current = next(cov_reader)
            pos_current = int(cov_current[1])
            
            in_reader = csv.reader(f_in, delimiter="\t")
            out_writer = csv.writer(f_out, delimiter="\t", lineterminator="\n")

            header = next(in_reader)
            out_writer.writerow(header + ["mean", "median", "over_1", "over_5", "over_10", "over_15" ,
                                          "over_20","over_25", "over_30", "over_50", "over_100"])
            for row in in_reader:
                chrom = row[2]
                pos = int(row[3])
                if chrom != cov_fname.split(".")[-1]:
                    print("changing chromosome ...")
                    cov_file.close()
                    cov_fname = "gnomad.exomes.lifted.coverage.summary.{}".format(chrom)
                    cov_file = open(os.path.join(cov_dir, cov_fname), "r")
                    cov_reader = csv.reader(cov_file, delimiter="\t")
                    cov_current = next(cov_reader)
                    pos_current = int(cov_current[1])
                while pos_current < pos:
                    try:
                        cov_current = next(cov_reader)
                        pos_current = int(cov_current[1])
                    except StopIteration:
                        break
                # if pos_current > pos:
                #     out_writer.writerow(row + ["NA"]*len(cov_current[2:]))
                if pos_current == pos:
                    out_writer.writerow(row + cov_current[2:])
                else:
                    out_writer.writerow(row + ["NA"]*len(cov_current[2:]))
                                                            
def chrom_num(chrom_str):
    try:
        return int(chrom_str[3:])
    except ValueError:
        print("Past numeric chromosomes")
        return None

def chrom_str(chrom_num):
    if "chr" in str(chrom_num):
        return chrom_num
    return "chr" + str(chrom_num)
                
def zoom_to(chrom, pos, variants, variant, window):
    if variant is None:
        print("starting as none..")
        return None
    curr_chrom = chrom_num(variant.chrom)
    # print(curr_chrom, variant.pos)
    if curr_chrom is None:
        return None
    if curr_chrom > chrom:
        return variant
    else:
        if (variant.pos > (pos + window)) and (curr_chrom == chrom):
            return variant
        while curr_chrom < chrom:
            # print("going to a new chrom! {}-{}".format(curr_chrom, chrom), end=" - ")
            try:
                variant = next(variants)
                curr_chrom = chrom_num(variant.chrom)
                if curr_chrom is None:
                    return None
            except StopIteration:
                print("end of file")
                return None
        while variant.pos < (pos - window) and (curr_chrom==chrom):
            try:
                variant = next(variants)
                curr_chrom = chrom_num(variant.chrom)
                if curr_chrom is None:
                    return None
                # print(curr_chrom, variant.pos)
            except StopIteration:
                print("end of file")
                return None
        return variant

def add_info(gnomad_info, variant, gnomad_info_fields=["AN", "AN_nfe", "AC", "AC_nfe"]):
    gnomad_info["pos"].append(variant.pos)
    gnomad_info["alt"].append(variant.alts[0])
    gnomad_info["filter"].append(";".join(variant.filter.keys()))

    for field in gnomad_info_fields:
        try:
            if type(variant.info[field]) is tuple:
                gnomad_info[field].append(variant.info[field][0])
            else:
                gnomad_info[field].append(variant.info[field])
        except KeyError:
            # print("missing {} at {}:{}".format(field, variant.chrom, variant.pos))
            gnomad_info[field].append("NA")

def remove_elements(annos, to_kill):
    return [anno for ii, anno in enumerate(annos) if ii not in to_kill]

def check_lengths(gnomad_info):
    all_same = True
    vals = gnomad_info.keys()
    pos_val = len(gnomad_info["pos"])
    len_dict = {}
    for val in vals:
        if val != "chrom":
            len_dict[val] = len(gnomad_info[val])
            if len_dict[val] != pos_val:
                all_same = False
    return all_same, len_dict

def update_gnomad(chrom, pos, gnomad_info, variants, variant, fname, window=500,
                  gnomad_info_fields=["AN", "AN_nfe", "AC", "AC_nfe"]):
    if chrom != gnomad_info["chrom"]:
        print("wrong_chrom {}:{}".format(chrom, gnomad_info["chrom"]))
        gnomad_info["chrom"] = chrom
        gnomad_info["pos"] = []
        gnomad_info["alt"] = []
        gnomad_info["filter"] = []
        for field in gnomad_info_fields:
            gnomad_info[field] = []
        vcf_in = pysam.VariantFile(fname)
        variants = vcf_in.fetch(chrom_str(chrom), start=max(0, pos - 2*window - 10), stop=None)
        try:
            variant = next(variants)
        except StopIteration:
            variant = None
    else:
        to_kill = [ii for ii in range(len(gnomad_info["pos"])) if
                   (gnomad_info["pos"][ii] < pos - window) or (gnomad_info["pos"][ii] > pos + window)]
        
        gnomad_info["pos"] = remove_elements(gnomad_info["pos"], to_kill)
        gnomad_info["alt"] = remove_elements(gnomad_info["alt"], to_kill)
        gnomad_info["filter"] = remove_elements(gnomad_info["filter"], to_kill)
        for field in gnomad_info_fields:
            gnomad_info[field] = remove_elements(gnomad_info[field], to_kill)

    # print(variant.pos, variant.chrom)
    variant = zoom_to(chrom, pos, variants, variant, window)
    # print("after zooming", variant.pos, variant.chrom)
    if variant is None:
        print("nowhere to zoom! {}:{}".format(chrom, pos))
        return None 
    var_pos = variant.pos
    var_chrom = chrom_num(variant.chrom)
    # print("adding info {}:{} ".format(chrom, pos), end=" , ")
    while ((var_pos >= (pos - window)) and (var_pos <= (pos + window)) and
           var_chrom==chrom):
        add_info(gnomad_info, variant, gnomad_info_fields=gnomad_info_fields)
        try:
            variant = next(variants)
            var_pos = variant.pos
            var_chrom = chrom_num(variant.chrom)
        except StopIteration:
            variant = None
            var_pos = None
            var_chrom = None
            return variant
    # print("{}:{}".format(var_chrom, var_pos))

    lengths_same, lengths = check_lengths(gnomad_info)
    assert lengths_same, (gnomad_info, lengths)
    
    return variant, variants

rule nearest_AN_splits:
    input:
        gnomad_file,
        # "/net/data/gnomAD/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz",
        os.path.join(data_dir, "region_splits/fake_{region_type}_variants_{{region_name}}_cov.tsv.gz".format(region_type=region_type))
    output:
        temp(os.path.join(data_dir,
                     "region_splits/fake_{region_type}_variants_{{region_name}}_gnomad.tsv.gz".format(region_type=region_type)))
    run:
        gnomad_info_fields = ["AN", "AN_nfe", "AC", "AC_nfe", "controls_AN", "controls_AC", "non_cancer_AN", "non_cancer_AC",
                              "BaseQRankSum", "ClippingRankSum", "DP", "FS", "InbreedingCoeff", "MQ", "MQRankSum",
                              "QD", "ReadPosRankSum", "SOR", "VQSLOD", "VQSR_NEGATIVE_TRAIN_SITE", "VQSR_POSITIVE_TRAIN_SITE",
                              "VQSR_culprit", "allele_type", "dp_hist_all_n_larger", "dp_hist_alt_n_larger",
                              "has_star", "lcr", "n_alt_alleles", "pab_max",
                              "rf_label", "rf_negative_label", "rf_positive_label", "rf_tp_probability_label"]
        gnomad_info = {}
        gnomad_info["pos"] = []
        gnomad_info["alt"] = []
        gnomad_info["filter"] = []
        for field in gnomad_info_fields:
            gnomad_info[field] = []
        
        window = 0

        with gzip.open(input[1], "rt") as f_in, gzip.open(output[0], "wt", newline="") as f_out:
            var_writer = csv.writer(f_out, delimiter="\t", lineterminator="\n")
            var_reader = csv.reader(f_in, delimiter="\t")
            header = next(var_reader)            
            var_writer.writerow(header + ["filter"] + gnomad_info_fields + ["interp_dist"])
            first = True
            for var in var_reader:
                chrom = int(var[2])
                pos = int(var[3])
                # print("current var: {} {}".format(chrom, pos))
                if first:
                    vcf_in = pysam.VariantFile(input[0])
                    variants = vcf_in.fetch(chrom_str(chrom), start=max(0, pos - 2*window - 10), stop=None)
                    try:
                        variant = next(variants)
                        gnomad_info["chrom"] = chrom_num(variant.chrom)
                    except StopIteration:
                        variant = None
                        gnomad_info["chrom"] = chrom
                    first = False
                    
                variant, variants = update_gnomad(chrom, pos, gnomad_info, variants, variant, input[0],
                                                  window=window, gnomad_info_fields=gnomad_info_fields)
                if len(gnomad_info["pos"]) == 0:
                    # print("nothing nearby")
                    var_writer.writerow(var + ["NA"]*(len(gnomad_info_fields)+1) + [-1])
                    continue
                if pos in gnomad_info["pos"]:
                    # print("exact position found")
                    matching = [ii for ii, pp in enumerate(gnomad_info["pos"]) if pp==pos]
                    alts = [gnomad_info["alt"][ii] for ii in matching]
                    alt = var[5]
                    if alt in alts:
                        ii = matching[alts.index(alt)]
                        var_writer.writerow(var + [gnomad_info["filter"][ii]] +
                                            [gnomad_info[field][ii] for
                                             field in gnomad_info_fields] + [0])
                    else:
                        dists = [abs(pos - pp) for pp in gnomad_info["pos"]]
                        min_d = min(dists)
                        ii = dists.index(min_d)
                        var_writer.writerow(var + [gnomad_info["filter"][ii]] +
                                            [gnomad_info[field][ii] for
                                             field in gnomad_info_fields] + [-1])

                else:
                    # print("distance interpolating")
                    dists = [abs(pos - pp) for pp in gnomad_info["pos"]]
                    min_d = min(dists)
                    ii = dists.index(min_d)
                    var_writer.writerow(var + [gnomad_info["filter"][ii]] +
                                        [gnomad_info[field][ii] for
                                         field in gnomad_info_fields] + [min_d])

rule add_GIAB:
    input:
        os.path.join(data_dir, "region_splits/fake_{region_type}_variants_{{region_name}}_gnomad.tsv.gz".format(region_type=region_type))
    output:
        os.path.join(data_dir, "region_splits/fake_{region_type}_variants_{{region_name}}_giab.tsv.gz".format(region_type=region_type))
    run:
        annotated_vars = synonAnnotate(input[0], GIAB_file)
        if annotated_vars is None:
            os.system("gzip -cd " + input[0] + " | sed '1  s/.*/&\tsite_in_genome_bottle/' | gzip -c > " + output[0])
        else:
            annotated_vars.to_csv(output[0], sep="\t", index=False, header=True, na_rep="NA")


##################################################Functions for Mutation Rate################################################

def get_Carlson_rate(chrom_df, mut_df):
    mut_df_copy = mut_df.copy()
    
    mut_df_copy = mut_df_copy[["Pos", "Carlson_rate"]]
    
    return chrom_df.join(mut_df.set_index('Pos'), on='Pos')
    
def chunks(l, n):
    n = max(1, n)
    return (l[i:i+n] for i in range(0, len(l), n))

############################################################################################################################

rule get_Carlson:
    input:
        os.path.join(data_dir,
                     "region_splits/fake_{region_type}_variants_{{region_name}}_giab.tsv.gz".format(region_type=region_type))
    output:
        os.path.join(data_dir,
                     "region_splits/fake_{region_type}_variants_{{region_name}}_carlson.tsv.gz".format(region_type=region_type))
    run:        
        var_filename = str(input[0])
        output_filename = str(output[0])

        var_dd = pd.read_csv(var_filename, sep = '\t', dtype={"VQSR_culprit":"object",
                                                             "allele_type":"object",
                                                             "cDNA_position":"object",
                                                             "filter":"object",
                                                              "rf_label":"object"}, low_memory=False)
        var_dd["Chrom"] = var_dd["Chrom"].astype(int)
        var_dd["ref_alt"] = var_dd["Allele_ref"] + var_dd["Allele"]
        
        first = 1
        
        chrom_list = list(var_dd["Chrom"].unique())
        
        header_list = ["Pos", "Carlson_rate"]

        for chrom in chrom_list:
            mut_rate_df_CT = dd.read_csv(carlson_rate_folder +
                                         str(chrom).strip() + "_C_T_v2", sep = ' ', header = None, names = header_list)
            mut_rate_df_CA = dd.read_csv(carlson_rate_folder +
                                         str(chrom).strip() + "_C_A_v2", sep = ' ', header = None, names = header_list)
            mut_rate_df_CG = dd.read_csv(carlson_rate_folder +
                                         str(chrom).strip() + "_C_G_v2", sep = ' ', header = None, names = header_list)
            
            mut_rate_df_TG = dd.read_csv(carlson_rate_folder +
                                         str(chrom).strip() + "_T_G_v2", sep = ' ', header = None, names = header_list)
            mut_rate_df_TC = dd.read_csv(carlson_rate_folder +
                                         str(chrom).strip() + "_T_C_v2", sep = ' ', header = None, names = header_list)
            mut_rate_df_TA = dd.read_csv(carlson_rate_folder +
                                         str(chrom).strip() + "_T_A_v2", sep = ' ', header = None, names = header_list)

            var_chrom_dd = var_dd.loc[var_dd["Chrom"] == chrom]
            pos_list_all = list(var_chrom_dd["Pos"].unique())

            for pos_list in list(chunks(pos_list_all, 100000)):
                mut_rate_df_filter_CT = mut_rate_df_CT.loc[mut_rate_df_CT["Pos"].isin(pos_list)]
                mut_rate_df_filter_CA = mut_rate_df_CA.loc[mut_rate_df_CA["Pos"].isin(pos_list)]
                mut_rate_df_filter_CG = mut_rate_df_CG.loc[mut_rate_df_CG["Pos"].isin(pos_list)]
                
                mut_rate_df_filter_TG = mut_rate_df_TG.loc[mut_rate_df_TG["Pos"].isin(pos_list)]
                mut_rate_df_filter_TC = mut_rate_df_TC.loc[mut_rate_df_TC["Pos"].isin(pos_list)]
                mut_rate_df_filter_TA = mut_rate_df_TA.loc[mut_rate_df_TA["Pos"].isin(pos_list)]

                mut_rate_df_filter_pd_CT = mut_rate_df_filter_CT.compute()
                mut_rate_df_filter_pd_CA = mut_rate_df_filter_CA.compute()
                mut_rate_df_filter_pd_CG = mut_rate_df_filter_CG.compute()
                
                mut_rate_df_filter_pd_TG = mut_rate_df_filter_TG.compute()
                mut_rate_df_filter_pd_TC = mut_rate_df_filter_TC.compute()
                mut_rate_df_filter_pd_TA = mut_rate_df_filter_TA.compute()

                var_chrom_dd_filter = var_chrom_dd.loc[var_chrom_dd["Pos"].isin(pos_list)]

                ## Break up current variants into mutation types
                var_chrom_dd_filter_CT = var_chrom_dd_filter[var_chrom_dd_filter["ref_alt"].isin(["CT", "GA"])]
                var_chrom_dd_filter_CA = var_chrom_dd_filter[var_chrom_dd_filter["ref_alt"].isin(["CA", "GT"])]
                var_chrom_dd_filter_CG = var_chrom_dd_filter[var_chrom_dd_filter["ref_alt"].isin(["CG", "GC"])]
                
                var_chrom_dd_filter_TG = var_chrom_dd_filter[var_chrom_dd_filter["ref_alt"].isin(["TG", "AC"])]
                var_chrom_dd_filter_TC = var_chrom_dd_filter[var_chrom_dd_filter["ref_alt"].isin(["TC", "AG"])]
                var_chrom_dd_filter_TA = var_chrom_dd_filter[var_chrom_dd_filter["ref_alt"].isin(["TA", "AT"])]

                var_chrom_dd_filter_CT = get_Carlson_rate(var_chrom_dd_filter_CT, mut_rate_df_filter_pd_CT)
                var_chrom_dd_filter_CA = get_Carlson_rate(var_chrom_dd_filter_CA, mut_rate_df_filter_pd_CA)
                var_chrom_dd_filter_CG = get_Carlson_rate(var_chrom_dd_filter_CG, mut_rate_df_filter_pd_CG)

                var_chrom_dd_filter_TG = get_Carlson_rate(var_chrom_dd_filter_TG, mut_rate_df_filter_pd_TG)
                var_chrom_dd_filter_TC = get_Carlson_rate(var_chrom_dd_filter_TC, mut_rate_df_filter_pd_TC)
                var_chrom_dd_filter_TA = get_Carlson_rate(var_chrom_dd_filter_TA, mut_rate_df_filter_pd_TA)

                var_chrom_dd_filter = pd.concat([var_chrom_dd_filter_CT, var_chrom_dd_filter_CA, var_chrom_dd_filter_CG,
                                                 var_chrom_dd_filter_TG, var_chrom_dd_filter_TC, var_chrom_dd_filter_TA])
                var_chrom_dd_filter.sort_values(by="Pos", inplace=True)

                if first == 1:
                    var_chrom_dd_filter.to_csv(output_filename, sep = '\t', index=False,
                                               compression="gzip", na_rep="NA")
                    first = 0
                else:
                    var_chrom_dd_filter.to_csv(output_filename, mode= 'a', header=False, sep = '\t',
                                               index=False, compression="gzip", na_rep="NA")


rule grab_syn:
    input:
        os.path.join(data_dir,
                     "region_splits/fake_{region_type}_variants_{{region_name}}_carlson.tsv.gz".format(region_type=region_type))
    output:
        os.path.join(data_dir,
                     "region_splits/fake_{region_type}_variants_{{region_name}}_syn.tsv.gz".format(region_type=region_type))
    run:
        use_cols = ["Gene", "Canonical_transcript", "Chrom", "Pos", "Allele_ref", "Allele", "Consequence",
                    "CDS_position", "Protein_position", "Amino_acids", "Extra", "mu", "mu_quality", "context",
                    "mean", "median", "over_100", "filter", "AN", "AN_nfe", "AC", "AC_nfe",
                    "BaseQRankSum", "ClippingRankSum", "ReadPosRankSum", "SOR", "allele_type", "n_alt_alleles",
                    "rf_label", "interp_dist", "site_in_genome_bottle", "Carlson_rate"]
        with gzip.open(input[0], "rt") as f_in, gzip.open(output[0], "wt", newline="") as f_out:
            var_writer = csv.writer(f_out, delimiter="\t", lineterminator="\n")
            var_reader = csv.reader(f_in, delimiter="\t")
            header = next(var_reader)
            col_inds = [header.index(col) for col in use_cols]
            var_writer.writerow([header[ii] for ii in col_inds])
            for var in var_reader:
                consequence = var[col_inds[use_cols.index("Consequence")]]
                if "synonymous" in consequence:
                    var_writer.writerow([var[ii] for ii in col_inds])
                
            

        
# rule merge_gpcr:
#     input:
#         [GPCR_fname + "_giab.tsv.gz" for GPCR_fname in GPCR_fnames]
#     output:
#         os.path.join(data_dir, "fake_transcript_variants_v4_GPCR.tsv.gz")
#     run:
#         with gzip.open(output[0], "wt", newline="") as f_out:
#             var_writer = csv.writer(f_out, delimiter="\t", lineterminator="\n")
#             with gzip.open(input[0], "rt") as f_in:
#                 var_reader = csv.reader(f_in, delimiter="\t")
#                 header = next(var_reader)
#                 var_writer.writerow(header)
#                 for row in var_reader:
#                     var_writer.writerow(row)
#             for ii in range(1, len(input)):
#                 with gzip.open(input[ii], "rt") as f_in:
#                     var_reader = csv.reader(f_in, delimiter="\t")
#                     header = next(var_reader)
#                     for row in var_reader:
#                         var_writer.writerow(row)

rule merge_all:
    input:
        [os.path.join(data_dir, "region_splits/fake_{region_type}_variants_".format(region_type=region_type) +
                      region_name + "_giab.tsv.gz")
         for region_name in region_names]
    output:
        os.path.join(data_dir, "fake_{region_type}_variants_all.tsv.gz".format(region_type=region_type))
    run:
        with gzip.open(output[0], "wt", newline="") as f_out:
            var_writer = csv.writer(f_out, delimiter="\t", lineterminator="\n")
            with gzip.open(input[0], "rt") as f_in:
                var_reader = csv.reader(f_in, delimiter="\t")
                header = next(var_reader)
                var_writer.writerow(header)
                for row in var_reader:
                    var_writer.writerow(row)
            for ii in range(1, len(input)):
                with gzip.open(input[ii], "rt") as f_in:
                    var_reader = csv.reader(f_in, delimiter="\t")
                    header = next(var_reader)
                    for row in var_reader:
                        var_writer.writerow(row)
    
            
rule merge_syn:
    input:
        syn_files = [os.path.join(data_dir, "region_splits/fake_{region_type}_variants_".format(region_type=region_type) +
                                  region_name + "_syn.tsv.gz")
                     for region_name in region_names]
    output:
        outfile = os.path.join(data_dir, "fake_{region_type}_variants_syn.tsv.gz".format(region_type=region_type))
    run:
        base_name = os.path.splitext(output.outfile)[0]
        print("gzip -cd " + input.syn_files[0] + " > " + base_name)
        os.system("gzip -cd " + input.syn_files[0] + " > " + base_name)
        for ii in range(1, len(input.syn_files)):
            print("gzip -cd " + input.syn_files[ii] + " | tail -n +2  >> " + base_name)
            os.system("gzip -cd " + input.syn_files[ii] + " | tail -n +2 >> " + base_name)
        os.system("sort -T /net/scratch/ -k 3,3 -k 4,4 -k 6,6 -n " + base_name + " | gzip -c > " + output.outfile)
        os.system("rm " + base_name)

# rule merge_ms:
#     input:
#         syn_files = [split_fname + "_ms_GIAB_v5.gz" for split_fname in split_fnames]
#     output:
#         outfile = os.path.join(data_dir, "fake_transcript_variants_sorted_all_ms_GIAB_v5.gz")
#     run:
#         base_name = os.path.splitext(output.outfile)[0]
#         print("gzip -cd " + input.syn_files[0] + " > " + base_name)
#         os.system("gzip -cd " + input.syn_files[0] + " > " + base_name)
#         for ii in range(1, len(input.syn_files)):
#             print("gzip -cd " + input.syn_files[ii] + " >> " + base_name)
#             os.system("gzip -cd " + input.syn_files[ii] + " >> " + base_name)
#         os.system("sort -k 1,1 -k 2,2 -n " + base_name + " | gzip -c > " + output.outfile)
#         os.system("rm " + base_name)
