import pandas as pd
import numpy as np
from tqdm import tqdm
import pysam
import os, sys
import subprocess as sp
import argparse
pd.options.mode.chained_assignment = None 
#Need to fix these hard-coded paths to user inputs
#also in addition change outputs to output directory instead of current directory. Lines 166,168,170

parser = argparse.ArgumentParser(description = 'identify DMD duplications, requires flye and samtools and bcftools in your environment path or supply the abs path to it. ')
parser.add_argument("-b", "--bamFn", help = "path for bam file")
parser.add_argument('-r', "--reference", help = "path for reference genome fasta file")
parser.add_argument("-sv", "--structure_variants", help = "path for delly or sniffles sv vcf")
parser.add_argument('-o', "--outdir", help = "path for output directory", default = "dmduper_output")
parser.add_argument('--force', action = "store_true", help = "force over-write, default = False")

parser.add_argument("-m", "--softclip_min_length",
                    dest="min_length",
                    help="Minimum length of softclip for DMD dup",
                    default=100)
parser.add_argument("-w", "--bin_length_cov",
                    dest="bin_length",
                    help="Bin size for coverage sliding window calculation",
                    default=100)
parser.add_argument("-u", "--spike_up_threshold",
                    dest="spike_up",
                    help="The fold for spike up detection",
                    default=1.2)
parser.add_argument("-l", "--spike_low_threshold",
                    dest="spike_low",
                    help="The fold for spike low detection",
                    default=0.2)
parser.add_argument("-g", "--gene_region_boundary",
                    dest="gene_region",
                    help="The fold for spike low detection",
                    default="chrX:31117222-33341388")
parser.add_argument("-t", "--cov_threshold",
                    dest="cvg_threshold",
                    help="The mnimum threshold for detection",
                    default=10)                   

args = parser.parse_args()

flye= "/research/bsi/tools/biotools/flye/2.9.4/bin/flye" #"/research/labs/dlmp/adl/ADL0074/SOW3/src/ADL0074_ONT/tools/flye"
seqkit="/research/bsi/tools8/biotools/seqkit/2.7.0/bin/seqkit"
minimap2="/research/labs/dlmp/adl/ADL0077/src/ADL0074_ONT/tools/minimap2"
reference=args.reference 

def subset_sv_file(delly_sv_bcf, dmd_sv, region):
    cmd = f'bcftools view {delly_sv_bcf} {region} > {dmd_sv}'
    #print(cmd)
    sp.call(cmd, shell =True)

def read_variants(DMD_SVFile):
    vcf_file_Stream = open(DMD_SVFile, "r")    
    deletions = {}
    for eachline in vcf_file_Stream:
        if eachline[0] == "#":
            continue
        eachline_split = eachline.split("\t")
        chrom = eachline_split[0]
        start = int(eachline_split[1])
        info_col = eachline_split[7].split(";")
        stop = 0
        if "SVTYPE=DEL" in info_col:
            #print(info_col)
            for eachcol in info_col:
                if "END" == eachcol.split("=")[0]:
                    stop = int(eachcol.split("=")[1])
            if chrom in deletions:
                deletions[chrom] = deletions[chrom] + ";" + chrom + ":" + str(start) + "-" + str(stop)
            else:
                deletions[chrom] = chrom + ":" + str(start) + "-" + str(stop)
    return deletions





###################################### parsing tandem calls #############################################
def parse_sam_file(sam_file, exclude_secondary = False):
    # Open the SAM file
    samfile = pysam.AlignmentFile(sam_file, "r")
    # Dictionary to store mapping information
    mapping_info = {}
    for read in samfile.fetch():
        if not read.is_unmapped:
            if exclude_secondary and read.is_secondary: continue
            query_name = read.query_name
            reference_name = read.reference_name
            reference_start = read.reference_start
            reference_end = read.reference_end
            # Calculate query positions
            query_start = read.query_alignment_start
            query_end = read.query_alignment_end
            
            if query_name not in mapping_info:
                mapping_info[query_name] = []
            mapping_info[query_name].append((reference_name, reference_start, reference_end, query_start, query_end))
    samfile.close()
    return mapping_info

def tandem_caller(ref_st, ref_ed, scale = 10):
    dmd_st, dmd_ed = 31072747, 33246985
    dmd_len = dmd_ed-dmd_st
    map_len = ref_ed - ref_st
    dmd_range = set(range(dmd_st, dmd_ed))
    plot_dist = 50000
    dmd_len_scale = dmd_len/scale
    map_len_scale = map_len/scale 
    outlist = []
    if len(set(range(ref_st,ref_ed)).intersection(dmd_range)) == 0:
        dist = ref_ed - dmd_ed
        if dist < 0 : 
            #### in upstream 
            real_dist = dmd_st - ref_ed
            print(f'Non-tandem duplication detected {abs(real_dist)}bp away upstream, start: {ref_st}, end: {ref_ed}  !!!')
            d_st, d_ed = plot_dist * -1, 0
            r_xmin, r_xmax = d_st - map_len_scale,  d_st + map_len_scale
            up_end_xmin, up_end_xmax = r_xmin - plot_dist, r_xmin
            dn_end_xmin, dn_end_xmax = dmd_len_scale, dmd_len_scale + plot_dist
            outlist = [d_st, d_ed, r_xmin, r_xmax, up_end_xmin, up_end_xmax, dn_end_xmin, dn_end_xmax, 'upstream', real_dist]
            #outlist = list(map(lambda x:x/scale, outlist))
            #padding_new = (split_padding + dmd_st)*-1
            ###### draw upstream genomic region  
        if dist > 0 :
            #### in downstream
            real_dist = ref_st - dmd_ed
            print(f'Non-tandem duplication detected {abs(real_dist)}bp away downstream, start: {ref_st}, end: {ref_ed}  !!!')
            d_st, d_ed = dmd_len_scale, plot_dist + dmd_len_scale  #dashed line indicate the out region distance 
            r_xmin, r_xmax = d_ed, d_ed + map_len_scale
            up_end_xmin, up_end_xmax = -plot_dist, 0
            dn_end_xmin, dn_end_xmax = r_xmax, r_xmax + plot_dist
            outlist = [d_st, d_ed, r_xmin, r_xmax, up_end_xmin, up_end_xmax, dn_end_xmin, dn_end_xmax,'downstream',real_dist]
            #outlist = list(map(lambda x:x/scale, outlist))

    else:
        outlist = []
        print('duplications is likely in tandem')
    return outlist

def draw_exon(geneBed, y, ax):
    import matplotlib.pyplot as plt
    import numpy as np
    texts = []
    label_pos_x = np.linspace(geneBed.st_scale.min()-2000, geneBed.ed_scale.max()+2000, geneBed.shape[0])
    geneBed['loc'] = label_pos_x
    for x in geneBed.index:
        exon_label, exon_loc = geneBed.loc[x, ['exon','loc']]
        #ax.hlines(xmin = st, xmax = ed, y = y, color = 'k', linewidth = 20)
        ax.text(exon_loc, y + .01, exon_label, rotation = 60, fontsize = 8 )
        
def label_distance(tandem_call, ax ):
    # Define the two points
    x1, y1 = tandem_call[0], 3
    x2, y2 = tandem_call[1], 3
    distance = str(tandem_call[-1]) + 'bp'
    # Define the control point for the Bezier curve to make it sharp
    control_x, control_y = (x1 + x2) / 2, 3.5  # Adjust the control_y for sharpness
    # Generate points for the Bezier curve
    t = np.linspace(0, 1, 100)
    curve_x = (1-t)**2 * x1 + 2 * (1-t) * t * control_x + t**2 * x2
    curve_y = (1-t)**2 * y1 + 2 * (1-t) * t * control_y + t**2 * y2
    # Plot the Bezier curve
    ax.plot(curve_x, curve_y, 'grey')
    # Label the distance at the tip of the curve
    ax.text(control_x, 3.3, distance, ha='center', va='bottom')
            
def draw_ref_nontandem(spike, dmd_bed_fn, tandem_call, ax, scale = 10, tandem=False):
    dmd_st, dmd_ed = 31072747, 33246985
    dmd_len = dmd_ed - dmd_st
    
    spike_st, spike_ed = list(map(int, spike.split(':')[1].split('-')))
    dmd_bed = pd.read_csv(dmd_bed_fn, sep = '\t',header = None)
    dmd_bed = dmd_bed.iloc[:, [0,1,2,4,3]]
    dmd_bed.columns = ['chrom','st','ed','exon','gene']
    dmd_bed.sort_values(by = 'st', ascending=False, inplace = True)
    spike_bed = dmd_bed.loc[(dmd_bed.st >= spike_st) & (dmd_bed.ed <= spike_ed)]
    
    spike_bed['st_scale'] = (spike_bed['st'] - dmd_st)/scale
    spike_bed['ed_scale'] = (spike_bed['ed'] - dmd_st)/scale
    
    spike_st_scale, spike_ed_scale = spike_st/scale, spike_ed/scale
    dmd_st_scale, dmd_ed_scale, dmd_len_scale = dmd_st/scale, dmd_ed/scale, dmd_len/scale
    

    ######## draw the dmd, coordinates from 0, so everythign shift by subtracting dmd_st
    ax.hlines(xmin = 0, xmax = dmd_len_scale, y = 3, linewidth = 10) ## draw gene region
    ax.hlines(xmin=(spike_st_scale - dmd_st_scale), xmax = (spike_ed_scale - dmd_st_scale), y = 3, color = 'red', linewidth = 10) ## raw spike region
    ax.text(1, 3.2, 'DMD')
    ######## draw exons ######
    draw_exon(spike_bed, 3.2, ax)
    
    if not tandem:
        ###################### draw out region ############
        
        dash_xmin, dash_xmax, r_xmin, r_xmax, up_end_xmin, up_end_xmax, dn_end_xmin, dn_end_xmax = tandem_call[:-2]
        ###### draw dash line
        ax.hlines(xmin = dash_xmin, xmax = dash_xmax, y = 3, linestyles='dashed', linewidth=3, color = 'grey')
        ###### draw out region
        ax.hlines(xmin = r_xmin, xmax = r_xmax, y = 3,  linewidth=10, color = 'grey')

        #ax.text(out_range_xmax + split_padding, 3.2, str(abs(up_stream_aln)) + 'bp', rotation = 30)

        ####### ends of the genomic region ####   
        ax.hlines(xmin = up_end_xmin, xmax = up_end_xmax, y = 3, linestyles='dashed', linewidth=3, color = 'grey')
        ax.hlines(xmin = dn_end_xmin, xmax = dn_end_xmax, y = 3, linestyles='dashed', linewidth=3, color = 'grey')


        ##### draw distance 
        label_distance(tandem_call,ax)
        


def draw_contig(assembly, contigID, ax):
    from Bio import SeqIO
    cDict = {}
    fa = SeqIO.parse(assembly, 'fasta')
    for rec in fa:
        if rec.id == contigID:
            contigLen = len(rec.seq)
            ax.hlines(xmin = 0, xmax = contigLen , y = 1, linewidth= 10, color = 'blue' )
            ax.text(contigLen/4, .7, contigID)

    
def draw_mapping_track(tracks, tandem_call, ax,  gene = False, scale = 10):
    
    dmd_st, dmd_ed = 31072747, 33246985
    dmd_len = dmd_ed - dmd_st
    dmd_st_scale, dmd_ed_scale, dmd_len_scale = dmd_st/scale, dmd_ed/scale, dmd_len/scale
    plot_dist = 50000
    y=[1.1,2.9]
    if gene:
        for track in tracks:
            qst, qed, rst, red = track
            rst_scale, red_scale = (rst - dmd_st)/scale, (red - dmd_st)/scale
            contig_x1 = [qst, rst_scale]
            contig_x2  =[qed, red_scale] 
            ax.plot(contig_x1, y, color = 'orange', linewidth=.1)
            ax.plot(contig_x2, y, color = 'orange', linewidth=.1)
            ax.fill_betweenx(y, contig_x1, contig_x2 ,color='orange', alpha=0.5)
    
    else:
        rst_scale, red_scale = tandem_call[2], tandem_call[3]  
        #print(shift_dist)
        for track in tracks:
            qst, qed, rst, red = track            
            #rst_shift, red_shift = (rst - dmd_st)/scale,  (red - dmd_st)/scale
            #r_len = red_shift - rst_shift
            
            #rst_scale, red_scale = shift_dist - r_len, shift_dist
            #print('scaled ref mapping is ', rst_scale, red_scale)
            contig_x1 = [qst, rst_scale]
            contig_x2 = [qed, red_scale]
            ax.plot(contig_x1, y, color = 'grey', linewidth=.1)
            ax.plot(contig_x2, y, color = 'grey', linewidth=.1)
            ### draw sections in contig that aligned to regions outside of dmd
            ax.hlines(xmin = qst, xmax = qed, y = 1, linewidth= 10, color = 'grey', alpha = .3 ) 
            ax.fill_betweenx(y, contig_x1, contig_x2 ,color='k', alpha=0.3)
            
def profile_tandem_calls(mapping_info):
    # Print the mapping information
    nontandem = False
    tandem_calls = {}

    for query_name, mappings in mapping_info.items():
        #print(f"Query: {query_name}")
        if not query_name in tandem_calls:
            tandem_calls[query_name] = {}
            tandem_calls[query_name]['tandem_call'] = []
            tandem_calls[query_name]['in_dmd'] = []
            tandem_calls[query_name]['out_dmd'] = []
            
        for ref_name, ref_start, ref_end, query_start, query_end in mappings:
            if not (ref_end - ref_start) >= 2000:continue
            #print(f"  Mapped to {ref_name} from {ref_start} to {ref_end} (Query: {query_start} to {query_end})")
            tandem_call = tandem_caller(ref_start, ref_end)

            if len(tandem_call) > 0 :
                tandem_calls[query_name]['out_dmd'].append([query_start, query_end, ref_start, ref_end])
                tandem_calls[query_name]['tandem_call'].append(tandem_call)
            else:
                tandem_calls[query_name]['in_dmd'].append([query_start, query_end, ref_start, ref_end])
            
    return tandem_calls        

def record_alignment_positions_with_softclipping(bam_file_path, deletions):
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")
    countdict = {}
    
    #with open(output_txt_path, 'w') as output:
        #output.write("Read_Name\tChromosome\tPosition_Start\tPosition_End\tIs_Supplementary\tSoft_Clip_Start\tSoft_Clip_End\n")
        
    softclip_dict = {}
    for read in bamfile:
        if read.is_unmapped:
            continue  # Skip unmapped reads
        
        if read.mapping_quality == 0:
            continue
        
        
        
        # Extract soft clipping details from the CIGAR string
        soft_clip_start = read.cigartuples[0][1] if read.cigartuples[0][0] == 4 else 0
        soft_clip_end = read.cigartuples[-1][1] if read.cigartuples[-1][0] == 4 else 0
        
        # Check for more than 100 bases of soft clipping on either end
        if soft_clip_start > 100 or soft_clip_end > 100:
            chrom = bamfile.get_reference_name(read.reference_id)
            start = read.reference_start + 1  # Convert to 1-based position
            end = read.reference_end
            #is_supplementary = 'yes' if read.is_supplementary else 'no'
            check_if_deletion_exists = 0
            if chrom in deletions:
                chrom_dels = deletions[chrom]
                for eachdel in chrom_dels.split(";"):
                    delstart = eachdel.split(":")[1].split("-")[0]
                    delstop = eachdel.split(":")[1].split("-")[1]
                    
                    if abs(int(start) - int(delstart)) <= 2 or abs(int(start) - int(delstop)) <=2:
                        #deletion is within 2 bases of spike start - signature of del and not a dup
                        check_if_deletion_exists = 1
                        break
                    if abs(int(end) - int(delstart)) <= 2 or abs(int(end) - int(delstop)) <=2:
                        #deletion is within 2 bases of spike start - signature of del and not a dup
                        check_if_deletion_exists = 1
                        break
                    
                if check_if_deletion_exists == 1: #If deletion exist coverage spike after deletion should not be considered dup
                    #if "DMD02" in bam_file_path:
                    #    print(start, end, delstart, delstop)
                    continue
            
            if start not in softclip_dict:
                softclip_dict[start] = {}
            if end not in softclip_dict:
                softclip_dict[end] = {}
            softclip_dict[start][read.query_name] = read.query_sequence
            softclip_dict[end][read.query_name] = read.query_sequence    

            if start not in countdict:
                countdict[start] = 1
            else:
                countdict[start] += 1
            if end not in countdict:
                countdict[end] = 1
            else:
                countdict[end] += 1
    
    
    #output.write(f"{read.query_name}\t{chrom}\t{start}\t{end}\t{is_supplementary}\t{soft_clip_start}\t{soft_clip_end}\n")

    bamfile.close()
    return countdict, softclip_dict


###################################### Run pipeline ##################
if not os.path.exists(args.outdir):
    os.mkdir(args.outdir)
dmd_sv = os.path.join(args.outdir, 'dmd_only_sv.vcf')
if not os.path.exists(dmd_sv) or args.force: 
    #print("subset")
    subset_sv_file(args.structure_variants, dmd_sv, args.gene_region) 

deletions = read_variants(dmd_sv)
print('detected deletions :', deletions)

dupcounter = 0

#result_df = pd.read_csv('spike.csv')

print("*** Step 2: Extract soft-clipped reads in DMD")
ignore_pos = [31465211,31471610]
countdict, softclip_dict = record_alignment_positions_with_softclipping(args.bamFn, deletions)
sorted_dict = dict(sorted(countdict.items(), key=lambda item: item[1]))
softclip_reads_file = open(args.outdir + "/softclip_reads.fa","w")
fasta_file_name = args.outdir + "/softclip_reads.fa"
regions = []
chr_req = args.gene_region.split(":")[0]
for eachpos in sorted_dict:
    if sorted_dict[eachpos] >= 5:
        if int(eachpos) not in ignore_pos:
            print(eachpos)
            for eachread in softclip_dict[eachpos]:
                softclip_reads_file.write(">" + eachread + "\n" + softclip_dict[eachpos][eachread] + "\n")
            regions.append(chr_req +":" +  str(int(eachpos)-1) + "-" + str(eachpos))
softclip_reads_file.close()    

potential_dups = open(args.outdir + "/dup_candidates.txt","w")
for eachregion in regions:
    potential_dups.write(eachregion)


dup_rm_fasta = args.outdir + "/softclip_reads_dup_rm.fa"
dup_rm_cmd = f"{seqkit} rmdup -s {fasta_file_name} -o {dup_rm_fasta}"
print(f'running assembly using fly, command: {dup_rm_cmd}')
sp.run(dup_rm_cmd, shell=True, capture_output=True, text=True)

print("*** Step 3: Denovo haplotype-aware assembly of soft-clip reads")
assembly_out = f"{args.outdir}/softclipped_region"
if not os.path.exists(assembly_out + '/assembly.fasta'):
    flye_cmd = f"{flye} --keep-haplotypes --nano-hq {dup_rm_fasta} -o {assembly_out}"
    print(f'running assembly using fly, command: {flye_cmd}')
    sp.run(flye_cmd, shell=True, capture_output=True, text=True)
else:
    print('assemly from softclipped reads already exists, skip this step...')
print("*** Step 4: Aligning assembled contigs to reference genome ****")
asm_aln_out = os.path.join(assembly_out, 'aln.sam')
minimap2_cmd = f"{minimap2} -ax asm5 -x map-ont -Y {reference} {assembly_out}/assembly.fasta > {asm_aln_out}"
if not os.path.exists(asm_aln_out):
    sp.call(minimap2_cmd, shell=True) 
else:
    print(f'assembly already aligned, skip aligning')


print("*** Step 5: Runing tamdem caller **********")

############## parsing minimap2 output and plot
for eachregion in regions:
    import matplotlib.pyplot as plt 
    print(f'processing spike identified in {eachregion}')
    dmd_bed_fn = '/research/labs/dlmp/adl/ADL0077/src/ADL0074_ONT/data/DMD_NM_00406.bed'
    #dmd_bed_fn = "/research/labs/dlmp/adl/ADL0077/dev/dmd/PMS2.bed"
    assembly = os.path.join(assembly_out, 'assembly.fasta') 

    #mapping_info = parse_sam_file(sam_file, exclude_secondary=True)
    mapping_info = parse_sam_file(asm_aln_out)
    tandem_calls = profile_tandem_calls(mapping_info)

    for query in tandem_calls:
        query_tandem_calls = tandem_calls[query]['tandem_call']
        #print(query, query_tandem_calls)
        fig, ax = plt.subplots(figsize=(12,4))
        ax.set_ylim(0,4)
        ax.set_xticks([])
        ax.set_yticks([1,3],['Assembly', 'Reference'])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        if len(query_tandem_calls) > 0:
            out_mapping = tandem_calls[query]['out_dmd']
            in_mapping = tandem_calls[query]['in_dmd']
            for call in query_tandem_calls:         
                draw_ref_nontandem(eachregion, dmd_bed_fn, call, ax)
                draw_contig(assembly, query, ax)
                draw_mapping_track(in_mapping, call, ax, gene=True)
                draw_mapping_track(out_mapping, call, ax)

        else:
            call = []
            in_mapping = tandem_calls[query]['in_dmd']
            draw_ref_nontandem(eachregion, dmd_bed_fn, call, ax, tandem=True)
            draw_contig(assembly, query, ax)
            draw_mapping_track(in_mapping, call, ax, gene=True)
        plt.savefig(f'{assembly_out}/{query}_tandem_call.png', facecolor = 'white', dpi = 600)
