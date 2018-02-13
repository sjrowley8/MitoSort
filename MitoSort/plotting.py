"""
Plotting module for mitosort.

@author Stephen Rowley

"""

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def plot_risk(risk_file, png_file):

    plt.style.use('ggplot')
    
    df = pd.read_csv(risk_file)
    df['Quality_Percentage'] = df['Quality_Score'] * 100
    risk_df = df[df['Quality_Category'] != '-']
    high_df = risk_df[risk_df['Quality_Category'] == 'high']
    med_df = risk_df[risk_df['Quality_Category'] == 'unknown']
    low_df = risk_df[risk_df['Quality_Category'] == 'low']
    
    unique_samples = np.unique(risk_df['Line'])
    fig, ax = plt.subplots(figsize=(10,5))
    #plt.figure(figsize=(10,5))
    
    plt.title('iPSC Risk')
   
    #dummy to build axes
    plt.scatter(unique_samples, np.linspace(0, 100, len(unique_samples)),alpha=0)    


    
    x_range = np.arange(len(unique_samples))
    plt.xticks(x_range, unique_samples, rotation=90)
    plt.ylabel('iPSC Risk Score')
    plt.xlabel('Subjects')
    for x, y in zip([high_df, med_df, low_df], ['r', 'g', 'b']):
        ax.plot(x['Line'], x['Quality_Percentage'], c=y, marker='o', linestyle='None')

    #build legend 
    v1 = plt.scatter([], [], c='r')
    v2 = plt.scatter([], [], c='g')
    v3 = plt.scatter([], [], c='b')
    leg = plt.legend([v1, v2, v3], ['High', 'Unknown', 'Low'],bbox_to_anchor=(1.16,0.65), title='iPSC Clones')
    plt.tight_layout(rect=[0,0,0.89,1])
    plt.savefig(png_file)
    plt.close()

def plotMutations(mut_file,png_file,group,gene_file):

    plt.style.use('ggplot')
    mut_df = pd.read_csv(mut_file,sep="\t")

        # dummy data to make room for gene labels
    dummy = {'Sample':' ', 'Sample_Type':'iPSC', 'Position':0, 'Variant':0, 
         'Variant_Type':'DEL', 'Heteroplasmy':0, 'Mut_Type':'globalPrivateMut'}

    mut_df = mut_df.append(dummy, ignore_index=True)
    sample_types = pd.unique(mut_df["Sample_Type"])

    genes = pd.read_csv(gene_file, sep="\t")
    genes['x_avg'] = (genes['Ending'] + genes['Starting']) / 2
    genes['xdiff'] = genes['Ending'] - genes['Starting']


    #slice data for graphing
    fb = mut_df[mut_df['Sample_Type'] == 'FB']
    ips = mut_df[mut_df['Sample_Type'] == 'iPSC']
    del_ips = ips[ips['Variant_Type']=='DEL']
    ins_ips = ips[ips['Variant_Type']=='INS']    
    snp_ips = ips[ips['Variant_Type']=='SNV']
    del_fb = fb[fb['Variant_Type']=='DEL']
    ins_fb = fb[fb['Variant_Type']=='INS']    
    snp_fb = fb[fb['Variant_Type']=='SNV']

    glob_del_ips = del_ips[del_ips['Mut_Type']=='globalPrivateMut']
    loc_del_ips = del_ips[del_ips['Mut_Type'] =='localPrivateMut']

    glob_ins_ips = ins_ips[ins_ips['Mut_Type'] == 'globalPrivateMut']
    loc_ins_ips = ins_ips[ins_ips['Mut_Type'] =='localPrivateMut']

    glob_snp_ips = snp_ips[snp_ips['Mut_Type'] == 'globalPrivateMut']
    loc_snp_ips = snp_ips[snp_ips['Mut_Type'] =='localPrivateMut']

    glob_del_fb = del_fb[del_fb['Mut_Type'] == 'globalPrivateMut']
    loc_del_fb = del_fb[del_fb['Mut_Type'] =='localPrivateMut']

    glob_ins_fb = ins_fb[ins_fb['Mut_Type'] == 'globalPrivateMut']
    loc_ins_fb = ins_fb[ins_fb['Mut_Type'] =='localPrivateMut']

    glob_snp_fb = snp_fb[snp_fb['Mut_Type'] == 'globalPrivateMut']
    loc_snp_fb = snp_fb[snp_fb['Mut_Type'] =='localPrivateMut']

        # dummy plots for building legends
    v1 = plt.scatter([], [], c='r') 
    v2 = plt.scatter([], [], c='g')
    v3 = plt.scatter([], [], c='b')
    mt1 = plt.scatter([], [], c='grey', marker='o') 
    mt2 = plt.scatter([], [], c='grey', marker='^')


    f1 = plt.scatter([], [], s=25, c='grey', edgecolors='none')
    f2 = plt.scatter([], [], s=50, c='grey', edgecolors='none')
    f3 = plt.scatter([], [], s=75, c='grey', edgecolors='none')
    f4 = plt.scatter([], [], s=100, c='grey', edgecolors='none')

    names =  ["DEL","INS","SNV"]
    colors = ["#F8766D","#00BA38","#619CFF"]
    num_samples = len(np.unique(mut_df['Sample']))

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True,
                                   gridspec_kw = \
                                   {'height_ratios':[1,num_samples]})

    fig.set_size_inches(13, 5)
    ax0.set_title(group + ' Private Mutations')
    ax1.set_xlabel('Position')
    ax0.margins(y=0.1)
    ax1.margins(y=0.1)
    fig.text(0.01, 0.5, 'Sample', ha='center', va='center', 
             rotation='vertical', size='smaller')


    # dummy to build axes
    ax0.scatter(0, np.unique(fb['Sample']), alpha=0)
    ax1.scatter(np.linspace(0, 16000, num_samples-1), np.unique(ips['Sample']), alpha=0)    

    # scatter plot w/ shapes, colors, and sizes            
    for x, y, z in zip(colors, [glob_del_fb, glob_ins_fb, glob_snp_fb],
                       [glob_del_ips, glob_ins_ips, glob_snp_ips]):
        ax0.scatter(y['Position'], y['Sample'], s=y['Heteroplasmy'] * 1.55, marker='o', c=x, alpha=0.66)
        ax1.scatter(z['Position'], z['Sample'], s=z['Heteroplasmy'] * 1.55, marker='o', c=x, alpha=0.66)
        

    for x, y, z in zip(colors, [loc_del_fb, loc_ins_fb, loc_snp_fb],
                       [loc_del_ips, loc_ins_ips, loc_snp_ips]):
        ax0.scatter(y['Position'], y['Sample'], s=y['Heteroplasmy'] * 1.55, marker='^', c=x, alpha=0.66)
        ax1.scatter(z['Position'], z['Sample'], s=z['Heteroplasmy'] * 1.55, marker='^', c=x, alpha=0.66)

    # scaling issues for text labels
    text_scale = -.1*(num_samples-1)/2

     # gene labels
    for i in range(len(genes)):
        ax1.add_patch(mpatches.Rectangle((genes['Starting'][i], 0), 
                                         genes['xdiff'][i], 0.02*num_samples))  
        ax1.text(x=genes['x_avg'][i], y=text_scale, 
                 s=genes['Shorthand'][i], fontsize=8, ha='center')

    # dummy patch for "second" legend 
    title_proxy = mpatches.Rectangle((0,0), 0, 0, color='w', alpha=0)

    leg = plt.legend([f1, f2, f3, f4, title_proxy, v1, v2, v3, title_proxy,
                      mt1, mt2], 
                     ['25', '50', '75', '100', 
                      ' ', "DEL", 'INS', 'SNV',
                      ' ', 'Global', 'Local'], 
                     bbox_to_anchor=(1,1)) 

    # axis limits
    y0_lim = ax0.get_ylim()
    y1_lim = ax1.get_ylim()
    x_lim = ax1.get_xlim()
    patch_start = x_lim[1] - x_lim[1]*0.015
    # add patches and labels to axes indicating sample type
    fb_patch = mpatches.Rectangle((patch_start, y0_lim[0]), 
                                  x_lim[1] - patch_start, 
                                  y0_lim[1] + (-1*y0_lim[0]), color='.75')
    ipsc_patch = mpatches.Rectangle((patch_start, y1_lim[0]), 
                                    x_lim[1] - patch_start, 
                                    y1_lim[1] + (-1*y1_lim[0]), 
                                    color='.75')
    ax0.add_patch(fb_patch)
    ax1.add_patch(ipsc_patch)

    ax0.text(patch_start, 0, 'FB', 
             color='black', rotation=270, size=7)
    ax1.text(patch_start, np.median(np.arange(num_samples-1)), 
             'iPSC', 
             color='black', rotation=270, size=7)
    fig.subplots_adjust(right=0.14)
    fig.tight_layout()
    plt.savefig(png_file)
    plt.close()

def plotCoverage(wig_file,png_file,sample,gene_file,log_scale=False):

	wig = pd.read_csv(wig_file)
	genes = pd.read_csv(gene_file, sep='\t')

	# parse out wig file
	span = int(list(wig)[0].split(' ')[-1].split('=')[1])
	col = list(wig)[0]
	xcoords = np.array([])
	for i in range(len(wig)):
		xcoords = np.append(xcoords, 
				    np.repeat(wig[col][i], span))

	# graph positions for gene labels
	genes['x_avg'] = (genes['Ending'] + genes['Starting']) / 2
	genes['xdiff'] = genes['Ending'] - genes['Starting']

	plt.style.use('ggplot')

	fig, ax = plt.subplots()
	fig.set_size_inches(13, 5)
	ax.set_xlabel('Position')
	ax.set_title('Mitochondria Coverage for ' + sample)
	ax.margins(y=0.1)


	if(log_scale):
		# filling to min value to make graph look better
		xcoords = np.log(xcoords)
		ymin = np.min(xcoords)
		ax.set_ylabel('Read Depth (log)')
	
	else: 
		ymin = 0
		ax.set_ylabel('Read Depth')
	
	ymax = np.max(xcoords)

	# place gene labels correctly
	graph_bottom = ymin*1.1 - (ymax*.1)
	middle = (graph_bottom + ymin) / 2
	top_middle = (middle + ymin) / 2
	rect_center = (top_middle + middle) / 2
	bot_middle = (middle + graph_bottom) /2
	rect_height = top_middle - middle
	text_anchor = (graph_bottom + bot_middle) / 2
	ax.fill_between(np.arange(len(xcoords)), ymin, xcoords)

	for i in range(len(genes)):
		ax.add_patch(mpatches.Rectangle((genes['Starting'][i], rect_center), 
						genes['xdiff'][i], top_middle - middle))  
	      	ax.text(x=genes['x_avg'][i], y=text_anchor, 
			s=genes['Shorthand'][i], fontsize=8, 
			ha='center')

	plt.savefig(png_file)

	plt.close()

if __name__ == '__main__':


    mut_file = '/data2/bsi/tertiary/Nelson_Timothy_m000917/s203481.mtDNA_tool_development/full_output/variants/plot/5H1.private_mutations.txt'
    png_file = '/home/m169420/test.group_mutations.png'
    group = '5H1'
    gene_file = '/data2/bsi/tertiary/Nelson_Timothy_m000917/s203481.mtDNA_tool_development/git/MitoSort/MitoSort/resources/mtdna_genes.plot.txt'
    plotMutations(mut_file, png_file, group, gene_file)

    png_file = '/home/m169420/test.MitoB.36.p7.coverage.png'
    gene_file = '/data2/bsi/tertiary/Nelson_Timothy_m000917/s203481.mtDNA_tool_development/git/MitoSort/MitoSort/resources/mtdna_genes.plot.txt'
    wig_file = '/home/m169420/output_for_testing/coverage/wig/MitoB.36.p7.wig'
    sample = 'MitoB.36.p7'

    plotCoverage(wig_file,png_file,sample,gene_file)

    risk_file = '/data2/bsi/tertiary/Nelson_Timothy_m000917/s200199.HLHS_mtDNA_sequencing/poster/out7/Healthy/iPSC_Risk.csv'
    png_file = '/home/m169420/test.Healthy_risk.png'

    plot_risk(risk_file, png_file)

	
