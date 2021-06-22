# Create figures showing donor site segments that make up the
# different sequences for a variable region.

# Currently specifically set up for choice V1 donor sites in SFig 5.
import os
import sys
import dna_features_viewer
import pandas as pd 
import numpy as np
from dna_features_viewer import GraphicFeature, GraphicRecord
from matplotlib import pyplot as plt
from matplotlib import patches
from collections import Counter

# Assigns certain colors to often found 4-nucleotide repeats in V1.
def find_color(seq):
	if seq == "GCAT":
		return "firebrick"

# Reads table from ds_names.R - "donorsites.csv."
df = pd.read_csv("V1_donorsites.csv",sep=",")

fig = plt.figure(figsize=(20,15))

# Loops through each unique sequence in a variable region
for num,fasta_seq in enumerate(df.qseqid.unique()):
	location = num + 1
	# Subfigures in rows, columns, format
	ax = fig.add_subplot(7,2,location)
	subset = df.loc[df['qseqid']==fasta_seq]
	features = []

	percentage = 0

	num_events = 0 
	qstart_list = []
	ds_list = []

	highlight_repeats=False
	# Loops through each donor site segment for each unique sequence
	for index,record in subset.iterrows():

		qstart = int(record['qstart'])
		qend = int(record['qend'])+1

		# Finds where donor site segment is in full donor site sequence
		dsstart = record['qseqseq'].find(record['sequence'])
		dsend = dsstart + (qend - qstart)

		features.append(
			GraphicFeature(start=record['qstart'],end=record['qend']+1,strand=+1,color=record['color'],label=record['ds_name']))

		# Defining variables for plotting
		sequence = record['qseqseq']
		aaseq = record['aaseq']
		#percentage = record['percentage']
		num_events += 1
		qstart_list.append(int(record['qstart']))
		ds_list.append(record['qseqseq'])

		if("_any" in record['ds_name']):
			highlight_repeats=True
	
	# Sort donor sites by order of position
	qstart_list,ds_list = zip(*sorted(zip(qstart_list,ds_list)))
	# For highlighting different repeat regions
	if(highlight_repeats):
		color = "cornflowerblue"
		ax.add_patch(patches.Rectangle((1.5,-0.85),4,0.35,facecolor=color,linewidth=0,clip_on=False))
		ax.add_patch(patches.Rectangle((6.5,-0.85),4,0.35,facecolor=color,linewidth=0,clip_on=False))
		ax.add_patch(patches.Rectangle((qend-4.5,-0.85),4,0.35,facecolor=color,linewidth=0,clip_on=False))
	
	# Plotting
	record = GraphicRecord(sequence_length=45,features=features,sequence=sequence,first_index = 1,labels_spacing=-500)
	record.plot(ax = ax,figure_width = 13,figure_height=30, annotate_inline=True)
	record.plot_sequence(ax)
	plot_title = aaseq
	ax.text(0,1.5,plot_title,fontsize=14,fontweight='bold')

fig = plt.tight_layout()
plt.subplots_adjust(hspace = -0.3)
plt.savefig("SFig 5.pdf")
