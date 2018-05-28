'''
This script is used to visualize high-quality circRNA potential circRNA transcripts
with output files produced by the KNIFE algorithm.

Run with -h arg to see input arguments
'''

import argparse

import matplotlib as mpl
mpl.use('Agg')
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
import matplotlib.patches as mpatches


def get_gene_from_gtf(attr_field, key='gene_name'):
	"""
	Extracts the gene name from the attribute field of a GTF
	To be used in vectorized/apply on a dataframe

	Returns a string- the gene name
	"""
	items = [x.strip() for x in attr_field.split(';')]
	kvp = {x.split(' ')[0]:x.split(' ')[1][1:-1] for x in items if len(x)>0}
	return kvp['gene_name']


def process_reads(row):
	"""
	Splits out the junction information into separate columns for easier parsing by other methods
	"""
	primary_junction = row['junction']
	primary_chrom, \
	primary_gene_and_pos1, \
	primary_gene_and_pos2, \
	primary_junction_style, \
	primary_strand = primary_junction.split('|')

	primary_gene1, primary_pos1 = primary_gene_and_pos1.split(':')
	primary_gene2, primary_pos2 = primary_gene_and_pos2.split(':')
	primary_pos1 = int(primary_pos1)
	primary_pos2 = int(primary_pos2)

	row['r1_junction_gene1'] = primary_gene1
	row['r1_junction_gene2'] = primary_gene2
	row['r1_junction_pos1'] = primary_pos1
	row['r1_junction_pos2'] = primary_pos2
	row['r1_junction_style'] = primary_junction_style
	row['r1_junction_strand'] = primary_strand

	second_junction = row['junctionR2']
	try:
		if len(second_junction.split('|')) > 1:
			secondary_chrom, \
			secondary_gene_and_pos1, \
			secondary_gene_and_pos2, \
			secondary_junction_style, \
			secondary_strand = second_junction.split('|')

			secondary_gene1, secondary_pos1 = secondary_gene_and_pos1.split(':')
			secondary_gene2, secondary_pos2 = secondary_gene_and_pos2.split(':')
			secondary_pos1 = int(secondary_pos1)
			secondary_pos2 = int(secondary_pos2)

			row['r2_junction_gene1'] = secondary_gene1
			row['r2_junction_gene2'] = secondary_gene2
			row['r2_junction_pos1'] = secondary_pos1
			row['r2_junction_pos2'] = secondary_pos2
			row['r2_junction_style'] = secondary_junction_style
			row['r2_junction_strand'] = secondary_strand
            
		else: # if something like 'chr10', we just ignore by inserting NaN
			row['r2_junction_gene1'] = np.nan
			row['r2_junction_gene1'] = np.nan
			row['r2_junction_pos1'] = np.nan
			row['r2_junction_pos2'] = np.nan
			row['r2_junction_style'] = np.nan
			row['r2_junction_strand'] = np.nan
	except:
		row['r2_junction_gene1'] = np.nan
		row['r2_junction_gene1'] = np.nan
		row['r2_junction_pos1'] = np.nan
		row['r2_junction_pos2'] = np.nan
		row['r2_junction_style'] = np.nan
		row['r2_junction_strand'] = np.nan

	all_same_gene = len(np.unique([
            row['r1_junction_gene1'],
            row['r1_junction_gene2'],
            row['r2_junction_gene1'],
            row['r2_junction_gene1']
        ])) == 1
	row['all_same_gene'] = all_same_gene
	return row


class Exon(object):
	"""
	A convenience object representing a single exon
	"""
	def __init__(self, chrom, start, end):
		self.chrom = chrom
		self.start = start
		self.end = end
		self.name = '%s:%d-%d' % (self.chrom, self.start, self.end)

	def __eq__(self, other):
		return self.__key() == other.__key()

	def __key(self):
		return (self.chrom, self.start, self.end)

	def __hash__(self):
		return hash(self.__key())

	def __repr__(self):
		return self.name 
        
        
def get_exon(gtf_df, pos):
	"""
	Given a GTF as a dataframe, and a position, search for exons that start or end with that base.
	Note that there can be multiple exons that share start or end positions.  We find and return
	them all.

	Note the +1 offset for the start position

	Returns a list of Exon objects.
	"""
	matches_start = gtf_df['start'] == (pos + 1)
	matches_end = gtf_df['end'] == pos 

	selected_rows = matches_start | matches_end
	subset_gtf = gtf_df.ix[selected_rows].drop_duplicates(subset=['chrom', 'start', 'end'])

	if subset_gtf.shape[0] >= 1:
		# sometimes can find multiple exons that start or end at the same position.
		matches = subset_gtf.apply(lambda x: Exon(x['chrom'], x['start'], x['end']), axis=1)
		return matches.values.tolist()
	else:
		gtf_df.to_csv('gtf_dump.tsv', sep='\t')
		raise Exception('No matching exons found')
      

def get_constituent_exons(row, gene_gtf):

	chrom = row['junction'].split('|')[0]
	r1_pos1 = row['r1_junction_pos1']
	r1_pos2 = row['r1_junction_pos2']
	r2_pos1 = row['r2_junction_pos1']
	r2_pos2 = row['r2_junction_pos2']
	strand = row['r1_junction_strand']

	exon_list = []

	# get exons corresponding to primary junction:
	es1 = get_exon(gene_gtf, r1_pos1)
	es2 = get_exon(gene_gtf, r1_pos2)
	exon_list.extend(es1)
	exon_list.extend(es2)

	# see if any r2 alignment:
	if (not np.isnan(r2_pos1)) and (not np.isnan(r2_pos2)):
		es1 = get_exon(gene_gtf, r2_pos1)
		es2 = get_exon(gene_gtf, r2_pos2)
		exon_list.extend(es1)
		exon_list.extend(es2)

	# keep only unique exons and then return those in a list
	exon_set = set(exon_list)
	return list(exon_set)


def orient_exons(reads_df, gene_gtf, primary_junction):
	"""
	Get exons explicitly involved in a circRNA construct and determine their orientation.
	"""	

	# get unique junctions in this subset of reads:
	unique_junction_rows = reads_df.drop_duplicates(subset=['junction','junctionR2'])

	# get a list of the exons which may or may not be involved in this construct.
	# Unsupported exons will be handled below
	exon_series = unique_junction_rows.apply(get_constituent_exons, args=(gene_gtf,), axis=1)

	# exon_series is a Series object where each entry is a list of Exon objects.
	# This turns it into a list of Exon objects
	exon_list = list(reduce(lambda x,y: set(x).union(set(y)), exon_series))

	# now order the exons by the start position.  Reverse if the transcript is on the minus strand
	exon_list.sort(key=lambda x: x.start)
	if primary_junction[-1] == '-':
		exon_list = exon_list[::-1]

	# get a set of all possible coordinates.  Some may be unused due to the overlapping exon issues
	# where multiple exons can share a start or end position
	coordinate_set = set()
	for x in exon_list:
		coordinate_set.add(x.start)
		coordinate_set.add(x.end)

	# now have to see how the exons link up to each other.  By using that information, we can possibly
	# remove exons that are unsupported.  That is, if two exons share the same end position, it's possible
	# that OTHER reads can span the junction and give us information about which of those two we should choose.
	#First, make some simplified junction strings which will let us identify junctions
	junction_strings = set()
	junction_strings.add(primary_junction)
	for i,row in unique_junction_rows.iterrows():
		secondary_junction = row['junctionR2']
		try:
			if len(secondary_junction.split('|')) > 1:
				junction_strings.add(secondary_junction)
		except AttributeError:
			# AttributeError is thrown if the read has NaN for a junction string
			pass

	# this block collects all the "used" starts/ends of the exons.  By "used" we mean a coordinate supported
	# by a read which mapped across a particular junction
	used_coordinates = set()
	if primary_junction[-1] == '-':
		used_coordinates = used_coordinates.union(unique_junction_rows['r1_junction_pos1'].dropna().unique()+1)
		used_coordinates = used_coordinates.union(unique_junction_rows['r1_junction_pos2'].dropna().unique())
		used_coordinates = used_coordinates.union(unique_junction_rows['r2_junction_pos1'].dropna().unique()+1)
		used_coordinates = used_coordinates.union(unique_junction_rows['r2_junction_pos2'].dropna().unique())
	else:
		used_coordinates = used_coordinates.union(unique_junction_rows['r1_junction_pos1'].dropna().unique())
		used_coordinates = used_coordinates.union(unique_junction_rows['r1_junction_pos2'].dropna().unique()+1)
		used_coordinates = used_coordinates.union(unique_junction_rows['r2_junction_pos1'].dropna().unique())
		used_coordinates = used_coordinates.union(unique_junction_rows['r2_junction_pos2'].dropna().unique()+1)

	# the set difference will give us those coordinates featured in our exon list that are not used by any reads
	# This can give us possible information to exclude some of the potential exons.
	unused_coordinates = coordinate_set.difference(used_coordinates)

	# now "group" exons that have the same start or end positions:
	# to do that, make a dataframe from the exons, so we can use pandas' groupby functionality:
	exon_df = pd.DataFrame()
	for i, exon in enumerate(exon_list):
		exon_df = pd.concat([exon_df,pd.Series([i, exon.chrom, exon.start, exon.end])], axis=1)
	exon_df = exon_df.T
	exon_df.columns = ['exon_idx','chrom','start','end']


	# check for those grouped by start position:
	keep_index1 = set()
	for i, sub_df in exon_df.groupby('start'):
		group_size = sub_df.shape[0]
		if group_size > 1:
			viable = sub_df.ix[sub_df['end'].isin(used_coordinates)]
			if viable.shape[0] != 1:
				# This covers two cases:
				#  - There are multiple supported junctions and we cannot reasonably decide which one to use.  Hence choose the first
				#  OR
				#  - neither of the exons in the sub_df dataframe were supported by any reads.  Of course, the OTHER
				# end of the exon is involved in a junction, but for this end we don't have any information to inform
				# which of the exons to choose.  Hence, just choose the first by default
				keep_index1.add(sub_df['exon_idx'].values[0])
			elif viable.shape[0] == 1:
				keep_index1.add(viable['exon_idx'].values[0])
		else: # equals 1
			keep_index1.add(sub_df['exon_idx'].values[0])

	keep_index2 = set()
	for i, sub_df in exon_df.groupby('end'):
		group_size = sub_df.shape[0]
		if group_size > 1:
			viable = sub_df.ix[sub_df['start'].isin(used_coordinates)]
			if viable.shape[0] != 1:
				# This covers two cases:
				#  - There are multiple supported junctions and we cannot reasonably decide which one to use.  Hence choose the first
				#  OR
				#  - neither of the exons in the sub_df dataframe were supported by any reads.  Of course, the OTHER
				# end of the exon is involved in a junction, but for this end we don't have any information to inform
				# which of the exons to choose.  Hence, just choose the first by default
				keep_index2.add(sub_df['exon_idx'].values[0])
			elif viable.shape[0] == 1:
				keep_index2.add(viable['exon_idx'].values[0])
		else: # equals 1
			keep_index2.add(sub_df['exon_idx'].values[0])
	
	# Get the intersection to see which of the exons we will keep, then filter them from the original list of exons 
	keep_index = keep_index1.intersection(keep_index2)
	exon_list = [exon_list[i] for i in range(len(exon_list)) if i in keep_index]

	# dictionary which maps the junction string (e.g. 'chr10|X:100|X:200|rev|+')
	# to a 2-length list which has the constitutive exons for said junction IN THE ORDER OF TRANSCRIPTION
	junction_to_exons_map = {}	
	for js in junction_strings:		
		chrom, p1, p2, t, strand = js.split('|')
		gene1,pos1 = p1.split(':')
		gene2,pos2 = p2.split(':')
		pos1 = int(pos1)
		pos2 = int(pos2)

		if strand == '-':
			e1s = [x for x in exon_list if x.start == (pos1+1)]
			e2s = [x for x in exon_list if x.end == pos2]			
		else:
			e1s = [x for x in exon_list if x.end == pos1]
			e2s = [x for x in exon_list if x.start == (pos2+1)]

		if (len(e1s)==1) and (len(e2s)==1):
			e1 = e1s[0]
			e2 = e2s[0]
			junction_to_exons_map[js] = [e1,e2]
		else:
			print exon_list
			raise Exception('Ambiguous junction')

	# It is possible to have "conflicting" transcripts-- i.e. those that support different exon junctions within the same 
	# circRNA construct.  Rather than enumerate the possibilities, just issue a warning
	unique_junctions_df = pd.DataFrame()
	inconsistent_junctions = set()
	for k,v in junction_to_exons_map.items():
		s = pd.Series({'j_id':k, 'e1':v[0].name, 'e2':v[1].name})
		unique_junctions_df = unique_junctions_df.append(s, ignore_index=True)

	for i,g in unique_junctions_df.groupby('e1'):
		if g.shape[0] > 1:
			inconsistent_junctions = inconsistent_junctions.union(g['j_id'].values)
	for i,g in unique_junctions_df.groupby('e2'):
		if g.shape[0] > 1:
			inconsistent_junctions = inconsistent_junctions.union(g['j_id'].values)
	inconsistent_junctions = list(inconsistent_junctions)
	if len(inconsistent_junctions) > 0:
		return None, None

	# Now that we know the exons and have them in order, this block of code
	# let's us locate any "gaps" in the construct.  That is, exons that are adjacent
	# in the ordering, but whose junction is NOT supported by any reads.  They COULD be
	# spliced, or there could be many exons in between.  For instance, if we had only alignments
	# across the primary junction, then we have no idea what is going on between the non-spliced end of those exons
	exon_list_with_gaps = []
	prior = None
	for i in range(len(exon_list)):
		e1 = exon_list[i]
		try:
			e2 = exon_list[i+1]
		except IndexError:
			# handles the circular nature by seeing if the last and first exon in the list are spliced
			e2 = exon_list[0]

		pairing = [e1, e2]
		if pairing in junction_to_exons_map.values():
			if prior is None:
				exon_list_with_gaps.extend(pairing)
			else:
				exon_list_with_gaps.append(e2)
			prior = 1 # any non-None value will suffice
		else:
			exon_list_with_gaps.append(None)
			prior = None
			
	return exon_list_with_gaps, junction_to_exons_map




def plot_circrna(reads, ordered_exons, junction_to_exons_map, primary_junction_id):

	fig, ax = plt.subplots(figsize=(15,15))

	# some constants for aesthetics:
	#colors = sns.xkcd_palette(['light blue', 'sage', 'silver', 'slate blue'])
	colors = sns.color_palette('deep', n_colors=10)
	
	total_arc_degrees = 45 # for the first read, how many degrees it covers
	total_arc_rad = total_arc_degrees * (np.pi/180.0)
	r0 = 5.0 # radius of exon circle
	s0 = r0*total_arc_rad # the arc length for the first read, which sets the arc length for all reads

	junction1_theta = -30 # the location of the first exon (0 is horizontal)
	exon_thickness = 1.0 # how fat to make the exons
	read_thickness = 0.1
	read_spacing = 0.15
	connector_thickness = read_thickness/5.0

	# since the junction index in KNIFE has j0=150bp up and downstream of the junction, 
	# the junction itself is located at position 150.  Thus, alignments to the junction will be relative to zero 
	# and alignments starting prior to the junction will begin at < j0.  To get the start of the read
	# relative to the junction, need this value
	j0 = 150 

	# if we happen to be making a contiguous transcript construct supported by these reads, the first and last exons will be the same
	try:
		contiguous_construct = ordered_exons[0] == ordered_exons[-1]
	except AttributeError:
		# thrown if either the first or last item in the list is None, which indicates a "gap" in the construct that we have no information about
		contiguous_construct = False
	if contiguous_construct:
		ordered_exons = ordered_exons[:-1]

	# setup the exon schematic
	N = len(ordered_exons)
	exon_boundaries = np.linspace(junction1_theta,360+junction1_theta,N+1)
	color_index = 0
	for i in range(N):
		if ordered_exons[i] is not None:	
			ax.add_patch(mpatches.Wedge((0,0), \
				r0,\
				exon_boundaries[i],\
				exon_boundaries[i+1], \
				ec="none", \
				fc=colors[color_index], \
				width=exon_thickness))
			color_index += 1
		else:
			ax.add_patch(mpatches.Wedge((0,0), \
				r0-0.5*exon_thickness, \
				exon_boundaries[i],\
				exon_boundaries[i+1], \
				ec="k", linestyle='dashed', fc="none", width=0.0))

	# add an arrow showing transcript direction:
	arrow_radius_factor = 0.7
	arrow_radius = 2*arrow_radius_factor*r0
	arrow_start_deg = junction1_theta-30
	arrow_end_deg = junction1_theta+30
	ax.add_patch(mpatches.Arc((0,0), arrow_radius, arrow_radius, theta1=arrow_start_deg, theta2 = arrow_end_deg, linewidth=3))
	arrow_end_x = arrow_radius_factor*r0*np.cos(arrow_end_deg*(np.pi/180.0))
	arrow_end_y = arrow_radius_factor*r0*np.sin(arrow_end_deg*(np.pi/180.0))

	ax.add_patch(                    #Create triangle as arrow head
		    mpatches.RegularPolygon(
		        (arrow_end_x, arrow_end_y),            # (x,y)
		        3,                       # number of vertices
		        r0/20,                # radius
		        (junction1_theta+30)*(np.pi/180.0),     # orientation
		        color='k'
		    )
		)

	exon_junction_thetas = {}
	exon_index = 0
	for i, exon in enumerate(ordered_exons):

		try:
			left_exon = exon
			right_exon = ordered_exons[i+1]

			for j, pairing in junction_to_exons_map.items():
				if pairing == [left_exon, right_exon]:
					exon_junction_thetas[j] = exon_boundaries[i+1]
		except:
			if contiguous_construct:
				exon_junction_thetas[primary_junction_id] = exon_boundaries[0]
			else:
				pass

		# number at the middle of each exon
		if exon is not None:
			label_number = str(exon_index+1)
			t = exon_boundaries[i] + 0.5*(exon_boundaries[i+1]-exon_boundaries[i])
			x1=(r0-0.5*exon_thickness)*np.cos(t*(np.pi/180.0))
			y1=(r0-0.5*exon_thickness)*np.sin(t*(np.pi/180.0))
			ax.add_patch(mpatches.Circle((x1,y1),radius=0.30,fc='w',ec='none'),)
			ax.text(x1,y1,label_number,horizontalalignment='center', verticalalignment='center', family='serif')
			exon_index += 1


	# a triangle marking the primary junction:
	location_ = (np.pi/180.0) * exon_junction_thetas[primary_junction_id]
	marker_x = (r0-1.25*exon_thickness)*np.cos(location_)
	marker_y = (r0-1.25*exon_thickness)*np.sin(location_)
	ax.add_patch(                    #Create triangle as arrow head
		    mpatches.RegularPolygon(
		        (marker_x, marker_y),            # (x,y)
		        3,                       # number of vertices
		        r0/20,                # radius
		        location_ - 0.5*np.pi,     # orientation
		        color='k'
		    )
		)

	arc_length_0 = r0 * total_arc_degrees*(np.pi/180.0)
	i = 0
	max_radius = r0
	for idx,row in reads.iterrows():
		if row['class'] != 'unmapped':
			start_pos = row['pos']
			read_length = row['readLen']
			junction = row['junction']

			# locate the junction (along the edge of the circle):
			primary_junction_theta = exon_junction_thetas[junction]

			# locate the start relative to the junction (NOT on the space of the circle)
			upstream_overlap = j0-start_pos-1
			upstream_proportion = upstream_overlap/float(read_length)

			radius = r0+(i+1)*read_spacing
			# given that radius, and knowing the arc length we want, calculate the coverage in degrees

			theta_1_rad = (upstream_proportion*s0)/radius
			theta_1_deg = theta_1_rad * (180.0/np.pi)

			start_pt_r1 = primary_junction_theta - theta_1_deg
			end_pt_r1 = start_pt_r1 + (s0/radius)*(180.0/np.pi)

			ax.add_patch(mpatches.Wedge((0,0),radius,start_pt_r1,end_pt_r1, ec="none", fc="k", width=read_thickness))

			# did the R2 map?
			junctionR2 = row['junctionR2'] # can be nan, or some string.
			try:
				junctionR2.isnull()
			except: # was a string
				if len(junctionR2.split('|')) > 1:
					start_pos = row['posR2']
					read_length = row['readLenR2']
					current_junction_theta = exon_junction_thetas[junctionR2]
					upstream_overlap = j0-start_pos-1
					upstream_proportion = upstream_overlap/float(read_length)
					theta_1_rad = (upstream_proportion*s0)/radius
					theta_1_deg = theta_1_rad * (180.0/np.pi)
					start_pt_r2 = current_junction_theta - theta_1_deg
					end_pt_r2 = start_pt_r2 + (s0/radius)*(180.0/np.pi)
					ax.add_patch(mpatches.Wedge((0,0),radius,start_pt_r2,end_pt_r2, ec="none", fc="k", width=read_thickness))
					ax.add_patch(mpatches.Wedge((0,0),radius-0.5*read_thickness+0.5*connector_thickness,current_junction_theta,primary_junction_theta, ec="none", fc="k", width=connector_thickness))
			i+=1
			max_radius = radius

	# add some radial lines 
	largest_x = r0
	largest_y = r0
	smallest_x = -r0
	smallest_y = -r0
	for j, t in exon_junction_thetas.items():
		t = t*(np.pi/180.0)
		start_pos_x = r0*np.cos(t)
		start_pos_y = r0*np.sin(t)
		end_pos_x = max_radius*np.cos(t)
		end_pos_y = max_radius*np.sin(t)
		ax.plot([start_pos_x, end_pos_x], [start_pos_y, end_pos_y], 'k:', linewidth=0.5)

	# add label text at bottom:
	label_x0 = 0.0
	label_offset = 0.005
	label_radius = 0.02 #
	exon_index = 0
	for exon in ordered_exons:
		if exon is not None:
			label_y0 = -label_offset-(2*exon_index+1)*(label_radius+label_offset)
			ax.add_patch(mpatches.Circle((label_x0,label_y0),radius=label_radius,fc='w',ec='k',transform=ax.transAxes, clip_on=False),)
			plt.text(label_x0, label_y0, str(exon_index+1), fontsize=14, transform=ax.transAxes, verticalalignment='center', horizontalalignment='center')
			label_x = label_x0 + 1.5*label_radius
			plt.text(label_x, label_y0, '%s:%s - %s' % (exon.chrom,"{:,}".format(int(exon.start)),"{:,}".format(int(exon.end))), verticalalignment='center', horizontalalignment='left', transform=ax.transAxes)
			exon_index += 1

	# get rid of the frame
	for spine in plt.gca().spines.values():
	    spine.set_visible(False)
	ax.set_xticks([])
	ax.set_yticks([])
	ax.autoscale_view()
	plt.axis('equal')
	fig.savefig('%s.png' % j, bbox_inches='tight')
	plt.close()


def parse_commandline_args():
	"""
	Reads the commandline args and prints help.
	Returns a dictionary of the parsed arguments
	"""

	parser = argparse.ArgumentParser()
	parser.add_argument('-j',
                        '--junctions',
                        required = True,
                        dest = 'junction_file',
						help = 'Path to the junction output file.  This file ends with "__circJuncProbs.txt"')

	parser.add_argument('-r',
                        '--reads',
                        required = True,
                        dest = 'reads_file',
						help = 'Path to the reads output file, which has information about the individual reads and their junction info.  This file ends with "__output.txt"')

	parser.add_argument('-g',
                        '--gtf',
                        required = True,
                        dest = 'gtf_file',
						help = 'Path to a GTF file (e.g. from UCSC or Ensembl)')

	parser.add_argument('-n',
                        '--num_reads',
                        required = True,
						type = int,
                        dest = 'read_count_threshold',
						help = 'The required minimum number of supporting reads')

	parser.add_argument('-p',
                        '--cdf',
                        required = True,
						type = float,
                        dest = 'junction_cdf_threshold',
						help = 'The required minimum threshold for the confidence as reported via the CDF.  Think of this like a p-value, except you want HIGH values.  The publication considers values greater than 0.9.')

	return vars(parser.parse_args())


if __name__ == '__main__':

	args = parse_commandline_args()
	junction_table_file = args['junction_file']
	reads_file = args['reads_file']
	gtf_file = args['gtf_file']
	read_count_threshold = args['read_count_threshold']
	junction_cdf_threshold = args['junction_cdf_threshold']

	junction_table = pd.read_table(junction_table_file)
	high_quality_junctions = junction_table.ix[
				(junction_table['numReads'] > read_count_threshold ) 
				& 
				(junction_table['junction_cdf.x'] > junction_cdf_threshold)]


	print "Number of high-quality junctions: %s" % high_quality_junctions.shape[0]
	# load the table of reads (truncated here):
	reads = pd.read_table(reads_file)
	#reads = pd.read_table('dup_output.tsv')
	print 'Total table size: %s' % reads.shape[0]
	reads = reads.ix[reads['junction'].isin(high_quality_junctions['junction'])]
	print 'Total table size (after subset): %s' % reads.shape[0]
	reads = reads.apply(process_reads,axis=1)
	#sys.exit(1)
	# load GTF, extract a column for the gene name, and subset to keep only exons:
	gtf = pd.read_table(gtf_file, names=['chrom','source','feature','start','end','score','strand', 'frame', 'attributes'])
	gtf['gene'] = gtf['attributes'].apply(get_gene_from_gtf)
	gtf = gtf.ix[gtf.feature == 'exon']

	for j in high_quality_junctions['junction']:

		print 'Image for %s' % j
		reads_subset = reads.ix[reads['junction'] == j]

		gene = j.split('|')[1].split(':')[0]
		gene_gtf = gtf.ix[gtf.gene == gene]

		ordered_exons, junction_to_exons_map = orient_exons(reads_subset, gene_gtf, j)
		if ordered_exons and junction_to_exons_map:
			plot_circrna(reads_subset, ordered_exons, junction_to_exons_map, j)
		else:
			with open('notes.txt', 'wa') as fout:
				fout.write('For junction %s, multiple transcripts were found with conflicting splice junctions.\n' % j)




















