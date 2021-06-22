# Takes a two-column csv file with internal repeats ("acl_internal_repeats.csv"), and 
# file with blast hits for each sequence ("hits3.csv"). 

# Outputs all donor site combinations for each sequence in table format and in visual format.

import re

def extend_sequence(index_seq, segment_list, qseqseq, repeats):
    """
    index_seq: a list of indices, corresponding to the entries in segment_list

    returns: a list of all valid one-step extensions of the index_seq
    """
    _, seg_start, seg_end, seg_string, ds_name = segment_list[index_seq[-1]]
    extensions = []
    if not(seg_string == qseqseq):
        for idx, start, end, string,ds_name in segment_list[index_seq[-1]+1:]:
            if (start <= seg_end + 4 and 
                #seg_start > start and
                end > seg_end and 
                    any(rep in qseqseq[max(start-4,0):seg_end+4] for rep in repeats)):
                extensions.append(idx)
            elif start > seg_end + 4:
                break

    return [index_seq + [idx] for idx in extensions]

def extend_sequence_list(sequence_list, segment_list, qseqseq, repeats):
    """
    sequence_list: a list of index_seq, as taken by extend_sequence

    returns: pair (extended_sequences, completed_sequences) where extended_sequences are the
        extensions of index_seq in sequence_list which are strictly longer, and completed_sequences
        are the index_seq with no valid extensions
    """
    extended_sequences = [] # list of sequences which were successfully extended this run
    completed_sequences = [] # sequences with no valid extension
    for index_seq in sequence_list:
        extensions = extend_sequence(index_seq, segment_list, qseqseq, repeats)
        if len(extensions) == 0:
            completed_sequences.append(index_seq)
        else:
            extended_sequences.extend(extensions)

    return extended_sequences, completed_sequences

def all_sequences_recurse(sequence_list, completed_sequences, segment_list, qseqseq, repeats, depth):
    """
    recursive helper function which repeatedly extends index_seq in sequence list, as far as
    possible, and stopping when depth = 0

    returns a list of completed sequences
    """
    if len(sequence_list) == 0 or depth == 0:
        return sequence_list + completed_sequences
    else:
        extended, completed = extend_sequence_list(sequence_list, segment_list, qseqseq, repeats)
        return all_sequences_recurse(extended, completed_sequences + completed, segment_list, qseqseq, repeats, depth-1)


def all_sequences(init, segment_list, qseqseq, repeats, max_length = 15):
    """
    init: which index in segment_list you want your sequences to start with
    max_length: maximum number of strings to combine

    returns: list of valid index_seq, where index_seq is a list of indices corresponding to
        elements in segment_list
    """
    return all_sequences_recurse([[init]], [], segment_list, qseqseq, repeats, max_length-1)


# compute the coverage for an index sequence
def coverage(index_seq, segment_list, qseqseq):
    """
    Coverage of an index_seq, which is just the ratio of the length of the concatenation to the
    length of the original string.
    """
    return (segment_list[index_seq[-1]][2] - segment_list[index_seq[0]][1]) / len(qseqseq)

# get rid of segments that are within other segments as adds computational time, is redundant, and 
# no way of knowing which one to use
def substring_sieve(df, vregion):
    segment_list = df["sequence"].tolist()
    segment_list.sort(key=lambda s: len(s), reverse=True)
    out = []
    new_df = pd.DataFrame()
    for s in segment_list:
        if not any([s in o for o in out]):
            # Only pull segments with correct variable region donor sites
            out.append(s)
            new_df = new_df.append(df.loc[(df["sequence"]==s) & (df["ds_name"].str.contains(vregion))])
    return new_df

if __name__ == "__main__":
    import pandas as pd
    
    full_repeats = pd.read_csv('acl_internal_repeats.txt', delimiter='\t')
    #full_repeats = pd.read_csv('added_internal_repeats.txt', delimiter='\t')

    repeats_dict = {}
    # Convert repeats dataframe into dictionary with vregion as key, repeats as values
    for x in range(len(full_repeats)):
        currentid = full_repeats.iloc[x,0]
        currentvalue = full_repeats.iloc[x,1]
        repeats_dict.setdefault(currentid, [])
        repeats_dict[currentid].append(currentvalue)

    
    full_df = pd.read_csv('hits3.csv', delimiter=',')

    output = open("acl_donorsite_combinations.csv","w+")
    output.write("region,qseqid,qseqseq,num_segments,percent_covered,seg1,ds1,seg2,ds2,seg3,ds3\n")
    visual_output = open("acl_donorsite_combinations_visual_output.txt", "w+")

    
    for qseqid in full_df.qseqid.unique():
        print(qseqid)
        df = full_df.loc[full_df['qseqid']==qseqid]
        qseqseq = df.qseqseq.unique()[0]
        region = df.region.unique()[0]

        repeats = repeats_dict[region]

        # NOTE: currently substring_sieve gets rid of redundant segments (segments that are within other segments)
        # to save computational time and because we don't know how to tiebreak. Remove this line for actual full 
        # possibilities list, which may take forever.
        df = substring_sieve(df, region)
        segment_list = ([(x[0]-1,x[1],x[2],x[3]) for x in df[['qstart','qend','sequence','ds_name']].values])

        # Add repeats into segment list as they can also overlap
        for repeat in repeats:
            repeat_indices = [m.start() for m in re.finditer(repeat, qseqseq)]
            for index in repeat_indices:
                repeat_tuple = (index,index+4,repeat,"repeat")
                segment_list.append(repeat_tuple)

        segment_list.sort() # note that this is sorting first by qstart then by qend - important!
        segment_list = [(idx,) + x for idx, x in enumerate(segment_list)]
        # segment_list is a list of (index, qstart, qend, string), sorted first by qstart then by qend

        full_length_found = False
        for segment in segment_list:
            if (segment[3] == qseqseq):
                visual_output.write("#" + qseqid+": "+"\n"+qseqseq+"\n")
                visual_output.write(qseqseq + "\n\n")
                output.write(region+","+qseqid+","+qseqseq+",1,1,"+qseqseq+","+segment[4]+"\n")
                full_length_found = True 

        if not full_length_found:
            # Visual output
            visual_output.write("#" + qseqid+": "+"\n"+qseqseq+"\n")
            for _, start, end, segment,ds_name in segment_list:
                visual_output.write(" "*start + segment+"\t\t\t\t (" + ds_name + ")\n")

            # compute sequences which include the first string in string_segments
            seqs = all_sequences(0, segment_list, qseqseq, repeats, 15)
            # for seq in seqs:
            #     for ds in seq:
            #         #print(segment_list[ds][4])
            possibilities_list = []

            # Remove segments that are within other segments
            for l in [(seq, coverage(seq, segment_list, qseqseq)) for seq in seqs]:
                pop_indexes = []
                for index, segment in enumerate(l[0][1:]):
                    current_start = (segment_list[segment][1])
                    current_end = segment_list[segment][2]
                    last_start = segment_list[l[0][index]][1]
                    last_end = segment_list[l[0][index]][2]
                    if(current_start <= last_start and current_end >= last_end):
                        pop_indexes.append(index)
                for index in sorted(pop_indexes, reverse=True):
                    del l[0][index]
                visual_output.write(' '.join(str(s) for s in l) + '\n')
                possibilities_list.append(l)
            visual_output.write("\n")

            #most_coverage = (max(possibilities_list, key=lambda i: i[1]))
            #num_segments = (len(most_coverage[0]))

            for possibility in possibilities_list:
                pos_coverage = possibility[1]
                num_segments = len(possibility[0])
                for index,segment in enumerate(possibility[0]):
                    segment_seq = (segment_list[segment][3])
                    segment_ds = (segment_list[segment][4])
                    if index == 0:
                        output.write(region+","+qseqid+","+qseqseq+","+str(num_segments)+","+str(pos_coverage)+","+segment_seq+","+segment_ds)
                    else: 
                        output.write("," + segment_seq + "," + segment_ds)

                output.write("\n")