# Taken from Lucas/Hans specI update



import sys
import itertools

clusterfile = sys.argv[1]
m8file = sys.argv[2]
rerep_m8file = sys.argv[3]

cluster_2_members = {}
with open(clusterfile, 'r') as handle:
    current_cluster = handle.readline().strip()[1:]
    members = []
    for line in handle:
        if line.startswith('>'):
            cluster_2_members[current_cluster] = members
            current_cluster = line.strip()[1:]
            members = []
        else:
            members.append(line.strip())
    cluster_2_members[current_cluster] = members
# '*' captures unaligned seq in vsearch, not sure it's needed here but let's not break anything
cluster_2_members['*'] = ['*']

def write_perfect_alignments(cluster, cluster_dict, outfile):
    """
    Permute all members of a 100% derep cluster and write
    the identitcal alignments to file
    """
    seq_len = int(cluster.split(';')[1].replace('length=', ''))
    for s1, s2 in [pair for pair in itertools.combinations(cluster_dict[cluster], 2)]:
        # write a mock identical alignment line in blast6out format
        # query, subject, percid, alignlength, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
        # evalue always -1 and bitscore 0
        newline = [s1, s2, 100.0, seq_len, 0, 0, 1, seq_len, 1, seq_len, -1, 0]
        newline = [str(i) for i in newline]
        outfile.write('\t'.join(newline) + '\n')

of = open(rerep_m8file, 'w')
# add identical alignments
for cluster in cluster_2_members.keys():
    if cluster != '*':
        write_perfect_alignments(cluster, cluster_2_members, of)
# propagate clustered alignments
with open(m8file, 'r') as handle:
    for line in handle:
        line = line.strip()
        splits = line.split('\t')
        query_cluster = splits[0]
        subject_cluster = splits[1]
        for query in cluster_2_members[query_cluster]:
            for subject in cluster_2_members[subject_cluster]:
                newline = [query, subject]
                newline.extend(splits[2:])
                of.write('\t'.join(newline))
                of.write('\n')
of.close()
