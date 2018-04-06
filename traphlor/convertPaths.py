#!/usr/bin/env python
"""
Call for one path file at a time, results are appended to the output file

Path files of format:
Path 1
Path 2
...

Author: aekuosma

"""

import optparse, pysam, sys, os

# returns True/False on if used id
# bamfile required to extract the strand
# Count the total expression in the caller
def convert_path_no_annotation(ofile, nodesfile, pathsfile, gene_id, bamfile, debug):
    if debug:
        print nodesfile, pathsfile
    if not os.path.exists(pathsfile):
        sys.stderr.write("convertPaths internal error: Paths file %s doesn't exist.\nSome transcripts might be missing from the results file.\n" % (pathsfile))
        return
        
    if not os.path.exists(nodesfile):
        sys.stderr.write("convertPaths internal error: Nodes file %s doesn't exist.\nSome transcripts might be missing from the results file.\n" % (nodesfile))
        return
        
    if not os.path.exists(bamfile):
        sys.stderr.write("convertPaths: BAM file doesn't exist.\nExiting.\n")
        sys.exit(2)

    out = open(ofile, "a")
    in_nodes = open(nodesfile, "r")
    in_paths = open(pathsfile, "r")
	
    nodes = list()
	
    while True:
        line = in_nodes.readline()
		
        if line == "":
            break
		
        parts = line.split("\t")
		
        nodes.append((parts[0], (long)(parts[1]), (long)(parts[2])))
		
    in_nodes.close()
	
    paths = list()

    while True:
        line = in_paths.readline()
        if line == "":
            break
        if "Running time" in line:
            break

        paths.append(line.strip("\n").strip(" ").split(" "))
	
    in_paths.close()
    
    if len(paths) == 0:
        out.close()
        return False

    # No accepting single exons
#    if len(paths) == 1 and len(paths[0]) == 1:
#       out.close()
#       return False

    paths = sorted(paths)
    
    transcript_count = 1
    
    bam_handle = pysam.Samfile(bamfile, "rb")
	
    for ind, path in enumerate(paths):

        if path[-1] == "":
            path.pop()

        source = "Traph"
        frame = "."
        score = 1
        
        strand = "."

        forward = False
        reverse = False

        # Iterate over the reads crossing exons on the path to check for strand tags
        for node in path:
            selected_exon = nodes[int(node)]
            it = bam_handle.fetch(selected_exon[0], selected_exon[1]-1, selected_exon[2]-1)
            try:
                read = it.next()
            
                while len(read.cigar) == 1:
                    read = it.next()

                for pair in read.tags:
                    if pair[0] == "XS":
                        if read.opt("XS") == "+":
                            forward = True
                        elif read.opt("XS") == "-":
                            reverse = True
                
            except StopIteration:
                pass

        if forward and not reverse:
            strand = "+"
        elif not forward and reverse:
            strand = "-"
            
        
        exons = list()
				
        last_end = -1
		
        for index, node in enumerate(path):
            current_node = nodes[(int)(node)]
            if index == 0:
                seqname = current_node[0]
                transcript_start = current_node[1]
            if index == len(path)-1:
                transcript_end = current_node[2]
            start = current_node[1]
            end = current_node[2]

            if last_end == start -1:
                last_exon = exons[-1].split("\t")
                exons.pop()
                exons.append("\t".join(last_exon[0:4]) + "\t" + str(end) + "\t" + "\t".join(last_exon[5:8]))
            else: 
                exons.append(seqname + "\t" + source  + "\t" + "exon" + "\t" + str(start) + "\t" + str(end) + "\t" + str(score) + "\t" + strand + "\t" + frame)
                
            last_end = end

        out.write('%s\t%s\ttranscript\t%s\t%s\t%s\t%s\t%s\tgene_id "TRAPH.%s"; transcript_id "TRAPH.%s.%s"\n' % (seqname, source, transcript_start, transcript_end, score, strand, frame, gene_id, gene_id, transcript_count))
        
        for index, exon in enumerate(exons):
            out.write('%s\tgene_id "TRAPH.%s"; transcript_id "TRAPH.%s.%s"; exon_number "%s"\n' % (exon, gene_id, gene_id, transcript_count, index+1))
            
        transcript_count = transcript_count + 1

    in_paths.close()
    out.close()
    bam_handle.close()
    return True


