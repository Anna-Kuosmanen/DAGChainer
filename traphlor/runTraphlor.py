#!/usr/bin/env python
"""
Script to do all the stuff required to run Traphlor.

Author:
aekuosma

"""

GRAPH_FILE_PREFIX = "gene.graph"
NODES_FILE_PREFIX = "gene.nodes"
PATHS_FILE_SUFFIX = ".sol"
ANNOTATION_PATHS_PREFIX = "annotation.paths"
TRANSCRIPT_FILE = "transcripts.gtf"
TEMP_GENE_FILE = "gene.tmp"

SAMTOOLS = "samtools"

DATE_FORMAT = "%Y.%m.%d %H:%M:%S"

import os, optparse, pysam, sys, subprocess, datetime

# Add the script location to python path that can import the other modules
sys.path.append(os.path.dirname(__file__))
    
import createGraphWithConstraints
import convertPaths

def create_graphs_no_annotation(alignment_file, odir, threshold, debug):

    file_tally = 0
    
    alignment_handle = pysam.Samfile(alignment_file, "rb")
    
    references = alignment_handle.references
    
    for reference in references:
        #print "Reference %s" % (reference)
        it = alignment_handle.fetch(reference)
        
        last_alignment_end = 0
        first_alignment_start = -1

        # Find the bounds within which all the alignments for this reference lie        
        for alignment in it:
         #   print alignment
            if first_alignment_start == -1 or first_alignment_start > alignment.pos:
                first_alignment_start = alignment.pos
            if last_alignment_end < alignment.aend:
                last_alignment_end = alignment.aend
        #print "Bounds: %s-%s" % (first_alignment_start, last_alignment_end)
        # No reads for this reference
        if first_alignment_start == -1:
            continue


        # Check base at a time till there's a gap, then process that region
        start = first_alignment_start
        end = -1

        prev_column = first_alignment_start-1        

        for pileupcolumn in alignment_handle.pileup(reference, first_alignment_start, last_alignment_end):
           # print "Column %s" % (pileupcolumn.pos)
            if pileupcolumn.pos != (prev_column + 1) and start != -1:
                end = prev_column
                #print "Creating graph for %s:%s-%s" % (reference, start, end)
                file_tally += createGraphWithConstraints.create_graph_without_annotation(alignment_file, reference, start, end, "%s/tmp/%s_%s" % (odir, GRAPH_FILE_PREFIX, file_tally + 1), "%s/tmp/%s_%s" % (odir, NODES_FILE_PREFIX, file_tally + 1), threshold, debug)
                # Start right after previous ended in case pileupcolumn misses some alignments starting (in some cases it does)
                start = prev_column + 2

            prev_column = pileupcolumn.pos

        end = last_alignment_end
        #print "Creating graph for %s:%s-%s" % (reference, start, end)
        file_tally += createGraphWithConstraints.create_graph_without_annotation(alignment_file, reference, start, end, "%s/tmp/%s_%s" % (odir, GRAPH_FILE_PREFIX, file_tally + 1), "%s/tmp/%s_%s" % (odir, NODES_FILE_PREFIX, file_tally + 1), threshold, debug)
 
    alignment_handle.close()
    return file_tally


def find_average_coverage(odir, no_of_files):
    coverages = list()

    for i in range(1, no_of_files+1):
        handle = open("%s/tmp/%s_%s" % (odir, GRAPH_FILE_PREFIX, i))
        node_number = int(handle.readline())
        for j in range(0, node_number):
            handle.readline()
        coverages = coverages + handle.readline().split(" ")
        handle.close()
    cov_sum = 0.0
    for cov in coverages:
       cov_sum += float(cov)
    return 1.0*cov_sum/len(coverages)

# Check that there's enough evidence that this graph isn't just created from mismappings
def check_graph(graph_file, nodes_file, avg_cov):

    handle = open(graph_file, "r")
    node_number = int(handle.readline())
    max_coverage = 0
    for i in range(0,node_number):
        handle.readline()
    coverages = handle.readline().split(" ")
    for coverage in coverages:
        if float(coverage) > max_coverage:
            max_coverage = float(coverage)
    handle.close()
    handle = open(nodes_file)
    max_length = 0
    for line in handle.readlines():
        parts = line.split("\t")
        max_length += (long(parts[2]) - long(parts[1]))
    handle.close()
    # If coverage of a graph is less than 1, it consists of a single read, scrap it
    if max_coverage > 1.0*avg_cov/10 and max_coverage > 1 :
        return True
    if node_number >= 3:
        return True
    if max_length > 1000:
        return True

    return False



# Called when there isn't enough evidence that the graph isn't just created from mismappings
def write_empty_solution_file(odir, i):
    handle = open("%s/tmp/%s_%s.sol" % (odir, GRAPH_FILE_PREFIX, i), "w")
    handle.close()


def run_traph(odir, no_of_files, opts):
    traph_exec = os.path.abspath(os.path.dirname(__file__)) + "/RNA_MPC_SC/traphlor"

    avg_cov = find_average_coverage(odir, no_of_files)

    for i in range(1,no_of_files+1):
	if check_graph("%s/tmp/%s_%s" % (odir, GRAPH_FILE_PREFIX, i), "%s/tmp/%s_%s" % (odir, NODES_FILE_PREFIX, i), avg_cov):
            cmd = "%s  %s/tmp/%s_%s %s/tmp/%s_%s" % (traph_exec, odir, GRAPH_FILE_PREFIX, i, odir, NODES_FILE_PREFIX, i)
            p = subprocess.Popen(cmd, shell = True, stdout = None, stderr = None)
            p.wait()
        else:
            write_empty_solution_file(odir, i)
         
 
def convert_paths_no_annotation(odir, no_of_files, ifile,  debug):

   # Dummy opening of output_file for writing to erase previous contents, as the conversion method calls append to file
   temp_handle = open("%s/%s" % (odir, TRANSCRIPT_FILE), "w")
   temp_handle.close()

   gene_id = 1
   
   for i in range(1,no_of_files+1):
        used_id = convertPaths.convert_path_no_annotation("%s/%s" % (odir, TRANSCRIPT_FILE), "%s/tmp/%s_%s" % (odir, NODES_FILE_PREFIX, i), "%s/tmp/%s_%s.sol" % (odir, GRAPH_FILE_PREFIX, i), gene_id, ifile, debug)
        
        if used_id:
            gene_id = gene_id + 1    

def index_bam(fn):
    ofn = "%s.bai" % (fn)
    cmd = "%s index %s %s" % (SAMTOOLS, fn, ofn)
    p = subprocess.Popen(cmd, shell = True, stdout = None, stderr = subprocess.STDOUT)
    p.wait()
    if p.returncode != 0:
        status("Failed to run SAMTOOLS", True)
        sys.exit(1)


def main(options, args):
    
    odir = options.output
    ifile = options.input
    
    if not os.path.exists(ifile):
        sys.stderr.write("Input file doesn't exist.\n")
        sys.exit(2)

    # If no index, create it (will cause program to exit if there's no samtools in path, but pysam would do that anyways)
    if not os.path.exists("%s.bai" % ifile):
        index_bam(ifile)
     
    # Might throw exception that directory exists
    try:
        os.mkdir(odir)
    except:
        pass
        
    try:
        os.mkdir("%s/tmp" % (odir))
    except:
        pass

    print "[%s]: Creating graphs..." % (datetime.datetime.now().strftime(DATE_FORMAT))
    no_of_files = create_graphs_no_annotation(ifile, odir, options.alternative_threshold, options.debug)
    if no_of_files == 0:
        sys.stderr.write("Fatal error, no graph files were created. Aborting.\n")
        sys.exit(2) 
    print "[%s]: Searching for paths..." % (datetime.datetime.now().strftime(DATE_FORMAT))
    run_traph(odir, no_of_files, options)
    print "[%s]: Converting paths..." % (datetime.datetime.now().strftime(DATE_FORMAT))
    convert_paths_no_annotation(odir, no_of_files, ifile, options.debug)

#    cmd = "rm -r %s/tmp/" % (odir)
#    p = subprocess.Popen(cmd, shell=True, stdout=None, stderr=None)
#    p.wait()

    print "[%s] Running Traphlor completed." % (datetime.datetime.now().strftime(DATE_FORMAT))

if __name__ == "__main__":
    parser = optparse.OptionParser(description="Traphlor - Traph with Long Reads")

    mandatory_group = optparse.OptionGroup(parser, "Required options")

    mandatory_group.add_option("--input", "-i", help="Input file")
    mandatory_group.add_option("--output", "-o", help="Output directory")

    parser.add_option_group(mandatory_group)
 
    other_group = optparse.OptionGroup(parser, "Other options")

    other_group.add_option("--alternative-threshold", type="float", help="Slope threshold for searching for alternative transcripts starts/ends. The smaller, the more sensitive, but too small value can cause false positives if coverage is very uneven. Default 0.3", default=0.3)


    other_group.add_option("--debug", "-d", dest="debug", default = False, action='store_true', help="Debug mode")

    parser.add_option_group(other_group)

    (options, args) = parser.parse_args()

    if options.input == None:
        sys.stderr.write("Specify the input file with -i\n")
        parser.print_help()
        sys.exit(2)    
        
    if options.output == None:
        sys.stderr.write("Specify the output directory with -o\n")
        parser.print_help()
        sys.exit(2)
    
    main(options, args)
