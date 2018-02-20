import glob
import os
import time

BIN_FOLDER = os.getenv("BIN_FOLDER", "../bin/")
OUTPUT_FOLDER = os.getenv("OUTPUT_FOLDER", "output/")

def generate_k_path_graph(k, n, m, output_file):
  os.system(BIN_FOLDER + "generator {0} {1} {2} {3}".format(output_file, k, n, m))

def generate_k_path_graph_alt(k, n, m, output_file):
  os.system(BIN_FOLDER + "generator {0} {1} {2} {3} --alt".format(output_file, k, n, m))

def generate_dag(n, prob, output_file):
  os.system(BIN_FOLDER + "generator {0} {1} {2}".format(output_file, n, prob))

def decompose(input_file, output_folder):
  os.system(BIN_FOLDER + "decompose {0} {1}".format(input_file, output_folder))

def create_clean_folders():
  # create necessary directories, if they don't exist
  os.system("mkdir -p {0}decomposition".format(OUTPUT_FOLDER))
  os.system("mkdir -p {0}results".format(OUTPUT_FOLDER))

  # cleanup possible previous decompositions & results
  os.system("rm {0}decomposition/*".format(OUTPUT_FOLDER))
  os.system("rm {0}results/*".format(OUTPUT_FOLDER))

def solve_with_decomposition(input_file):

  os.system("mkdir -p {0}results".format(OUTPUT_FOLDER))
  os.system("{0}mc-mpc {1} {2} --decomp".format(BIN_FOLDER, input_file, OUTPUT_FOLDER+"results/"+"result"))

  arc_sum = 0
  f = open(OUTPUT_FOLDER+"results/result", 'r')
  for line in f:
    line = line.split()
    amount = int(line[1])
    weight = int(line[2])
    arc_sum += amount*weight
  f.close()
  
  return arc_sum #arc_sum

def solve_without_decomposition(input_file):

  os.system("mkdir -p {0}results".format(OUTPUT_FOLDER))
  os.system("{0}mc-mpc {1} {2}".format(BIN_FOLDER, input_file, OUTPUT_FOLDER+"results/"+"result"))

  arc_sum = 0
  f = open(OUTPUT_FOLDER+"results/result", 'r')
  for line in f:
    line = line.split()
    amount = int(line[1])
    weight = int(line[2])
    arc_sum += amount*weight
  f.close()

  return arc_sum
