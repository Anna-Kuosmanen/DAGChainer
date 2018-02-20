import matplotlib.pyplot as plt
import math
import time
import solvers
import timeit

TEST_GRAPH_NAME = "test_graph"

def plot(k_values, n_values , m_values, title, filename):
  plt.clf()
  m = 0.5
  results_decomp = [0]*len(values)
  results_normal = [0]*len(values)
  index = 0
  for n in n_values:
    solvers.generate_k_path_graph(k_values[index], n, m_values[index], TEST_GRAPH_NAME)
    #solvers.generate_dag(n, 0.9, TEST_GRAPH_NAME)
    print "GRAPH OF SIZE: " +str(n)
    start = time.time()
    decomposed_sum = solvers.solve_with_decomposition(TEST_GRAPH_NAME)
    print decomposed_sum
    decomp_end = time.time()
    normal_sum = solvers.solve_without_decomposition(TEST_GRAPH_NAME)
    print normal_sum
    end = time.time()
    results_decomp[index] = decomp_end - start
    print "decomp: " + str(results_decomp[index])
    results_normal[index] = end - decomp_end
    print "normal: " + str(results_normal[index])
    index += 1
  plt.plot(n_values, results_decomp, 'r', label="with decomposition")
  plt.plot(n_values, results_normal, 'b', label="without decomposition")
  plt.ylabel("time in seconds")
  plt.xlabel("n (path length)")
  plt.title(title)
  plt.legend()
  plt.savefig(filename)

def plot_alt(k_values, n_values , m_values, title, filename):
  plt.clf()
  m = 0.5
  results_decomp = [0]*len(values)
  results_normal = [0]*len(values)
  index = 0
  for n in n_values:
    solvers.generate_k_path_graph_alt(k_values[index], n, m_values[index], TEST_GRAPH_NAME)
    #solvers.generate_dag(n, 0.9, TEST_GRAPH_NAME)
    print "GRAPH OF SIZE: " +str(n)
    start = time.time()
    decomposed_sum = solvers.solve_with_decomposition(TEST_GRAPH_NAME)
    print decomposed_sum
    decomp_end = time.time()
    normal_sum = solvers.solve_without_decomposition(TEST_GRAPH_NAME)
    print normal_sum
    end = time.time()
    results_decomp[index] = decomp_end - start
    print "decomp: " + str(results_decomp[index])
    results_normal[index] = end - decomp_end
    print "normal: " + str(results_normal[index])
    index += 1
  plt.plot(n_values, results_decomp, 'r', label="with decomposition")
  plt.plot(n_values, results_normal, 'b', label="without decomposition")
  plt.ylabel("time in seconds")
  plt.xlabel("n (path length)")
  plt.title(title)
  plt.legend()
  plt.savefig(filename)




values = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000]
#values = [10, 100]
def k_10():
  k_values = len(values)*[10]
  plot(k_values, values, values, "k = 10, m = n", "k10.png")

def k_logn():
  k_values = map(lambda n: max(10, math.ceil(math.log(n))), values)
  plot(k_values, values, values, "k = log(n), m = n", "klogn.png")

def k_sqrt():
  k_values = map(lambda n: max(10, math.ceil(math.sqrt(n))), values)
  plot(k_values, values, values, "k = sqrt(n), m = n", "ksqrtn.png")

def k_div10():
  k_values = map(lambda n: n/10, values)
  plot(k_values, values, values, "k = n/10, m = n", "kdiv10.png")

def k_10_alt():
  k_values = len(values)*[10]
  plot_alt(k_values, values, values, "k = 10, m = n", "k10_alt.png")

def k_sqrt_alt():
  k_values = map(lambda n: max(10, math.ceil(math.sqrt(n))), values)
  plot_alt(k_values, values, values, "k = sqrt(n), m = n", "ksqrtn_alt.png")

def k_div10_alt():
  k_values = map(lambda n: n/10, values)
  plot_alt(k_values, values, values, "k = n/10, m = n", "k_div10_alt.png")

k_10()
k_sqrt()
k_10_alt()
k_sqrt_alt()
k_div10()
k_div10_alt()