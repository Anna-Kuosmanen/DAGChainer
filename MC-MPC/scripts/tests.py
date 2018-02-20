import solvers

TEST_GRAPH_NAME = "test_graph"

def test_sum(k, n, m):

  solvers.generate_k_path_graph(k, n , m, TEST_GRAPH_NAME)

  decomposed_sum = solvers.solve_with_decomposition(TEST_GRAPH_NAME)
  normal_sum = solvers.solve_without_decomposition(TEST_GRAPH_NAME)
  print decomposed_sum
  print normal_sum
  if( decomposed_sum != normal_sum):
    print "Test failed"
    return 0

  print "Test passed."
  return 1


#for i in range(1, 10):
test_sum(2, 3, 10)