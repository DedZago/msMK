from inst.hilbertUtils import hilbertOrder, decodeNorm, splitNodeBounds, idx_sh, comparePts
import csv


# Generate a sequence of Hilbert integers.

def hilbertSplit(d, smax):
  """
  Compute upper and lower bounds for the p-dimensional hypercube dyadic
  partition and save them in a .csv file.

  Parameters
  ----------
  p : integer,
      Dimension of the hypercube to partition.

  smax : integer,
      Maximum depth of the binary tree.

  Returns
  -------
  None

  Side effects
  ------------
  Creates a file called `bounds.csv` where the bounds are stored, in order
  to be later loaded by the R functions.

  """
  thrs = {i:[] for i in range(smax+1)}

  l = [0.0 for i in range(d)]
  u = [1.0 for i in range(d)]
  bounds = [l, u]

  s = 0
  thrs[s].append(bounds)

  for s in range(smax):
    for h in range(1, 2**s + 1):
      bounds = thrs[s][h-1]
      k = hilbertOrder(s+1, d)                    # Compute Hilbert order
      idx = idx_sh(s, h, d, type = "python")      # Get next indices to compare
      y1 = decodeNorm(idx[0], d, k)[0]            # First vector
      y2 = decodeNorm(idx[1], d, k)[0]            # Second vector
      j = comparePts(y1, y2)                      # Find unequal dimension
      newBounds = splitNodeBounds(bounds, j, y1)  # Split current nodes
      thrs[s+1] = thrs[s+1] + newBounds

  with open('bounds.csv', 'w') as f:
    # using csv.writer method from CSV package
    write = csv.writer(f)
    for s in range(smax+1):
      for h in range(1, 2**s + 1):
        bounds = thrs[s][h-1]
        write.writerow(bounds[0])
        write.writerow(bounds[1])



