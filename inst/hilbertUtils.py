import numpy as np
import math

from hilbert import decode

def hilbertOrder(s, d):
  """
  Compute order of the d-dimensional Hilbert curve which corresponds to
  scale s.

  Parameters
  ----------
  s : int
      Current scale of the binary tree
  d : int
      Dimension of the Hilbert curve
  Returns
  -------
  out : int

  """
  return math.ceil(s / d)

def nHilbert(s, d):
  """
  Number of Hilbert points in each rectangle at scale s for dimension d
  """
  return abs((d-s)%d)


def idx_sh(s, h, d, type = "R"):
  """
  Indices of elements to compare at node (s,h) for dimension d
  """


  if type == "R":
    return ( (h-1) * 2**(nHilbert(s+1, d) + 1) + 2**(nHilbert(s+1,d)),
             (h-1) * 2**(nHilbert(s+1, d) + 1) + 2**(nHilbert(s+1,d)) + 1
            )

  elif type == "python":
    return ( (h-1) * 2**(nHilbert(s+1, d) + 1) + 2**(nHilbert(s+1,d))-1,
             (h-1) * 2**(nHilbert(s+1, d) + 1) + 2**(nHilbert(s+1,d))
            )


def idx_s(s, d, type = "R"):
  """
  List of elements to compare at scale s for dimension d
  """
  if type == "R":
    return [idx_sh(s, i, d, type="R") for i in range(1, 2**s + 1)]

  elif type == "python":
    return [idx_sh(s, i, d, type="python") for i in range(1, 2**s + 1)]


def decodeNorm(pts, d, k):
  """
  Compute Hilbert points on a [0,1]^d hypercube

  Parameters
  ----------
  pts : array-like, integers
  d : int, dimension of the hypercube
  k : int, order of the Hilbert curve

  Returns
  -------
  out : ndarray of dimension (len(pts), d), floats
        Coordinates of Hilbert points on the [0,1]^d hypercube

  """
  norm = 2**k - 1                       # Normalize to unit square
  transl = 1 / (2**(k + 1))             # Center to Hilbert's original
  scale = (1 - 2*transl)                # Scale to Hilbert's original
  out = decode(pts, d, k) / norm * scale + np.array([transl for _ in range(d)])
  return out



def comparePts(x, y):
  """
  Compare two points on the hilbert curve and return the index of the dimension
  which is not shared by the points.

  Parameters
  ----------
  x : array-like of dimension d, floats
  y : array-like of dimension d, floats

  Returns
  -------
  idx : int,
        index of dimension where the two points are different
  """
  idx = np.where(x != y)
  return idx[0][0]



def splitNodeBounds(bounds, j, y):
  """
  Split the subcube {x : l <= x <= u} along the j-th dimension, for 0 <= j < d.

  Parameters
  ----------
  bounds : tuple of array-like
            Contains (l, u) where l are current node lower bounds and u are
            current node upper bounds.
  s : int
      Current scale at which to split the cube
  h : int
      Current node index
  y : ndarray,
      First of the two points being compared on the hilbert curve

  Returns
  -------
  out : list of tuple containing new bounds
        [(l1, u1), (l2, u2)]
        where
          l1, u1: bounds for node (s+1, 2h - 1)
          l2, u2: bounds for node (s+1, 2h)
  """
  l = bounds[0]
  u = bounds[1]
  d = len(l)
  length = u[j] - l[j]                  # Length of segment which joins u[j], l[j]
  l1 = [x for x in l]
  l2 = [x for x in l]
  u1 = [x for x in u]
  u2 = [x for x in u]

  # WRONG, must check in which subcube the Hilbert point is to preserve
  # ordering of the subcubes. Use coordinates of Hilbert points to
  # determine order of return
  u1[j] = u[j] - length/2               # Update lower bound
  l2[j] = u[j] - length/2               # Update upper bound

  #--------------------------------------------------#
  # Check where the first Hilbert point lies and return accordingly

  if np.all(y >= l1) and np.all(y <= u1):
    return [(l1, u1), (l2, u2)]

  else:
    return [(l2, u2), (l1, u1)]


