def downscale(A, factor, axis):
	indexer = [slice(None) for _i in range(len(A.shape))]
	O = None
	for i in range(factor):
		indexer[axis] = slice(i, None, factor)
		iter_indexer = tuple(indexer)
		if O is None:
			O = A[iter_indexer]
		else:
			O += A[iter_indexer]
	return O / factor

def downscale_3d(A, factor):
	A = downscale(A, factor, 2)
	A = downscale(A, factor, 1)
	A = downscale(A, factor, 0)
	return A