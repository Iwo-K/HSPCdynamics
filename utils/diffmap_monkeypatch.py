import numpy as np
import scipy

def compute_eigen_fix(
    self,
    n_comps: int = 15,
    sym = None,
    sort = 'decrease',
    random_state = 0
):
    """\
    Compute eigen decomposition of transition matrix.

    Parameters
    ----------
    n_comps
        Number of eigenvalues/vectors to be computed, set `n_comps = 0` if
        you need all eigenvectors.
    sym
        Instead of computing the eigendecomposition of the assymetric
        transition matrix, computed the eigendecomposition of the symmetric
        Ktilde matrix.
    random_state
        A numpy random seed

    Returns
    -------
    Writes the following attributes.

    eigen_values : numpy.ndarray
        Eigenvalues of transition matrix.
    eigen_basis : numpy.ndarray
         Matrix of eigenvectors (stored in columns).  `.eigen_basis` is
         projection of data matrix on right eigenvectors, that is, the
         projection on the diffusion components.  these are simply the
         components of the right eigenvectors and can directly be used for
         plotting.
    """
    np.set_printoptions(precision=10)
    if self._transitions_sym is None:
        raise ValueError('Run `.compute_transitions` first.')
    matrix = self._transitions_sym
    # compute the spectrum
    if n_comps == 0:
        evals, evecs = scipy.linalg.eigh(matrix)
    else:
        n_comps = min(matrix.shape[0] - 1, n_comps)
        # ncv = max(2 * n_comps + 1, int(np.sqrt(matrix.shape[0])))
        ncv = None
        which = 'LM' if sort == 'decrease' else 'SM'
        # it pays off to increase the stability with a bit more precision
        matrix = matrix.astype(np.float64)

        # Setting the random initial vector
#         random_state = check_random_state(random_state)
        np.random.seed(random_state)
        print('running patched function!')
        v0 = np.random.randn((matrix.shape[0]))

        evals, evecs = scipy.sparse.linalg.eigsh(
            matrix, k=n_comps, which=which, ncv=ncv, v0=v0
        )
        evals, evecs = evals.astype(np.float32), evecs.astype(np.float32)
    if sort == 'decrease':
        evals = evals[::-1]
        evecs = evecs[:, ::-1]
#     logg.info(
#         '    eigenvalues of transition matrix\n'
#         '    {}'.format(str(evals).replace('\n', '\n    '))
#     )
#     if self._number_connected_components > len(evals) / 2:
#         logg.warning('Transition matrix has many disconnected components!')
    self._eigen_values = evals
    self._eigen_basis = evecs
