def unclog_umap_caching(x, file='umap_unclog_jhsd7.temp'):
    '''
    This is an ugly fix to allow fit_transform into an UMAP object,
    it gets rids of "TypeError: cannot pickle 'weakref' object" error
    Provide a UMAP object as x.
    IMPORTANT: this function will break the UMAP object, recompute
    afterwards.
    '''
    import pickle
    import tempfile
    #  umap_ref._knn_search_index._visited = np.zeros((umap_ref._knn_search_index._raw_data.shape[0] // 8) + 1, dtype=np.uint8, order="C")
    for i in range(3):
        try:
            with tempfile.TemporaryFile() as f:
                pickle.dump(x, f)
        except:
            pass