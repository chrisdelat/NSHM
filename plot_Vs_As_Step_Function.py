import numpy as np

def get_depth_versus_Vs_for_step_plot(layerThicks, layerVs, startDepth=0):

    #get the depth to bottom of layers
    depth_BoL = np.cumsum(layerThicks)

    depth4plot = np.concatenate([[x, x] for x in depth_BoL])
    depth4plot = np.insert(depth4plot, 0, 0)
    depth4plot = np.delete(depth4plot, -1)

    depth4plot += startDepth

    # make array of Vs for plotting stepwise Vs

    Vs4plot = np.concatenate([[x, x] for x in layerVs])

    return Vs4plot, depth4plot