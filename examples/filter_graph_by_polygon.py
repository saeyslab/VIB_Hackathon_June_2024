
# %%
import spatialdata as sd
import squidpy as sq

sdata = sd.read_zarr("/Users/giovanni.palla/Downloads/WT_trajectory.zarr")

# %%
sq.gr.spatial_neighbors(
    sdata,
    elements_to_coordinate_systems={"segmentation_mask_inverted_maize_model_diameter_150_boundaries":"global"},
    table_key="table",
    coord_type="generic",
)
graph_original = sdata["table"].obsp["spatial_connectivities"].copy()
# %%
sq.gr._build.mask_graph(
    sdata,
    "table",
    "trajectory",
)
# %%
graph_filtered = sdata["table"].obsp["spatial_connectivities"].copy()

# %%
import matplotlib.pyplot as plt
import spatialdata_plot
from networkx import Graph
from networkx.drawing import draw_networkx_edges

edges_width=1
edges_color="grey"

fig, ax = plt.subplots(1,1,figsize = (8, 8))
sdata.pl.render_shapes("segmentation_mask_inverted_maize_model_diameter_150_boundaries",color="GA2ox3").pl.show(ax=ax)
g = Graph(graph_filtered)
coords = sdata.tables["table"].obsm["spatial"]
edge_collection = draw_networkx_edges(
    g, coords, width=edges_width, edge_color=edges_color, arrows=False, ax=ax
)
ax.add_collection(edge_collection)