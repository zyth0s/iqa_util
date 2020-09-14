# iqa_util
Some scripts to analyze IQA results. IQAJSON, we shall call it, is an approach to store
all IQA information using JSON format.

Surface files store, for each basin, the intersection of rays starting at the
nuclei with the basin boundary, in a standard HDF5 format.

### `iqanalysis.py`
Start IPython and execute `%run iqanalysis.py`. A window will open
requesting you to select the IQAJSON file. Then, you will return back to the
interpreter. A new variable `m` will be defined, that contains all extracted
information. `m.TAB` will show you all accessible properties or methods.

### `get_coords.py`
`xyz = get_coords()` and select the DGRID WFN file.

### `measure_error_surf.py`

Compare two surface files using one as a reference.

### `surfconv.py`
Converts basin surface files between Promolden (modular) and
ChemInt format.

### `surf_resolution.py`
Finds the angular resolution of a given surface file. Start IPython,
`%run surf_resolution.py` and select the desired surface file.

### `ray2pointscloud.py`
Convert each surface file to a point cloud of surface points (each basin will have
a separate file) with the nuclear coordinate as origin. `%run ray2pointscloud.py`
and select the associated DGRID WFN file. The output (`*.csv files) can be
visualized with Avizo.

### `equivalent_surfaces.py`
Generate symmetry equivalent surfaces. Needs the cell parameters,
symmetry operations, which is the asymmetric unit and their surfaces, nuclei coordinates, and filename.

### `nearest_neighbors.py`
Analyze interatomic distances.

### `gen_chmnt_input.py`
Generate ChemInt input from a template.

### `required_resources.jl`
Estimate required resources to do IQA decomposition.

### `plot_aom.py`
Plot Atomic Overlap Matrix (AOM) as a 2D image.

## Dependencies

* json package
* matplotlib
* numpy
* PyQt5 (not strictly necessary but better to have it)
* h5py
* sympy
