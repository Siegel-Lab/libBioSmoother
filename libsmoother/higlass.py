import ipywidgets as widgets


def display_interface(filename):
    def map_q_val_change(change):
        print(change)

    def normalization_val_change(change):
        print(change)

    def ddd_val_change(change):
        print(change)

    mapping_quality = widgets.IntSlider(max=256, continuous_update=False)
    mapping_quality.observe(map_q_val_change, names="value")
    mapping_quality_box = widgets.HBox(
        [widgets.Label("Mapping Quality Minimum"), mapping_quality]
    )

    normalization = widgets.Dropdown(
        options=[("None", 1), ("Iterative Correction", 2)],
        value=1,
    )
    normalization.observe(normalization_val_change, names="value")
    normalization_box = widgets.HBox([widgets.Label("Normalization"), normalization])

    ddd = widgets.Checkbox(
        value=False,
        description="Remove Distance Dependant Decay",
        disabled=False,
        indent=False,
    )
    ddd.observe(ddd_val_change, names="value")

    display(mapping_quality_box, normalization_box, ddd)
