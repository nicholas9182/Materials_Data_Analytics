import _plotly_utils.basevalidators


class XValidator(_plotly_utils.basevalidators.NumberValidator):
    def __init__(
        self, plotly_name="x", parent_name="layout.scene.aspectratio", **kwargs
    ):
        super(XValidator, self).__init__(
            plotly_name=plotly_name,
            parent_name=parent_name,
            edit_type=kwargs.pop("edit_type", "plot"),
            implied_edits=kwargs.pop("implied_edits", {"^aspectmode": "manual"}),
            min=kwargs.pop("min", 0),
            **kwargs,
        )