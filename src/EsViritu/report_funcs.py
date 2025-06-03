from dash import Dash, html
import dash_ag_grid as dag

import json
import os

app = Dash()

with open("dash_docs/assets/example-data/data.json") as json_file:
    data = json.load(json_file)

columnDefs = [
    {"field": "symbol", "maxWidth": 120},
    {"field": "name", "minWidth": 250},
    {
        "field": "change",
        "cellRenderer": "agSparklineCellRenderer",
    },
    {
        "field": "volume",
        "type": "numericColumn",
        "maxWidth": 140,
    },
]

app.layout = html.Div(
    [
        dag.AgGrid(
            id="sparklines-basic-example",
            enableEnterpriseModules=True,
            licenseKey = os.environ['AGGRID_ENTERPRISE'],
            columnDefs=columnDefs,
            rowData=data,
            dashGridOptions={"animateRows": False}
        )
    ]
)


if __name__ == "__main__":
    app.run(debug=True)