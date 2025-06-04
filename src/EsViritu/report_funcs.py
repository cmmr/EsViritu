from dash import Dash, html
import dash_ag_grid as dag

import json
import os
import polars as pl
import tempfile
import numpy as np

## make json, load json, make report with sparklines
def make_dash_aggrid_html_report(windows: 'pl.DataFrame', meta: 'pl.DataFrame', output_html: str = 'report.html') -> str:
    """
    Create an HTML file with Dash from two polars DataFrames, rendering sparklines using the Sparks-Bar-Medium font.
    - windows: DataFrame with columns including Accession and average_coverage
    - meta: DataFrame with columns Accession, mean_coverage, read_count, family, genus, species, RPKMF
    - output_html: Path to save the HTML file
    """
    

    # Helper: map normalized values to Sparks-Bar-Medium glyphs (8 bars: U+2581 to U+2588)
    def sparkline(vals):
        if not vals:
            return ''
        # Unicode block chars: 0x2581 (lowest) to 0x2588 (highest)
        bars = [chr(code) for code in range(0x2581, 0x2589)]
        arr = np.array(vals, dtype=float)
        if np.all(arr == arr[0]):  # flat
            idxs = [0 for _ in arr]
        else:
            arr_norm = (arr - arr.min()) / (arr.max() - arr.min())
            idxs = (arr_norm * (len(bars) - 1)).round().astype(int)
        return ''.join([bars[i] for i in idxs])

    # Group windows by Accession, flatten average_coverage to cov_list
    grouped = (
        windows.groupby('Accession')
        .agg([
            pl.col('average_coverage').list().alias('cov_list')
        ])
    )
    # Join with meta on Accession
    joined = grouped.join(meta, on='Accession', how='inner')
    # Select and order columns
    columns = [
        'Accession', 'mean_coverage', 'read_count', 'family', 'genus', 'species', 'RPKMF', 'cov_list'
    ]
    joined = joined.select([col for col in columns if col in joined.columns])
    # Save to temp JSON file
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.json', delete=False) as tmp:
        json_path = tmp.name
        joined.write_json(json_path, row_oriented=True)
    # Load JSON
    with open(json_path) as f:
        data = json.load(f)

    # Build Dash app with a pure HTML table
    from dash import Dash, html
    app = Dash(__name__, assets_folder='assets')
    # Table header
    table_header = html.Thead([
        html.Tr([
            html.Th(col) for col in columns
        ])
    ])
    # Table rows
    table_rows = []
    for row in data:
        row_cells = []
        for col in columns:
            val = row.get(col, "")
            if col == "cov_list":
                spark = sparkline(val)
                cell = html.Td(
                    spark,
                    style={"fontFamily": "Sparks-Bar-Medium", "fontSize": "18px", "letterSpacing": "1px"},
                    className="bar-medium"
                )
            else:
                cell = html.Td(val)
            row_cells.append(cell)
        table_rows.append(html.Tr(row_cells))
    table_body = html.Tbody(table_rows)

    app.layout = html.Div([
        html.H2("Coverage Report"),
        html.Link(rel="stylesheet", href="/assets/sparks-fonts.css"),
        html.Table([
            table_header,
            table_body
        ], style={"width": "100%", "borderCollapse": "collapse"})
    ])

    # Save as HTML (without launching server)
    from dash._utils import AttributeDict
    app._setup_server()
    app._validate_layout()
    rendered = app._generate_html(
        AttributeDict({"request": None}),
        app._layout_value(),
        app._get_assets_url(),
        app._collect_and_register_resources(),
    )
    with open(output_html, 'w') as f:
        f.write(rendered)
    print(f"HTML report saved to {output_html}")
    # Remove temp file
    os.remove(json_path)
    return output_html


