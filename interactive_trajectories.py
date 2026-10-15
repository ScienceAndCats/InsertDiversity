"""Standalone trajectory HTML with sample selection and offline Plotly assets."""

from __future__ import annotations

import html
from pathlib import Path


SAMPLE_CONTROLS = r"""
const graph = document.getElementById('{plot_id}');
const panel = document.createElement('section');
panel.className = 'sample-controls';
panel.setAttribute('aria-label', 'Trajectory sample selection');
const heading = document.createElement('h2');
heading.textContent = 'Samples';
panel.appendChild(heading);
const note = document.createElement('p');
note.textContent = 'Uncheck samples to gray out their points and connecting segments. Counts, frequencies, and barcode selection stay unchanged. All plotted barcodes are real in the first sample only; later real flags do not restrict plots.';
panel.appendChild(note);
const boxes = [];
const sampleNames = graph.layout.meta.sample_names;
const coloredIndices = [];
const originalY = [];
graph.data.forEach((trace, index) => {
    if (index % 2 === 1) {
        coloredIndices.push(index);
        originalY.push(Array.from(trace.y));
    }
});
function updateSamples() {
    const enabled = boxes.map(box => box.checked);
    const values = originalY.map(row => row.map((value, i) => enabled[i] ? value : null));
    return Plotly.restyle(graph, {y: values}, coloredIndices);
}
['Select all', 'Gray out all'].forEach((text, index) => {
    const button = document.createElement('button');
    button.type = 'button';
    button.textContent = text;
    button.addEventListener('click', () => {
        boxes.forEach(box => { box.checked = index === 0; });
        updateSamples();
    });
    panel.appendChild(button);
});
const list = document.createElement('div');
list.className = 'sample-list';
sampleNames.forEach(sample => {
    const label = document.createElement('label');
    const box = document.createElement('input');
    box.type = 'checkbox';
    box.checked = true;
    box.addEventListener('change', updateSamples);
    label.appendChild(box);
    label.appendChild(document.createTextNode(sample));
    list.appendChild(label);
    boxes.push(box);
});
panel.appendChild(list);
graph.parentNode.insertBefore(panel, graph);
"""


def save_interactive_trajectory_plot(barcodes, freqs, outpath: Path, title: str, cfg: dict):
    import plotly.graph_objects as go
    import plotly.io as pio
    from plotly.colors import qualitative

    samples = list(freqs.columns)
    pseudo = float(cfg["pseudocount_for_log_plots"])
    fig = go.Figure()
    for i, barcode in enumerate(barcodes):
        frequencies = freqs.loc[barcode].astype(float).tolist()
        values = [value + pseudo for value in frequencies]
        hover = [html.escape(f"{barcode} | {sample}") for sample in samples]
        for color in ("#b8b8b8", qualitative.Plotly[i % len(qualitative.Plotly)]):
            fig.add_trace(go.Scattergl(
                x=list(range(len(samples))), y=values, mode="lines+markers",
                line={"color": color, "width": 1}, marker={"color": color, "size": 5},
                name=html.escape(barcode), text=hover, customdata=frequencies,
                hovertemplate="%{text}<br>Frequency: %{customdata:.6g}<extra></extra>",
                showlegend=False, connectgaps=False,
            ))
    fig.update_layout(
        title=html.escape(title), template="plotly_white", height=650,
        xaxis={"tickmode": "array", "tickvals": list(range(len(samples))),
               "ticktext": [html.escape(sample) for sample in samples], "tickangle": -90},
        yaxis={"title": "Barcode frequency", "type": "log" if cfg.get("log_scale_y_axis", True) else "linear"},
        meta={"sample_names": samples, "real_barcode_sample": "first"},
        uirevision="sample-selection",
    )
    # Embed the JavaScript bundle: these files can be opened with no server/network.
    fragment = pio.to_html(
        fig, full_html=False, include_plotlyjs=True,
        div_id="barcode-trajectories", post_script=SAMPLE_CONTROLS,
        config={"responsive": True, "displaylogo": False},
    )
    document = """<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Barcode trajectories</title><style>
body { font: 15px system-ui, sans-serif; margin: 24px; color: #222; }
.sample-controls { padding: 16px; background: #f5f5f5; border-radius: 8px; }
.sample-controls h2 { margin: 0; }
.sample-controls button { margin: 0 8px 12px 0; padding: 6px 12px; }
.sample-list { display: flex; flex-wrap: wrap; gap: 12px; }
.sample-list label { display: inline-flex; align-items: center; gap: 4px; }
</style></head><body>""" + fragment + "</body></html>"
    outpath.write_text(document, encoding="utf-8")
