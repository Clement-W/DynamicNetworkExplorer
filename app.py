from pathlib import Path
from typing import Callable, Mapping
import dash
from dash import html, dcc, Input, Output, State
from dash import ctx
import dash_cytoscape as cyto
import pandas as pd
import networkx as nx
import argparse

# argument parser for the tool
parser = argparse.ArgumentParser(description="Dynamic Gene Interaction Network Explorer")
parser.add_argument("--topn", type=int, default=5000,
                    help="Number of top interactions to display per time-point (default: 1000)")
parser.add_argument("--path", type=str, default="data",
                    help="Path to the directory containing the CSV files (default: data)")
parser.add_argument("--pattern", type=str, default="ranked_interactions_*.csv",
                    help="File pattern to match CSV files (default: ranked_interactions_*.csv)")
parser.add_argument("--scale", type=int, default=1200,
                    help="Scale factor for node layout (default: 1200). Increase/Decrease if you increase/decrease topn.")
parser.add_argument("--organism", type=str, default="mouse",
                    help="Organism for TF highlighting (default: mouse). Options: mouse, human, none")
parser.add_argument(
    "--subnetwork", type=str,  default="",
    help=(
        "Comma-separated list of source genes (gene1) to restrict the network to, "
        "e.g. --subnetwork Tcf21,Zeb2. If empty, the full network is used." ))

args = parser.parse_args()
print(f"Arguments: {args}")

if args.subnetwork:
    SUBNETWORK_GENES = {
        g.strip() for g in args.subnetwork.replace(";", ",").split(",") if g.strip()
    }
    print(f"Subnetwork mode ON. Source genes: {sorted(SUBNETWORK_GENES)}")
else:
    SUBNETWORK_GENES = set()
    print("Subnetwork mode OFF. Using full network.")

#  CONFIG 
TOP_EDGES = args.topn
DATA_DIR = Path(args.path)
FILE_PATTERN = args.pattern
SCALE = args.scale
ORGANISM = args.organism.lower()
LAYOUT_ALGO   = "kamada_kawai" 
SPRING_K      = 0.6
SEED          = 42

# colours
NODE_BASE     = "#dbdbdb"
NODE_HL       = "#FFD700"
EDGE_POS      = "#2CA02C"
EDGE_NEG      = "#D62728"
EDGE_POS_HL = "#66ff33"  # lighter of #2CA02C
EDGE_NEG_HL = "#ff3300"  # lighter of #D62728

#  LOAD CSVs 
print("Loading csv data…")
frames: dict[str, pd.DataFrame] = {}

for fp in sorted(DATA_DIR.glob(FILE_PATTERN)):
    label = fp.stem.split("_", 2)[-1]

    if SUBNETWORK_GENES:
        # Subnetwork mode:
        #   1) Read the full file 
        #   2) Keep only edges with gene1 in SUBNETWORK_GENES
        #   3) Take top TOP_EDGES within THAT filtered set
        df_full = pd.read_csv(
            fp,
            usecols=["gene1", "gene2", "strength", "sign"],
            dtype={ "gene1": "string",   "gene2": "string", "strength": "float32", "sign": "int8"},
            engine="c",
            memory_map=True,
        )

        df_pass1 = df_full[df_full["gene1"].isin(SUBNETWORK_GENES)].groupby("gene1", group_keys=False).head(TOP_EDGES)

        # Node set = seeds + all targets from pass 1
        node_set = set(SUBNETWORK_GENES) | set(df_pass1["gene2"].dropna().tolist())

        # For any gene1 in node_set, keep top K edges where gene2 is also in node_set
        df = (
            df_full[
                df_full["gene1"].isin(node_set)
                & df_full["gene2"].isin(node_set)
            ]
            .groupby("gene1", group_keys=False)
            .head(TOP_EDGES)
        )

        del df_full
        
        print(f"  {label}: {len(df)} edges (subnetwork-filtered)")
    else:
        # Read only the first N rows; assumes CSV is already ranked by |strength| desc
        df = pd.read_csv(
            fp,
            nrows=TOP_EDGES,  
            usecols=["gene1", "gene2", "strength", "sign"],
            dtype={"gene1": "string", "gene2": "string", "strength": "float32", "sign": "int8"}, 
            engine="c",
            memory_map=True
        )
        print(f"  {label}: {len(df)} edges")
    
    frames[label] = df


if SUBNETWORK_GENES:
    total_edges = sum(len(df) for df in frames.values())
    if total_edges == 0:
        raise RuntimeError(
            f"No edges found where gene1 is in {sorted(SUBNETWORK_GENES)}. "
            "Check your gene names or the input files."
        )
    
if not frames:
    raise RuntimeError(f"No CSVs like {FILE_PATTERN} in {DATA_DIR}")

# Load TF list for highlighting
# tf set for mouse downloaded on 9 July 2025 on https://guolab.wchscu.cn/AnimalTFDB4/#/
# tf set for human downloaded on 9 July 2025 on https://humantfs.ccbr.utoronto.ca/download.php
if ORGANISM in ("mouse", "human"):
    # get folder of app.py as the tf files are next to it in assets folder
    file_name = "Mus_musculus_TF" if ORGANISM == "mouse" else "human_TF_names_v_1.01.txt"
    tf_fp = Path(__file__).parent / "assets" / file_name
    if ORGANISM == "human":
        with open(tf_fp, 'r') as fh:
            # strip whitespace, ignore empty lines, store as a set for fast lookup
            TF_SET = {line.strip().upper() for line in fh if line.strip()}
    else:
        df_tf = pd.read_csv(tf_fp,sep="\t")["Symbol"]
        TF_SET = set(df_tf.str.upper().tolist())
else:
    TF_SET = set()
print(f"Loaded {len(TF_SET)} transcription factors for {ORGANISM}")

def is_tf(gene: str) -> bool:
    return gene.upper() in TF_SET

timepoints = list(frames)

#  Build union graph + layout 
print("Building graph layout…")
G = nx.DiGraph()
for df in frames.values(): #may take time
    G.add_edges_from(zip(df.gene1, df.gene2))

#make it general if we want to add more layout algorithms later
pos_fn: Mapping[str, Callable] = {
    "spring": lambda g: nx.spring_layout(g, k=SPRING_K, seed=SEED),
    "kamada_kawai": nx.kamada_kawai_layout,
}

# Compute positions on an undirected view to avoid KK degeneracy
H = G.to_undirected()         

#adapt scale in case there are too much nodes so we can manage how close they are. Increase scale if too many nodes
pos = pos_fn[LAYOUT_ALGO](H)
for p in pos.values():
    p[0] *= SCALE; p[1] *= SCALE

print(f"Graph has {len(G)} nodes and {G.number_of_edges()} unique edges")
#ranked by degree node
deg    = dict(G.degree())
maxdeg = max(deg.values()) or 1
neigh  = {n: list(G.neighbors(n)) for n in G.nodes()}
all_genes  = sorted(G.nodes())
gene_options = [{"label": g, "value": g} for g in all_genes]

nodes_template = [
    {"data": {"id": n, "label": n, "degree": deg[n], "is_tf": is_tf(n)},
     "position": {"x": float(x), "y": float(y)},
     "classes": "tf" if is_tf(n) else ""}
    for n, (x, y) in pos.items()
]

# edges per time‑point
def edges_for(df: pd.DataFrame):
    return [
        {"data": {"id": f"{r.gene1}-{r.gene2}",
                   "source": r.gene1, "target": r.gene2,
                   "strength": float(r.strength), "sign": int(r.sign)},
         "classes": ""}
        for r in df.itertuples(index=False)
    ]
edges_by_tp = {tp: edges_for(df) for tp, df in frames.items()}
print("Data loading complete.")

#  html stylesheet 
base_styles = [
    # nodes
    {"selector": "node", "style": {
        "content": "data(label)", "font-size": 12,
        "text-valign": "center", "color": "#222",
        "background-color": NODE_BASE,
        "border-width": 1, "border-color": "#555",
        "width": f"mapData(degree, 0, {maxdeg}, 15, 60)",
        "height": f"mapData(degree, 0, {maxdeg}, 15, 60)",
    }},
    # edges
    {"selector": "edge", "style": {
        "curve-style": "bezier", "opacity": .75,
        "width": "mapData(strength, 0, 10, 1, 6)", 
        "target-arrow-shape": "triangle",
    }},
    {"selector": "edge[sign = 1]", "style": {"line-color": EDGE_POS, 
                                             "target-arrow-color": EDGE_POS
                                             }},
    {"selector": "edge[sign = -1]", "style": {"line-color": EDGE_NEG,
                                              "target-arrow-color": EDGE_NEG
                                             }},

    
    {"selector": "node.tf", "style": {
        "shape": "diamond",
        "background-color": "#9bc9e4",
        "border-width": 2,
        "border-color": "#08519C",
    }},

    {"selector": "node.tf.hl", "style": {
    "background-color": NODE_HL,
    "border-color": "#000",
    "border-width": 3,
    "z-index": 10000
    }},

    # highlighted nodes & edges and bring to front
    {"selector": ".hl", "style": {
        "background-color": NODE_HL,
        "border-width": 3, "border-color": "#000", "z-index": 9999
    }},

    {"selector": ".hlEdge[sign = 1]", "style": {
        "curve-style": "bezier",               # <-- arrows need non-haystack
        "target-arrow-shape": "triangle",
        "line-color": EDGE_POS_HL,
        "opacity": 1,
        "width": 6,
        "z-index": 999,
        # arrow color
        "target-arrow-color": EDGE_POS_HL,
    }},
    {"selector": ".hlEdge[sign = -1]", "style": {
        "curve-style": "bezier",               # <-- arrows need non-haystack
        "target-arrow-shape": "triangle",
        "line-color": EDGE_NEG_HL,
        "opacity": 1,
        "width": 6,
        "z-index": 999,
        # arrow color
        "target-arrow-color": EDGE_NEG_HL,
    }},
]

print("Starting Dash app…")
#  Dash layout 
app = dash.Dash(__name__)
app.title = "Gene Interaction Explorer"
app.layout = html.Div([
    html.H3("Dynamic gene network explorer"),
    html.Div(id="tp-label", style={"margin": "0.6rem 0", "fontWeight": "bold"}),
    dcc.Slider(id="slider", min=0, max=len(timepoints)-1, step=None,
               marks={i: tp for i, tp in enumerate(timepoints)}, value=0,
               updatemode="drag"),
    html.Div([
        html.Label("Find gene"),
        dcc.Dropdown(id="gene-search", options=gene_options, placeholder="Type a gene…", clearable=True, style={"width":"420px"}),
    ], style={"marginBottom": "0.8rem"}),

   
    cyto.Cytoscape(id="net",
                   elements=nodes_template + edges_by_tp[timepoints[0]],
                   layout={"name": "preset"},
                   style={"width": "100%", "height": "750px"},
                   stylesheet=base_styles,
                   userZoomingEnabled=True,
                   userPanningEnabled=True,
                   responsive=True,),
    html.Div(
        f"Top {TOP_EDGES} interactions per time-point. Click or drag-select a node to highlight its neighbourhood; click blank space to clear.",
        style={"marginTop": "0.6rem", "fontSize": "0.85rem", "color": "#555"}),
    dcc.Store(id="store", data=None),  # holds current highlight id
])

# callback 
@app.callback(
    Output("net", "elements"), Output("tp-label", "children"), Output("store", "data"),
    Input("slider", "value"),        # time change
    Input("net", "selectedNodeData"), # fires on click & when selection cleared
    Input("gene-search", "value"),
    State("net", "elements"), State("store", "data"))

def update(tp_idx, selected, gene_query,cur_elems, stored_id):
    """Rebuild the visible sub-graph.

    * `hl_id` = currently highlighted node (may be None)
    * Highlight **only neighbours that exist in the
      *current* time-point’s edge list**.
    """
    #  decide what triggered the callback
    trig = ctx.triggered_id            # "slider" or "net"

    if trig == "net" and selected:     # a node was clicked
        hl_id = selected[0]["id"]
    elif trig == "net":                # click on blank canvas so clear highlight
        hl_id = None
    elif trig == "gene-search" and gene_query:  # search box changed
        hl_id = gene_query              # highlight the searched gene
        if hl_id not in G:              # if it doesn't exist, clear highlight
            hl_id = None
    else:                              # slider moved so previous highlight
        hl_id = stored_id

    # pick the edge list of the current time‑point
    tp    = timepoints[tp_idx]
    edges = edges_by_tp[tp]            # list[dict] for this frame only

    #  build the neighbour set for this time‑point 
    if hl_id:
        # neighbours = {
        #     e["data"]["target"] if e["data"]["source"] == hl_id else e["data"]["source"]
        #     for e in edges
        #     if hl_id in (e["data"]["source"], e["data"]["target"])
        # }
        neighbours = {e["data"]["target"] for e in edges if e["data"]["source"] == hl_id} # just consider outgoing edges for highlighting
    else:
        neighbours = set()

    # rebuild nodes, preserving positions 
    pos_map = {el["data"]["id"]: el.get("position", {})
               for el in cur_elems if "source" not in el["data"]}

    nodes = []
    for tmpl in nodes_template:
        n_id = tmpl["data"]["id"]
        n    = tmpl.copy()
        n["position"] = pos_map.get(n_id, n["position"])
  
        base_classes = tmpl.get("classes", "")
        hl_class = "hl" if hl_id and (n_id == hl_id or n_id in neighbours) else ""
        n["classes"] = " ".join(c for c in [base_classes, hl_class] if c)
        nodes.append(n)

    # add highlight class to edges for this frame 
    new_edges = []  
    for e in edges:
        ed = e.copy()
        s, t = ed["data"]["source"], ed["data"]["target"]
        #ed["classes"] = "hlEdge" if hl_id and (s == hl_id or t == hl_id) else ""
        ed["classes"] = "hlEdge" if hl_id and s == hl_id else "" # just consider outgoing edges for highlighting
        new_edges.append(ed)

    return nodes + new_edges, f"Current display: Network underlying the transition from time-point {tp} to {int(tp)+1} ", hl_id

print("App setup complete.")

if __name__ == "__main__":
    app.run(debug=True, port=8050,use_reloader=False) 