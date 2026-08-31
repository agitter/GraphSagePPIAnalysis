"""Step 0: verify the mirrored PPI files match the published GraphSAGE release.
Makes no assumptions about dtype, binarity, node ordering, or graph count."""
import json, hashlib, numpy as np, networkx as nx
from networkx.readwrite import json_graph

D = "/home/claude/dl/PPI-Inductive/ppi/"

for f in ["ppi-G.json", "ppi-feats.npy", "ppi-class_map.json", "ppi-id_map.json"]:
    h = hashlib.sha256(open(D + f, "rb").read()).hexdigest()
    print(f"sha256 {f:22s} {h}")

feats = np.load(D + "ppi-feats.npy")
print("\nfeats shape", feats.shape, "dtype", feats.dtype)

G_data = json.load(open(D + "ppi-G.json"))
print("top-level keys in G.json:", list(G_data.keys()))
G = json_graph.node_link_graph(G_data, edges="links")
print("G nodes", G.number_of_nodes(), "edges", G.number_of_edges(),
      "directed", G.is_directed(), "selfloops", nx.number_of_selfloops(G))

id_map = json.load(open(D + "ppi-id_map.json"))
class_map = json.load(open(D + "ppi-class_map.json"))
print("id_map entries", len(id_map), "class_map entries", len(class_map))

# What do the id_map keys/values look like? Are they gene identifiers?
ks = list(id_map.items())[:5]
print("id_map sample", ks)
print("id_map keys all numeric-string?",
      all(str(k).lstrip("-").isdigit() for k in id_map))
print("id_map is identity permutation?",
      all(int(k) == v for k, v in id_map.items()))

lab_lens = {len(v) for v in class_map.values()}
print("label vector lengths", lab_lens)

# split flags carried on nodes
n_val = sum(1 for n, d in G.nodes(data=True) if d.get("val"))
n_test = sum(1 for n, d in G.nodes(data=True) if d.get("test"))
n_train = G.number_of_nodes() - n_val - n_test
print(f"split from node flags: train={n_train} val={n_val} test={n_test}")
print("node attribute keys seen:",
      set(k for _, d in list(G.nodes(data=True))[:50] for k in d))

ncc = nx.number_connected_components(G)
print("connected components", ncc)
