import argparse
import pandas as pd
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
from Bio import ExPASy
from Bio import SwissProt
matplotlib.use('Agg')  # Use the Agg backend
import ssl
import urllib.request

# Specify the path to the CA certificate bundle
ca_cert_path = '/etc/ssl/certs/ca-bundle.crt'

# Set the CA certificate file path for SSL verification
ssl._create_default_https_context = ssl.create_default_context(cafile=ca_cert_path)

def get_unique_nodes(df1, df2):
    df1_interactions = df1[["complex_name", "ranking_confidence"]]
    df2_interactions = df2[["complex_name", "ranking_confidence"]]
    # Round the values in the "ranking_confidence" column to 3 decimal places
    df1_interactions["ranking_confidence"] = df1_interactions["ranking_confidence"].round(3)
    df2_interactions["ranking_confidence"] = df2_interactions["ranking_confidence"].round(3)

    unique_nodes = []
    for i in df1_interactions["complex_name"]:
        protein1, protein2 = i.split("_")
        unique_nodes.append(protein1)
        unique_nodes.append(protein2)
    for i in df2_interactions["complex_name"]:
        protein1, protein2 = i.split("_")
        unique_nodes.append(protein1)
        unique_nodes.append(protein2)
    unique_nodes = list(set(unique_nodes))
    for i in range(len(unique_nodes)):
        unique_nodes[i] = get_protein_name(unique_nodes[i])
    return unique_nodes, df1_interactions, df2_interactions

def generate_hex_color(value):
    # Ensure the input value is within the specified range and normalize it
    normalized_value = max(0.0, min(1.0, (value - 0.7) / 0.3))
    
    # Calculate the RGB values for the gradient from blue to red
    blue = int((1 - normalized_value) * 255)
    red = int(normalized_value * 255)
    green = 0
    
    # Convert RGB values to hexadecimal format
    hex_color = "#{:02X}{:02X}{:02X}".format(red, green, blue)
    
    return hex_color

def get_protein_name(uniprot_id):
    try:
        url = f'https://www.uniprot.org/uniprot/{uniprot_id}.txt'
        req = urllib.request.Request(url)
        
        with urllib.request.urlopen(req, cafile=ca_cert_path) as f:
            record_lines = f.read().decode('utf-8')

        ID = record_lines.splitlines()[0].split()[1]

        
        if not ID:
            print(f"Gene name not found for UniProt ID {uniprot_id}. Using entry name instead.")
            return uniprot_id
        else:
            return ID
    except Exception as e:
        print(f"Error fetching data for UniProt ID {uniprot_id}: {str(e)}")
        return None

def generate_graph(unique_nodes, df1_interactions, output_path):
    G = nx.Graph()

    # Add nodes to the graph
    G.add_nodes_from(unique_nodes)

    # Add edges with weights based on ranking_confidence from df1
    for _, row in df1_interactions.iterrows():
        protein1, protein2 = row["complex_name"].split("_")
        protein1_name = get_protein_name(protein1)
        protein2_name = get_protein_name(protein2)
        weight = round(row["ranking_confidence"], 3)
        G.add_edge(protein1_name, protein2_name, weight=weight)

    # Draw and save the graph
    pos = nx.spring_layout(G, k=10, seed=42, iterations=100)
    edge_labels = nx.get_edge_attributes(G, 'weight')

    plt.figure(figsize=(10, 6))
    
    # Draw nodes with larger size
    nx.draw_networkx_nodes(G, pos, node_size=800, node_color="skyblue", label=None)

    # Draw edges with larger width for self-loops
    for edge, weight in edge_labels.items():
        nx.draw_networkx_edges(G, pos, edgelist=[edge], edge_color=generate_hex_color(weight), width=5)

    # Draw labels for all nodes
    nx.draw_networkx_labels(G, pos, font_size=10, font_color="black", font_weight="bold", labels={node: node for node in G}, verticalalignment="center", horizontalalignment="center")
    # Draw other edges and their labels
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=10, font_color="DarkSlateGray", font_weight="bold")

    for edge, weight in edge_labels.items():
        if edge[0] == edge[1]:
            x, y = pos[edge[0]]
            # Shift the label position outside the node
            label_pos = (x - 0.12, y + 0.12)
            pos[edge[0]] = label_pos
            nx.draw_networkx_edge_labels(G, pos, edge_labels={edge: edge_labels[(edge[0],edge[1])]}, font_size=10, bbox=dict(alpha=0), font_color="DarkSlateGray", font_weight="bold")

    plt.title("Protein-Protein Interactions")

    # Create a custom colormap for the legend (blue to purple to red)
    custom_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'custom_cmap', [(0, 0, 1), (0.5, 0, 0.5), (1, 0, 0)], N=256
    )

    # Create colorimetric scale legend with the custom colormap
    sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=plt.Normalize(vmin=0.7, vmax=1.0))
    sm._A = []  # empty array for the data range
    cbar = plt.colorbar(sm, orientation='vertical', fraction=0.03, pad=0.1)
    cbar.set_label('Ranking Confidence (0.8*ipTM + 0.2*pTM)')

    # Automatically adjust layout to prevent overlap
    plt.tight_layout()

    plt.savefig(output_path)
    print(f"Graph saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description='Generate a graph from two CSV files.')
    parser.add_argument('csv1', type=str, help='Path to the first CSV file')
    parser.add_argument('csv2', type=str, help='Path to the second CSV file')
    parser.add_argument('output', type=str, help='Path to the output PNG file')

    args = parser.parse_args()

    df1 = pd.read_csv(args.csv1)
    df2 = pd.read_csv(args.csv2)

    unique_nodes, df1_interactions, df2_interactions = get_unique_nodes(df1, df2)

    generate_graph(unique_nodes, df1_interactions, args.output)

if __name__ == "__main__":
    main()
