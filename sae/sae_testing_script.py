import torch
import pandas as pd
from esm.models import ESMC_6B # Load the ESMC base model
from esm.sdk.api import SAEHead

# 1. Load the model and its pre-trained Sparse Autoencoder head
print("Loading ESMC-6B and SAE head...")
model = ESMC_6B.from_pretrained()
sae_head = SAEHead.from_pretrained(layer=30) # SAEs map specific hidden layers

# 2. Define your incoming metagenomic query sequence
query_sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQF"

# 3. Step 1: Identify the top features activated by the sequence
with torch.no_grad():
    # Pass sequence through the pLM to get hidden states
    hidden_states = model.compute_hidden_states(query_sequence)
    
    # Pass hidden states through the SAE head to get sparse activations
    sae_output = sae_head(hidden_states)
    
    # Extract the top-K features ranked by normalized activation
    # Each protein activates a small subset of the 16,384 total feature space
    top_features = torch.topk(sae_output.activations, k=5)
    active_feature_ids = top_features.indices.tolist()[0]
    print(f"Top Activated SAE Features for query: {active_feature_ids}")

# 4. Step 2: Query the local SAE index tables to find the cluster
# Load your downloaded representative index map
print("Querying index tables...")
df_representatives = pd.read_parquet("representative_proteins.parquet")

# Find clusters in the Atlas that share these exact dominant SAE features
# (Matches the "Nearest Neighbor Fingerprint Voting" pattern)
matched_clusters = df_representatives[
    df_representatives['top_sae_features'].apply(lambda x: any(f in x for f in active_feature_ids))
]

print(f"Found {len(matched_clusters)} candidate cluster matches in the ESM Atlas!")
print(matched_clusters[['cluster_id', 'top_pfam_domains', 'taxonomy']].head())
