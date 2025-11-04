"""
============================================================
Photo2Struct — Structure Recognition & Similarity App
 Author: Ricardo M. Borges | GitHub: RicardoMBorges
 
 This Streamlit app allows users to:
 1. Upload a chemical structure image or SMILES string.
 2. Extract molecular structure and descriptors.
 3. Upload a custom database (CSV with SMILES column).
 4. Suggest experimental methods based on molecular similarity.
 5. Visualize a similarity network (nodes = molecules, edges = similarity).

Core features:
 - Molecular structure rendering
 - RDKit-based descriptor calculation
 - Tanimoto similarity search with threshold filtering
 - Custom column selection for method suggestion and hover labels
 - Interactive network plot using Plotly and NetworkX

Requirements:
 - streamlit
 - rdkit-pypi
 - pandas, numpy, plotly, networkx
 - pillow

 Run with: streamlit run app.py
 ============================================================
"""

import streamlit as st
import subprocess
import tempfile
import os
import pandas as pd
from pathlib import Path
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, Crippen, Lipinski, AllChem, rdMolTransforms
import py3Dmol
from typing import Optional

st.set_page_config(layout="wide")
st.title("DBsimilarity-unique – Chemical Structure Recognition")

# Load logo
STATIC_DIR = Path(__file__).parent / "static"
LOGO_PATH = STATIC_DIR / "LAABio.png"

try:
    logo = Image.open(LOGO_PATH)
    st.sidebar.image(logo, use_container_width=True)
except FileNotFoundError:
    st.sidebar.warning("Logo not found at static/LAABio.png")

# Load logo
STATIC_DIR = Path(__file__).parent / "static"
LOGO_DBsim_PATH = STATIC_DIR / "DBsimilarity_unique.png"

try:
    logo = Image.open(LOGO_DBsim_PATH)
    st.sidebar.image(logo, use_container_width=True)
except FileNotFoundError:
    st.sidebar.warning("Logo not found at static/LAABio.png")




st.sidebar.markdown("""---""")

# Input from sidebar
smiles_input = st.sidebar.text_input("📋 Paste a SMILES string (optional):", placeholder="e.g., C1=CC=CC=C1")
uploaded_file = st.sidebar.file_uploader("📷 Or upload an image of the structure", type=["png", "jpg", "jpeg"])

st.sidebar.markdown("""---""")
db_file = st.sidebar.file_uploader("📁 Upload CSV database file (must contain 'SMILES' column)", type=["csv"])

mol = None
source = ""

def osra_image_to_smiles(image_bytes):
    with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as tmp:
        tmp.write(image_bytes)
        tmp.close()
        wsl_path = "/mnt/" + tmp.name[0].lower() + tmp.name[2:].replace('\\', '/').replace('\\\\', '/')
    try:
        result = subprocess.run(["wsl", "osra", wsl_path], capture_output=True, text=True, timeout=15)
        smiles = result.stdout.strip()
        if smiles == "":
            return None, "⚠️ No SMILES recognized."
        return smiles, None
    except Exception as e:
        return None, f"Error running OSRA: {e}"

def mol_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES."
    Chem.rdDepictor.Compute2DCoords(mol)
    return mol, None

# --- Flexible SMILES column finder ---
def get_smiles_column(df: pd.DataFrame) -> str | None:
    # 1) exact (case-insensitive)
    for c in df.columns:
        if c.strip().lower() == "smiles":
            return c

    # 2) common aliases (case-insensitive)
    aliases = [
        "canonical_smiles", "isomeric_smiles",
        "smiles_string", "smile", "structure_smiles"
    ]
    for alias in aliases:
        for c in df.columns:
            if c.strip().lower() == alias:
                return c

    # 3) any column that contains "smiles"
    for c in df.columns:
        if "smiles" in c.strip().lower():
            return c

    return None

# --- Flexible SMILES column finder ---
def get_smiles_column(df: pd.DataFrame) -> Optional[str]:
    # 1) exact (case-insensitive)
    for c in df.columns:
        if c.strip().lower() == "smiles":
            return c

    # 2) common aliases (case-insensitive)
    aliases = [
        "canonical_smiles", "isomeric_smiles",
        "smiles_string", "smile", "structure_smiles"
    ]
    for alias in aliases:
        for c in df.columns:
            if c.strip().lower() == alias:
                return c

    # 3) any column that contains "smiles"
    for c in df.columns:
        if "smiles" in c.strip().lower():
            return c

    return None

def suggest_method(mol, similarity_threshold, db_file):
    try:
        if db_file is None:
            st.warning("⚠️ Please upload a database CSV file in the sidebar.")
            return

        database_df = pd.read_csv(db_file, sep=";")
        smiles_col = get_smiles_column(database_df)
        if smiles_col is None:
            st.error("❌ Your database must include a SMILES column (e.g., SMILES, smiles, canonical_smiles).")
            st.stop()

        st.markdown("### 📊 Suggested Method Parameters")

        # Let the user choose which column to use for method and comments
        available_columns = list(database_df.columns)
        method_col = st.selectbox(
            "Column for Descriptive Method",
            available_columns,
            index=available_columns.index("Descriptive_Method") if "Descriptive_Method" in available_columns else 0
        )
        comment_col = st.selectbox(
            "Column for Additional Comments",
            available_columns,
            index=available_columns.index("Additional_comments") if "Additional_comments" in available_columns else 0
        )

        from rdkit import DataStructs
        fp_query = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        best_sim = 0.0
        best_method = None
        best_comment = None

        for _, row in database_df.iterrows():
            try:
                s = row.get(smiles_col, None)
                if pd.isna(s) or not isinstance(s, str):
                    continue
                db_mol = Chem.MolFromSmiles(s)
                if db_mol is None:
                    continue
                fp_db = AllChem.GetMorganFingerprintAsBitVect(db_mol, 2, 2048)
                sim = DataStructs.TanimotoSimilarity(fp_query, fp_db)
                if sim > best_sim:
                    best_sim = sim
                    best_method = row.get(method_col, "(no method provided)")
                    best_comment = row.get(comment_col, None)
            except Exception:
                continue

        if best_sim >= similarity_threshold:
            st.success(f"Suggested method (similarity {best_sim:.2f}): **{best_method}**")
            if best_comment:
                st.info(f"🗒️ Additional comments: {best_comment}")
        elif best_sim > 0:
            st.info(f"Most similar compound has similarity {best_sim:.2f}, but below threshold ({similarity_threshold}).")
        else:
            st.warning("No similar compound found in the database.")
    except Exception as e:
        st.error(f"❌ Could not perform similarity search: {e}")


if smiles_input:
    mol = Chem.MolFromSmiles(smiles_input)
    if mol:
        st.success("✅ SMILES parsed successfully!")
        st.image(Draw.MolToImage(mol, size=(250, 250)), caption="Structure from SMILES")
        source = "SMILES"
    else:
        st.error("❌ Invalid SMILES string.")

elif uploaded_file is not None:
    st.info("🔍 Processing uploaded image...")
    image_bytes = uploaded_file.read()
    smiles, error = osra_image_to_smiles(image_bytes)

    if error:
        st.error(error)
    else:
        st.success("✅ SMILES recognized from image:")
        st.code(smiles, language="none")
        mol, mol_error = mol_from_smiles(smiles)
        if mol_error:
            st.error(mol_error)
        else:
            #st.image(Draw.MolToImage(mol, size=(250, 250)), caption="Structure from Image")
            source = "Image"

def calculate_descriptors(mol):
    try:
        inchi = Chem.MolToInchi(mol)
    except:
        inchi = "Unavailable"

    try:
        inchikey = Chem.MolToInchiKey(mol)
    except:
        inchikey = "Unavailable"

    data = {
        'Molecular Formula': rdMolDescriptors.CalcMolFormula(mol),
        'Exact Mass': round(Descriptors.ExactMolWt(mol), 5),
        'Molecular Weight (g/mol)': round(Descriptors.MolWt(mol), 2),
        'LogP (Crippen)': round(Crippen.MolLogP(mol), 2),
        'TPSA (Polar Surface Area)': round(rdMolDescriptors.CalcTPSA(mol), 2),
        'HBD (H-Donors)': Lipinski.NumHDonors(mol),
        'HBA (H-Acceptors)': Lipinski.NumHAcceptors(mol),
        'Rotatable Bonds': Lipinski.NumRotatableBonds(mol),
        'Heavy Atom Count': rdMolDescriptors.CalcNumHeavyAtoms(mol),
        'Fraction Csp3': round(rdMolDescriptors.CalcFractionCSP3(mol), 2),
        'Aliphatic Ring Count': rdMolDescriptors.CalcNumAliphaticRings(mol),
        'Aromatic Ring Count': rdMolDescriptors.CalcNumAromaticRings(mol),
        'InChI': inchi,
        'InChIKey': inchikey
    }
    return pd.DataFrame(data.items(), columns=["Descriptor", "Value"])

def generate_3d_model(mol):
    mol3d = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol3d) != 0:
        return None
    AllChem.UFFOptimizeMolecule(mol3d)
    return mol3d

def render_3d(mol3d):
    mb = Chem.MolToMolBlock(mol3d)
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(mb, "mol")
    viewer.setStyle({"stick": {}})
    viewer.zoomTo()
    return viewer

def suggest_method(mol, similarity_threshold, db_file):
    try:
        if db_file is None:
            st.warning("⚠️ Please upload a database CSV file in the sidebar.")
            return

        database_df = pd.read_csv(db_file, sep=";")  # as you already do
        smiles_col = get_smiles_column(database_df)
        if smiles_col is None:
            st.error("❌ Your database must include a SMILES column (e.g., SMILES, smiles, canonical_smiles).")
            st.stop()

        # use smiles_col everywhere instead of hardcoding "SMILES"
        db_mol = Chem.MolFromSmiles(row[smiles_col])

        st.markdown("### 📊 Suggested Method Parameters")

        # Let the user choose which column to use for method and comments
        available_columns = list(database_df.columns)
        method_col = st.selectbox("Column for Descriptive Method", available_columns, index=available_columns.index("Descriptive_Method") if "Descriptive_Method" in available_columns else 0)
        comment_col = st.selectbox("Column for Additional Comments", available_columns, index=available_columns.index("Additional_comments") if "Additional_comments" in available_columns else 0)

        from rdkit import DataStructs
        fp_query = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        best_sim = 0
        best_method = None
        best_comment = None

        for _, row in database_df.iterrows():
            try:
                db_mol = Chem.MolFromSmiles(row["SMILES"])
                if db_mol:
                    fp_db = AllChem.GetMorganFingerprintAsBitVect(db_mol, 2, 2048)
                    sim = DataStructs.TanimotoSimilarity(fp_query, fp_db)
                    if sim > best_sim:
                        best_sim = sim
                        best_method = row.get(method_col, "(no method provided)")
                        best_comment = row.get(comment_col, None)
            except:
                continue

        if best_sim >= similarity_threshold:
            st.success(f"Suggested method (similarity {best_sim:.2f}): **{best_method}**")
            if best_comment:
                st.info(f"🗒️ Additional comments: {best_comment}")
        elif best_sim > 0:
            st.info(f"Most similar compound has similarity {best_sim:.2f}, but below threshold ({similarity_threshold}).")
        else:
            st.warning("No similar compound found in the database.")
    except Exception as e:
        st.error(f"❌ Could not perform similarity search: {e}")


if mol:
    st.subheader("2D Structure")
    st.image(Draw.MolToImage(mol, size=(250, 250)))

    similarity_threshold = st.slider("Select similarity threshold", min_value=0.0, max_value=1.0, value=0.7, step=0.01)

    st.subheader("Molecular Descriptors")
    desc_df = calculate_descriptors(mol)
    def render_html_table(df):
        html = df.to_html(index=False, escape=False, border=0)
        styled_html = f"""
        <div style="overflow-x:auto;">
            <style>th, td {{ padding: 8px 12px; text-align: left; }}</style>
            {html}
        </div>
        """
        st.markdown(styled_html, unsafe_allow_html=True)

    render_html_table(desc_df.astype(str))

    if st.button("🔍 Suggest method based on chemical similarity"):
        suggest_method(mol, similarity_threshold, db_file)

    st.subheader("Interactive 3D Structure")
    mol3d = generate_3d_model(mol)
    if mol3d:
        viewer = render_3d(mol3d)
        st.components.v1.html(viewer._make_html(), height=400)

        with st.expander("Measure Bond and Dihedral Angles"):
            st.subheader("Bond Angle (3 atoms)")
            col1, col2, col3 = st.columns(3)
            atom_v1 = col1.number_input("Atom 1", 0, mol.GetNumAtoms() - 1)
            atom_v2 = col2.number_input("Atom 2 (central)", 0, mol.GetNumAtoms() - 1)
            atom_v3 = col3.number_input("Atom 3", 0, mol.GetNumAtoms() - 1)

            if st.button("Calculate Bond Angle"):
                conf = mol3d.GetConformer()
                angle_deg = rdMolTransforms.GetAngleDeg(conf, int(atom_v1), int(atom_v2), int(atom_v3))
                st.success(f"Angle: {angle_deg:.2f}° between atoms {atom_v1}-{atom_v2}-{atom_v3}")

            st.subheader("Dihedral Angle (4 atoms)")
            col4, col5, col6, col7 = st.columns(4)
            atom_d1 = col4.number_input("Atom 1", 0, mol.GetNumAtoms() - 1, key="d1")
            atom_d2 = col5.number_input("Atom 2", 0, mol.GetNumAtoms() - 1, key="d2")
            atom_d3 = col6.number_input("Atom 3", 0, mol.GetNumAtoms() - 1, key="d3")
            atom_d4 = col7.number_input("Atom 4", 0, mol.GetNumAtoms() - 1, key="d4")

            if st.button("Calculate Dihedral Angle"):
                conf = mol3d.GetConformer()
                angle_dihedral = rdMolTransforms.GetDihedralDeg(conf, int(atom_d1), int(atom_d2), int(atom_d3), int(atom_d4))
                st.success(f"Dihedral: {angle_dihedral:.2f}° between atoms {atom_d1}-{atom_d2}-{atom_d3}-{atom_d4}")
    else:
        st.warning("⚠️ Failed to generate 3D structure.")

import networkx as nx
import plotly.graph_objects as go

def generate_similarity_network(mol, db_df, threshold, hover_col=None):
    from rdkit import DataStructs

    fp_query = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
    G = nx.Graph()
    G.add_node("Query")

    node_data = {}

    for i, row in db_df.iterrows():
        try:
            db_mol = Chem.MolFromSmiles(row["SMILES"])
            if db_mol:
                fp_db = AllChem.GetMorganFingerprintAsBitVect(db_mol, 2, 2048)
                sim = DataStructs.TanimotoSimilarity(fp_query, fp_db)
                if sim >= threshold:
                    label = f"Mol_{i}"
                    hover_info = f"Similarity: {sim:.2f}"

                    if hover_col and hover_col in row and pd.notna(row[hover_col]):
                        hover_info += f"<br>{hover_col}: {row[hover_col]}"
                    
                    G.add_node(label)
                    G.add_edge("Query", label, weight=sim)
                    node_data[label] = hover_info
        except:
            continue

    if G.number_of_edges() == 0:
        st.warning("⚠️ No connections found above the selected similarity threshold.")
        return

    pos = nx.spring_layout(G, seed=42)
    edge_x, edge_y = [], []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    node_x, node_y, node_text = [], [], []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        if node == "Query":
            node_text.append("Query molecule")
        else:
            node_text.append(node_data.get(node, node))

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines", line=dict(width=1, color="gray"), hoverinfo="none"))
    fig.add_trace(go.Scatter(
        x=node_x, y=node_y,
        mode="markers+text",
        text=node_text,
        hoverinfo="text",
        textposition="top center",
        marker=dict(size=20, color="skyblue")
    ))

    fig.update_layout(
        title="Similarity Network (Filtered by Threshold)",
        margin=dict(l=20, r=20, t=40, b=20),
        showlegend=False
    )

    st.plotly_chart(fig, use_container_width=True)

    # Add download button
    with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmpfile:
        fig.write_html(tmpfile.name)
        tmpfile.flush()
        with open(tmpfile.name, "rb") as f:
            st.download_button(
                label="📥 Download similarity network as HTML",
                data=f.read(),
                file_name="similarity_network.html",
                mime="text/html"
            )



# === Call the function at the end ===
# === Call the function at the end ===
if mol and db_file is not None:
    try:
        db_df = pd.read_csv(db_file, sep=";", on_bad_lines="skip")
        smiles_col = get_smiles_column(db_df)
        if smiles_col:
            with st.expander("Similarity Network", expanded=False):
                hover_col = st.selectbox("Choose column to display on nodes", db_df.columns.tolist(), index=0)
                # rename so the network code can keep using "SMILES"
                generate_similarity_network(
                    mol,
                    db_df.rename(columns={smiles_col: "SMILES"}),
                    similarity_threshold,
                    hover_col
                )
        else:
            st.error("❌ Your database must contain a SMILES-like column (SMILES/smiles/canonical_smiles...).")
    except Exception as e:
        st.error(f"❌ Error reading database: {e}")
