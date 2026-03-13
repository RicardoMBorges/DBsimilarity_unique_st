"""
============================================================
DBsimilarity-unique — Structure Recognition & Similarity App
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
import plotly.io as pio
import io

from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, Crippen, Lipinski, AllChem, rdMolTransforms
import py3Dmol
from typing import Optional
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import MolsToGridImage

import numpy as np
import plotly.express as px
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

# Mordred (2D/3D descriptors)
try:
    from mordred import Calculator, descriptors
    MORDRED_AVAILABLE = True
except Exception as e:
    import sys, traceback
    MORDRED_AVAILABLE = False
    print("==== Mordred import failed in DBsimilarity-unique ====")
    print("Python executable:", sys.executable)
    traceback.print_exc()




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

def read_csv_safe(uploaded_file, sep=";"):
    raw = uploaded_file.getvalue()  # bytes, stable
    for enc in ["utf-8", "utf-8-sig", "cp1252", "latin1", "ISO-8859-1"]:
        try:
            return pd.read_csv(io.BytesIO(raw), sep=sep, encoding=enc, engine="python")
        except UnicodeDecodeError:
            continue
    # last-resort: replace invalid chars
    return pd.read_csv(io.BytesIO(raw), sep=sep, encoding="latin1", engine="python", encoding_errors="replace")

def suggest_method(mol, similarity_threshold, db_file, fp_kind: str):
    try:
        if db_file is None:
            st.warning("⚠️ Please upload a database CSV file in the sidebar.")
            return

        reset_uploaded_file(db_file)
        database_df = read_csv_safe(db_file, sep=";")
        smiles_col = get_smiles_column(database_df)
        if smiles_col is None:
            st.error("❌ Your database must include a SMILES-like column (e.g., SMILES/canonical_smiles).")
            return

        st.markdown("### 📊 Suggested Method Parameters")

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

        fp_query = get_fingerprint(mol, fp_kind)
        if fp_query is None:
            st.error("Could not compute query fingerprint.")
            return

        best_sim = -1.0
        best_method = None
        best_comment = None

        for _, row in database_df.iterrows():
            s = row.get(smiles_col, None)
            if pd.isna(s) or not isinstance(s, str):
                continue
            db_mol = Chem.MolFromSmiles(s)
            if db_mol is None:
                continue

            fp_db = get_fingerprint(db_mol, fp_kind)
            if fp_db is None:
                continue

            sim = DataStructs.TanimotoSimilarity(fp_query, fp_db)
            if sim > best_sim:
                best_sim = sim
                best_method = row.get(method_col, "(no method provided)")
                best_comment = row.get(comment_col, None)

        if best_sim >= similarity_threshold:
            st.success(f"Suggested method (similarity {best_sim:.2f}): **{best_method}**")
            if best_comment:
                st.info(f"🗒️ Additional comments: {best_comment}")
        elif best_sim >= 0:
            st.info(f"Most similar compound has similarity {best_sim:.2f}, below threshold ({similarity_threshold}).")
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

def reset_uploaded_file(f):
    """Reset file pointer for a Streamlit UploadedFile so it can be read again."""
    try:
        f.seek(0)
    except Exception:
        pass
    return f


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

from rdkit import DataStructs

def _bemis_murcko_scaffold(mol: Chem.Mol, generic: bool = False) -> Chem.Mol:
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        if scaf is None or scaf.GetNumAtoms() == 0:
            return mol
        if generic:
            scaf = MurckoScaffold.MakeScaffoldGeneric(scaf)
        return scaf
    except Exception:
        return mol

def get_fingerprint(mol: Chem.Mol, kind: str):
    """Return an RDKit fingerprint object for the given molecule and type."""
    if mol is None:
        return None
    try:
        if kind == "Morgan (ECFP4)":
            return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        elif kind == "Bemis–Murcko (ECFP4 on scaffold)":
            scaf = _bemis_murcko_scaffold(mol)
            return AllChem.GetMorganFingerprintAsBitVect(scaf, radius=2, nBits=2048)
        elif kind == "RDKit":
            return Chem.RDKFingerprint(mol)
        elif kind == "MACCS":
            from rdkit.Chem import MACCSkeys
            return MACCSkeys.GenMACCSKeys(mol)
        elif kind == "AtomPair":
            from rdkit.Chem.AtomPairs import Pairs
            return Pairs.GetAtomPairFingerprintAsBitVect(mol)
        elif kind == "TopologicalTorsion":
            from rdkit.Chem.AtomPairs import Torsions
            return Torsions.GetTopologicalTorsionFingerprintAsIntVect(mol)
        else:
            # default to Morgan if unknown
            return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    except Exception:
        return None


def compute_similarity_table(
    query_mol: Chem.Mol,
    db_df: pd.DataFrame,
    smiles_col: str,
    use_generic_scaffold: bool = False,
    include_mcs: bool = True,
    top_k_images: int = 0  # number of top rows to render as PNG panel
) -> pd.DataFrame:
    from rdkit import DataStructs
    rows = []

    # Precompute query fingerprints
    fp_query_morgan = get_fingerprint(query_mol, "Morgan (ECFP4)")
    fp_query_murcko, q_scaf = scaffold_fp(query_mol, generic=use_generic_scaffold)

    for i, row in db_df.iterrows():
        s = row.get(smiles_col, None)
        if not isinstance(s, str) or pd.isna(s):
            continue
        db_mol = Chem.MolFromSmiles(s)
        if db_mol is None:
            continue

        # Morgan similarity on full molecules
        fp_db_morgan = get_fingerprint(db_mol, "Morgan (ECFP4)")
        sim_morgan = DataStructs.TanimotoSimilarity(fp_query_morgan, fp_db_morgan) if fp_db_morgan is not None else None

        # Scaffold similarity
        fp_db_murcko, d_scaf = scaffold_fp(db_mol, generic=use_generic_scaffold)
        sim_murcko = DataStructs.TanimotoSimilarity(fp_query_murcko, fp_db_murcko) if fp_db_murcko is not None else None

        mcs_smarts, mcs_na, mcs_nb = "", None, None
        if include_mcs and q_scaf is not None and d_scaf is not None:
            mcs_res, mcs_mol = mcs_between_scaffolds(q_scaf, d_scaf)
            if mcs_res is not None and mcs_mol is not None:
                mcs_smarts = mcs_res.smartsString
                mcs_na = mcs_res.numAtoms
                mcs_nb = mcs_res.numBonds

        rows.append({
            "row_idx": i,
            "SMILES": s,
            "Similarity_Morgan": round(sim_morgan, 4) if sim_morgan is not None else None,
            "Similarity_BemisMurcko": round(sim_murcko, 4) if sim_murcko is not None else None,
            "Query_Scaffold_SMILES": mol_to_smiles(q_scaf),
            "DB_Scaffold_SMILES": mol_to_smiles(d_scaf),
            "MCS_SMARTS": mcs_smarts,
            "MCS_NumAtoms": mcs_na,
            "MCS_NumBonds": mcs_nb,
        })

    out = pd.DataFrame(rows)

    # Optionally render PNGs for top K by scaffold similarity
    if top_k_images and not out.empty:
        out = out.sort_values("Similarity_BemisMurcko", ascending=False).reset_index(drop=True)
        png_bytes_list = []
        for j in range(min(top_k_images, len(out))):
            idx = int(out.loc[j, "row_idx"])
            db_smiles = out.loc[j, "SMILES"]
            db_mol = Chem.MolFromSmiles(db_smiles)
            png = None
            if isinstance(out.loc[j, "MCS_SMARTS"], str) and out.loc[j, "MCS_SMARTS"]:
                png = draw_mcs_panel(query_mol, db_mol, out.loc[j, "MCS_SMARTS"], "Query", f"DB idx {idx}")
            png_bytes_list.append(png)
        out["MCS_Image_PNG"] = png_bytes_list
    return out



if mol:
    st.subheader("2D Structure")
    st.image(Draw.MolToImage(mol, size=(250, 250)))

    # --- Descriptors ---
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

    # --- 3D structure + angles ---
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

def generate_similarity_network(mol, db_df, threshold, hover_col=None, fp_kind="Morgan (ECFP4)"):
    fp_query = get_fingerprint(mol, fp_kind)
    if fp_query is None:
        st.error("Could not compute query fingerprint.")
        return

    G = nx.Graph()
    G.add_node("Query")
    node_data = {}

    for i, row in db_df.iterrows():
        try:
            s = row.get("SMILES", None)
            if pd.isna(s) or not isinstance(s, str):
                continue
            db_mol = Chem.MolFromSmiles(s)
            if db_mol is None:
                continue

            fp_db = get_fingerprint(db_mol, fp_kind)
            if fp_db is None:
                continue

            sim = DataStructs.TanimotoSimilarity(fp_query, fp_db)
            if sim >= threshold:
                label = f"Mol_{i}"
                hover_info = f"Similarity: {sim:.2f}"
                if hover_col and hover_col in row and pd.notna(row[hover_col]):
                    hover_info += f"<br>{hover_col}: {row[hover_col]}"
                G.add_node(label)
                G.add_edge("Query", label, weight=sim)
                node_data[label] = hover_info
        except Exception:
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
    # Separate DB nodes and Query node
    db_x, db_y, db_text = [], [], []
    q_x, q_y, q_text = [], [], []

    for node in G.nodes():
        x, y = pos[node]
        if node == "Query":
            q_x.append(x)
            q_y.append(y)
            q_text.append("Query molecule")
        else:
            db_x.append(x)
            db_y.append(y)
            db_text.append(node_data.get(node, node))

    fig = go.Figure()

    # Edges
    fig.add_trace(go.Scatter(
        x=edge_x,
        y=edge_y,
        mode="lines",
        line=dict(width=1, color="gray"),
        hoverinfo="none"
    ))

    # DB nodes (skyblue)
    fig.add_trace(go.Scatter(
        x=db_x,
        y=db_y,
        mode="markers+text",
        name="Database molecules",
        text=db_text,
        hoverinfo="text",
        textposition="top center",
        marker=dict(size=18, color="skyblue")
    ))

    # Query node (highlighted red star)
    fig.add_trace(go.Scatter(
        x=q_x,
        y=q_y,
        mode="markers+text",
        name="Query",
        text=q_text,
        hoverinfo="text",
        textposition="bottom center",
        marker=dict(
            size=26,
            #symbol="star",
            color="red",
            line=dict(width=1, color="black")
        )
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


def mol_to_smiles(m: Chem.Mol) -> str:
    try:
        return Chem.MolToSmiles(m, isomericSmiles=True) if m is not None else ""
    except Exception:
        return ""


def get_scaffold(m: Chem.Mol, generic: bool = False) -> Chem.Mol:
    scaf = _bemis_murcko_scaffold(m, generic=generic)
    # sanitize to avoid depiction issues
    try:
        Chem.SanitizeMol(scaf)
    except Exception:
        pass
    return scaf


def scaffold_fp(m: Chem.Mol, generic: bool = False):
    scaf = get_scaffold(m, generic=generic)
    return AllChem.GetMorganFingerprintAsBitVect(scaf, radius=2, nBits=2048), scaf


def mcs_between_scaffolds(scaf1: Chem.Mol, scaf2: Chem.Mol):
    """Return MCS object + molecule from SMARTS; None-safe."""
    if scaf1 is None or scaf2 is None:
        return None, None
    try:
        res = rdFMCS.FindMCS(
            [scaf1, scaf2],
            ringMatchesRingOnly=True,
            completeRingsOnly=True,
            timeout=5
        )
        if res.canceled:
            return None, None
        mcs_mol = Chem.MolFromSmarts(res.smartsString)
        return res, mcs_mol
    except Exception:
        return None, None


def draw_mcs_panel(query_mol: Chem.Mol, db_mol: Chem.Mol, mcs_smarts: str, legend_q="Query", legend_d="DB"):
    """Return a PNG bytes image showing both molecules with MCS highlighted."""
    try:
        patt = Chem.MolFromSmarts(mcs_smarts)
        q_hits = query_mol.GetSubstructMatches(patt) if patt is not None else []
        d_hits = db_mol.GetSubstructMatches(patt) if patt is not None else []
        q_atoms = list(q_hits[0]) if q_hits else []
        d_atoms = list(d_hits[0]) if d_hits else []
        imgs = MolsToGridImage(
            [query_mol, db_mol],
            molsPerRow=2,
            highlightAtomLists=[q_atoms, d_atoms],
            legends=[legend_q, legend_d],
            subImgSize=(350, 300),
            useSVG=False
        )
        # RDKit returns PIL Image
        import io
        buf = io.BytesIO()
        imgs.save(buf, format="PNG")
        return buf.getvalue()
    except Exception:
        return None

def sanitize_matrix_for_modeling(df_vals: pd.DataFrame) -> pd.DataFrame:
    """
    Keep numeric, finite values only and drop constant-variance columns.
    This is a lightweight version of what you already use in DBsimilarity.
    """
    X = df_vals.select_dtypes(include=[np.number]).copy()
    if X.empty:
        return X
    # Replace inf with NaN, then fill with 0.0
    X = X.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    # Drop constant columns
    variances = X.var(axis=0)
    X = X.loc[:, variances > 0.0]
    return X

def get_similar_compounds_df(mol, db_df, smiles_col, threshold, fp_kind="Morgan (ECFP4)"):
    fp_query = get_fingerprint(mol, fp_kind)
    if fp_query is None:
        return pd.DataFrame()

    rows = []
    for i, row in db_df.iterrows():
        try:
            s = row.get(smiles_col, None)
            if pd.isna(s) or not isinstance(s, str):
                continue

            db_mol = Chem.MolFromSmiles(s)
            if db_mol is None:
                continue

            fp_db = get_fingerprint(db_mol, fp_kind)
            if fp_db is None:
                continue

            sim = DataStructs.TanimotoSimilarity(fp_query, fp_db)
            if sim >= threshold:
                row_dict = row.to_dict()
                row_dict["Similarity"] = round(float(sim), 4)
                row_dict["Row_Index"] = i
                rows.append(row_dict)
        except Exception:
            continue

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(rows).sort_values("Similarity", ascending=False).reset_index(drop=True)

# === Call the function ===
if mol and db_file is not None:
    try:
        #reset_uploaded_file(db_file)
        db_df = read_csv_safe(db_file, sep=";")
        smiles_col = get_smiles_column(db_df)

        if smiles_col:
            st.subheader("Similarity-based Method Suggestion & Network")

            # Controls close to the network
            col_thr, col_fp = st.columns(2)
            with col_thr:
                similarity_threshold = st.slider(
                    "Select similarity threshold",
                    min_value=0.0, max_value=1.0, value=0.9, step=0.01
                )
            with col_fp:
                fp_kind = st.selectbox(
                    "Fingerprint",
                    ["Morgan (ECFP4)", "Bemis–Murcko (ECFP4 on scaffold)", "RDKit", "MACCS", "AtomPair", "TopologicalTorsion"],
                    index=0
                )

            if st.button("🔍 Suggest method based on chemical similarity"):
                suggest_method(mol, similarity_threshold, db_file, fp_kind)

            # --- Similarity network ---
            with st.expander("Similarity Network", expanded=False):
                hover_col = st.selectbox(
                    "Choose column to display on nodes",
                    db_df.columns.tolist(),
                    index=0,
                    key="hover_col_similarity_network"
                )
                st.session_state["hover_col_selected"] = hover_col

                generate_similarity_network(
                    mol,
                    db_df.rename(columns={smiles_col: "SMILES"}),
                    similarity_threshold,
                    hover_col,
                    fp_kind=fp_kind
                )

                similar_df = get_similar_compounds_df(
                    mol=mol,
                    db_df=db_df,
                    smiles_col=smiles_col,
                    threshold=similarity_threshold,
                    fp_kind=fp_kind
                )

                st.markdown("### Wordcloud from similar compounds")

                if similar_df.empty:
                    st.info("No similar compounds available above the selected threshold.")
                else:
                    text_candidate_cols = [
                        c for c in similar_df.columns
                        if similar_df[c].dtype == "object" and c != smiles_col
                    ]

                    if text_candidate_cols:
                        wordcloud_col = st.selectbox(
                            "Choose column to build the wordcloud",
                            text_candidate_cols,
                            index=0,
                            key="wordcloud_column_similarity_network"
                        )

                        st.caption(f"{len(similar_df)} compounds matched the threshold.")

                        text_data = (
                            similar_df[wordcloud_col]
                            .dropna()
                            .astype(str)
                            .str.strip()
                        )
                        text_data = text_data[text_data != ""]

                        if text_data.empty:
                            st.warning("Selected column has no text values among the similar compounds.")
                        else:
                            try:
                                from wordcloud import WordCloud
                                import matplotlib.pyplot as plt

                                wc_text = " ".join(text_data.tolist())
                                wc = WordCloud(
                                    width=1200,
                                    height=500,
                                    background_color="white"
                                ).generate(wc_text)

                                fig_wc, ax_wc = plt.subplots(figsize=(12, 5))
                                ax_wc.imshow(wc, interpolation="bilinear")
                                ax_wc.axis("off")
                                st.pyplot(fig_wc)

                            except ModuleNotFoundError:
                                st.warning("Package 'wordcloud' is not installed. Showing term frequencies instead.")

                                term_counts = (
                                    text_data.str.split(r"[;,/\|\s]+")
                                    .explode()
                                    .dropna()
                                    .astype(str)
                                    .str.strip()
                                )
                                term_counts = term_counts[term_counts != ""]
                                freq_df = term_counts.value_counts().reset_index()
                                freq_df.columns = ["Term", "Count"]

                                if not freq_df.empty:
                                    fig_bar = px.bar(freq_df.head(20), x="Term", y="Count", title="Top term frequencies")
                                    st.plotly_chart(fig_bar, use_container_width=True)
                                else:
                                    st.info("No terms available to plot.")

                            except Exception as e:
                                st.warning(f"Could not render wordcloud: {e}")

                        with st.expander("Preview similar compounds used for wordcloud", expanded=False):
                            st.dataframe(similar_df.head(30), use_container_width=True)
                    else:
                        st.info("No text-like columns available for wordcloud generation.")

            # --- Export (Morgan & Bemis–Murcko) ---
            with st.expander("Export similarities (Morgan & Bemis–Murcko + common scaffold)", expanded=False):
                use_generic = st.checkbox(
                    "Use generic Bemis–Murcko scaffold",
                    value=False,
                    help="Normalize atom types in the scaffold to emphasize topology."
                )
                id_col = st.selectbox(
                    "Optional identifier column to include",
                    ["(none)"] + db_df.columns.tolist(), index=0
                )
                top_k_imgs = st.number_input(
                    "Render PNG panels for top-K by scaffold similarity (0 = none)",
                    min_value=0, max_value=50, value=0, step=1
                )

                if st.button("⚙️ Compute scaffold similarities + MCS"):
                    sim_df = compute_similarity_table(
                        mol, db_df, smiles_col,
                        use_generic_scaffold=use_generic,
                        include_mcs=True,
                        top_k_images=top_k_imgs
                    )

                    # Attach optional columns from DB
                    extra_cols = []
                    if id_col != "(none)":
                        extra_cols.append(id_col)
                    for maybe in ["Descriptive_Method", "Additional_comments"]:
                        if maybe in db_df.columns:
                            extra_cols.append(maybe)

                    if extra_cols:
                        sim_df = sim_df.merge(
                            db_df[extra_cols].reset_index().rename(columns={"index": "row_idx"}),
                            on="row_idx",
                            how="left"
                        )

                    st.dataframe(
                        sim_df.drop(columns=["MCS_Image_PNG"], errors="ignore").head(30),
                        use_container_width=True
                    )

                    # CSV download
                    csv_bytes = sim_df.drop(
                        columns=["row_idx", "MCS_Image_PNG"], errors="ignore"
                    ).to_csv(index=False).encode("utf-8")
                    st.download_button(
                        "📥 Download CSV (Morgan + Bemis–Murcko + MCS)",
                        data=csv_bytes,
                        file_name="similarities_scaffolds_with_MCS.csv",
                        mime="text/csv"
                    )

                    # Optional: zip the PNGs (top-K rows only)
                    if top_k_imgs and "MCS_Image_PNG" in sim_df.columns and sim_df["MCS_Image_PNG"].notna().any():
                        import io, zipfile
                        buf = io.BytesIO()
                        with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
                            for k, row in sim_df.head(top_k_imgs).iterrows():
                                png = row.get("MCS_Image_PNG", None)
                                if png:
                                    zf.writestr(
                                        f"mcs_panel_rank{k+1}_row{int(row['row_idx'])}.png",
                                        png
                                    )
                        st.download_button(
                            "🖼️ Download MCS PNGs (ZIP)",
                            data=buf.getvalue(),
                            file_name="mcs_panels.zip",
                            mime="application/zip"
                        )
        else:
            st.error("❌ Your database must contain a SMILES-like column (SMILES/smiles/canonical_smiles...).")

    except Exception as e:
        st.error(f"❌ Error reading database: {e}")

#################
# ======================================================================
# 3D / 2D Mordred descriptors + PCA / t-SNE on the uploaded database
# ======================================================================
with st.expander("3D descriptors + PCA / t-SNE on database", expanded=False):
    st.markdown(
        """
        This section computes **Mordred descriptors** (2D and optional 3D) for the
        molecules in your uploaded CSV database (`db_file`), then runs:
        
        - **PCA** (2D scores plot)  
        - **t-SNE** (2D embedding, with optional cluster coloring later if you want)
        
        ⚠️ This can be computationally heavy — use the *Max molecules* limit below
        to keep things responsive.
        """
    )

    if not MORDRED_AVAILABLE:
        st.warning(
            "Mordred is not available in this environment. "
            "Install `mordred` and redeploy to enable this section."
        )
    else:
        if db_file is None:
            st.info("Upload a CSV database file in the sidebar to enable this section.")
        else:
            try:
                reset_uploaded_file(db_file)
                db_desc_df = read_csv_safe(db_file, sep=";")
            except Exception as e:
                st.error(f"Could not read database file for descriptors: {e}")
                db_desc_df = pd.DataFrame()

            if db_desc_df.empty:
                st.info("Database appears empty or could not be parsed.")
            else:
                smiles_col = get_smiles_column(db_desc_df)
                if smiles_col is None:
                    st.error(
                        "Your database must contain a SMILES-like column "
                        "(e.g., SMILES/canonical_smiles/isomeric_smiles)."
                    )
                else:
                    # Get the column selected in the Similarity Network (if any)
                    hover_col_sel = st.session_state.get("hover_col_selected", None)
                    if hover_col_sel and hover_col_sel not in db_desc_df.columns:
                        hover_col_sel = None

                    # --- Controls for this analysis ---
                    c1, c2, c3 = st.columns(3)

                    c1, c2, c3 = st.columns(3)
                    with c1:
                        max_mols = st.number_input(
                            "Max molecules to process",
                            min_value=10, max_value=2000, value=300, step=10,
                            help="For performance. Molecules are truncated to the first N rows."
                        )
                    with c2:
                        use_3d_desc = st.checkbox(
                            "Use 3D descriptors (generate 3D conformers)",
                            value=True,
                            help="If checked, RDKit will generate 3D conformers before computing Mordred descriptors."
                        )
                    with c3:
                        tsne_perp = st.slider(
                            "t-SNE perplexity",
                            min_value=5, max_value=50, value=20, step=1,
                            help="Roughly, the effective number of neighbors used by t-SNE."
                        )

                    tsne_seed = st.number_input(
                        "t-SNE random_state (for reproducibility)",
                        min_value=0, max_value=9999, value=0, step=1
                    )

                    # Prepare SMILES list
                    work_df = (
                        db_desc_df
                        .dropna(subset=[smiles_col])
                        .drop_duplicates(subset=[smiles_col])
                        .head(int(max_mols))
                        .copy()
                    )
                    smiles_list = work_df[smiles_col].astype(str).tolist()

                    st.caption(f"Using up to {len(smiles_list)} database molecules for descriptor calculation.")

                    # --- Build molecules (with optional 3D), including the query ---
                    mols = []
                    meta = []  # track which row is Query vs DB
                    failed = 0

                    # 0) Add query molecule first (if available)
                    if mol is not None:
                        # Define hover text for the query
                        if hover_col_sel:
                            query_hover = f"Query (no {hover_col_sel} in DB)"
                        else:
                            query_hover = "Query molecule"

                        try:
                            q_m = mol
                            if use_3d_desc:
                                q3d = Chem.AddHs(q_m)
                                status = AllChem.EmbedMolecule(q3d, randomSeed=0xF00D)
                                if status == 0:
                                    AllChem.UFFOptimizeMolecule(q3d)
                                    q_m = q3d
                            mols.append(q_m)
                            meta.append({"role": "Query", "label": "Query", "hover": query_hover})
                        except Exception:
                            # if 3D fails, still keep 2D version
                            mols.append(mol)
                            meta.append({"role": "Query", "label": "Query", "hover": query_hover})


                    # 1) Add DB molecules
                    # 1) Add DB molecules
                    for idx, s in zip(work_df.index, smiles_list):
                        m = Chem.MolFromSmiles(s)
                        if m is None:
                            failed += 1
                            continue
                        try:
                            if use_3d_desc:
                                m3d = Chem.AddHs(m)
                                status = AllChem.EmbedMolecule(m3d, randomSeed=0xF00D)
                                if status == 0:
                                    AllChem.UFFOptimizeMolecule(m3d)
                                    m = m3d
                            mols.append(m)

                            # Choose a meaningful label if possible
                            if "Name" in work_df.columns:
                                label = str(work_df.loc[idx, "Name"])
                            elif "ID" in work_df.columns:
                                label = str(work_df.loc[idx, "ID"])
                            else:
                                label = f"DB_{idx}"

                            # Hover text uses the chosen column if available
                            if hover_col_sel:
                                hover_val = str(work_df.loc[idx, hover_col_sel])
                            else:
                                hover_val = label

                            meta.append({"role": "DB", "label": label, "hover": hover_val})
                        except Exception:
                            failed += 1
                            continue


                    if not mols:
                        st.error("No valid molecules could be generated from SMILES.")
                    else:
                        if failed:
                            st.warning(f"{failed} SMILES could not be parsed and were skipped.")

                        # --- Mordred descriptors ---
                        with st.spinner("Calculating Mordred descriptors (may take a while)…"):
                            try:
                                calc = Calculator(descriptors, ignore_3D=not use_3d_desc)
                                df_desc = calc.pandas(mols)
                            except Exception as e:
                                st.error(f"Mordred failed: {e}")
                                df_desc = pd.DataFrame()
                        meta_df = pd.DataFrame(meta)

                        if df_desc.empty:
                            st.info("Descriptor table is empty. Skipping PCA / t-SNE.")
                        else:
                            # Clean numeric matrix
                            df_vals_raw = df_desc.select_dtypes(include=[np.number])
                            X_df = sanitize_matrix_for_modeling(df_vals_raw)

                            if X_df.empty or X_df.shape[0] < 2:
                                st.info("Not enough data for PCA/t-SNE after cleaning.")
                            else:
                                # --- Align meta_df and X_df row counts ---
                                n_meta = len(meta_df)
                                n_X = len(X_df)

                                if n_meta > n_X:
                                    # More metadata entries than descriptor rows -> trim meta
                                    meta_df = meta_df.iloc[:n_X].reset_index(drop=True)
                                elif n_meta < n_X:
                                    # More descriptor rows than metadata entries -> trim X_df
                                    X_df = X_df.iloc[:n_meta, :].reset_index(drop=True)

                                # Now they MUST be equal
                                n_samples, n_features = X_df.shape
                                st.caption(
                                    f"Descriptors after cleaning: {n_samples} samples × {n_features} features."
                                )

                                # Optional: cap number of features by variance
                                max_features = 1000
                                if n_features > max_features:
                                    topk = X_df.var().sort_values(ascending=False).index[:max_features]
                                    X_df = X_df[topk]
                                    n_samples, n_features = X_df.shape

                                # Standardization (zero mean, unit variance)
                                try:
                                    X = StandardScaler(with_mean=True, with_std=True).fit_transform(X_df.values)
                                except Exception as e:
                                    st.warning(f"Standardization failed ({e}); using raw values.")
                                    X = X_df.values


                                # Standardization (zero mean, unit variance)
                                try:
                                    X = StandardScaler(with_mean=True, with_std=True).fit_transform(X_df.values)
                                except Exception as e:
                                    st.warning(f"Standardization failed ({e}); using raw values.")
                                    X = X_df.values

                                # IDs for hover (try InChIKey if present, else row index)
                                if "InChIKey" in db_desc_df.columns:
                                    ids = (
                                        db_desc_df.dropna(subset=[smiles_col])
                                        .drop_duplicates(subset=[smiles_col])
                                        .head(len(X_df))
                                        ["InChIKey"]
                                        .astype(str)
                                        .tolist()
                                    )
                                else:
                                    ids = [f"mol_{i+1}" for i in range(len(X_df))]

                                # --- PCA ---
                                st.markdown("### PCA (2D)")
                                try:
                                    ncomp = int(min(2, n_features, n_samples))
                                    if ncomp >= 1:
                                        pca = PCA(n_components=ncomp, svd_solver="auto")
                                        scores = pca.fit_transform(X)

                                        pc1 = scores[:, 0]
                                        pc2 = scores[:, 1] if ncomp > 1 else np.zeros_like(pc1)

                                        coords_pca = pd.DataFrame({
                                            "PC1": pc1,
                                            "PC2": pc2,
                                            "role": meta_df["role"],
                                            "label": meta_df["label"],
                                            "hover": meta_df["hover"],
                                        })


                                        fig_pca = go.Figure()

                                        # DB molecules
                                        db_mask = coords_pca["role"] == "DB"
                                        fig_pca.add_trace(go.Scatter(
                                            x=coords_pca.loc[db_mask, "PC1"],
                                            y=coords_pca.loc[db_mask, "PC2"],
                                            mode="markers",
                                            name="Database",
                                            hovertext=coords_pca.loc[db_mask, "hover"],
                                            hoverinfo="text",
                                            marker=dict(size=8, opacity=0.9)
                                        ))


                                        # Query molecule (highlighted)
                                        q_mask = coords_pca["role"] == "Query"
                                        if q_mask.any():
                                            fig_pca.add_trace(go.Scatter(
                                                x=coords_pca.loc[q_mask, "PC1"],
                                                y=coords_pca.loc[q_mask, "PC2"],
                                                mode="markers+text",
                                                name="Query",
                                                text=["Query"],
                                                textposition="top center",
                                                hovertext=coords_pca.loc[q_mask, "label"],
                                                hoverinfo="text",
                                                marker=dict(size=14, symbol="star", line=dict(width=1),color="red")
                                            ))

                                        fig_pca.update_layout(
                                            title="PCA on Mordred descriptors",
                                            xaxis_title="PC1",
                                            yaxis_title=("PC2" if ncomp > 1 else "PC2 (zero)"),
                                            legend_title="Molecule type",
                                            margin=dict(l=20, r=20, t=40, b=20),
                                        )

                                        st.plotly_chart(fig_pca, use_container_width=True)

                                        # Downloads
                                        pca_df = coords_pca[["PC1", "PC2", "role", "label"]]
                                        st.download_button(
                                            "Download PCA scores (CSV)",
                                            data=pca_df.to_csv(index=False).encode("utf-8"),
                                            file_name="pca_scores_mordred.csv",
                                            mime="text/csv"
                                        )
                                        st.download_button(
                                            "📥 Download PCA plot (HTML)",
                                            data=pio.to_html(fig_pca, include_plotlyjs='cdn', full_html=True).encode("utf-8"),
                                            file_name="pca_mordred.html",
                                            mime="text/html"
                                        )
                                    else:
                                        st.info("PCA skipped: fewer than 1 usable component.")
                                except Exception as e:
                                    st.warning(f"PCA failed: {e}")


                                # --- t-SNE ---
                                st.markdown("### t-SNE (2D)")
                                try:
                                    if n_samples < 3:
                                        st.info("t-SNE skipped: need at least 3 samples.")
                                    else:
                                        max_perp = max(
                                            5,
                                            min(int(tsne_perp), n_samples - 1, max(5, (n_samples - 1) // 3))
                                        )
                                        tsne = TSNE(
                                            n_components=2,
                                            learning_rate="auto",
                                            random_state=int(tsne_seed),
                                            init="random",
                                            perplexity=max_perp,
                                        )
                                        ts = tsne.fit_transform(X)

                                        coords_ts = pd.DataFrame({
                                            "tSNE1": ts[:, 0],
                                            "tSNE2": ts[:, 1],
                                            "role": meta_df["role"],
                                            "label": meta_df["label"],
                                            "hover": meta_df["hover"],
                                        })


                                        fig_ts = go.Figure()

                                        # DB molecules
                                        db_mask = coords_ts["role"] == "DB"
                                        if q_mask.any():
                                            fig_ts.add_trace(go.Scatter(
                                                x=coords_ts.loc[q_mask, "tSNE1"],
                                                y=coords_ts.loc[q_mask, "tSNE2"],
                                                mode="markers+text",
                                                name="Query",
                                                text=["Query"],
                                                textposition="top center",
                                                hovertext=coords_ts.loc[q_mask, "hover"],
                                                hoverinfo="text",
                                                marker=dict(size=14, symbol="star", line=dict(width=1), color="red")
                                            ))



                                        # Query molecule (highlighted)
                                        q_mask = coords_ts["role"] == "Query"
                                        if q_mask.any():
                                            fig_ts.add_trace(go.Scatter(
                                                x=coords_ts.loc[db_mask, "tSNE1"],
                                                y=coords_ts.loc[db_mask, "tSNE2"],
                                                mode="markers",
                                                name="Database",
                                                hovertext=coords_ts.loc[db_mask, "hover"],
                                                hoverinfo="text",
                                                marker=dict(size=8, opacity=0.9)
                                            ))


                                        fig_ts.update_layout(
                                            title=f"t-SNE on Mordred descriptors (perplexity={max_perp})",
                                            xaxis_title="tSNE-1",
                                            yaxis_title="tSNE-2",
                                            legend_title="Molecule type",
                                            margin=dict(l=20, r=20, t=40, b=20),
                                        )

                                        st.plotly_chart(fig_ts, use_container_width=True)

                                        tsne_df = coords_ts[["tSNE1", "tSNE2", "role", "label"]]
                                        st.download_button(
                                            "Download t-SNE scores (CSV)",
                                            data=tsne_df.to_csv(index=False).encode("utf-8"),
                                            file_name="tsne_scores_mordred.csv",
                                            mime="text/csv"
                                        )
                                        st.download_button(
                                            "📥 Download t-SNE plot (HTML)",
                                            data=pio.to_html(fig_ts, include_plotlyjs='cdn', full_html=True).encode("utf-8"),
                                            file_name="tsne_mordred.html",
                                            mime="text/html"
                                        )
                                except Exception as e:
                                    st.warning(f"t-SNE failed: {e}")


                                # --- Optionally export full descriptor table ---
                                with st.expander("Raw Mordred descriptor table (download)", expanded=False):
                                    st.dataframe(df_desc.head(10), use_container_width=True)
                                    st.download_button(
                                        "📥 Download Mordred descriptors (CSV)",
                                        data=df_desc.to_csv(index=False).encode("utf-8"),
                                        file_name="mordred_descriptors_full.csv",
                                        mime="text/csv"
                                    )


