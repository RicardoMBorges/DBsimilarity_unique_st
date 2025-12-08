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
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors, Crippen, Lipinski, AllChem, rdMolTransforms
import py3Dmol
from typing import Optional
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import MolsToGridImage



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

def suggest_method(mol, similarity_threshold, db_file, fp_kind: str):
    try:
        if db_file is None:
            st.warning("⚠️ Please upload a database CSV file in the sidebar.")
            return

        database_df = pd.read_csv(db_file, sep=";")
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

from rdkit import DataStructs

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


# === Call the function ===
if mol and db_file is not None:
    try:
        db_df = pd.read_csv(db_file, sep=";", on_bad_lines="skip")
        smiles_col = get_smiles_column(db_df)

        if smiles_col:
            st.subheader("Similarity-based Method Suggestion & Network")

            # Controls close to the network
            col_thr, col_fp = st.columns(2)
            with col_thr:
                similarity_threshold = st.slider(
                    "Select similarity threshold",
                    min_value=0.0, max_value=1.0, value=0.7, step=0.01
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
                    index=0
                )
                generate_similarity_network(
                    mol,
                    db_df.rename(columns={smiles_col: "SMILES"}),  # network expects "SMILES"
                    similarity_threshold,
                    hover_col,
                    fp_kind=fp_kind
                )

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
