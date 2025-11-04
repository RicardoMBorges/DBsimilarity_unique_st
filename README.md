Here’s a complete **`README.md` tutorial** for your project, ready to be used on **GitHub** or **Streamlit Cloud** (like Streamlit Community Cloud or Hugging Face Spaces):

---

# photo2Struct: Chemical Structure Recognition & Similarity Analysis App

This Streamlit web application allows you to **extract chemical structure information from images or SMILES**, compute descriptors, compare to a user-uploaded database, and **visualize structural similarity as a network**.

## Features

* Upload a **chemical structure image** or enter a **SMILES** string
* Calculate molecular descriptors
* Perform similarity search against a custom database
* Suggest experimental methods based on most similar known structures
* Adjust similarity threshold with a slider
* Visualize **similarity networks** with node hover metadata from any database column

---

## 🧰 Requirements

Install dependencies (preferably in a virtual environment):

```bash
pip install streamlit rdkit-pypi pandas numpy plotly networkx pillow
```

---

## ▶️ How to Run

From the project folder, run the app using:

```bash
streamlit run app.py
```

The app will open in your default browser.

---

## 📂 Input Files

### 1. Structure Input

* Upload a **chemical structure image** (`.png`, `.jpg`) to convert it into SMILES
  *(Requires internet connection — uses `img2mol` behind the scenes)*
* OR paste a **SMILES string** manually

### 2. Database Input

Upload a `.csv` file containing at least:

* `SMILES` (required) – structure for similarity comparison
* Optional: `Descriptive_Method`, `Additional_comments`, or any custom metadata

ℹ️ Other columns will be used for network hover tooltips or method suggestion (you can select them).

---

## Interface Overview

### Sidebar

* Upload **structure image** or enter SMILES
* Upload your **database CSV**
* Select **columns** for:

  * Suggested method
  * Additional comments
  * Node hover metadata

### Main Tabs

* **Molecular Descriptors** – Show structure, descriptors, 3D info
* **Method Suggestion** – Suggest closest method if similarity ≥ threshold
* **Similarity Network** – Visual network of all structures above threshold

---

## Similarity Network

* Nodes = molecules (structure from SMILES)
* Edges = Tanimoto similarity above threshold
* Hover shows user-selected metadata from database
* Interactive plot using **Plotly**

---

## Example Database Format (CSV)

```csv
SMILES,Descriptive_Method,Additional_comments,Group
CCO,Method A,Fast elution,Alcohols
CCCC,Method B,High retention,Alkanes
...
```

---

## Future Ideas

* Export network to HTML/PNG
* Use molecular fingerprints other than Morgan
* Multi-image input or batch processing
* Predict class labels from structure
* SMILES standardization / tautomer handling

---

## 🧑‍Author

Developed by **Ricardo M. Borges**
Contact: [GitHub](https://github.com/RicardoMBorges) • [LinkedIn](https://www.linkedin.com/in/ricardomborges/) • [UFRJ](http://www.ippn.ufrj.br)

---

## 📄 License

MIT License

