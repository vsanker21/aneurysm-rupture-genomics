# GTEx File Requirements Checklist

## Required Vessel Tissues (3 files minimum)

### ✅ Required:
1. **Artery - Aorta** (Artery_Aorta)
2. **Artery - Coronary** (Artery_Coronary)
3. **Artery - Tibial** (Artery_Tibial)

### 📋 Optional but Recommended:
4. **Sample Metadata** - Information about samples (age, sex, etc.)
5. **Gene Annotation** - Gene IDs, names, biotypes
6. **Subject Phenotypes** - Donor information

## Expected File Formats

### Expression Data:
- **GCT format** (Gene Cluster Text) - Most common
  - Extension: `.gct` or `.gct.gz`
  - Contains: Gene IDs, sample IDs, expression values (TPM or counts)
  
- **TSV/CSV format**
  - Extension: `.txt`, `.tsv`, `.csv` or compressed `.gz`
  - Contains: Tab/comma-separated expression matrix

### Metadata:
- **Text/TSV format**
  - Sample attributes, subject phenotypes, etc.

## File Naming Patterns

GTEx files typically contain:
- Tissue name: `Aorta`, `Coronary`, `Tibial`
- Data type: `expression`, `TPM`, `counts`, `metadata`
- Version: `v8`, `v9`
- Format: `GCT`, `TSV`, `CSV`

## What You Should Have (6 files)

If you downloaded 6 files, they likely include:

1. **Artery_Aorta expression** (GCT or TSV)
2. **Artery_Coronary expression** (GCT or TSV)
3. **Artery_Tibial expression** (GCT or TSV)
4. **Sample attributes/metadata** (TSV or TXT)
5. **Gene annotation** (TSV or TXT)
6. **Subject phenotypes** (TSV or TXT)

OR

1-3. Three vessel tissue expression files
4-6. Three additional vessel tissues or other data types

## Verification Steps

1. **Check file names** - Do they contain "Aorta", "Coronary", "Tibial"?
2. **Check file sizes** - Expression files are typically 10-100 MB each
3. **Check file formats** - GCT files start with `#1.2`
4. **Check location** - Files should be in `data/gtex/` directory

## If Files Are Missing

If you're missing any of the 3 required vessel tissues:
- Return to GTEx Portal: https://gtexportal.org/home/datasets
- Download the missing tissue(s)
- Ensure you download the same format (TPM or counts) for all tissues

## Next Steps

Once you have the 3 required vessel tissue files:
1. Move them to `data/gtex/` directory
2. Verify they can be opened/read
3. Check that they contain expression data (not just metadata)
4. Ensure consistent format across all files






