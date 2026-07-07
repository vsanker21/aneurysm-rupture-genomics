@echo off
REM Batch script to run all individual generation scripts in Command Prompt
REM This avoids PowerShell output buffering issues

echo ================================================================================
echo RUNNING ALL INDIVIDUAL OUTPUT GENERATION SCRIPTS
echo ================================================================================
echo.

cd /d "E:\Aneurysms Genomics Project"

echo [1/4] Generating risk stratification plot...
python generate_risk_stratification.py
if errorlevel 1 (
    echo ERROR: Risk stratification script failed
) else (
    echo SUCCESS: Risk stratification complete
)
echo.

echo [2/4] Generating gene lists...
python generate_gene_lists.py
if errorlevel 1 (
    echo ERROR: Gene lists script failed
) else (
    echo SUCCESS: Gene lists complete
)
echo.

echo [3/4] Generating feature importance plot...
python generate_feature_importance.py
if errorlevel 1 (
    echo ERROR: Feature importance script failed
) else (
    echo SUCCESS: Feature importance complete
)
echo.

echo [4/4] Generating pathway enrichment plots...
python generate_pathway_plots.py
if errorlevel 1 (
    echo ERROR: Pathway plots script failed
) else (
    echo SUCCESS: Pathway plots complete
)
echo.

echo ================================================================================
echo ALL SCRIPTS COMPLETED
echo ================================================================================
pause

