## **[中文版本](https://www.misaraty.com/2026-08-04_caia/)**

## Chemistry Artificial Intelligence Assistant (CAIA)

CAIA is a demonstration platform for teaching physical chemistry theory and laboratory courses. Each module is provided as an independent Python script that can be run and modified separately. The experimental data, assignment texts, learning records, and knowledge-base content included in the project are for demonstration purposes.

## Module Overview

### 1. Sucrose Hydrolysis Kinetics Analysis

```bash
python 1_sucrose_hydrolysis/caia_sucrose_hydrolysis.py
```

The script includes 11 sets of time and optical rotation data collected at 30.0 °C. It compares linearized fitting, nonlinear fitting, and symbolic regression. The outputs include the rate constant, half-life, R2, RMSE, MAE, SSE, maximum absolute residual, point-by-point predictions, residuals, and figures.

By default, the symbolic-regression result uses an explicit expression from an existing laboratory report to ensure stable reproduction. Set `USE_PYSR` to `True` to run a new expression search.

### 2. Ethyl Acetate Saponification Kinetics Analysis

```bash
python 2_saponification/caia_saponification.py
```

The script includes 10 sets of time and conductivity data at both 28.0 °C and 38.0 °C, with an initial concentration of 0.0100 mol/L. It outputs the rate constants obtained by the three methods at both temperatures, fitting metrics on the original conductivity scale, and the activation energy calculated by the two-temperature method.

By default, the symbolic-regression result uses an explicit expression from an existing laboratory report to ensure stable reproduction. Set `USE_PYSR` to `True` to run a new expression search.

### 3. RAG Course Knowledge Base

```bash
python 3_rag_course_assistant/caia_rag_course_assistant.py
python 3_rag_course_assistant/caia_rag_course_assistant.py --question "What is the relationship between the half-life of a first-order reaction and the initial concentration?"
```

The script retrieves evidence using character-level TF-IDF, BM25, and reaction-order keyword correction. It outputs the vector score, BM25 score, base hybrid score, topic-correction value, final score, selected evidence, and generated answer.

In offline mode, the answer is organized only from the retrieved course snippets. When DeepSeek is enabled, the program sends the filtered evidence and the question to the API, and the model then reformulates the response. Numerical calculations and retrieval scores are still produced locally.

### 4. Formula Derivation and Exercises

```bash
python 4_formula_derivation/caia_formula_derivation.py
```

The script uses SymPy to verify the integrated rate equations for zero-order, first-order, and equal-initial-concentration second-order reactions. It also generates formula explanations, exercises, and illustrative curves.

### 5. Analysis of Theory Assignments and Laboratory Reports

```bash
python 5_homework_report_feedback/caia_homework_report_feedback.py
```

The script includes a demonstration assignment, a demonstration laboratory report, and a keyword-based rubric. Scores are calculated from keyword coverage and are intended only to demonstrate rule-based scoring and formative-feedback workflows. The program cannot determine the scientific correctness of a student's statements solely from keyword occurrence, nor can it replace teacher evaluation.

DeepSeek is used only for optional language polishing. Deterministic scores and the original rule-based feedback are retained in the JSON output. When the API is enabled, student text is sent to an external service. Names, student IDs, and other identifiable information must therefore be removed before formal use.

### 6. Learning Diagnosis

```bash
python 6_learning_diagnosis/caia_learning_diagnosis.py
```

The script outputs pre-test mastery, post-test mastery, average gain, average normalized gain, weak knowledge points, and tiered recommendations at both the student and knowledge-point levels. The built-in records are intended only to demonstrate data aggregation and figure generation.

### 7. Knowledge Games

Demonstration mode:

```bash
python 7_knowledge_games/caia_knowledge_games.py
```

Interactive command-line mode:

```bash
python 7_knowledge_games/caia_knowledge_games.py --interactive
```

By default, the script uses preset demonstration answers and outputs the overall accuracy, category-level accuracy, item-by-item results, and a figure.

### 8. Integrated Web Platform

```bash
python 8_integrated_web_platform/caia_integrated_web_platform.py
```

Open `http://127.0.0.1:7860` in a browser. This page is a simplified single-file demonstration that includes nonlinear kinetic fitting, local course-knowledge retrieval, rule-based feedback, learning-analysis figures, and knowledge games. It does not fully reuse all algorithms from the first seven modules and does not call DeepSeek.

For local-area-network access, set the following environment variable.

Windows PowerShell:

```powershell
$env:CAIA_SHARE="1"
python 8_integrated_web_platform/caia_integrated_web_platform.py
```

Linux or macOS:

```bash
export CAIA_SHARE=1
python 8_integrated_web_platform/caia_integrated_web_platform.py
```

In local-area-network mode, the program listens on all network interfaces. Basic access authentication can be enabled by setting both `CAIA_USERNAME` and `CAIA_PASSWORD`. Do not enter real student information or other sensitive content when authentication is not configured.

## Citation

To be added after the paper is officially published.
