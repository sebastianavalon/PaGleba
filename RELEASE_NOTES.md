# PaGleba 1.0.0

Portable Linux/WSL ONT amplicon workflow, with English documentation and output labels.

- CLI flags plus project JSON and explicit sample/metadata/task tables support arbitrary paths, barcode selections, technical runs and sample counts.
- EUKARYOME SSU and ITS are separate marker-specific references; plants remain by default. Generic SSU and ITS profiles require explicit assay length limits. The `18s-f566-r1200` profile supplies the historical example window.
- Hash-verified classification chunks support stopping and resuming. Functional adapters, ecological statistics, per-task networks/PLS-PM and vector figures expose evidence and design limitations.
- Final installed wheel: 27 tests passed, none skipped in the full validation environment. CLI and source/wheel checks passed.

Start with `README.md`, then `docs/USER_GUIDE.txt`. Actual real ITS controls/pilot validation and a clean environment solve on another computer remain operator validation steps. No automatic GitHub publication occurs.
