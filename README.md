# Zeus-256 — Supplementary Materials

This repository contains the supplementary data, quantum circuits, analysis code, and figures accompanying:

**Experimental Benchmarking of a 256-Qubit Room-Temperature NV-Diamond Quantum Processor with an Optically Programmable Biomolecular Quantum Bus**

**Adams Ivanov, Daniela Herrmann, and Erick Giovani Sperandio Nascimento**  
Dynex Quantum Technologies AG, Switzerland; University of Surrey, Great Britain  
2026

Repository: [github.com/DynexQuantum/zeus_paper](https://github.com/DynexQuantum/zeus_paper)

## Contents

| Material | Description | Manuscript section |
|---|---|---|
| Shor QASM circuits | Complete circuit files for all 11 Shor-family factorization instances, covering 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, and 84-bit integers and 28–256 processor-level logical qubits. The implementation uses Gidney–Ekerå-inspired Ekerå–Håstad circuits with windowed exact arithmetic. | 6.6 |
| Shor measurement summary | Exact test integers, qubit and gate counts, attempt counts, reported execution times, and factor-recovery flags. | 6.6, Table 6b |
| GHZ-256 measurement data and acquisition note | Computational-basis counts, all 257 phase-resolved parity settings, shot allocations, full-precision phases, fidelity results, and the acquisition summary. | 6.3.1 |
| Benchmark data and analysis | Supporting circuit and timing records, Bell–CHSH and Grover measurement results and analyses, and supplementary figure-generation materials. | Sections 5–6 |

The Shor staircase records 10 attempts per instance and factor recovery for all 11 test integers. Its largest instance is the 84-bit integer `16503807113155763618449753`, using 256 qubits and 40,250,342 reported gates, with a reported execution time of 551.694 ms. The manuscript defines the scope and interpretation of these measurements.

## Using the materials

Consult the manuscript for experimental procedures, circuit descriptions, statistical estimators, and timing definitions. QASM files provide the circuit definitions; execution on a particular backend requires the corresponding compilation and control environment.

For GHZ numerical reproduction, use the full-precision phase values in the JSON data rather than rounded values from tables or the acquisition note. When reporting analyses based on this repository, record the commit hash or release used so that the exact material can be identified.

## Cite the paper

If you use these data, circuits, code, or figures, please cite the accompanying manuscript:

> Ivanov, A., Herrmann, D., & Nascimento, E. G. S. (2026). *Experimental Benchmarking of a 256-Qubit Room-Temperature NV-Diamond Quantum Processor with an Optically Programmable Biomolecular Quantum Bus*. Manuscript.

```bibtex
@unpublished{ivanov2026zeus256,
  author = {Ivanov, Adams and Herrmann, Daniela and
            Nascimento, Erick Giovani Sperandio},
  title  = {Experimental Benchmarking of a 256-Qubit Room-Temperature
            {NV}-Diamond Quantum Processor with an Optically Programmable
            Biomolecular Quantum Bus},
  year   = {2026},
  note   = {Manuscript}
}
```

This citation identifies the manuscript version. Use the final publication citation and DOI when available.

## Cite the supplementary materials

Please also identify the repository and the commit or release used:

> Dynex (2026). *Supplementary materials for “Experimental Benchmarking of a 256-Qubit Room-Temperature NV-Diamond Quantum Processor with an Optically Programmable Biomolecular Quantum Bus”*. GitHub. https://github.com/DynexQuantum/zeus_paper

## Contact

Scientific correspondence: [office@dynex.ch](mailto:office@dynex.ch).
