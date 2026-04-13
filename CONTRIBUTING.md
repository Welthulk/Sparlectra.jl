# Contributing to Sparlectra.jl (Strict Governance Model + QA Requirements)

Thank you for your interest in contributing to **Sparlectra.jl**.

This project is primarily a **research and educational effort** with limited maintainer capacity.  
To keep the codebase consistent and maintainable, contributions follow a stricter process.

---

## ⚠️ Mandatory Proposal First Policy

**All non-trivial contributions require prior discussion.**

Before writing code, you MUST:

1. Open an Issue or Discussion
2. Describe your idea clearly
3. Wait for maintainer feedback and alignment

➡️ Pull Requests without prior agreement are typically not reviewed.

---

## Contribution Scope

Sparlectra prioritizes:

- numerical robustness  
- deterministic behavior  
- conceptual clarity  
- minimal complexity  

---

## Development Workflow

1. Fork the repository  
2. Create a feature branch  
3. Implement agreed changes  
4. Open a Pull Request  

---

## Coding Guidelines

- clear and structured code  
- descriptive naming  
- minimal dependencies  
- deterministic algorithms  

All comments and docstrings must be in **English**.

---

## Numerical Code Requirements

- preserve numerical stability  
- avoid performance regressions  
- document mathematical changes clearly  
- justify algorithmic decisions  

---

## Testing Requirements (Mandatory)

- Every change MUST include **dedicated unit tests**
- Tests must be **deterministic**
- Edge cases must be covered where relevant
- Tests must be part of the repository (not only described in PR)

---

## Coverage Requirements

- New code should include **test coverage**
- Contributions should not significantly reduce overall coverage
- Critical numerical paths must be explicitly tested

(Optional but recommended)
- Contributors may include coverage reports (e.g. via CI)

---

## Benchmark Requirements (for Solver / Numerical Changes)

For any change affecting:

- solver behavior  
- numerical routines  
- performance-critical code  

the contributor MUST provide:

- a **before/after benchmark**
- description of the test case (network size, type)
- comparison of:
  - runtime  
  - iteration count  
  - convergence behavior  

Benchmarks must be:

- reproducible  
- documented  
- runnable by the maintainer  

---

## Test & Validation Documentation

Each contribution must include:

- how to run tests  
- expected results  
- assumptions / limitations  

PR-only descriptions are NOT sufficient.

---

## Pull Request Acceptance Policy

A PR may be rejected if:

- no prior discussion  
- missing tests  
- missing benchmarks (if required)  
- insufficient documentation  
- increased complexity without benefit  

---

## Maintainer Availability

- Reviews are best-effort  
- Response times may vary  
- Not all PRs will be accepted  

---

## Summary

- Discuss first  
- Provide tests  
- Provide benchmarks (if relevant)  
- Document clearly  
- Keep changes minimal  

This ensures long-term maintainability.

Thank you for your understanding.
