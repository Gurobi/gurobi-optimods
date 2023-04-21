---
name: New Mod Proposal
about: Suggest a new Mod to add
title: ''
labels: 'mod contribution'
assignees: ''

---

**Why this Mod?**
Justify why this Mod is a useful addition. Criteria:

- Self-contained: the user doesn't need to interact directly with gurobipy objects when using the function/class exposed by gurobi-optimods
- Solve model, get solution (no follow up analysis, e.g. IIS)
- Stateless: data-in data-out, no follow-up optimization steps
  - Model can be kept in memory for the purpose of returning solutions
- Covers a well-known concept from a non-optimization field
  - In particular, we must be able to explain what the Mod does using domain-specific terminology
- Clear purpose / single focus
  - There may be natural parameters to switch on/off various side-constraints, but these should be kept to a minimum

**What will the API be?**
Add a code snippet for the proposed API. Does this link cleanly to:

- numpy
- scipy
- pandas
- Other (please specify)

to provide input data to the Mod?

**Additional context**
Any other details that should be communicated in the Mod writeup to help users/readers understanding?
