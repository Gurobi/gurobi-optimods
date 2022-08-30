---
name: New Nup Proposal
about: Suggest a new Nup to add to the Stup
title: ''
labels: ''
assignees: ''

---

**Why this Nup?**
Justify why this Nup is a useful addition. Criteria:

- Self-contained: the user doesn't need to interact directly with gurobipy objects when using the function/class exposed by nupstup
- Solve model, get solution (no follow up analysis, e.g. IIS)
- Stateless: data-in data-out, no follow-up optimization steps
  - Model can be kept in memory for the purpose of returning solutions
- Covers a well-known concept from a non-optimization field
  - In particular, we must be able to explain what the Nup does using domain-specific terminology
- Clear purpose / single focus
  - There may be natural parameters to switch on/off various side-constraints, but these should be kept to a minimum

**Does it fall under an existing category?**

- Machine Learning
- Classical OR
- Graphs
- Other (please specify)

**What will the API be?**
Add a code snippet for the proposed API. Does this link cleanly to:

- numpy
- scipy
- pandas
- Other (please specify)

to provide input data to the Nup?

**Additional context**
Any other details that should be communicated in the Nup writeup to help users/readers understanding?
