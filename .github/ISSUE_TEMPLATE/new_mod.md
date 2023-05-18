---
name: New Mod Proposal
about: Suggest a new Mod to add
title: ''
labels: 'mod proposal'
assignees: ''

---

**Why this Mod?**

Justify why this Mod is a useful addition.

**Design requirements**

Verify the mod will meet the following criteria:

- [ ] Self-contained: the user doesn't need to interact directly with gurobipy
  objects when using the function or class exposed by the mod.
- [ ] Stateless: data-in data-out, no follow-up optimization steps, and gurobipy
  objects are disposed as soon as solutions are retrieved.
- [ ] Covers a well-known concept from a non-optimization field. In particular,
  we must be able to explain what the Mod does using domain-specific terminology
  without reference to the mathematical model.
- [ ] Clear purpose / single focus. There may be natural parameters to switch on/off
  various side-constraints, but these should be kept to a minimum.

**What will the API be?**

Add a code snippet for the proposed API. Does this link cleanly to:

- [ ] numpy
- [ ] scipy
- [ ] pandas
- [ ] Other (please specify)

to provide input data to the Mod?

**Additional context**

Any other details that should be communicated in the Mod writeup to help
users/readers understanding?
