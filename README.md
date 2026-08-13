# Sewage-to-Signal: Early Warning for Immune-Escape Mutations

**NIAID-BRCs AI Codeathon 2.0** · September 16–18, 2026 · Argonne National Laboratory

An AI-assisted wastewater surveillance workflow that detects pathogens of concern, tracks their prevalence, and flags emerging immune-escape and therapeutic-resistance mutations.

Project page: https://niaid-brc-codeathons.github.io/projects/sewage-to-signal/

---

> **This is a draft pitch, not a plan.**
>
> What follows is a one-slide proposal from the organizing team. It exists
> to seed a team, not to constrain one. Scope, methods, target organism,
> and success criteria are all still open — expect them to change
> substantially. Turning this into a real plan is the team's first job, and
> it lands in the project charter due August 28, 2026.

---

## Goal (proposed)

Develop an AI-assisted wastewater surveillance workflow to detect pathogens of concern, track their prevalence, identify emerging mutations, and assess potential immune-escape or therapeutic impact.

## Three-Day MVP (proposed)

Starting from target-capture wastewater sequencing, classify reads and identify NIAID pathogens of concern. Call protein mutations, track their geographic, temporal, and lineage distribution using public sequence repositories, and use literature RAG plus curated databases to identify mutations associated with immune escape, therapeutic resistance, or altered pathogenicity.

Generate a provenance-linked early-warning report highlighting high-priority mutations.

## Evaluation (proposed)

Taxonomic classification accuracy; concordance of mutation calls with standard pipelines; accuracy of geographic/temporal context; precision of literature-derived mutation–phenotype associations; and ability to recover known escape or resistance mutations as high-priority signals.

## Leads

- Alexander Taepper
- Andrew Warren

Team assignments are still being finalized. Participants can review their project, and request a reassignment, in the participant spreadsheet circulated by the organizing team.

## Working here

This repository is the team's working space for the codeathon — code, notebooks, data pointers, and notes. Replace this README with the real thing once the charter is written. Team members get access through the [NIAID-BRC-Codeathons](https://github.com/NIAID-BRC-Codeathons) organization; accept the invitation if you have not already.
