# Cryptographic Status

Status: **DESIGN ONLY / NOT IMPLEMENTED**

The repository contains the historical design document:

`docs/CRYPTOGRAPHIC_BOUNDARY_V0_9.md`

Its `v0.9` label predates the 2026-10-05 repository reconciliation and collides with the Origins-DNA research-stage v0.9 namespace.

## Current interpretation

The historical document specifies a candidate ChaCha20-Poly1305 authenticated-encryption layer and a DNA transport encoding, but current main does **not** contain:

- an `origins_dna.crypto` implementation;
- seal/open APIs;
- key-generation code;
- AEAD envelope serialization;
- dedicated cryptographic tests.

Therefore the document is a design freeze candidate, not a software release or validated security claim.

## Namespace rule

Future cryptographic work must use a separate profile namespace such as:

`AEAD-PROFILE-0.1`

Research versions `v0.x` remain reserved for the Origins-DNA representation/benchmark lineage.

## Security boundary

When implemented, cryptographic strength must come from the selected established AEAD primitive and its key/nonce contract. DNA representation, PhaseNav, DHSE, RDSS or On-Primes mathematics must not be presented as a substitute for conventional cryptographic security.

This file supersedes only the **current-status interpretation** of the legacy filename. It does not delete or rewrite the original design document.
