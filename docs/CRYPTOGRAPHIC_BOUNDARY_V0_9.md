# Origins-DNA v0.9 — Cryptographic Boundary

Status: design freeze candidate
Branch: `feat/dna-aead-v0.9`

## Objective

Add an authenticated-encryption layer for DNA payloads while preserving the existing reversible 2-bit DNA representation and the repository's append-only research discipline.

The security target is a ciphertext representation for which a strong structural analyser, including a future PhaseNav adversarial harness, obtains no reproducible plaintext-recovery or key-recovery advantage over declared controls.

## Separation of responsibilities

The v0.9 design keeps three layers explicit:

1. **DNA representation** — canonical A/C/G/T normalization and reversible 2-bit packing already provided by `origins_dna.encoding`.
2. **Cryptographic protection** — an established authenticated-encryption primitive with a 256-bit secret key and explicit nonce handling.
3. **Adversarial validation** — conventional statistical tests plus an optional external PhaseNav attacker adapter that receives ciphertext only.

DNA structure may be used for transport, storage and controlled experiments. Cryptographic strength is assigned to the authenticated-encryption primitive and key-management contract.

## Initial primitive profile

The first implementation candidate uses ChaCha20-Poly1305 with:

- 256-bit uniformly random key;
- 96-bit nonce generated from a cryptographically secure random source;
- authenticated associated data binding protocol version and representation metadata;
- 128-bit authentication tag supplied by the AEAD primitive;
- explicit envelope versioning for future migrations.

The implementation should expose key bytes directly rather than inventing a password-derived key scheme in this version. Password-based key derivation, if required later, receives a separate frozen protocol and test suite.

## Binary-to-DNA transport

Arbitrary ciphertext bytes can be represented exactly as DNA by a fixed reversible mapping of every two bits to one base:

- `00 -> A`
- `01 -> C`
- `10 -> G`
- `11 -> T`

This transport mapping is deterministic and reversible. Envelope length and version metadata are authenticated so parsing cannot silently reinterpret ciphertext.

## Envelope contract

A sealed record should bind at least:

- schema/version identifier;
- AEAD algorithm identifier;
- nonce;
- ciphertext plus authentication tag;
- encoded plaintext nucleotide length;
- optional caller-supplied associated data.

Parsing and authentication errors fail closed.

## PhaseNav adversarial gate

A later validation harness should treat PhaseNav as an attacker model rather than as part of encryption.

For frozen plaintext families and independently generated keys/nonces, compare PhaseNav against matched controls on tasks such as:

- distinguishing ciphertexts produced from different plaintext families;
- recovering plaintext-derived labels above chance;
- ranking candidate plaintexts;
- exploiting nonce, length or representation metadata;
- detecting stable cross-message structure that survives key and nonce changes.

Primary cryptanalytic validation target:

`attacker_advantage = observed_success - matched_null_success`

Promotion requires the confidence interval for attacker advantage to remain inside a preregistered negligible band across declared message families, lengths and independent keys.

## Required implementation tests

The first implementation must include:

- exact encrypt/decrypt round-trip for empty, short and long DNA payloads;
- exact binary-to-DNA transport round-trip;
- rejection of malformed keys and nonces;
- authentication failure after any ciphertext mutation;
- authentication failure after associated-data mutation;
- distinct ciphertexts under fresh nonces for the same key/plaintext;
- deterministic parsing of a frozen envelope;
- no secret-key material in serialized records, logs or exceptions.

## Repository discipline

Implementation changes remain on `feat/dna-aead-v0.9` until explicitly promoted. Existing scientific results, legacy records and frozen benchmark artifacts remain byte-preserved. Cryptographic software verification and biological claims remain separate evidence tracks.

## Next implementation unit

Create a narrow `origins_dna.crypto` API providing:

- secure key generation;
- DNA sealing;
- authenticated opening;
- binary ciphertext <-> DNA transport conversion;
- versioned envelope serialization;
- fail-closed validation errors.

A dedicated test module should be added before any benchmark or performance claim is recorded.
