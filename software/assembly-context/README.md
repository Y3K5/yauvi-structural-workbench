# AssemblyContext

AssemblyContext measures contacts and surface burial in a supplied biological
assembly. It consumes a checksum-bound StructQC manifest and refuses to infer an
assembly from a structure title. PDB and mmCIF inputs are supported; an mmCIF
asymmetric unit that declares unapplied assembly operators is blocked until an
expanded biological assembly is supplied.

```bash
assembly-context run \
  --manifest STRUCTURE_EVIDENCE.json \
  --isolated query.pdb --assembly assembly1.cif \
  --subject-chain A --relationship exact_protein --out out
```

Contact geometry, buried solvent-accessible surface area, and occlusion are
evidence from one coordinate state. They are not native exposure, binding, or
functional claims.
