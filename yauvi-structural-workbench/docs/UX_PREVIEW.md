# A clearer local workbench

This private development preview organizes each analysis into **Inputs → Process → Findings → Evidence**. Start with **New**, or choose **Try an offline example → Load example**. Find previous work by subject, question, workflow or case ID in the library.

**Inputs** separates required files from additional evidence. Each file explains why it matters and exposes its content checksum on demand. Parameters remain available behind a disclosure. Check readiness before running; missing inputs and unavailable tools explain why execution is blocked. Saving parameters preserves the earlier case revision.

**Process** shows the actual worker stages, elapsed time, queue position and cancellation. These are execution events, not an estimated completion percentage. Events persist across restarts. Older jobs without timestamps are labeled accordingly. A run that stops unexpectedly remains interrupted; restarting does not silently resubmit it.

**Findings** keeps the interpretation limit beside the measurements. StructQC summary cards link to their exact engine document. All workflows expose their method documents and measurement tables; table previews disclose their limits and full files remain in the evidence bundle. The local molecular viewer supports coordinate selection, chain and residue selection, three display styles, non-polymer visibility and view reset.

**Evidence** traces the selected run through inputs, method, measurements, interpretation and qualification. The run selector reopens historical results without replacing the editable case. Results from earlier input or parameter revisions show a warning. Report, data, full evidence, checksums and reproduction-manifest downloads are individually described.

Background progress updates preserve unsaved parameters. Case selection, upload, reference adoption and run selection retain their intended case identity. Public retrieval stays off by default and the interface explains how to opt in. A source that only supports local evidence says so.

## Verification on 6 September 2026

The browser checks exercised:

- All six workflow choices, creation forms, required input roles and available parameter controls.
- Offline example loading, readiness checks and successful StructQC execution.
- A manual coordinate upload with checksum confirmation.
- Browser cancellation, preserved cancellation evidence and a deliberate retry that correctly ended scientifically incomplete.
- Saved parameter revisions, warning on earlier evidence, historical-run selection and protection of unsaved edits during navigation.
- Chain/residue selection, sticks, non-polymer visibility, view reset, structure loading and evidence navigation.
- Search, activity, dialogs, download activation and an actual repaired template download.
- Desktop and compact-window visual review. The observed compact viewport was 694 × 784 pixels with no page-width overflow. This does not establish coverage of every mobile browser or screen size.

The endpoint check verified served assets against canonical source bytes, five report downloads, and every template role (20 checks total). It discovered and repaired a real defect: pre-encoded JSON template bytes were mistakenly serialized a second time, returning HTTP 400. A regression test now covers all template definitions. The archive screen found no matches for its declared local-path, username, private-project or secret-shaped patterns; it is not an exhaustive privacy certification.

See the version-bound [verification record](../implementation-evidence/2026-09-06/UX_VERIFICATION.json) and [preview artifacts](../artifacts/local-preview-20260906-ux/).

## Remaining acceptance work

All six workflows have accessible controls; end-to-end scientific execution through the revised browser was exercised for StructQC only. Live public-source retrieval was deliberately not enabled or tested. The broader six-workflow equivalence matrix, external usability pilots, independent reproduction and scientific qualification remain open. Membrane orientation stays experimental; SF-CSA qualification still requires recovery of its full archived reference universe.

The visual and software checks do not establish biological validity, publication readiness or usability acceptance by independent researchers. No artifact is authorized for publication.
