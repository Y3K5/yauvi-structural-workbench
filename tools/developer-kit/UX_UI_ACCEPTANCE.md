# UX/UI acceptance review

Status starts **Not tested** for every reviewer task. Record Pass, Fail, Blocked, or Not tested; include OS, browser/version, viewport, input case and evidence. Automated checks do not substitute for these tasks.

## Priority journeys

| Task | Expected behavior |
|---|---|
| First visit | A newcomer can explain the product and find a next action without reading commands. |
| Open guided showcase | Welcome-page button reaches its embedded tour; no missing sibling file or blank frame. Test both file opening and the local preview server. |
| Six scientific questions | All six panels load, selection matches the URL, arrows/Home/End operate the tabs, and layout does not jump unexpectedly. |
| Conformational state | Region/residue selection, slider, play/pause, reset and geometry download work. Endpoint colors are explained; interpolation is explicitly not a physical simulation. |
| Two comparative studies | Both fold-family examples identify species, source models, research citations and limits. Similar fold is not presented as proof of identical function. |
| Complete synthetic case | Create example → check readiness → run → findings → evidence → save report. Each step explains what happened and what to do next. |
| Missing-evidence case | Omitting validation visibly reports the gap. It remains distinguishable from a software crash and from validated science. |
| Input mistakes | Missing file, wrong format, unsupported input, stale revision and unavailable runtime provide actionable messages, preserve inputs and avoid false success. |
| Evidence and export | Measurements trace to sources. Exported report reopens and preserves limitations. No private filesystem paths should appear in a shareable artifact without review. |
| Feedback | Welcome form downloads locally, preserves entered text, and sends nothing online. Unsaved-feedback behavior is clear. |

## Visual and accessibility matrix

- Desktop: 1440 × 900 and 1280 × 720. Narrow view: 390 × 844. Text zoom: 200%.
- Chrome, Firefox, Safari and Edge as available; record untested combinations explicitly.
- Keyboard-only navigation, visible focus, logical order, no dialog focus trap failures, Escape/close behavior, labels and useful screen-reader announcements.
- Check text contrast (WCAG AA: 4.5:1 normal text; 3:1 large text) and non-text control contrast. Do not rely on color alone for state.
- Long filenames, long scientific names, empty library, incomplete evidence, loading, failure, disabled controls and repeated clicks.
- Reduced-motion preference; readable tables without overlapping controls; usable molecule controls on touch devices.
- Offline showcase: confirm no required network requests. External research citations may require internet when explicitly opened.
- Fresh install, reopening, port collision, missing Python/dependency, stopped server and outdated asset behavior. Verify Windows and Linux launchers on real systems.

## Release decision

A guided pilot may proceed with explicitly recorded minor issues. Block release for inaccessible primary actions, broken navigation, lost inputs, incorrect evidence labels, private-data exposure or false biological claims. Full scientific readiness across arbitrary inputs requires separate workflow qualification; UX sign-off cannot establish it.

## Suggested evaluation

Give a newcomer 20 minutes: find an example; explain one result and one limitation; inspect a moving residue; find one research source; save feedback. Record success without help, success with help, or blocked. This is a proposed protocol, not evidence of successful user testing.
