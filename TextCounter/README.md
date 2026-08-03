# Text Counter

Single-file HTML app: paste or type text, get live counts.

## Counts

Main tiles (update on every keystroke):

- **Words** — whitespace-separated tokens (`\S+`).
- **Characters (including spaces)** — all Unicode code points except
  line breaks (matches the Microsoft Word "Characters (with spaces)"
  convention, since journal/proposal limits are usually checked
  against Word).
- **Characters (excluding spaces)** — all non-whitespace code points.

Secondary row: sentences (terminal `.!?…` followed by space/end),
paragraphs (blank-line separated), lines, and reading time (~200 wpm).

Selecting part of the text counts only the selection (a note appears
while a selection is active).

## Notes

- Counting is done in Unicode code points, so an emoji or `Ω` counts
  as one character.
- "Paste from clipboard" uses the async Clipboard API; if the browser
  blocks it, a hint suggests Ctrl/Cmd+V (which always works).
- Self-contained, works offline, no external requests.
