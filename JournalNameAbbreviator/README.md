# Journal Name Abbreviator

Single-file, offline web app that converts a full journal title to its
standard (ISO 4) abbreviation.

## Use

Open `JournalAbbreviator.html` in any browser (double-click). No internet
needed. Type or paste a journal name, pick a suggestion, click **Copy**.

- **Periods** checkbox toggles `J. Chromatogr. A` vs `J Chromatogr A`.
- **Auto-copy on select** copies the abbreviation as soon as you pick a match.
- Blue badge "from database" = exact entry from the embedded list
  (40,527 titles, merged from the JabRef abbreviation lists — ACS, general,
  life science, mathematics, mechanical, IEEE, PubMed/Entrez, etc.).
- Orange badge "rule-based" = the journal was not in the list and the
  abbreviation was generated from LTWA-derived word rules — double-check it.
- Green badge "custom entry" = one of your own additions or edits.

## Adding and editing entries

- **Edit** a result to correct or override its abbreviation (works on
  database entries too — your version takes priority).
- When a journal isn't found, use **＋ Save as entry** (on a rule-based
  result) or **＋ Add it as a new entry** to add it.
- **My entries (n)** at the bottom lists your custom entries, with edit and
  remove buttons. Removing an override restores the original database entry.
- Custom entries are saved in the browser (localStorage), so they persist on
  your machine but do not travel with the file.
- **⬇ Download shareable copy** produces a new single-file copy of the app
  with your custom entries built in — use this to share your curated version
  with colleagues or move it to another computer.

## Sharing

The HTML file is fully self-contained (~2.8 MB). Email it or share the
Drive link; recipients just open it in a browser.

## Rebuilding the data

Data was merged from https://github.com/JabRef/abbrv.jabref.org (journals/
CSV lists) with dotless MEDLINE abbreviations converted to dotted ISO 4
style by aligning abbreviation tokens with the full title. To refresh,
re-download the CSVs and re-run the merge (see project chat history,
July 2026).
