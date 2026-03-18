from pathlib import Path
path = Path(r"c:\Users\cw1245\Documents\PhD\Callum_Scripts\LUTS.R")
lines = path.read_text(encoding='utf-8').splitlines(keepends=True)
start_lines = [i for i,l in enumerate(lines) if l.startswith('# NEW FUNCTION: logreg_table')]
if len(start_lines) < 2:
    raise SystemExit(f"Expected >=2 logreg_table definitions, found {len(start_lines)}")
first = start_lines[0]
second = start_lines[1]
del lines[first:second]
path.write_text(''.join(lines), encoding='utf-8')
print('Removed duplicate logreg_table definition; kept the later one.')
