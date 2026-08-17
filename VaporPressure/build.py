"""build.py - inline vp_data.json into the template to produce VaporPressure.html."""
import pathlib
import sys

here = pathlib.Path(__file__).parent
tpl = (here / "vp_template.html").read_text(encoding="utf-8")
data = (here / "vp_data.json").read_text(encoding="utf-8")

start = tpl.index("/*__VP_DATA__*/")
end = tpl.index("/*__END__*/") + len("/*__END__*/")
out = tpl[:start] + data + tpl[end:]

dest = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else here / "VaporPressure.html"
dest.write_text(out, encoding="utf-8")
print(f"{dest}  ({len(out)/1024:.0f} kB)")
