"""Rebuild UI around frozen public showcase data; no scientific recomputation."""
from pathlib import Path
import html

SOURCE = Path(__file__).resolve().parents[1]
ROOT = SOURCE.parent
s = SOURCE / 'showcase'
template = (s/'template.html').read_text()
style = (s/'style.css').read_text()
app = (s/'app.js').read_text()
viewer = (SOURCE/'software/structural-workbench/src/yauvi_structural_workbench/ui/vendor/3Dmol-min.js').read_text().replace('</script', '<\\/script')
standalone = template.replace('__SHOWCASE_STYLE__', style).replace('__SHOWCASE_VIEWER__', '<script>'+viewer+'</script>').replace('__SHOWCASE_APP__', '<script>\n'+app+'</script>')
out = ROOT/'preview/showcase'
(out/'index.html').write_text(standalone)
workbench = template.replace('__SHOWCASE_STYLE__', style).replace('__SHOWCASE_VIEWER__', '<script src="/vendor/3Dmol-min.js"></script>').replace('__SHOWCASE_APP__', '<script src="/showcase/showcase.js"></script>')
workbench = workbench.replace('<header>', '<header><p><a href="/">← Back to your workbench</a></p>', 1)
(out/'workbench/index.html').write_text(workbench)
(out/'workbench/showcase.js').write_text(app)
(ROOT/'preview/START_HERE.html').write_text((s/'welcome-template.html').read_text().replace('__EMBEDDED_SHOWCASE__',html.escape(standalone, quote=True)))
print('Rebuilt standalone, embedded and workbench showcases from frozen public examples.')
print('Preview files changed; original package checksums no longer describe these edited previews.')
