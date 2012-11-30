# -*- mode: python -*-

# A specification file for PyInstaller which creates a
# distribution directory with an executable for Windows.
#
# Run with:
#
#    python /path/to/pyinstaller.py qt_win.spec
#
# This will produce a `dist' directory which contains all the files.

a = Analysis(
    ['cpg_islands/qt/main.py'],
    pathex=['.'],
    hiddenimports=[],
    hookspath=None)
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    exclude_binaries=1,
    name=os.path.join('build', 'cpg_islands.exe'),
    debug=False,
    strip=None,
    upx=True,
    console=False)
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=None,
    upx=True,
    name='dist')
