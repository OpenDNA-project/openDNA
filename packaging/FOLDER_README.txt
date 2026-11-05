OpenGenEdit, folder version
==========================

Run OpenGenEdit.exe in this folder. Keep the folder together: the program needs the _internal folder next to it.
Free software (AGPL-3.0). Download it only from https://github.com/Ant-Champ/GenomeEdit/releases . If you paid for it, you were scammed.

Why this version exists
-----------------------
OpenGenEdit uses Qt (through PySide6), which is under the LGPL-3.0. The LGPL lets you run the program with your own build of Qt.
The single-file OpenGenEdit.exe packs Qt inside itself; in this folder version Qt is made of plain files, so you can swap them.

How to use another Qt (another PySide6 version)
-----------------------------------------------
1. Get the wheels of the version you want (they are zip files): PySide6-Essentials and shiboken6, for Windows 64-bit, Python 3.13
   (the version this program was built with: see THIRD_PARTY_NOTICES.md). For example:
       pip download PySide6-Essentials==6.10.1 shiboken6==6.10.1 --no-deps -d wheels
   or take your own build of Qt for Python.
2. Close OpenGenEdit. In _internal, rename the folders  PySide6  and  shiboken6  (for example to PySide6.orig) so you can go back.
3. Open each wheel as a zip file (rename .whl to .zip if needed) and copy its  PySide6  folder and its  shiboken6  folder into _internal.
4. Start OpenGenEdit.exe. Help > About shows the Qt version in use.
If it does not start, put the original folders back (or unzip the download again). A version far from the one it was built with may not work.

This folder version does not update itself. Download the new OpenGenEdit-folder.zip from the releases page when there is one.

Licences and the source of everything in here: THIRD_PARTY_NOTICES.md and LICENSE (next to this file). Source code: https://github.com/Ant-Champ/GenomeEdit
